#!/usr/bin/env python
#
# Original filename: centeroflight.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: April 2011
# 
# Summary:  Find the center of light of an input image, fit an atmospheric
# dispersion model using unsaturated frames.  Output the fitted centroids.
# 

import numpy as np
import parallel

def fit_centroids(filesetup, flux, method='moffat', storeall=True):

    """
    Function fit_centroids runs the centeroflight centroiding routine
    in parallel on a list of HiCIAO frames and fits the unsaturated
    frames to the model delta alpha \propto tan(alpha), where alpha
    is the zenith angle and is retrieved from the fits headers.
    fit_centroids takes two arguments:

    1.  A list of HiCIAO frames (as filenames)
    2.  An array with all of the destriped HiCIAO data
    """

    print '\nFitting Centroids in Parallel'
    nframes = len(filesetup.framelist)

    tan_zd = np.ndarray(nframes, np.float)
    t = np.ndarray(nframes, np.float)
    newcenters = np.ndarray((2, nframes), np.float)

    ##################################################################
    # Take a first guess at the centroid from the center of light,
    # then fit a Moffat profile.
    ##################################################################

    if method == 'moffat':
        flux_avg = np.mean(flux, axis=0)
        dimy, dimx = flux_avg.shape
        x, y = np.meshgrid(np.arange(dimx), np.arange(dimy))
        x = np.sum(x * flux_avg) / np.sum(flux_avg)
        y = np.sum(y * flux_avg) / np.sum(flux_avg)
        #print y, x
        
        centerguess = np.ndarray((nframes, 2), np.float)
        centerguess[:, 0] = y
        centerguess[:, 1] = x
        y, x = parallel._moffat_centroid(filesetup.framelist, flux,
                                         centerguess=centerguess)
        mu_y, mu_x = [np.mean(y), np.mean(x)]
        sigy, sigx = [1000, 1000]
        for i in range(4):
            sigy = np.std(y[np.where((y - mu_y)**2 < 10 * sigy**2)])
            sigx = np.std(x[np.where((x - mu_x)**2 < 10 * sigx**2)])
            mu_y = np.mean(y[np.where((y - mu_y)**2 < 10 * sigy**2)])
            mu_x = np.mean(x[np.where((x - mu_x)**2 < 10 * sigx**2)])
        #np.copyto(y, mu_y, where=(y - mu_y)**2 > 10 * sigy**2)
        #np.copyto(x, mu_x, where=(x - mu_x)**2 > 10 * sigx**2)
        np.putmask(y, (y - mu_y)**2 > 10 * sigy**2, mu_y)
        np.putmask(x, (x - mu_x)**2 > 10 * sigx**2, mu_x)

        newcenters[:, :] = [y, x]
            
        #for i in range(nframes):
        #    print newcenters[:, i]
        return newcenters.T
    
    ##################################################################
    # Centroid the frames by cross-correlating against a library
    # of PSF templates
    ##################################################################

    elif method == 'crosscorr':
        y, x = parallel._crosscorr_centroid(filesetup.framelist, flux)
        newcenters[:, :] = [y, x]
        return newcenters.T

    else:
        print '\nCentroiding method ' + method + ' not recognized in fit_centroids.'
        sys.exit()

