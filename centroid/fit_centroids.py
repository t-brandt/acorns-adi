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
import re, pickle
import pyfits as pyf
import parallel
import centroid

def fit_centroids(filesetup, flux, pa, method='moffat', storeall=True,
                  objname='Unknown_Object', psf_dir='psfref', side=None,
                  ref_psf=None):

    """
    Function fit_centroids is a wrapper for routines with the full
    centroid calculations.
    fit_centroids takes two arguments:

    1.  A list of HiCIAO frames (as filenames)
    2.  An array with all of the destriped HiCIAO data
    """

    print '\nFitting Centroids in Parallel'
    nframes = len(filesetup.framelist)

    newcenters = np.ndarray((2, nframes), np.float)
    dr_rms = None

    ##################################################################
    # Take a first guess at the centroid from the center of light,
    # then fit a Moffat profile.
    ##################################################################

    if method == 'moffat':
        if flux is not None:
            flux_avg = np.mean(flux, axis=0)
        else:
            filename = re.sub(".fits", "_ds.fits",
                              filesetup.framelist[nframes // 2])
            filename = re.sub(".*/", filesetup.reduce_dir + "/", filename)
            flux_avg = pyf.open(filename)[-1].data

        if side is not None:
            if re.search('[Ll]eft', side):
                flux_avg[:, dimx // 2:] = 0
            elif re.search('[Rr]ight', side):
                flux_avg[:, :dimx // 2] = 0
            
        dimy, dimx = flux_avg.shape
        x, y = np.meshgrid(np.arange(dimx), np.arange(dimy))
        x = np.sum(x * flux_avg) / np.sum(flux_avg)
        y = np.sum(y * flux_avg) / np.sum(flux_avg)

        centerguess = np.ndarray((nframes, 2), np.float)
        centerguess[:, 0] = y
        centerguess[:, 1] = x
        y, x = parallel._moffat_centroid(filesetup.framelist, flux,
                                         centerguess=centerguess)
        mu_y, mu_x = [np.mean(y), np.mean(x)]
        sigy, sigx = [1000, 1000]

        ##################################################################
        # Use sigma-reject to mask discrepant frames.
        ##################################################################
        
        for i in range(4):
            sigy = np.std(y[np.where((y - mu_y)**2 < 10 * sigy**2)])
            sigx = np.std(x[np.where((x - mu_x)**2 < 10 * sigx**2)])
            mu_y = np.mean(y[np.where((y - mu_y)**2 < 10 * sigy**2)])
            mu_x = np.mean(x[np.where((x - mu_x)**2 < 10 * sigx**2)])
            
        np.putmask(y, (y - mu_y)**2 > 10 * sigy**2, mu_y)
        np.putmask(x, (x - mu_x)**2 > 10 * sigx**2, mu_x)

        newcenters[:, :] = [y, x]
        centers = newcenters.T
    
    ####################################################################
    # Main centroiding algorithm, discussed in Brandt+ 2012
    # Centroid a map of chi2 residuals after subtracting PCA components
    # The success array is false where the algorithm fails.
    ####################################################################

    elif method == 'crosscorr':
        success, centers, dr_rms = parallel._cc_centroid(filesetup.framelist,
                                                         flux, psf_dir=psf_dir,
                                                         side=side,
                                                         ref_psf=ref_psf)
        igood = np.extract(success, np.arange(nframes))
        for i in np.arange(nframes - 1, -1, -1):
            if not success[i]:
                del filesetup.framelist[i]
        pickle.dump(filesetup, open('./dirinfo', 'w'))
        centers = centers[igood]
        #print centers
        pa = pa[igood]
        nframes = len(igood)
        if flux is not None:
            flux = flux[igood]

    else:
        print '\nCentroiding method ' + method + ' not recognized in fit_centroids.'
        sys.exit()


    ####################################################################
    # Interactively set the absolute centroid.  
    ####################################################################

    fluxcen = parallel._rotate_recenter(filesetup, flux, storeall=True,
                                        centers=centers, newdimen=201,
                                        write_files=False)
    fluxcen = np.median(fluxcen, axis=0)
    yc, xc = centroid.finecenter(fluxcen, objname, filesetup.output_dir)
    centers[:, 0] += yc
    centers[:, 1] += xc

    return centers, dr_rms

