#!/usr/bin/env python
#
# Original filename: flat_coadd.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: 7 June 2011
# 
# Summary:  Reduce HiCIAO ADI Data
# 
"""
%prog [options]
"""

import optparse, sys, re, os
import pyfits as pyf
import numpy as np
from scipy import signal
from adiparam import *
import centroid
import transform
import parallel
import combine
import loci
import pca
import utils
import pickle
import addsource
import locitools

def main():

    """
    Main program for ADI data reduction, configured with a call to
    adiparam.GetConfig(), which brings up a GUI to set parameters.

    The pipeline is currently designed for SEEDS data taken without
    an occulting mask.  
    
    You must have scipy, numpy, pyephem, multiprocessing, and matplotlib
    installed to use this pipeline.
    """

    exec_path = os.path.dirname(os.path.realpath(__file__))
    filesetup, adipar, locipar = GetConfig()

    nframes = len(filesetup.framelist)
    flat = pyf.open(filesetup.flat)
    if filesetup.pixmask is not None:
        hotpix = pyf.open(filesetup.pixmask)
    else:
        hotpix = None

    dimy, dimx = pyf.open(filesetup.framelist[0])[0].data.shape
    mem, ncpus, storeall = utils.config(nframes, dimy * dimx)
    pa = np.asarray([transform.get_pa(frame) * -1 * np.pi / 180
                     for frame in filesetup.framelist])
    
    fullframe = re.sub("-C.*fits", ".fits", filesetup.framelist[0])
    objname = pyf.open(fullframe)[0].header['OBJECT']
    np.savetxt(filesetup.output_dir + '/' + objname + '_palist.dat', pa)
    dr_rms = None

    ####################################################################
    # Default save/resume points: destriping, recentering, final files
    ####################################################################
    
    if np.all(utils.check_files(filesetup, ext="_r")):
        print "\nResuming reduction from recentered files."
        flux = utils.read_files(filesetup, ext="_r") 
    else:
        if np.all(utils.check_files(filesetup, ext="_ds")):
            flux = utils.read_files(filesetup, ext="_ds")
        else:
            flux = parallel._destripe(filesetup, flat, hotpix, mem, adipar,
                                      write_files=True, storeall=storeall)
        flux = parallel._dewarp(filesetup, mem, flux=flux, storeall=storeall)

        ####################################################################
        # Main centroiding algorithm, discussed in Brandt+ 2012
        # Centroid a map of chi2 residuals after subtracting PCA components
        # The success array is false where the algorithm fails.
        ####################################################################
        
        if adipar.center == 'crosscorr':
            success, centers, dr_rms = parallel._cc_centroid(filesetup.framelist, flux, psf_dir=exec_path + '/psfref')
        
            igood = np.extract(success, np.arange(nframes))
            for i in np.arange(nframes - 1, -1, -1):
                if not success[i]:
                    del filesetup.framelist[i]
            pickle.dump(filesetup, open('./dirinfo', 'w'))
            centers = centers[igood]
            pa = pa[igood]
            nframes = len(igood)
            if flux is not None:
                flux = flux[igood]

        ####################################################################
        # Backup centroiding algorithm: fit a Moffat profile
        ####################################################################
        
        else:
            centers = centroid.fit_centroids(filesetup, flux, method='moffat')
            
        np.savetxt(filesetup.output_dir + '/' + objname + '_centers.dat',
                       centers)

        ####################################################################
        # Interactively set the absolute centroid.  
        ####################################################################

        mask_id = pyf.open(fullframe)[0].header['P_FMID']
        fluxcen = parallel._rotate_recenter(filesetup, flux, storeall=True,
                                            centers=centers, newdimen=201,
                                            write_files=False)
        fluxcen = np.sum(fluxcen, axis=0)
        yc, xc = centroid.finecenter(fluxcen, objname, filesetup.output_dir)
        centers[:, 0] += yc
        centers[:, 1] += xc

        ####################################################################
        # Recenter the data onto a square array of the largest dimension
        # such that the entire array has data
        ####################################################################

        mindim = min(dimy - centers[:, 0].max(), centers[:, 0].min(),
                     dimx - centers[:, 1].max(), centers[:, 1].min())
        mindim = int(mindim) * 2 - 1
        flux = parallel._rotate_recenter(filesetup, flux, storeall=storeall,
                                         centers=centers, newdimen=mindim,
                                         write_files=True)

    ####################################################################
    # Perform scaled PCA on the flux array; alternatively, read in an
    # array of principal components.  Neither is currently used.
    ####################################################################
    
    if False:
        pcapath = '/scr/wakusei1/users/tbrandt'
        flux, pca_arr = pca.pca(flux, ncomp=0, nread=5, dosub=True,
                                pcadir=pcapath + '/psfref')
        if dr_rms is None:
            dr_rms = 20
    elif False:
        pca_dir = '.'
        npca = 40
        pca_arr = np.zeros((npca, flux.shape[1], flux.shape[2]), np.float32)
        for i in range(npca):
            tmp = pyf.open(pca_dir + '/pcacomp_' + str(i) + '.fits')[0].data
            dy, dx = [tmp.shape[0] // 2, tmp.shape[1] // 2]
            pca_arr[i, yc - dy:yc + dy + 1, xc - dx:xc + dx + 1] = tmp
    else:
        pca_arr = None

    ####################################################################
    # Find the n closest matches to each frame.  Not currently used.
    ####################################################################

    if False:
        corr = pca.allcorr(range(int(locipar.rmax)), flux, n=80)    
    else:
        corr = None
        
    ####################################################################
    # Subtract a radial profile from each frame.  Not currently used.
    ####################################################################

    if False:
        flux = parallel._radialsub(filesetup, flux, mode='median', 
                                   center=None, rmax=None, smoothwidth=0)

    ####################################################################
    # Run LOCI if that ADI reduction method is chosen
    ####################################################################

    partial_sub = None
    x = np.arange(flux.shape[1]) - flux.shape[1] // 2
    x, y = np.meshgrid(x, x)
    r = np.sqrt(x**2 + y**2)
    
    if adipar.adi == 'LOCI':

        ####################################################################
        # Set the maximum radius at which to perform LOCI
        ####################################################################
        
        deltar = np.sqrt(np.pi * locipar.fwhm**2 / 4 * locipar.npsf)
        rmax = int(flux.shape[1] // 2 - deltar - 50)
        locipar.rmax = min(locipar.rmax, rmax)

        if dr_rms is None:
            nf, dy, dx = flux.shape
            fluxmed = np.median(flux, axis=0)[dy // 2 - 100:dy // 2 + 101,
                                              dx // 2 - 100:dx // 2 + 101]
            sat = fluxmed > 0.7 * fluxmed.max()
            r2 = r[dy // 2 - 100:dy // 2 + 101, dx // 2 - 100:dx // 2 + 101]**2
            dr_rms = np.sqrt(np.sum(r2 * sat) / np.sum(sat))

        ####################################################################
        # This is regular LOCI
        ####################################################################
        
        if locipar.feedback == 0:
            partial_sub = loci.loci(flux, pa, locipar, mem, mode='LOCI',
                                    pca_arr=pca_arr, r_ex=dr_rms, corr=corr,
                                    method='matrix', do_partial_sub=True,
                                    sub_dir=exec_path)
            
        ####################################################################
        # The next block runs LOCI once, de-rotates, takes the median,
        # and re-rotates to each frame's position angle.  It then runs
        # LOCI again to over-correct the result.  Not recommended for
        # SEEDS data with AO188.
        ####################################################################
        
        else:
            fluxref = np.ndarray(flux.shape, np.float32)
            fluxref[:] = flux
            
            loci.loci(fluxref, pca_arr, pa, locipar, mem, mode='LOCI',
                      r_ex=dr_rms, pca_arr=pca_arr,
                      corr=corr, method='matrix', do_partial_sub=False)
            
            for i in range(nframes):
                np.putmask(fluxref[i], r > locipar.rmax - 1, 0)
                np.putmask(fluxref[i], r < dr_rms + 1, 0)
            locipar.rmax -= 100
            fluxref = parallel._rotate_recenter(filesetup, fluxref, theta=pa)
            
            for i in range(nframes):
                np.putmask(fluxref[i], r > locipar.rmax - 1, 0)
                np.putmask(fluxref[i], r < dr_rms + 1, 0)
            locipar.rmax -= 100
            fluxmed = np.median(fluxref, axis=0)
            for i in range(nframes):
                fluxref[i] = fluxmed * locipar.feedback
            fluxref = parallel._rotate_recenter(filesetup, fluxref, theta=-pa)
            
            loci.loci(flux, pa, locipar, mem, mode='refine', fluxref=fluxref,
                      pca_arr=pca_arr, rmin=dr_rms, r_ex=dr_rms)

        ####################################################################
        # Mask saturated areas (< dr_rms), do median subtraction at radii
        # beyond the limit of the LOCI reduction
        ####################################################################

        fluxmed = np.median(flux, axis=0)
        for i in range(nframes):
            np.putmask(flux[i], r < dr_rms + 2, 0)
            np.putmask(flux[i], r > locipar.rmax - 1, flux[i] - fluxmed)
            
    ####################################################################
    # Alternative to LOCI: median PSF subtraction
    ####################################################################

    elif adipar.adi == 'median':
        medpsf = np.median(flux, axis=0)
        for i in range(nframes):
            flux[i] -= medpsf

    else:
        print "Error:  ADI reduction method " + adipar.adi + " not recognized."
        #sys.exit(1)

    ####################################################################
    # Derotate, combine flux array using mean/median hybrid (see
    # Brandt+ 2012), measure standard deviation at each radius
    ####################################################################

    newhead = utils.makeheader(flux[0], pyf.open(fullframe)[0].header,
                               filesetup, adipar, locipar)
    
    flux = parallel._rotate_recenter(filesetup, flux, theta=pa)        
    fluxbest, noise = combine.meanmed(flux)
    
    x, y = np.meshgrid(np.arange(7) - 3, np.arange(7) - 3)
    window = (x**2 + y**2 < 2.51**2) * 1.0
    window /= np.sum(window)
    fluxbest = signal.convolve2d(fluxbest, window, mode='same')
    noise = combine.radprof(fluxbest, mode='std', smoothwidth=2, sigrej=4.5)[0]
    
    np.putmask(fluxbest, r < dr_rms + 3, np.nan)
    np.putmask(fluxbest, r > locipar.rmax - 2, np.nan)
    
    np.savetxt(filesetup.output_dir + '/' + objname + '_noiseprofile.dat',
               noise[noise.shape[0] // 2, noise.shape[1] // 2:].T)
    fluxsnr = (fluxbest / noise).astype(np.float32)

    ####################################################################
    # 5-sigma sensitivity maps--just multiply by the scaled aperture
    # photometry of the central star
    ####################################################################

    if partial_sub is not None:
        sensitivity = np.ones(fluxsnr.shape, np.float32)
        dx = (sensitivity.shape[0] - partial_sub.shape[0]) // 2
        if dx > 0:
            sensitivity[dx:-dx, dx:-dx] = 1 / partial_sub
        else:
            sensitivity[:] = 1 / partial_sub
        sensitivity *= noise * 5
        
    ####################################################################
    # Write the output fits files.  Rescale all arrays to 2001x2001
    # so that the center is pixel number (1000, 1000) indexed from 0.
    # Use NaN to pad arrays.
    ####################################################################

    dimy, dimx = fluxbest.shape
    
    if dimy > 2001:
        di = (dimy - 2001) // 2
        snr = pyf.HDUList(pyf.PrimaryHDU(fluxsnr[di:-di, di:-di], newhead))
        final = pyf.HDUList(pyf.PrimaryHDU(fluxbest[di:-di, di:-di].astype(np.float32), newhead))
        if partial_sub is not None:
            contrast = pyf.HDUList(pyf.PrimaryHDU(sensitivity[di:-di, di:-di],
                                                  newhead))
    elif dimy < 2001:
        di = (2001 - dimy) // 2
        snr_pad = np.zeros((2001, 2001), np.float32) + np.nan
        fluxbest_pad = np.zeros((2001, 2001), np.float32) + np.nan
        contrast_pad = np.zeros((2001, 2001), np.float32) + np.nan

        snr_pad[di:-di, di:-di] = fluxsnr
        fluxbest_pad[di:-di, di:-di] = fluxbest
        
        snr = pyf.HDUList(pyf.PrimaryHDU(snr_pad, newhead))
        final = pyf.HDUList(pyf.PrimaryHDU(fluxbest_pad, newhead))

        if partial_sub is not None:
            contrast_pad[di:-di, di:-di] = sensitivity
            contrast = pyf.HDUList(pyf.PrimaryHDU(contrast_pad, newhead))
    else:
        snr = pyf.HDUList(pyf.PrimaryHDU(fluxsnr, newhead))
        final = pyf.HDUList(pyf.PrimaryHDU(fluxbest.astype(np.float32), newhead))
        if partial_sub is not None:
            contrast = pyf.HDUList(pyf.PrimaryHDU(sensitivity, newhead))

    name_base = filesetup.output_dir + '/' + objname
    snr.writeto(name_base + '_snr.fits', clobber=True)
    final.writeto(name_base + '_final.fits', clobber=True)
    if partial_sub is not None:
        contrast.writeto(name_base + '_5sigma_sensitivity.fits', clobber=True)

#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()
