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
import photometry

def main():

    """
    Main program for ADI data reduction, configured with a call to
    adiparam.GetConfig(), which brings up a GUI to set parameters.

    The pipeline is currently designed for SEEDS data taken without
    an occulting mask.  
    
    You must have scipy, numpy, pyephem, multiprocessing, and matplotlib
    installed to use this pipeline.
    """

    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-p", "--prefix", dest="prefix", default="HICA",
                      help="Specify raw file name prefix (default=%default)")
    opts, args = parser.parse_args()

    exec_path = os.path.dirname(os.path.realpath(__file__))
    filesetup, adipar, locipar = GetConfig(prefix=opts.prefix)

    nframes = len(filesetup.framelist)
    ngroup = 1 + int((nframes - 1) / locipar.max_n)
    flat = pyf.open(filesetup.flat)
    if filesetup.pixmask is not None:
        hotpix = pyf.open(filesetup.pixmask)
    else:
        hotpix = None

    dimy, dimx = pyf.open(filesetup.framelist[0])[-1].data.shape
    mem, ncpus, storeall = utils.config(nframes, dimy * dimx)
    
    if filesetup.scale_phot:
        x, y = np.meshgrid(np.arange(7) - 3, np.arange(7) - 3)
        window = (x**2 + y**2 < 2.51**2) * 1.0
        window /= np.sum(window)
        ref_phot, ref_psf = photometry.calc_phot(filesetup, adipar, flat,
                                                 hotpix, mem, window)
    else:
        ref_psf = None
        ref_phot = None
    
    ################################################################
    # WCS coordinates are not reliable in HiCIAO data with the image
    # rotator off.  Compute parallactic angle.  Otherwise, trust the
    # WCS coordinates.
    ################################################################

    if 'HICA' in filesetup.framelist[0]:
        pa = np.asarray([transform.get_pa(frame) * -1 * np.pi / 180
                         for frame in filesetup.framelist])
    else:
        pa = np.ones(len(filesetup.framelist))
        for i in range(len(filesetup.framelist)):
            cd2_1 = pyf.open(filesetup.framelist[i])[0].header['cd2_1']
            cd2_2 = pyf.open(filesetup.framelist[i])[0].header['cd2_2']
            pa[i] = -np.arctan2(cd2_1, cd2_2)
            
    fullframe = re.sub("-C.*fits", ".fits", filesetup.framelist[0])
    try:
        objname = pyf.open(fullframe)[0].header['OBJECT']
    except:
        objname = "Unknown_Object"
    objname = re.sub(' ', '_', objname)
    np.savetxt(filesetup.output_dir + '/' + objname + '_palist.dat', pa)
    dr_rms = None

    ####################################################################
    # Default save/resume points: destriping, recentering, final files
    # Configuration gives the option to skip the destriping step (only
    # performing a flat-field), the dewarping, and the centering.
    ####################################################################
    
    if np.all(utils.check_files(filesetup, ext="_r")):
        print "\nResuming reduction from recentered files."
        if ngroup == 1:
            flux = utils.read_files(filesetup, ext="_r")
        else:
            flux = utils.read_files(filesetup, ext="_r")
    else:
        if storeall and np.all(utils.check_files(filesetup, ext="_ds")):
            flux = utils.read_files(filesetup, ext="_ds")
        elif not np.all(utils.check_files(filesetup, ext="_ds")):
            flux = parallel._destripe(filesetup, flat, hotpix, mem, adipar,
                                      write_files=True, storeall=storeall,
                                      full_destripe=adipar.full_destripe,
                                      do_horiz=adipar.full_destripe)
        else:
            flux = None
            
        if adipar.dewarp:
            flux = parallel._dewarp(filesetup, mem, flux=flux, storeall=storeall)

        if adipar.do_centroid:
            centers, dr_rms = centroid.fit_centroids(filesetup, flux, pa,
                                                     storeall=storeall,
                                                     objname=objname,
                                                     method=adipar.center,
                                                     psf_dir=exec_path+'/psfref', ref_psf=ref_psf)
            #centers = np.ndarray((nframes, 2))
            #centers[:, 0] = 1026 - 128
            #centers[:, 1] = 949 + 60

            #dr_rms = 30
            np.savetxt(filesetup.output_dir + '/' + objname +
                       '_centers.dat', centers)

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
            nframes = len(filesetup.framelist)

    ####################################################################
    # Perform scaled PCA on the flux array; alternatively, read in an
    # array of principal components.  Neither is currently used.
    ####################################################################
    
    if False:
        pcapath = '/scr/wakusei1/users/tbrandt'
        flux, pca_arr = pca.pca(flux, ncomp=20, nread=2, dosub=True,
                                pcadir=pcapath + '/psfref')
        for i in range(nframes):
            out = pyf.HDUList(pyf.PrimaryHDU(flux[i].astype(np.float32),
                                             pyf.open(filesetup.framelist[i])[0].header))
            rootfile = re.sub('.*/', '', filesetup.framelist[i])
            out.writeto(filesetup.reduce_dir + '/' + re.sub('.fits', '_r.fits', rootfile), clobber=True)
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
        ngroup = 1
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
    full_pa = pa.copy()
    full_framelist = [frame for frame in filesetup.framelist]
    for igroup in range(ngroup):

        if ngroup > 1:
            filesetup.framelist = full_framelist[igroup::ngroup]
            if np.all(utils.check_files(filesetup, ext="_r")):
                flux = utils.read_files(filesetup, ext="_r")
            else:
                print "Unable to read recentered files for LOCI."
                sys.exit()
            pa = full_pa[igroup::ngroup]
        
        x = np.arange(flux.shape[1]) - flux.shape[1] // 2
        x, y = np.meshgrid(x, x)
        r = np.sqrt(x**2 + y**2)
        
        if adipar.adi == 'LOCI':

            ################################################################
            # Set the maximum radius at which to perform LOCI
            ################################################################
        
            deltar = np.sqrt(np.pi * locipar.fwhm**2 / 4 * locipar.npsf)
            rmax = int(flux.shape[1] // 2 - deltar - 50)
            locipar.rmax = min(locipar.rmax, rmax)
                        
            if dr_rms is None:
                nf, dy, dx = flux.shape
                fluxmed = np.median(flux, axis=0)[dy // 2 - 100:dy // 2 + 101,
                                                  dx // 2 - 100:dx // 2 + 101]
                sat = fluxmed > 0.7 * fluxmed.max()
                r2 = r[dy//2 - 100:dy//2 + 101, dx//2 - 100:dx//2 + 101]**2
                dr_rms = np.sqrt(np.sum(r2 * sat) / np.sum(sat))

            ################################################################
            # This is regular LOCI
            ################################################################
        
            if locipar.feedback == 0:
                partial_sub = loci.loci(flux, pa, locipar, mem, mode='LOCI',
                                        pca_arr=None, r_ex=dr_rms, corr=corr,
                                        method='matrix', do_partial_sub=True,
                                        sub_dir=exec_path)
                
            ################################################################
            # The next block runs LOCI once, de-rotates, takes the median,
            # and re-rotates to each frame's position angle.  It then runs
            # LOCI again to over-correct the result.  Not recommended for
            # SEEDS data with AO188.
            ################################################################
        
            else:
                fluxref = np.ndarray(flux.shape, np.float32)
                fluxref[:] = flux
            
                loci.loci(fluxref, pca_arr, pa, locipar, mem, mode='LOCI',
                          r_ex=dr_rms, pca_arr=pca_arr,
                          corr=corr, method='matrix', do_partial_sub=False)
                
                for i in range(flux.shape[0]):
                    np.putmask(fluxref[i], r > locipar.rmax - 1, 0)
                    np.putmask(fluxref[i], r < dr_rms + 1, 0)
                locipar.rmax -= 100
                fluxref = parallel._rotate_recenter(filesetup, fluxref, theta=pa)
            
                for i in range(flux.shape[0]):
                    np.putmask(fluxref[i], r > locipar.rmax - 1, 0)
                    np.putmask(fluxref[i], r < dr_rms + 1, 0)
                locipar.rmax -= 100
                fluxmed = np.median(fluxref, axis=0)
                for i in range(flux.shape[0]):
                    fluxref[i] = fluxmed * locipar.feedback
                fluxref = parallel._rotate_recenter(filesetup, fluxref, theta=-pa)
            
                loci.loci(flux, pa, locipar, mem, mode='refine', fluxref=fluxref,
                          pca_arr=pca_arr, rmin=dr_rms, r_ex=dr_rms)

            ################################################################
            # Mask saturated areas (< dr_rms), do median subtraction at radii
            # beyond the limit of the LOCI reduction
            ################################################################

            fluxmed = np.median(flux, axis=0)
            for i in range(flux.shape[0]):
                np.putmask(flux[i], r < dr_rms + 2, 0)
                np.putmask(flux[i], r > locipar.rmax - 1, flux[i] - fluxmed)
                
         ####################################################################
         # Alternative to LOCI: median PSF subtraction
         ####################################################################

        elif adipar.adi == 'median':
            medpsf = np.median(flux, axis=0)
            for i in range(flux.shape[0]):
                flux[i] -= medpsf

        else:
            print "Error:  ADI reduction method " + adipar.adi + " not recognized."
            #sys.exit(1)

        ####################################################################
        # Derotate, combine flux array using mean/median hybrid (see
        # Brandt+ 2012), measure standard deviation at each radius
        ####################################################################

        if igroup == 0:
            newhead = utils.makeheader(flux[0], pyf.open(fullframe)[0].header,
                                       full_framelist, adipar, locipar)
            
            flux = parallel._rotate_recenter(filesetup, flux, theta=pa)
            fluxtmp, noise = combine.meanmed(flux)
            fluxbest = fluxtmp / ngroup
            if partial_sub is not None:
                partial_sub_tot = partial_sub / ngroup
        else:
            flux = parallel._rotate_recenter(filesetup, flux, theta=pa)
            fluxtmp, noise = combine.meanmed(flux)
            fluxbest += fluxtmp / ngroup
            if partial_sub is not None:
                partial_sub_tot += partial_sub / ngroup
            
    filesetup.framelist = full_framelist
    if partial_sub is not None:
        partial_sub = partial_sub_tot            
    
    ####################################################################
    # Rescale all arrays to 2001x2001 so that the center is pixel number
    # (1000, 1000) indexed from 0.  Use NaN to pad arrays.
    ####################################################################

    fluxbest = utils.arr_resize(fluxbest)
    if partial_sub is not None:
        partial_sub = utils.arr_resize(partial_sub, newdim=fluxbest.shape[0]).astype(np.float32)
        fluxbest /= partial_sub
        out = pyf.HDUList(pyf.PrimaryHDU(partial_sub))
        out.writeto('partial_sub2.fits', clobber=True)
        
    x, y = np.meshgrid(np.arange(7) - 3, np.arange(7) - 3)
    window = (x**2 + y**2 < 2.51**2) * 1.0
    window /= np.sum(window)
    fluxbest = signal.convolve2d(fluxbest, window, mode='same')
    noise = combine.radprof(fluxbest, mode='std', smoothwidth=2, sigrej=4.5)[0]

    r = utils.arr_resize(r)
    if dr_rms is not None:
        np.putmask(fluxbest, r < dr_rms + 3, np.nan)
    np.putmask(fluxbest, r > locipar.rmax - 2, np.nan)
    
    fluxsnr = (fluxbest / noise).astype(np.float32)

    ####################################################################
    # 5-sigma sensitivity maps--just multiply by the scaled aperture
    # photometry of the central star
    ####################################################################
    
    if partial_sub is not None:
        sensitivity = noise * 5 / partial_sub

        ####################################################################
        # Photometry of the central star
        ####################################################################

        if filesetup.scale_phot:
            #ref_phot = photometry.calc_phot(filesetup, adipar, flat,
            #                                    hotpix, mem, window)[0]
            sensitivity /= ref_phot
            fluxbest /= ref_phot
            noise /= ref_phot
        
        sig_sens = combine.radprof(sensitivity, mode='std', smoothwidth=0)[0]
        outfile = open(filesetup.output_dir + '/' + objname +
                       '_5sigma_sensitivity.dat', 'w')
        for i in range(sig_sens.shape[0] // 2, sig_sens.shape[0]):
            iy = sig_sens.shape[0] // 2
            if np.isfinite(sensitivity[iy, i]):
                outfile.write('%8d  %12.5e  %12.5e  %12e\n' %
                              (i - iy, sensitivity[iy, i], sig_sens[iy, i],
                               partial_sub[iy, i]))
        outfile.close()
        
    else:
        np.savetxt(filesetup.output_dir + '/' + objname + '_noiseprofile.dat',
                   noise[noise.shape[0] // 2, noise.shape[1] // 2:].T)
            
    ####################################################################
    # Write the output fits files. 
    ####################################################################

    snr = pyf.HDUList(pyf.PrimaryHDU(fluxsnr.astype(np.float32), newhead))
    final = pyf.HDUList(pyf.PrimaryHDU(fluxbest.astype(np.float32), newhead))
    if partial_sub is not None:
        contrast = pyf.HDUList(pyf.PrimaryHDU(sensitivity.astype(np.float32), newhead))

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
