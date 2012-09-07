#!/usr/bin/env python
#
# Original filename: loci.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: July 2011
# 
# Summary:  Run LOCI on the input frame sequence.
# 

import sys, re, gc
import numpy as np
from scipy import interpolate
import multiprocessing
from progressbar import ProgressBar
from locitools import *
from parallel import *

def loci(flux, pa, locipar, mem, mode='LOCI', fluxref=None, 
         pca_arr=None, rmin=None, r_ex=None, niter=0,
         corr=None, method='matrix', do_partial_sub=True,
         sub_dir='.'):

    """

    Function loci performs Locally Optimized Combination of Images
    on an input dataset (argument 1).  It takes at least 4 arguments:
    
    1.  A three-dimensional flux array.  The first dimension should run
        over frames.
    2.  The position angle of the image rotator (in radians) in each frame
    3.  The LOCI configuration parameters.  This is an object defined in
        adiparam.lociparam
    4.  The available memory (in bytes)

    Optional arguments:
    
    5.  'mode': The mode in which to run the routine.  Currently
        accepted values are 'LOCI' (normal LOCI) and 'refine', in
        which case the LOCI coefficients minimize the difference
        between a frame and and a linear combination of the other
        frames with fluxref (argument 6) added.  Default 'LOCI'
    6.  'fluxref': A reference array, to be computed using PCA, to fit
        the PSF locally in each frame.  This is treated as an extra
        set of reference frames with infinite angular displacement.
        Ignored if None.  Default None.
    7.  'pca_arr': The flux array used to compute the LOCI subtraction
        coefficients.  This can be different from the array (argument
        1) on which those subtraction coefficients are used.  Ignored
        if None, mandatory if mode=='refine'.  Default None.
    8.  'rmin': The minimum radius at which to perform LOCI.  If None,
        this is the smaller of the saturation radius and minimum
        radius at which there are always frames to satisfy LOCI's
        angular displacement criterion.  Default None.
    9.  'r_ex': Exclusion/saturation radius.  Ignored if None; should
        be calculated using the centroiding algorithm.  Default None.
    10. 'niter': Number of steps in which to split execution to limit
        memory usage.  If 0, this is calculated within the routine.
        Default 0.
    11. 'corr': Cross-correlation coefficients between frames, to use
        only good matches in LOCI.  Possibly useful for very large
        datasets, but it can increase run time significantly.  Ignored
        if None.  Default None.
    12. 'method': Method for solving the LOCI matrix equation.
        Default is 'matrix', in which case the problem matrix is
        formed and solved using LU- decomposition.  Other choices:
        'lstsq' (SVD), 'eqcoef' (all coefficients are equal and sum to
        1--similar to mean subtraction).  Default: 'matrix'.  'lstsq'
        is more robust but will be significantly slower in most cases.
    13. 'do_partial_sub': Compute the fractional flux suppression
        using the method described in Brandt+2012.  This requires
        argument 'method' to be 'matrix', as it uses the LU
        decomposition of the problem matrix, and it requires argument
        'mode' to be 'LOCI'.  Default True

    Function loci returns the flux suppression map if argument
    'do_partial_sub' is True; otherwise, it returns None.
    
    """

    np.seterr(all='ignore')
    
    ######################################################################
    # Check arrays and sizes.
    ######################################################################
    
    if mode == 'refine' and fluxref is None:
        print "Cannot apply feedback to LOCI without reference output."
        sys.exit(1)
    if fluxref is not None:
        assert fluxref.shape == flux.shape, \
               "Flux and reference arrays must have the same shape."
    assert len(flux.shape) == 3, "Flux array must have three dimensions."
    assert niter >= 0 and niter < flux.shape[0], \
           "Argument niter to loci must be at least 0 and less than nframes"
    assert mode == 'refine' or mode == 'LOCI', \
           "Argument mode = " + mode + " to routine loci is not understood." \
           + "\nPlease use either 'refine' or 'LOCI'."
    
    ######################################################################
    # Create and populate array fluxsmooth to calculate LOCI coefficients
    ######################################################################
    
    nframes = flux.shape[0]
    ngroup = 1 + int((nframes - 1) / locipar.max_n)
    oldshape = flux.shape
    fluxsmooth = np.ndarray(flux.shape, np.float32)

    if mode == 'LOCI' and fluxref is not None:
        fluxsmooth[:] = fluxref[i]
    elif mode == 'refine' and fluxref is not None:
        fluxsmooth[:] = flux + fluxref
    elif mode == 'refine':
        print "Cannot call loci in 'refine' mode without a reference array."
        sys.exit()
    else:
        fluxsmooth[:] = flux

    fluxsmooth = np.reshape(fluxsmooth, (nframes, -1))
    flux = np.reshape(flux, (nframes, -1))
    if pca_arr is not None:
        pca_arr = np.reshape(pca_arr, (pca_arr.shape[0], -1))

    ######################################################################
    # Create arrays sorted by radius to facilitate division of frames
    # into optimization and subtraction regions.
    ######################################################################

    print "\nSorting radius array to compute indices at annuli."
    rindex, rsort, thetasort = sorted_arrays(oldshape)

    optreg = int(np.pi * locipar.fwhm**2 / 4 * locipar.npsf)
    minsep = locipar.nfwhm * locipar.fwhm
    if rmin is None:
        rmin = int(minsep * 2 / (pa.max() - pa.min())) + 1
        rmin = max(r_ex - r_ex % 5, rmin - rmin % 5) + 5
    if r_ex is None:
        r_ex = rmin

    ######################################################################
    # Calculate the total number of subtraction regions and the 
    # (r, theta) boundaries of those regions.  
    ######################################################################
    
    totn, r1, r2, theta1, theta2 = calczones(rmin, locipar.rmax, locipar.dr0, 
                                             optreg)

    ######################################################################
    # Coefficients for flux suppression by a displaced source
    # They are pre-computed for an azimuthally averaged HiCIAO PSF.
    ######################################################################

    if do_partial_sub:
        try:
            sub_coefs = np.zeros((totn, nframes), np.float32)
            x, y = np.loadtxt(sub_dir + '/partial_sub.dat').T
            partial_sub = interpolate.UnivariateSpline(x, y, k=3, s=0)
        except:
            print "Unable to read flux loss data from file partial_sub.dat."
            print "Cannot calculate a fractional flux loss map in loci."
            do_partial_sub = False
            
    ######################################################################
    # Calculate where to add mock sources for flux loss analysis.
    # We will add a positive source (one at a time) at the center of
    # each subtraction region, and negative sources azimuthally
    # displaced from it.  See Appendix of Brandt+2012 for explanation.
    ######################################################################
    
    r_avg = 0.5 * (r1 + r2)
    phi = 0.5 * (theta1 + theta2)
    dphi = minsep * 1.6 / r_avg
    if do_partial_sub:
        np.putmask(dphi, dphi < np.std(pa), np.std(pa))

        ######################################################################
        # Size of PSF regions will be 2 * n + 1.  We will compute them
        # later; just store r**2 for now.
        ######################################################################
        
        n = 7
        x, y = np.meshgrid(np.arange(2 * n + 1) - n, np.arange(2 * n + 1) - n)
        r_sq_psf = (x**2 + y**2) * 1.
        dy, dx = [oldshape[1] // 2 + 0.5, oldshape[2] // 2 + 0.5]

        xp0 = (-r_avg * np.cos(phi) + dx).astype(int) - n
        yp0 = (-r_avg * np.sin(phi) + dy).astype(int) - n
        xp1 = (-r_avg * np.cos(phi - dphi) + dx).astype(int) - n
        yp1 = (-r_avg * np.sin(phi - dphi) + dy).astype(int) - n
        xp2 = (-r_avg * np.cos(phi + dphi) + dx).astype(int) - n
        yp2 = (-r_avg * np.sin(phi + dphi) + dy).astype(int) - n
        xp3 = (-r_avg * np.cos(phi - 2 * dphi) + dx).astype(int) - n
        yp3 = (-r_avg * np.sin(phi - 2 * dphi) + dy).astype(int) - n
        xp4 = (-r_avg * np.cos(phi + 2 * dphi) + dx).astype(int) - n
        yp4 = (-r_avg * np.sin(phi + 2 * dphi) + dy).astype(int) - n

        ######################################################################
        # Two references to the same array, so that we can modify it
        # with two indices and access it with one
        ######################################################################
        
        sub_arr_full = np.zeros((oldshape[1], oldshape[2]))
        sub_arr = np.reshape(sub_arr_full, -1)

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    
    if mode == 'LOCI':
        print "Running LOCI on {0} regions in {1} frames out to {2} pixels in radius.".format(totn, pa.size, int(locipar.rmax))
    elif mode == 'refine':
        print "Refining LOCI on {0} regions in {1} frames out to {2} pixels in radius.".format(totn, pa.size, int(locipar.rmax))
    else:
        print "Argument mode=" + mode + " to loci.py not recognized."
        sys.exit()
        
    ######################################################################
    # Run LOCI on each optimization section.  Do the sections in parallel,
    # with all of the frames processed in each job.  Split into multiple
    # execution blocks due to the large memory demands of the index
    # arrays.
    ######################################################################

    ######################################################################
    # Automatically set niter to use < 10% of the total memory.
    # In the 'refine' mode, this will be doubled, and incomplete freeing
    # seems to add a further factor of 2.
    ######################################################################

    if niter == 0:
        mem_iter = optreg * 4. * nframes / mem
        if mode == 'refine':
            niter = max(int(totn * mem_iter * 60), 1)
        else:
            niter = max(int(totn * mem_iter * 30), 1)

    for i in range(niter):

        tasks = multiprocessing.Queue()
        results = multiprocessing.Queue()
        ncpus = multiprocessing.cpu_count()
        consumers = [ Consumer(tasks, results)
                      for j in range(ncpus) ]
        for w in consumers:
            w.start()

        i1 = totn // niter * i
        i2 = totn // niter * (i + 1)
        if i == niter - 1:
            i2 = totn

        ##################################################################
        # Compute the indices in the original flux array corresponding
        # to the optimization and subtraction regions, using the 
        # sorted radius and angle arrays to make initial guesses.
        ##################################################################

        iopt, isub, nopt, nsub = annulus_indices(
            r_ex, r1[i1:i2], r2[i1:i2], theta1[i1:i2], theta2[i1:i2], 
            optreg, rindex, rsort, thetasort, p, i1, totn,
            locipar.fwhm, innerfrac=locipar.innerfrac) 
        
        for j in range(i1, i2):
            k = j - i1

            ##################################################################
            # Add pseudo-test sources to compute the fractional flux loss.
            # We need to do this independently for each region, so we will
            # remove these same sources after submitting the jobs.
            ##################################################################

            if do_partial_sub:
                n = r_sq_psf.shape[0]
                
                sig = np.sort([2.5, (dphi[j] * r_avg[j] / 3), 4])[1]
                psf = np.exp(-r_sq_psf / 10) + 0.05 * np.exp(-r_sq_psf / 30)
                psf1 = np.exp(-r_sq_psf / (2 * sig**2))
                psf2 = np.exp(-r_sq_psf / (2 * 2 * sig**2))
                psf1 *= np.sum(psf) / np.sum(psf1) / 3
                psf2 *= np.sum(psf) / np.sum(psf2) / 6
                
                sub_arr_full[yp0[j]:yp0[j] + n, xp0[j]:xp0[j] + n] += psf
                sub_arr_full[yp1[j]:yp1[j] + n, xp1[j]:xp1[j] + n] -= psf1
                sub_arr_full[yp2[j]:yp2[j] + n, xp2[j]:xp2[j] + n] -= psf1
                sub_arr_full[yp3[j]:yp3[j] + n, xp3[j]:xp3[j] + n] -= psf2
                sub_arr_full[yp4[j]:yp4[j] + n, xp4[j]:xp4[j] + n] -= psf2
                
                sub = sub_arr[iopt[k, :nopt[k]]].copy()
                
                sub_arr_full[yp0[j]:yp0[j] + n, xp0[j]:xp0[j] + n] = 0
                sub_arr_full[yp1[j]:yp1[j] + n, xp1[j]:xp1[j] + n] = 0
                sub_arr_full[yp2[j]:yp2[j] + n, xp2[j]:xp2[j] + n] = 0
                sub_arr_full[yp3[j]:yp3[j] + n, xp3[j]:xp3[j] + n] = 0
                sub_arr_full[yp4[j]:yp4[j] + n, xp4[j]:xp4[j] + n] = 0
            else:
                sub = None
            
            if corr is not None:
                corr_rad = corr[int(r1[j])]
            else:
                corr_rad = None
            if pca_arr is not None:
                pcaopt = pca_arr[:, iopt[k, 0:nopt[k]]]
                pcasub = pca_arr[:, isub[k, 0:nsub[k]]]
            else:
                pcaopt = None
                pcasub = None
                
            tasks.put(Task(j, loci_sub, (fluxsmooth[:, iopt[k, 0:nopt[k]]],
                                         flux[:, isub[k, 0:nsub[k]]],
                                         pa, minsep, r_avg[j], pcaopt, pcasub,
                                         partial_sub, ngroup,  method,
                                         corr_rad, sub))) 
    
        ##################################################################
        # "Poison Pills" to break execution
        ##################################################################
        
        for j in range(ncpus):
            tasks.put(None)

        for j in range(i1, i2):
            p.render([i2 * 100 / totn, (j + 1) * 100 / totn],
                     ['Calculating Indices for Region {0}'.format(i2),
                      'Processing Region {0} at Radius {1}'.format(j + 1, int(r1[j]))])
            index, result = results.get()
            if mode == 'LOCI':
                flux[:, isub[index - i1, 0:nsub[index - i1]]], \
                        sub_coefs[index] = result
            elif mode == 'refine':
                flux[:, isub[index - i1, 0:nsub[index - i1]]] = result
            
        del iopt, isub, nopt, nsub
        del tasks, results, consumers
        gc.collect()

    flux = np.reshape(flux, oldshape)
    if pca_arr is not None:
        pca_arr = np.reshape(pca_arr, (pca_arr.shape[0], oldshape[1], oldshape[2]))

    ##################################################################
    # Go from fractional flux loss points to a smoothed map
    ##################################################################
        
    if do_partial_sub:
        print "Calculating map of fractional flux loss."
        return partial_sub_map(r1, r2, theta1, theta2, sub_coefs, pa)
    else:
        return None
