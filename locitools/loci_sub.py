#!/usr/bin/env python
#
# Original filename: loci_sub.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: July 2011
# 
# Summary:  Calculate the optimal linear combination of frames for
# PSF subtraction in a subregion in LOCI, perform subtraction,
# calculate fractional flux loss.  
# 

import sys
import numpy as np
from scipy import linalg
import warnings

def loci_sub(annulus_opt, annulus_sub, pa, minsep, r0, pcaopt=None,
             pcasub=None, partial_sub=None, ngroup=1, method='matrix',
             corr=None, sub_arr=None):
    
    """
    
    Function loci_sub calculates the optimal linear combination of PSF
    reference frames for each saturated frame and subtracts that
    combination.  Run loci_sub on each region of an input image to
    implement the LOCI algorithm.

    Input arguments:
    1.  A 2D array of optimization regions.  Shape:  (nframes, npixels_opt)
    2.  A 2D array of subtraction regions.  (nframes, npixels_sub)
    3.  The parallactic angles to calculate angular separation.  (nframes)
    4.  Minimum angular separation for subtraction in pixels
    5.  Radius at which to check angular separation

    Optional arguments:
    6.  'pcaopt': a 2D array of comparison regions from PSF templates.
        Ignored if None.  Default None.
    7.  'pcasub': a 2D array of subtraction regions from PSF templates.
        Ignored if None.  Default None.
    8.  'partial_sub': a UnivariateSpline set of ticks for interpolation
        of the fractional flux suppression.  Called with the separation in
        pixels.  Ignored if None.  Default None.
    9.  'ngroup': Number of independent groups of frames to reduce.
        There will be nframes / ngroup frames in each LOCI reduction.
        Default 1.
    10. 'method': Method to solve for the LOCI coefficients.  'matrix':
        construct normal equations, solve with LU decomposition.  'lstsq':
        solve with SVD.  'eqcoef': all coefficients are equal; similar to
        a mean PSF subtraction.  Default 'matrix'.
    11. 'corr': matrix of cross-correlation coefficients to select best
        matches in frames for LOCI references.  Ignored if None.
        Default None.
    12. 'sub_arr': array of extra sources to compute the fractional flux
        from LOCI.  Ignored if None.  Default None.    
    
    """
    ######################################################################
    # Duplicate the subtraction regions, append reference frames
    # Cast optimization regions to double precision floating point
    # Append reference frames, set reference parallactic angles to inf
    ######################################################################
    
    np.seterr(all='ignore')
    nframes = pa.shape[0]    
    subshape = np.asarray(annulus_sub.shape)
    optshape = np.asarray(annulus_opt.shape)
    sub_coefs = np.ones(nframes, np.float32)

    if pcasub is not None:
        subshape[0] += pcasub.shape[0]
        optshape[0] += pcaopt.shape[0]
        
    fluxsub = np.ndarray(tuple(subshape), np.float64)
    fluxopt = np.ndarray(tuple(optshape), np.float64)
    pa_full = np.zeros(optshape[0])

    fluxsub[:annulus_sub.shape[0]] = annulus_sub
    fluxopt[:annulus_opt.shape[0]] = annulus_opt
    pa_full[:pa.shape[0]] = pa
    
    if pcaopt is not None:
        fluxsub[annulus_sub.shape[0]:] = pcasub
        fluxopt[annulus_opt.shape[0]:] = pcaopt
        pa_full[pa.shape[0]:] = np.inf
        
    ######################################################################
    # Cannot use the normal equations redundantly if the correlation
    # matrix is used.  Set solver to use SVD.
    ######################################################################

    if ngroup > 1 and corr is not None:
        ngroup = 1

    for i in range(ngroup):
        indx = range(i, nframes, ngroup)
        if pcaopt is not None:
            indx += range(nframes + 1, optshape[0])
        indx = np.asarray(indx)
        n = len(indx)
        
        ######################################################################
        # Multiply out the problem matrix, and the trial source (if not None)
        ######################################################################

        if method == 'matrix':
            bigmat = np.dot(fluxopt[indx], fluxopt[indx].T)
            if sub_arr is not None:
                submat = np.ndarray(n) 
                all_y = np.ndarray(n)
                for j in range(n):
                    submat[j] = np.sum(sub_arr * fluxopt[indx[j]])
                    all_y[j] = np.sum(fluxopt[indx[j]])
            else:
                submat = None

        for iframe in range(i, nframes, ngroup):

            ##################################################################
            # Fetch the indices of the comparison frames that satisfy the
            # angular displacement criterion
            ##################################################################

            if corr is None:
                rot_ok = np.abs(pa_full[indx] - pa_full[iframe]) * r0 > minsep
                optframes = np.extract(rot_ok, np.arange(n))
                opt_ref = np.extract(rot_ok, indx)
                padiff = np.extract(rot_ok, pa_full[indx] - pa_full[iframe])
                sep = r0 * np.sqrt(2 - 2 * np.cos(padiff))
            else:
                opt_ref = np.extract(np.abs(pa_full[corr[iframe]] -
                                            pa_full[iframe]) * r0 > minsep,
                                     corr[iframe])
                
            ##################################################################
            # Solve the linear system using the user-input method, either
            # the normal equations, SVD, or equal coefficients (mean
            # subtraction).  If LU decomposition fails (matrix is singular),
            # warn the user and use SVD.
            ##################################################################

            if opt_ref.size > 0:
                if method == 'matrix':
                    A = bigmat[optframes[:], :][:, optframes[:]]
                    b = bigmat[optframes[:], :][:, iframe // ngroup] 
                    if submat is not None:
                        sub = submat[optframes[:]]
                        y = all_y[optframes[:]]

                    with warnings.catch_warnings():
                        warnings.simplefilter('error')
                        try:
                            lu, piv = linalg.lu_factor(A)
                            coef = linalg.lu_solve((lu, piv), b)
                            if submat is not None:
                                sub_coef = linalg.lu_solve((lu, piv), sub)
                        
                        except:
                            #print err
                            A = fluxopt[opt_ref[:]]
                            b = fluxopt[iframe]
                            coef = linalg.lstsq(A.T, b)[0]
                            if sub_arr is not None:
                                sub_coef = linalg.lstsq(A.T, sub_arr)[0]
                                
                    
                elif method == 'lstsq':
                    A = fluxopt[opt_ref[:]]
                    b = fluxopt[iframe]
                    coef = linalg.lstsq(A.T, b)[0]
                    if sub_arr is not None:
                        sub_coef = linalg.lstsq(A.T, sub_arr)[0]
                    
                elif method == 'eqcoef':
                    coef = np.ones(opt_ref.shape)
                    coef /= np.sum(coef)
                    sub_coefs = np.zeros(opt_ref.shape)
                else:
                    print "Method " + method + "not recognized in " + \
                        "subroutine locitools.loci_sub ."
                    sys.exit(1)

                if partial_sub is not None:
                    sub_coefs[iframe] -= np.sum(coef * partial_sub(sep))
                    if submat is not None:
                        
                        ######################################################
                        # Compute the fractional loss in flux within
                        # 0.5*max isophot, multiply by the flux loss from
                        # angularly displaced copies of the source
                        ######################################################
                        
                        loss = sub_coef[0] * fluxopt[opt_ref[0]]
                        for j in range(1, sub_coef.size):
                            loss += sub_coef[j] * fluxopt[opt_ref[j]]
                            
                        mask = sub_arr >= 0.5 * sub_arr.max()
                        loss = np.sum(loss * mask) / np.sum(mask * sub_arr)
                        
                        sub_coefs[iframe] *= 1 - loss

                for j in range(coef.size):
                    annulus_sub[iframe] -= coef[j] * fluxsub[opt_ref[j]]
            else:
                annulus_sub[iframe] = np.float('nan')

    return annulus_sub, sub_coefs

