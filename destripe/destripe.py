#!/usr/bin/env python
#
# Original filename: destripe.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: 21 Dec 2010
#
# Summary: A set of routines for bias-subtracting, flat-fielding, and
# hot pixel masking of H2RG images

import re
import sys
import numpy as np
import pyfits as pyf
import warnings
from scipy.signal import medfilt2d
from destripe_utils import *

def destripe(frame, flat, hotpix, write_files, output_dir, bias_only,
             clean=True, storeall=True, r_ex=0, extraclean=True,
             full_destripe=True, do_horiz=True, PDI=False):

    """
    Function destripe takes two arguments:
    1.  A (usually) 2048 x 2048 array of flux values
    2.  A (usually) 2048 x 2048 flatfield image
    3.  The coordinates of the hot pixels to be masked
    4.  Write destriped data to files? 
    5.  Directory to write data (ignored if write_files=False)
    6.  Use only reference pixels?

    Optional arguments:
    7.  Interpolate over hot pixels?  default True
    8.  Store all frames in memory?  default True
    9.  Radial exclusion region for vertical destriping, default 0.
        Ignored if using only reference pixels.
    10. Mask deviant pixels by smoothing with a large median filter
        and looking for discrepancies?  default True
    11. Perform the full H2RG analysis? default True
    12. Calibrate the zero point in readout channels?  default True
    13. Use separate left and right channels, as for HiCIAO's PDI
        mode?  default False

    This function returns the destriped data.  It uses verticalmed,
    verticalref, horizontal, and interpbadpix from destripe_utils.
    """
        
    np.seterr(all='ignore')
    if not (storeall or write_files):
        print "Error: attempting to run destripe without saving files to either disk or memory"

    ncoadd = 1
    try:
        fluxfits = pyf.open(frame, "readonly")        
        header = fluxfits[0].header
        try:
            ncoadd = header['COADD']
        except:
            try:
                ncoadd = header['COADDS']
            except:
                ncoadd = 1
        flux = fluxfits[-1].data.astype(np.float32)
        dimy, dimx = flux.shape
        if hotpix is not None:
            flux[hotpix] = np.nan
    except:
        print "Error reading file " + frame
        exit()

    ##############################################################
    # reference voltage scaled by a number less than one provides
    # the best estimate of the vertical pattern, 0.87 in my tests.
    ##############################################################

    if do_horiz:
        try:
            for stripe in range(32):      
                horizontal(flux, stripe)
        except:
            print "Horizontal destriping failed on frame " + frame
            exit()
            
    ##############################################################
    # Calculate and subtract the vertical pattern.
    ##############################################################

    if full_destripe:
        if bias_only:
            sub_coef = 0.87
        else:
            sub_coef = 1
            
        try:
            if bias_only:
                oddstripe = verticalref(flux, 1)
                evenstripe = oddstripe[::-1, :]
            else:
                oddstripe, evenstripe = verticalmed(flux, flat, r_ex=r_ex,
                                                    PDI=PDI)
        except:
            print "Vertical destriping failed on frame " + frame
            #exit()
        
        for i in range(1, 33, 2):
            flux[64 * i:64 * i + 64] -= oddstripe * sub_coef
            flux[64 * i - 64:64 * i] -= evenstripe * sub_coef

    ##############################################################
    # Four rows on each edge are reference pixels--don't
    # flatfield them
    ##############################################################

        flux[4:-4, 4:-4] /= flat[4:-4, 4:-4]

        np.putmask(flux, flux < -1000, 0)
        np.putmask(flux, flux > 5e4 * ncoadd, np.nan)
    else:
        flux[4:-4, 4:-4] /= flat[4:-4, 4:-4]

    try:
        if clean:
            if extraclean:
                
                #############################################################
                # We'll be taking a median, so make half the bad pixels
                # inf and the other half ninf
                #############################################################
                
                np.putmask(flux[::2, :], np.isnan(flux[::2, :]), np.NINF)
                np.putmask(flux[1::2, :], np.isnan(flux[1::2, :]), np.inf)
                resid = medfilt2d(flux, 11)
                
                fluxresid = np.abs(flux - resid)
                sigval = medfilt2d(fluxresid, 9)
                
                #############################################################
                # Mask everything deviant by at least 3.5 'sigma'.  Since
                # sigval is a median, for Gaussian errors, this is
                # 3.5 * sqrt(2*ln(2)) ~= 4.1 sigma.
                #############################################################
                
                mask = fluxresid > 4.5 * sigval
                mask[:10] = 0
                mask[-10:] = 0
                mask[:, :10] = 0
                mask[:, -10:] = 0

                np.putmask(flux, mask, np.nan)
                np.putmask(flux, np.isinf(flux), np.nan)
                
            interpbadpix(flux, n=6)
    except:
        print "Cleaning bad pixels failed on frame " + frame
        sys.exit(1)
        
    ##############################################################
    # We don't want any NaNs or infs in the returned data
    ##############################################################

    np.putmask(flux, np.logical_not(np.isfinite(flux)), 0)

    if write_files:
        try:
            fluxout = pyf.HDUList()
            flux_hdu = pyf.PrimaryHDU(flux, header)
            fluxout.append(flux_hdu)
    
            outname = re.sub(".fits", "_ds.fits", frame)
            outname = re.sub(".*/", output_dir + "/", outname)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fluxout.writeto(outname, clobber=True)
                fluxout.close()
        except IOError, err:
            print err
            sys.exit(1)

    if storeall:
        return flux
    else:
        return 


