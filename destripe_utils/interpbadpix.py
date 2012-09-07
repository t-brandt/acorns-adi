#!/usr/bin/env python
# Original filename: interpbadpix.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: Jan 2012
# 
# Summary:  Replace bad pixels with the median of the surrounding good pixels
# 

import numpy as np

def interpbadpix(flux, n=5):

    """
    Function interpbadpix takes one argument:
    1.  A 2048 x 2048 array of flux values, with bad pixels flagged by NaN

    Optional argument:
    2.  Box size from which to select value (default 5, for a total 
        of 25 pixels)

    interpbadpix 'cleans' bad pixels using the median of surrounding
    good pixels.  It uses an n x n box (default 5x5) to select good
    pixels.  Bad pixels must be flagged by NaN.
    """
    
    dimy, dimx = flux.shape
    flux = np.reshape(flux, -1)

    ##################################################################
    # Map the bad pixels, create and populate an array to store
    # their neighbors
    ##################################################################

    badpix = np.logical_not(np.isfinite(flux))

    fluxmed = np.ndarray((n**2 - 1, np.sum(badpix)), np.float32)
    indx = np.extract(badpix, np.arange(flux.shape[0]))
    
    for i in range(0, n):
        if i <= n // 2 and i > 0:
            di = -i * dimy
        else:
            di = (i - n // 2) * dimy

        if i > 0:
            fluxmed[i - 1] = flux[indx + di]
            
        for j in range(1, n):
            if j <= n // 2:
                dj = -j
            else:
                dj = j - n // 2
            fluxmed[(i + 1) * (n - 1) + j - 1] = flux[indx + di + dj]

    ##################################################################
    # This is significantly faster than scipy.stats.nanmedian (I don't
    # know how that works).  Sort arrays, NaN are at the end.  Then
    # count the non-NaN values for each pixel and take their median.
    ##################################################################

    fluxmed = np.sort(fluxmed, axis=0)
    imax = np.sum(np.logical_not(np.isnan(fluxmed)), axis=0) 

    for i in range(1, np.amax(imax) + 1):
        indxnew = np.where(imax == i)
        if len(indxnew) > 0:
            flux[indx[indxnew]] = np.median(fluxmed[:i, indxnew], axis=0)

    flux = np.reshape(flux, (dimy, dimx))
    
    return

    

