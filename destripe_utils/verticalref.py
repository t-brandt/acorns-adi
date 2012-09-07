#!/usr/bin/env python
# Original filename: verticalref.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: 21 Dec 2010
# 
# Summary:  Construct a reference vertical striping pattern.
# 

import numpy as np

def verticalref(flux, refstripe=None, smoothwidth=300, bias_only=True):

    """
    Function verticalref takes two arguments:
    1.  A 2048 x 2048 array of flux values

    Optional arguments:
    2.  The stripe to smooth, constructing a reference stripe.  Ignored if
        bias_only=True, required otherwise.
    3.  The width of the Gaussian filter (default 300) -- see paper
    4.  Use only the reference pixels at the top and bottom of the detector
        (default True)

    verticalref convolves the input stripe with a Gaussian to suppress
    high frequency pickup.  It returns the smoothed stripe.
    """
    
    nsig = 5
    returntime = 3
    dimy, dimx = flux.shape
    if not bias_only:
        refrow = refstripe * 64

    ##################################################################
    # Extract the stripe as a 1-D array in the order it was read out.  
    # Take care to leave extra blank pixels for the start of a new row.
    ##################################################################

    pixels = np.ndarray((64 + returntime, dimx), np.float32)

    if refstripe % 2 == 0 and not bias_only:
        pixels[:] = flux[refrow:refrow + 64 + returntime]
    elif not bias_only:
        pixels[:] = flux[refrow + 63:refrow - 1 - returntime:-1]
    else:
        pixels[:] = np.nan
        pixels[:4] = 0.5 * (flux[:4] + flux[dimy - 1:dimy - 5:-1])
    
    pixels[64:64 + returntime, 4:dimx - 4] = np.nan

    ##################################################################
    # Mask deviant pixels, blank ('nan') pixels
    ##################################################################

    if bias_only:
        imax = 4
    else:
        imax = pixels.shape[0]

    ##################################################################
    # Sigma-reject to get the useful pixels.
    ##################################################################

    mu = 0
    sig = 1e4
    niter = 4

    for i in range(niter):
        pts = np.extract(np.abs(pixels[:imax] - mu) < nsig * sig, 
                         pixels[:imax])
        mu = np.mean(pts)
        sig = np.std(pts)
            
    pixmask = np.zeros(pixels.shape, np.float32)
    np.putmask(pixmask[:imax], np.abs(pixels[:imax] - mu) < nsig * sig, 1)
    np.putmask(pixels[:imax], pixmask[:imax] < 1, 0)

    pixelrow = np.reshape(pixels, -1, order='F')
    pixmask = np.reshape(pixmask, -1, order='F')
            
    ##################################################################
    # Convolve the pixel mask with a normalized Gaussian
    ##################################################################

    window = 3 * smoothwidth - np.arange(6 * smoothwidth + 1)
    window = np.exp(-(window * 1.0)**2 / (2 * smoothwidth**2))
        
    if bias_only:
        pixsmooth = np.zeros(pixelrow.shape, np.float32)
        pixnorm = np.zeros(pixelrow.shape, np.float32)
        for i in range(dimx):
            j = i * (64 + returntime)
            dj = window.shape[0] // 2 + 2
            i1 = max(0, j - dj)
            i2 = min(pixelrow.shape[0], j + dj)
            j1 = max(0, dj - j)
            j2 = window.shape[0] + 3 - max(0, j + dj - pixelrow.shape[0])
            
            pixsmooth[i1:i2] += np.convolve(pixelrow[j:j + 4], window)[j1:j2]
            pixnorm[i1:i2] += np.convolve(pixmask[j:j + 4], window)[j1:j2]
    else:
        pixsmooth = np.convolve(pixelrow, window, 'same')
        pixnorm = np.convolve(pixmask, window, 'same')
        
    pixsmooth /= pixnorm
    
    ##################################################################
    # Reshape and return the smoothed array
    ##################################################################

    pixsmooth = np.reshape(pixsmooth, (64 + returntime, -1), order = 'F')
    if refstripe % 2 == 1 or refstripe is None and bias_only:
        vstripe = pixsmooth[63::-1]
    else:
        vstripe = pixsmooth[0:64]
        
    return vstripe
    


