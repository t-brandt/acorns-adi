# Original filename: verticalmed.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: Jan 2012
# 
# Summary:  Construct a reference vertical striping pattern.
# 

import numpy as np

def verticalmed(flux, flat, r_ex=0, PDI=False):

    """
    Function verticalmed takes two arguments:
    1.  A 2048 x 2048 array of flux values
    2.  A 2048 x 2048 flat-field array

    Optional arguments:
    3.  Exclusion radius for calculting the median of the horizontal stripes
          (default zero, recommended values from 0 to 800)
          See Kandori-san's IDL routine for the equivalent.
    4.  Use separate left and right channels, as for HiCIAO's PDI
        mode?  default False
        
    verticalmed takes the median of the horizontal stripes to calculate a
    vertical template, as in Kandori-san's routine.  The routine ignores a
    circular region about the array center if r_ex > 0, and also ignores
    all pixels flagged with NaN.
    """

    ###############################################################
    # Construct radius array, copy flux to mask
    ###############################################################

    dimy, dimx = flux.shape
    x = np.arange(dimx)
    y = np.arange(dimy)
    x, y = np.meshgrid(x, y)

    if r_ex > 0:
        if not PDI:
            r_ok = ((x - dimx / 2)**2 + (y - dimy / 2)**2) > r_ex**2
        else:
            r_ok = ((x - dimx / 4)**2 + (y - dimy / 2)**2) > r_ex**2
            r_ok *= ((x - 3 * dimx / 4)**2 + (y - dimy / 2)**2) > r_ex**2
            
        flux2 = np.ndarray(flux.shape, np.float32)
        flux2[:] = flux
        np.putmask(flux2, np.logical_not(r_ok), np.nan)
    else:
        flux2 = flux

    ###############################################################
    # Estimate background level
    ###############################################################

    backgnd = np.ndarray(flux2.shape)
    backgnd[:] = flux2 / flat
    backgnd = np.sort(np.reshape(backgnd, -1))
    ngood = np.sum(np.isfinite(backgnd))
    level = np.median(backgnd[:ngood])
    flux2 -= level * flat

    ###############################################################
    # Sort the flux values.  NaN values will be at the end for
    # numpy versions >= 1.4.0; otherwise this routine may fail
    ###############################################################
    
    tmp = np.ndarray((32, dimy // 32, dimx), np.float32)
    for i in range(1, 33, 2):
        tmp[i] = flux2[64 * i:64 * i + 64]
        tmp[i - 1] = flux2[64 * i - 64:64 * i]
        
    tmp = np.sort(tmp, axis=0)
    oldshape = tmp[0].shape
    tmp = np.reshape(tmp, (tmp.shape[0], -1))
    
    oddstripe = np.zeros(tmp[0].shape, np.float32)

    ###############################################################
    # imax = number of valid (!= NaN) references for each pixel.
    # Calculate the median using the correct number of elements,
    # doing it only once for each pixel.
    ###############################################################
    
    imax = np.sum(np.logical_not(np.isnan(tmp)), axis=0)
    for i in range(np.amin(imax), np.amax(imax) + 1):
        indxnew = np.where(imax == i)
        if len(indxnew) > 0:
            oddstripe[indxnew] = np.median(tmp[:i, indxnew], axis=0)

    ###############################################################
    # Set the median of the pattern to be subtracted equal to the 
    # median of the difference between the science and reference
    # pixels.
    ###############################################################

    oddstripe -= np.median(oddstripe)
    oddstripe += 0.5 * (np.median(flux[:4]) + np.median(flux[-4:]))
    
    oddstripe = np.reshape(oddstripe, oldshape)
    evenstripe = oddstripe[::-1]

    return [oddstripe, evenstripe]

