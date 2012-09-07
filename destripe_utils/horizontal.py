#!/usr/bin/env python
# Original filename: horizontal.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: January 2011
# 
# Summary:  Compute the best-fit offset between two stripes, subtract 
# from the first stripe supplied.  
# 

import numpy as np
import sys

def horizontal(flux, stripe, minx=0, xdist=4):

    """
    Function horizontal takes three arguments:
    1.  A 2048 x 2048 array of flux values
    2.  The stripe to correct (0 - 31)

    Optional arguments:
    3.  x pixel from which to start counting (default = 0)
    4.  Width of x pixel range to be used (default = 4)

    horizontal finds the best-fit difference (bias) between two stripes.  
    It then subtracts that bias from the first stripe supplied.
    """
    
    dimy, dimx = flux.shape
    row = stripe * 64

    ##################################################################
    # Even/odd rows are read in opposite directions
    ##################################################################

    diffdist_l = flux[row:row + 64, minx:minx + xdist]
    diffdist_r = flux[row:row + 64, dimx - xdist - minx:dimx - minx]

    diffdist = np.hstack((diffdist_l, diffdist_r))

    ##################################################################
    # Gaussian data ==> use mean.  Exclude > nsig sigma outliers as bad
    # pixels.  Loop over the data niter times to get a good guess.
    ##################################################################
    
    mu = 0
    nsig = 3
    niter = 4
    sig = 1e4
    for i in range(niter):
        pts = np.extract(np.abs(diffdist - mu) < nsig * sig, diffdist)
        mu = np.mean(pts)
        sig = np.std(pts)
    
    for i in range(row, row + 64):
        flux[i] -= mu
    
    return 

