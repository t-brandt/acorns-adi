#!/usr/bin/env python
#
# Original filename: sorted_arrays.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: August 2011
# 
# Summary:  
# 
import numpy as np
import sys

def sorted_arrays(arrshape):

    """
    Function sorted_arrays takes a single argument:
    1.  arrshape:  tuple, second and third values are the shape
        of the the arrays to sort and return

    It returns a list of three arrays: [rindex, rsort, thetasort]
    rindex:  the indices of an array of shape arrshape, which has been
        sorted by distance from the center
    rsort:  the (sorted) distances from the array center.  Reshaped to 1D.
    thetasort:  the angles corresponding to rindex. 
    
    """

    ######################################################################
    # Create radius, polar angle arrays, make them one-dimensional.
    # Force theta to run from 0 to 2 pi (default is -pi to pi).
    ######################################################################
        
    x = np.arange(arrshape[2] * 1.) - arrshape[2] // 2
    y = np.arange(arrshape[1] * 1.) - arrshape[1] // 2
    x, y = np.meshgrid(x, y)

    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x) + np.pi
    r = np.reshape(r, (-1))
    theta = np.reshape(theta, (-1))
    
    ######################################################################
    # Sort everything by radius.  rindex maps sorted arrays to
    # unsorted arrays, which makes it easy to compute annulus indices.
    # zip returns tuples, which are converted to arrays.
    ######################################################################
        
    rindx = np.arange(r.shape[0])
    rsort, rindx = zip(*sorted(zip(r, rindx)))
    rindex = np.arange(r.shape[0])
    rindex[:] = rindx[:]
    rsort = r[rindex[:]]
    thetasort = theta[rindex[:]]

    return [rindex, rsort, thetasort]



