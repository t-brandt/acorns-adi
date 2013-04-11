#!/usr/bin/env python
#
# Original filename: arr_resize.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: Dec 2012
#
# Summary: Resize a square array.  Crop or pad with NaN as necessary.
# If padding, create and return a new array.
#

import numpy as np

def arr_resize(arr, newdim=None, padval=np.nan):

    """
    Function arr_resize takes a 2D array and returns a square array
    with the same data, either cropped to the given new dimension, or
    centered in a newly created array.  It takes one argument:

    1. arr, a 2D numpy ndarray
    
    Two optional arguments:

    2. newdim: new dimension.  Forced to be a positive, odd integer.
       Default: the smallest multiple of 500 plus 1 at least as large
       as the largest dimension of arr
    3. padval: floating point value with which to pad the new array.
       Default: np.nan

    """

    assert len(np.asarray(arr).shape) == 2, "Function arr_resize expects a 2D array"

    arr = np.asarray(arr)
    dimy, dimx = arr.shape

    ####################################################################
    # Default new size is at least as large as old size, and is a
    # multiple of 500 plus 1.  Must be odd.
    ####################################################################

    if newdim is None:
        newdim = ((np.max([dimy, dimx]) - 2) // 500 + 1) * 500 + 1
    else:
        try:
            newdim = np.absolute(newdim // 2) * 2 + 1
        except:
            print "Error: new dimension in arr_resize must be an integer."
            return None
    try:
        padval = float(padval)
    except:
        print "Error: padding value in arr_resize must be a float."
        return None
    
    ####################################################################
    # Return a subarray of the original array is possible, otherwise
    # make a larger array and copy into it.
    ####################################################################

    dx1 = dimx // 2
    dy1 = dimy // 2
    dx2 = dy2 = newdim // 2

    if newdim < dimy and newdim < dimx:
        return arr[dy1 - dy2:dy1 + dy2 + 1,
                   dx1 - dx2:dx1 + dx2 + 1]
    else:
        newarr = np.ones((newdim, newdim), dtype=arr.dtype) * padval
        if newdim > dimy and newdim > dimx:
            newarr[dy2 - dy1:dy2 + dy1 + 1,
                   dx2 - dx1:dx2 + dx1 + 1] = arr
        elif newdim > dimy:
            newarr[dy2 - dy1:dy2 + dy1 + 1, :] = \
                       arr[:, dx1 - dx2:dx1 + dx2 + 1]
        else:
            newarr[:, dx2 - dx1:dx2 + dx1 + 1] = \
                       arr[dy1 - dy2:dy1 + dy2 + 1]
        return newarr
