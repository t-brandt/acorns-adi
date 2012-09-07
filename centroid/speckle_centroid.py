#!/usr/bin/env python
#
# Original filename: speckle_centroid.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: April 2011
# 
# Summary:  Estimate the centroid of an image by matching pairs of
# "speckles," or peaks in the intensity.
# 

import numpy as np
from scipy import ndimage, signal, optimize
import pyfits as pyf
import re

def chi2speckle(p, x, y):
    chi2 = np.zeros(x.shape)
    for i in range(x.shape[0]):
        chi2[i] = min((2 * p[1] - x[i] - x)**2 + (2 * p[0] - y[i] - y)**2)
        
    chi2 = sorted(chi2)
    neff = max(x.shape[0] // 2, 15)
    neff = max(1, neff)

    return np.sqrt(chi2[0:neff])

def speckle_centroid(frame, image=None, center=None):

    """
    
    Function speckle_centroid estimates the centroid of a (usually)
    saturated image by finding local maxima in intensity and
    attempting to identify symmetric pairs of maxima.  The function
    takes one argument:

    1.  The original filename.  If no flux array is supplied, it will
        be read from this file.

    Optional arguments:
    2.  A 2D flux array to centroid.
    3.  The first guess for the centroid.  Default is the center of
        the image.

    Function speckle_centroid returns the centroid as a list [yc, xc]
    
    """

    if image is None:
        frame_dw = re.sub(".fits", "_dw.fits", frame)
        image = pyf.open(frame_dw)[0].data

    if center is None:
        center = [image.shape[0] // 2, image.shape[1] // 2]
    yc = center[0] + 0.0
    xc = center[1] + 0.0
    
    #################################################################
    # Pull out the center of the image.  This requires a good
    # initial guess. 
    #################################################################

    x = np.linspace(-60., 60, 601) + xc
    y = np.linspace(-60., 60, 601) + yc
    x, y = np.meshgrid(x, y)

    imcenter = ndimage.map_coordinates(image, [y, x], order=1)

    #################################################################
    # Convolve with a 2-D Gaussian to look for local maxima.
    #################################################################

    sig = 2
    x2 = np.linspace(-2, 2, 21)
    y2 = np.linspace(-2, 2, 21)
    x2, y2 = np.meshgrid(x2, y2)
    r = np.sqrt(x2**2 + y2**2)
    psf = r < 1.8

    speckles = signal.convolve(imcenter, psf, mode='same')
    
    #################################################################
    # Pull out the coordinates of the local maxima.
    #################################################################

    nmax = 9
    domain = np.ones((nmax, nmax))
    specklemax = signal.order_filter(speckles, domain, nmax**2 - 1)
    
    maxima = np.where(np.all([specklemax == speckles, specklemax > 0], axis=0))
    
    x = np.reshape(x[maxima], (-1))
    y = np.reshape(y[maxima], (-1))

    newarr = np.zeros((2, x.shape[0]))
    newarr[0] = y
    newarr[1] = x
    
    #################################################################
    # Minimize the squared offsets of the best-fit points.  
    #################################################################

    p0 = [yc, xc]
    p1, success = optimize.leastsq(chi2speckle, p0[:], args=(x, y))

    return [p1[0], p1[1]] 

    

