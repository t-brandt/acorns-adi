#!/usr/bin/env python
#
# Original filename: addsource.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: Feb 2012
# 
# Summary: inject a fake source into an ADI sequence
#

import numpy as np
import random
from scipy import ndimage
from parallel import *

def addsource(flux, pos, pa, psf, norm=1, jitter=0):

    """
    addsource assumes the frames are centered.  It will add a source
    at position is [y, x] relative to the center.

    addsource takes four arguments:
    1.  A 2D or 3D flux array.  If 3D, the first index should be the
        frame number.
    2.  The position [y, x], relative to the center, at which to
        add the source.  This can be a list of positions, one for
        each frame.
    3.  The position angle (in radians) at which to add the source.
        If there is more than one frame, this should be a list or array.
    4.  The source PSF, as a 2D array

    Optional arguments:
    5.  The normalization of the additional source's PSF
    6.  The random positional jitter, in pixels (both horizontal
        and vertical) of the additional source

    addsource returns the 

    """

    ####################################################################
    # Rotate to the appropriate position angle
    ####################################################################

    x = pos[1] * np.cos(pa) - pos[0] * np.sin(pa) + flux.shape[-1] // 2
    y = pos[1] * np.sin(pa) + pos[0] * np.cos(pa) + flux.shape[-2] // 2
    normpsf = psf * norm

    try:
        
        ################################################################
        # Case 1: image sequence (flux is a 3D array)
        ################################################################

        for i in range(x.size):
            x[i] += random.gauss(0, jitter)
            y[i] += random.gauss(0, jitter)   
        
        ################################################################
        # x, y > 0.  Decompose into integer, fractional parts.
        ################################################################
        
        xdiff = x % 1
        x = x.astype(int)
        ydiff = y % 1
        y = y.astype(int)

        x1 = np.zeros(x.shape, int)
        x2 = np.zeros(x.shape, int)
        y1 = np.zeros(x.shape, int)
        y2 = np.zeros(x.shape, int)
        z1 = np.zeros(x.shape, int)
        z2 = np.zeros(x.shape, int)
        w1 = np.zeros(x.shape, int)
        w2 = np.zeros(x.shape, int)
        
        xref = np.arange(psf.shape[1])
        yref = np.arange(psf.shape[0])
        xref, yref = np.meshgrid(xref - 0., yref - 0.)

        ################################################################
        # Interpolate the template PSF onto the fractional part, add
        # it to the correct part of the array.
        ################################################################

        for i in range(x.size):
            x1[i] = x[i] - min([x[i], psf.shape[1] // 2])
            x2[i] = x[i] + min([flux.shape[-1] - x[i], psf.shape[1] // 2])
            y1[i] = y[i] - min([y[i], psf.shape[0] // 2])
            y2[i] = y[i] + min([flux.shape[-2] - y[i], psf.shape[0] // 2])
            z1[i] = x1[i] + psf.shape[1] // 2 - x[i]
            z2[i] = x2[i] + psf.shape[1] // 2 - x[i]
            w1[i] = y1[i] + psf.shape[0] // 2 - y[i]
            w2[i] = y2[i] + psf.shape[0] // 2 - y[i]

            newpsf = ndimage.map_coordinates(normpsf, [yref - ydiff[i], xref - xdiff[i]], order=3)

            flux[i, y1[i]:y2[i], x1[i]:x2[i]] += newpsf[w1[i]:w2[i], z1[i]:z2[i]]
        return flux
            
    except:
        
        ################################################################
        # Same algorithm for a single image.
        ################################################################

        xdiff = x % 1
        x = x.astype(int)
        ydiff = y % 1
        y = y.astype(int)

        x1 = x - min([x, psf.shape[1] // 2])
        x2 = x + min([flux.shape[-1] - x, psf.shape[1] // 2])
        y1 = y - min([y, psf.shape[0] // 2])
        y2 = y + min([flux.shape[-2] - y, psf.shape[0] // 2])
        z1 = x1 + psf.shape[1] // 2 - x
        z2 = x2 + psf.shape[1] // 2 - x
        w1 = y1 + psf.shape[0] // 2 - y
        w2 = y2 + psf.shape[0] // 2 - y
        
        xref = np.arange(psf.shape[1])
        yref = np.arange(psf.shape[0])
        xref, yref = np.meshgrid(xref - xdiff, yref - ydiff)

        newpsf = ndimage.map_coordinates(normpsf, [yref, xref], order=3)
        flux[y1:y2, x1:x2] += newpsf[w1:w2, z1:z2]
        

    return flux




    
