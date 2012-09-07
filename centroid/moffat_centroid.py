#!/usr/bin/env python
#
# Original filename: moffat_centroid.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: April 2011
# 
# Summary:  Find the center of light of an input image.
# 

from scipy import optimize, stats
import numpy as np

def r2(p, x, y):
    return p[4] * (x - p[3])**2 + p[5] * (y - p[2])**2

def moffat(p, rsq, satlevel):
    return np.clip(p[0] + p[1] / np.abs(1 + rsq)**p[6], 0, satlevel)

def errfunc(p, x, y, f, satlevel):
    rsq = r2(p, x, y)   
    return (moffat(p, rsq, satlevel) - f) 
    
def moffat_centroid(frame, flux=None, lastcen=None):

    """

    Function moffat_centroid fits a Moffat profile to an image.  It
    returns the centroid of the best-fit profile as a list [y, x].
    Arguments:

    1.  The original filename.  If no flux array is supplied, it will
        be read from this file.

    Optional arguments:
    2.  A 2D flux array to centroid.
    3.  The first guess for the centroid.  Default is the center of
        the image.

    """
    
    if flux == None:
        frame_dw = re.sub(".fits", "_dw.fits", frame)
        flux = pyfits.open(frame_dw)[0].data

    if lastcen == None:
        lastcen = [flux.shape[0] // 2, flux.shape[1] // 2]
    lastcen[0] = int(lastcen[0])
    lastcen[1] = int(lastcen[1])

    #################################################################
    # Pull out the central part of the image.  Reduce the dimensions
    # by half, for 1/4 of the pixels.  
    #################################################################

    dim = flux.shape[0] // 4
    if np.mod(dim, 2) == 0:
        dim += 1

    flux2 = np.ndarray((dim, dim))
    flux2[:, :] = flux[lastcen[0] - dim // 2:lastcen[0] + dim // 2 + 1,
                       lastcen[1] - dim // 2:lastcen[1] + dim // 2 + 1]
                       
    x = np.linspace(0, dim - 1., dim) - dim // 2
    y = np.linspace(0, dim - 1., dim) - dim // 2
    x, y = np.meshgrid(x, y)

    #################################################################
    # Estimate the errors as photon+read noise, fit a Moffat profile.
    # Apply an intensity cut at the estimated saturation level;
    # otherwise, the method will fail.
    #################################################################

    flux2 = np.reshape(flux2, (-1))
    f_ierr = np.ndarray(flux2.shape, np.float)
    f_ierr[:] = 1 / np.sqrt(flux2**2 + 250)
    x = np.reshape(x, (-1))
    y = np.reshape(y, (-1))
    satlevel = stats.scoreatpercentile(flux2, 99)

    p0 = [0., 1e6, 0., 0., 5e-4, 5e-4, 1.8]
    
    p1, success = optimize.leastsq(errfunc, p0[:],
                                   args=(x, y, flux2, satlevel))

    center = [p1[2] + lastcen[0], p1[3] + lastcen[1]]
    
    return center

