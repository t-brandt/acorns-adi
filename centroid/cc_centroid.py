#!/usr/bin/env python
#
# Original filename: cc_centroid.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: April 2012
# 
# Summary:  Find the centroid of a (usually) saturated frame
# 

import numpy as np
from scipy import linalg, optimize
import pyfits as pyf
import warnings
import re
import sys


def errorfunc(p, y, x, chi2):
    return chi2 - (p[0] + p[1] * (y - p[3])**2 + p[2] * (x - p[4])**2)

def cc_centroid(refimage, image=None, frame=None, usemask=True, side=None):

    """
    function cc_centroid(refimage, image=None, frame=None)
    
    refimage should be a 2D or 3D numpy.ndarray.  If a 3D array,
    the first index runs over the template images.

    Must supply either image, a 2D numpy.ndarray to be centroided,
    or frame, the filename from which that image may be loaded.

    The function returns the centroid [yc, xc] if successful, None
    otherwise.

    Description:
    cc_centroid finds the centroid of the input image using the
    following algorithm:
    1. Flag saturated pixels, centroid the greatest concentration of
    such pixels to compute a provisional center.
    2. Mask pixels near the provisional center, compute a variance for
    all other pixels.  Variance = shot noise + read noise.
    3. Fit the PSF templates using \chi^2 at a grid of offsets.
    4. Centroid the map of \chi^2 merit statistics.

    """

    np.seterr(all='ignore')

    ####################################################################
    # Locate data if no frame supplied, load data
    ####################################################################

    if image is None and frame is None:
        print "Error: must supply either data or a filename to crosscorr_centroid."
        sys.exit(1)
    elif image == None:
        if not "_dw.fits" in frame:
            frame_dw = re.sub(".fits", "_dw.fits", frame)
        else:
            frame_dw = frame
        try:
            image = pyf.open(frame_dw)[-1].data
        except:
            frame_ds = re.sub(".fits", "_ds.fits", frame)
            try:
                image = pyf.open(frame_ds)[-1].data
            except:
                print "Error, cannot read data from " + frame_ds
                sys.exit(1)


    ####################################################################
    # Add the capability to only search the left or right half of the
    # image 
    ####################################################################

    image_save = image.copy()
    dimy, dimx = image.shape
    if side is not None:
        if re.search('[Ll]eft', side):
            image[:, dimx // 2:] = 0
        elif re.search('[Rr]ight', side):
            image[:, :dimx // 2] = 0                

    ####################################################################
    # Find approximate centroid by flagging (near-)saturated pixels
    # and locating the greatest concentration of them
    ####################################################################

    sat = min(image.max() * 0.7, 1e5)
    x = np.arange(image.shape[1])
    y = np.arange(image.shape[0])
    x, y = np.meshgrid(x, y)
    satpts = image > 0.8 * sat
    image = image_save
    
    maxpts = 0
    imax, jmax = [0, 0]
    for i in range(100, image.shape[0] - 100, 100):
        for j in range(100, image.shape[1] - 100, 100):
            npts = np.sum(satpts[i - 100:i + 100, j - 100:j + 100])
            if npts > maxpts:
                maxpts = npts
                imax, jmax = [i, j]

    ####################################################################
    # Check to see that this guess is in the central half of the FOV.
    # Then refine the estimate by calculating the mean position of the
    # (near-)saturated pixels in the neighborhood of the guess.
    # Do this iteratively, with the final estimate computed from a
    # 100x100 pixel region.
    ####################################################################

    di, dj = [image.shape[0] // 2, image.shape[1] // 2]
    if side is None and (np.abs(imax - di) > di / 2 or np.abs(jmax - dj) > dj / 2):
        return None # failure

    for di in range(100, 70, -10):
        npts = 1. * np.sum(satpts[imax - di:imax + di, jmax - di:jmax + di])
        yc = np.sum(satpts[imax - di:imax + di, jmax - di:jmax + di] *
                    y[imax - di:imax + di, jmax - di:jmax + di]) / npts
        xc = np.sum(satpts[imax - di:imax + di, jmax - di:jmax + di] *
                    x[imax - di:imax + di, jmax - di:jmax + di]) / npts
        try:
            imax, jmax = [int(yc), int(xc)]
        except:
            return None # failure

    ####################################################################
    # Calculate the typical saturation radius; cap at 700 mas
    ####################################################################

    dr_rms = np.sum(satpts[imax - di:imax + di, jmax - di:jmax + di] *
                    (y[imax - di:imax + di, jmax - di:jmax + di] - yc)**2)
    dr_rms += np.sum(satpts[imax - di:imax + di, jmax - di:jmax + di] *
                    (x[imax - di:imax + di, jmax - di:jmax + di] - xc)**2)
    dr_rms = np.sqrt(dr_rms / npts)
    dr_rms = min(dr_rms, 70)

    center = [imax, jmax]

    ####################################################################
    # Verify shape of reference PSF
    ####################################################################
        
    if len(refimage.shape) == 2:
        dimy, dimx = refimage.shape
        nref = 1
    elif len(refimage.shape) == 3:
        nref, dimy, dimx = refimage.shape
    else:
        print "Reference image must be a single 2D image or an array of 2D images."
        sys.exit(1)
    if dimy % 2 == 0 or dimx % 2 == 0 or dimy != dimx:
        print "Reference image to crosscorr_centroid must be square and\nhave an odd dimension."
        sys.exit(1)

    ####################################################################
    # Mask questionable data in the image, reshape arrays
    ####################################################################

    di = dimy // 2
    r_im = np.sqrt((x[imax - di:imax + di, jmax - di:jmax + di] - jmax)**2 +
                   (y[imax - di:imax + di, jmax - di:jmax + di] - imax)**2)
    mask = np.all([image < 0.5 * sat, image > 0], axis=0)

    baddata = np.all([image[imax-di:imax+di, jmax-di:jmax+di] < 0.2 * sat,
                      r_im < 2 * dr_rms], axis=0)

    if usemask:
        np.putmask(mask[imax - di:imax + di, jmax - di:jmax + di], 
                   r_im < 1.5 * dr_rms, 0)
        np.putmask(mask[imax - di:imax + di, jmax - di:jmax + di], 
                   baddata, 0)
    refimage2 = np.reshape(refimage, (nref, -1))

    sub_istd = np.ndarray(refimage2.shape)
    if usemask:
        istd = np.sqrt(mask / (np.abs(image) + 200))
    else:
        istd = np.sqrt(1 / (np.abs(image) + 200))

    ####################################################################
    # Produce an nxn map of chi2 as a function of offset.
    # Use SVD to do the fitting at each offset.
    ####################################################################  

    chi2_best = np.inf
    n = 21
    x = np.arange(n) - n // 2
    x, y = np.meshgrid(x, x)
    
    chi2 = np.zeros((n, n))
    ybest, xbest = [0, 0]
    for i in range(n):
        for j in range(n):
            y1 = center[0] + y[i, j] - dimy // 2
            x1 = center[1] + x[i, j] - dimx // 2

            subarr = np.reshape(image[y1:y1 + dimy, x1:x1 + dimx], -1)
            for k in range(nref):
                sub_istd[k] = np.reshape(istd[y1:y1 + dimy, x1:x1 + dimx], -1)
            A = sub_istd * refimage2
            b = sub_istd[0] * subarr
            coef = linalg.lstsq(A.T, b)[0]

            # Compute residuals, sum to get chi2
            resid = subarr - coef[0] * refimage2[0]
            for k in range(1, nref):
                resid -= coef[k] * refimage2[k]
            chi2[i, j] = np.sum((resid * sub_istd[0])**2)
                
            if chi2[i, j] < chi2_best:
                chi2_best = chi2[i, j]
                ibest, jbest = [i, j]

    ####################################################################
    # Take a 5x5 map around the best chi2, centroid this.
    # If that 5x5 map would be off the grid, return the initial guess.
    # If the centroiding fails (result falls outside the 5x5 grid),
    # return the [y, x] with the best chi2.
    ####################################################################      

    ybest0, xbest0 = [y[ibest, jbest], x[ibest, jbest]]
    p0 = [chi2_best, 2., 2., ybest0, xbest0]
    if ibest < 2 or ibest >= n - 3 or jbest < 2 or jbest >= n - 3:
        return None #failure

    x = np.reshape(x[ibest - 2:ibest + 3, jbest - 2:jbest + 3], -1)
    y = np.reshape(y[ibest - 2:ibest + 3, jbest - 2:jbest + 3], -1)
    chi2 = np.reshape(chi2[ibest - 2:ibest + 3, jbest - 2:jbest + 3], -1)
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p1, success = optimize.leastsq(errorfunc, p0[:], args=(y, x, chi2))
    ybest, xbest = [p1[3], p1[4]]
    
    if ybest > y.min() and ybest < y.max() and xbest > x.min() and xbest < x.max():
        return [center[0] + ybest, center[1] + xbest, dr_rms]
    else:
        return None # failure
    

