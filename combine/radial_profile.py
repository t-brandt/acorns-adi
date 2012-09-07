#!/usr/bin/env python
#
# Original filename: radial_profile.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: June 2011
# 
# Summary:  Construct a radial profile of a given quantity, implement an
# optimized trimmed mean
# 

import sys
import numpy as np
from progressbar import ProgressBar

def radprof(flux, mode='median', rsort=None, rindex=None, r_ref=None,
                       center=None, rmax=None, smoothwidth=2, sigrej=None):

    np.seterr(all='ignore')
    assert smoothwidth >= 0 and smoothwidth < 10, \
           "Smoothing width must be between 0 and 10 pixels."
    assert len(flux.shape) == 2, "Input array must be two-dimensional."

    dimy, dimx = flux.shape
    flux = np.reshape(flux, -1)
    if center is None:
        center = [dimy // 2, dimx // 2]

    ##################################################################
    # Construct, sort radial profile arrays.
    ##################################################################

    if rsort is None or rindex is None or r_ref is None:
        x = np.arange(dimx) - dimx // 2
        y = np.arange(dimy) - dimy // 2
        x, y = np.meshgrid(x, y)
        r_ref = np.sqrt(x**2 + y**2)
        r = np.reshape(np.sqrt(x**2 + y**2).astype(int), -1)
        rindx = np.arange(r.shape[0])
        rsort, rindx = zip(*sorted(zip(r, rindx)))
        
        rindex = np.arange(r.shape[0])
        rindex[:] = rindx
        rsort = r[rindex]

    rmax1 = min(dimx // 2, dimy // 2)
    rmax2 = int(np.sqrt((dimx // 2)**2 + (dimy // 2)**2))
    
    radprof = np.zeros(dimy + dimx, np.float32)
    radprof2d = np.zeros((dimy, dimx), np.float32)

    for i in range(1, rmax2):
        if i < rmax1:
            i1 = int(np.pi * (i - 1)**2)
            i2 = int(np.pi * (i + 1)**2)
        else:
            area = 0
            if i >= dimx // 2:
                theta = np.arccos((dimx // 2) * 1. / i)
                area += 4 * i * theta - 2 * (dimx // 2) * i * np.sin(theta)
            if i >= dimy // 2:
                phi = np.arccos((dimy // 2) * 1. / i)
                area += 4 * i * phi - 2 * (dimy // 2) * i * np.sin(phi)
            i1 = int(np.pi * (i - 1.5)**2 - area)
            i2 = int(np.pi * (i + 1.5)**2 - area)

        if sigrej is not None:
            try:
                sig = 1e5
                for j in range(4):
                    indx = np.extract(np.abs(rsort[i1:i2] - i) <= 0.5,
                                     rindex[i1:i2])
                    indices = indx[np.where(np.abs(flux[indx]) < sigrej * sig)]
                    sig = np.std(flux[indices])
            except:
                indices = np.extract(np.abs(rsort[i1:i2] - i) <= 0.5,
                                     rindex[i1:i2])
        else:
            indices = np.extract(np.abs(rsort[i1:i2] - i) <= 0.5, rindex[i1:i2])
        if mode == 'median':
            radprof[i] = np.median(flux[indices])
        elif mode == 'mean':
            radprof[i] = np.mean(flux[indices])
        elif mode == 'std':
            radprof[i] = np.std(flux[indices])
        elif mode == 'var':
            radprof[i] = np.var(flux[indices])
        else:
            print "Mode '" + mode + "' not recognized in routine combine.radprof"
            sys.exit(1)
            
    if smoothwidth > 0:
        window = np.exp(-(3*smoothwidth - np.arange(6*smoothwidth + 1))**2 /
                           (2.0*smoothwidth**2))
        window /= np.sum(window)
        radprof_save = np.zeros(radprof.shape, np.float32)
        radprof_save[:] = radprof  
        mask = np.zeros(radprof.shape, np.int)
    
        radprof = np.convolve(radprof, window, mode="same")
        np.putmask(mask, np.abs(radprof) < float('inf'), 1)
        np.putmask(radprof, mask == 0, radprof_save)
        np.putmask(radprof, radprof_save == 0, radprof_save)

    radprof2d = radprof[(r_ref + 0.5).astype(int)]
    flux = np.reshape(flux, (dimy, dimx))
 
    return radprof2d, rsort, rindex, r_ref




def meanmed(flux, nmed=10, smoothwidth=2):

    nframes, dimy, dimx = flux.shape
    nmed = min(nmed, nframes // 2)

    fluxsort = np.sort(flux, axis=0)
    meanarr = np.zeros((nmed, dimy, dimx), np.float32)
    n1 = np.zeros(nmed, int)
    n2 = np.zeros(nmed, int)

    # Median first.  This case handles either odd or even n
    if nframes % 2 == 0:
        n1[0] = nframes // 2 - 1
        n2[0] = n1[0] + 2
    else:
        n1[0] = nframes // 2
        n2[0] = n1[0] + 1

    meanarr[0] = np.sum(fluxsort[n1[0]:n2[0]], axis=0)
    for i in range(1, nmed):
        # 1 <= dn, minimum index referenced >= 1
        dn = int(i / (nmed - 0.5) * ((nframes - 1) // 2))
        n1[i] = n1[0] - dn
        n2[i] = n2[0] + dn
        meanarr[i] = meanarr[i - 1] + np.sum(fluxsort[n1[i]:n1[i - 1]], axis=0)
        meanarr[i] += np.sum(fluxsort[n2[i - 1]:n2[i]], axis=0)

    for i in range(nmed):
        meanarr[i] /= n2[i] - n1[i]

    noise = np.zeros(meanarr.shape, np.float32)
    
    noise[0], rsort, rindex, r_ref = radprof(meanarr[0], mode='std')
    for i in range(1, nmed):
        noise[i] = radprof(meanarr[i], mode='std', rsort=rsort,
                                      rindex=rindex, r_ref=r_ref)[0]

    noisebest = np.nanmin(noise, axis=0)
    im_best = np.zeros((dimy, dimx), np.float32)
    for i in range(nmed - 1, -1, -1):
        np.putmask(im_best, noise[i] == noisebest, meanarr[i])

    return im_best, noisebest.astype(np.float32)
