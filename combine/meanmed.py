#!/usr/bin/env python
#
# original filename: meanmed.py
#
# Author: Tim Brandt
# email: tbrandt@astro.princeton.edu
# Date: July 2012
#
# Description: Implement the mean/median hybrid estimator described in
# Brandt+ (2012).  Combine an image sequence using the optimal estimator
# at each separation from the central star.
#

import numpy as np

def meanmed(flux, smoothwidth=2):

    """


    
    """

    np.seterr(all='ignore')
    nframes, dimy, dimx = flux.shape

    center = [dimy // 2, dimx // 2]

    ##################################################################
    # Construct, sort radial profile arrays.
    # Default max radius is half the width of the domain.
    ##################################################################

    x = np.arange(dimx) - dimx // 2
    y = np.arange(dimy) - dimy // 2
    x, y = np.meshgrid(x, y)
    r = np.reshape(np.sqrt(x**2 + y**2).astype(int), -1)
    rindx = np.arange(r.shape[0])
    rsort, rindx = zip(*sorted(zip(r, indx)))

    rindex = np.arange(r.shape[0])
    rindex[:] = rindx
    rsort = r[rindex]

    rmax1 = min(dimx // 2, dimy // 2)
    rmax2 = int(np.sqrt((dimx // 2)**2 + (dimy // 2)**2))

    for i in range(1, rmax2):
        if i < rmax1:
            i1 = int(math.pi * (i - 1)**2)
            i2 = int(math.pi * (i + 1)**2)
        else:
            i1 = int(math.pi * (rmax1 - 1)**2)
            i2 = None

        indices = np.extract(np.abs(rsort - i) <= 0.5, rindex)
    
    rmax = dimx // 2
    
    fluxsort = np.extract(r < rmax, fluxsort)
    r_tmp = np.extract(r < rmax, r)

    r_tmp, fluxsort = zip(*sorted(zip(r_tmp, fluxsort)))



    x_tmp = np.arange(dimx * 1.) - center[1]
    y_tmp = np.arange(dimy * 1.) - center[0]
    x, y = np.meshgrid(x_tmp, y_tmp)

    r = np.ndarray(flux.shape, np.int)
    radprof = np.zeros(dimy + dimx, np.float32)
    radprof2d = np.zeros((dimy, dimx), np.float32)
    fluxsort = np.ndarray(dimy, dimx), np.float32)
    fluxsort[:, :] = flux
    r[:, :] = np.sqrt(x**2 + y**2)

    if rmax is None:
        rmax = dimx // 2
    fluxsort = np.extract(r < rmax, fluxsort)
    r_tmp = np.extract(r < rmax, r)

    r_tmp, fluxsort = zip(*sorted(zip(r_tmp, fluxsort)))

    ##################################################################
    # Measure the properties of the flux array according to the
    # method given in 'mode' at each radius
    ##################################################################

    for i in range(rmax):
        n1 = int(math.pi * i**2)
        dn = int(2 * i * math.pi)
        n2 = min(n1 + dn, len(fluxsort))
        
        if mode == 'median':
            radprof[i] = np.median(fluxsort[n1:n1 + dn])
        elif mode == 'mean':
            radprof[i] = np.mean(fluxsort[n1:n1 + dn])
        elif mode == 'std':
            radprof[i] = np.std(fluxsort[n1:n1 + dn])
        elif mode == 'var':
            radprof[i] = np.var(fluxsort[n1:n1 + dn])
        else:
            print "Mode '" + mode + "' not recognized in routine combine.radial_profile"
            sys.exit(1)

    ##################################################################
    # Smooth the radial profile with a Gaussian filter.
    # If the smoothing width is zero, skip this step.
    ##################################################################

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

    radprof2d = radprof[r]
    
    return radprof2d




    
    






