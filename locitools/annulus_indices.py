#!/usr/bin/env python
#
# Original filename: annulus_indices.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: August 2011
# 
# Summary:  Compute the indices for each LOCI sub-annular region
#

import numpy as np
import multiprocessing
from parallel import *
from locitools import *

def annulus_indices(rmin, r1, r2, theta1, theta2, optreg, 
                    rindex, rsort, thetasort, 
                    progress, i1, itot, fwhm=6,
                    nfwhm_ex=0, innerfrac=0.15):

    """

    Function annulus_indices computes the indices for each LOCI
    sub-annular region.  It takes at least 12 arguments:
    1.  rmin:  float, the minimum allowable radius
    2.  r1:  1D array, the minimum radius of each subtraction region
    3.  r2:  1D array, the maximum radius of each subtraction region
    4.  theta1:  1D array, the minimum angle of each subtraction
        region.  In the range [0, 2\pi).
    5.  theta2:  1D array, the maximum angle of each subtraction
        region.  In the range [0, 2\pi).
    6.  optreg:  integer, the size (in pixels) of the optimization regions
    7.  rindex:  1D array, the indices of the radius array sorted by
        radius.  Output from sorted_arrays.
    8.  rsort:  1D array, the radii at rindex.  Output from sorted_arrays.
    9.  thetasort:  1D array, the angles at rindex.  Output from
        sorted_arrays.
    10. progress:  class ProgressBar
    11. i1:  integer, tracks progress
    12. itot:  integer, total number of calcluations (to track progress)

    Optional arguments:
    13. fwhm:  float, PSF FWHM in pixels.  Only used if optimization
        regions exclude an area around the PSF.  Default 6.
    14. nfwhm_ex:  float, buffer region around subtraction zone to be
        ignored when calculating LOCI coefficients
    15. innerfrac:  float, amount that optimization regions extend
        interior to subtraction regions.  Default 15%; i.e., 15% of the
        radial extent of the optimization regions is inward.

    Returns a list of four variables: [iopt, isub, nopt, nsub]
    iopt, isub are 2D arrays.  First dimension runs over the zone ID,
    while the second dimension runs over the pixel indices.
    nopt, nsub are 1D arrays with the number of elements in each
    optimization and subtraction region, respectively.
    
    """
    
    ######################################################################
    # Compute maximum subtraction, optimization region areas, then
    # create large arrays iopt and isub to hold indices for all of the
    # regions.  1-D arrays nsub and nopt give the number of elements in
    # each region.
    ######################################################################

    rmax = rsort[-1]    
    totn = r1.size
    maxarea = 0
    for i in range(totn):
        maxarea = max(maxarea, (r2[i]**2  - r1[i]**2) * 0.5 * 
                      (theta2[i] - theta1[i]))
    maxarea += 100 + np.sqrt(maxarea)
        
    iopt = np.ndarray((totn, int(optreg + np.sqrt(optreg) + 100)), np.int)
    isub = np.ndarray((totn, int(maxarea)), np.int)
    
    nsub = np.ndarray(totn, np.int)
    nopt = np.ndarray(totn, np.int)

    ######################################################################
    # Compute the indices for each subtraction and optimization region.
    # The arrays sorted by radius provide initial guesses and greatly
    # reduce the amount of time spent searching.
    # Copy and pass the region (in radius) we need to search.
    ######################################################################
        
    j = np.arange(0, rsort.shape[0], 100)

    ######################################################################
    # Calculate index of minimum radius needed
    ######################################################################

    dtheta_opt1 = min(np.sqrt(optreg) / r2[0], np.pi / 2)
    dtheta_opt1 = max(dtheta_opt1, theta2[0] - theta1[0])
    dr_opt1 = np.sqrt(2 * optreg / dtheta_opt1 + r1[0]**2) - r1[0]
    r_opt1 = max(rmin, r1[0] - innerfrac * dr_opt1)
    n1 = len(np.extract(rsort[j] < r_opt1, j))
    j1 = j[max(n1 - 1, 0)]
    
    ######################################################################
    # Calculate index of maximum radius needed
    ######################################################################

    dtheta_opt2 = min(np.sqrt(optreg) / r2[-1], np.pi / 2)
    dtheta_opt2 = max(dtheta_opt2, theta2[-1] - theta1[-1])
    dr_opt2 = np.sqrt(2 * optreg / dtheta_opt2 + r1[-1]**2) - r1[-1]
    r_opt2_min = max(rmin, r1[-1] - innerfrac * dr_opt2)
    r_opt2 = np.sqrt(2 * optreg / dtheta_opt2 + r_opt2_min**2)
    n2 = len(np.extract(rsort[j] < r_opt2, j))
    j2 = min(j[min(n2 + 1, j.shape[0] - 1)], rsort.shape[0])

    ######################################################################
    # Copy the relevant parts of the arrays for passing
    ######################################################################

    theta = thetasort[j1:j2].copy()
    r = rsort[j1:j2].copy()
    rindx = rindex[j1:j2].copy()

    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for i in range(ncpus) ]
    for w in consumers:
        w.start()

    for i in range(totn):
        tasks.put(Task(i, indxregion, (rmin, rmax, r1[i], r2[i], theta1[i],
                                       theta2[i], optreg, rindx, r, theta, 
                                       fwhm, nfwhm_ex, innerfrac)))
    for i in range(ncpus):
        tasks.put(None)

    for i in range(totn):
        progress.render([(i1 + i + 1) * 100 / itot, i1 * 100 / itot], 
                        ['Calculating Indices for Region {0}'.format(i1 + i + 1),
                         'Processing Region {0} at Radius {1}'.format(i1, int(r1[0]))])
                         
        index, result = results.get()
        iopt_tmp, isub_tmp, nopt[index], nsub[index] = result
        iopt[index, 0:nopt[index]] = iopt_tmp[:]
        isub[index, 0:nsub[index]] = isub_tmp[:]

    return [iopt, isub, nopt, nsub]

