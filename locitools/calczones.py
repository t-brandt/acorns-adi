#!/usr/bin/env python
#
# Original filename: loci.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: July 2011
# 
# Summary:  Calculate the boundaries of the subtraction zones.
# 

import numpy as np

def calczones(rmin, rmax, dr, optreg):

    """
    Function calczones computes the boundaries of the subtraction
    zones for LOCI.  It takes four inputs:
    1.  rmin:  float, minimum radius at which to perform LOCI
    2.  rmax:  float, maximum radius at which to perform LOCI
    3.  dr:  float, the minimum increment in radius
    4.  optreg:  the size of the optimization regions (in pixels)

    Function calczones returns five values: [totzones, r1, r2, theta1, theta2]
    totzones:  integer, the total number of subtraction regions
    r1:  1D array, the minimum radii of the subtraction zones
    r2:  1D array, the maximum radii of the subtraction zones
    theta1:  1D array, the maximum angle of the subtraction zones.
       Angles run from 0 to 2pi.
    theta2:  1D array, the maximum angle of the subtraction zones
    
    """

    totzones = 0
    dr0 = int(dr)
    deltar = dr0
    r = rmin

    ################################################################
    # Tile twice as finely in azimuthal angle as in the original
    # LOCI algorithm.
    # Beginning at a radius of 200 pixels (about 2'' for HiCIAO),
    # increase the radial stepsize to a maximum of dr0 + 10.  
    ################################################################

    theta1 = []
    theta2 = []
    r1 = []
    r2 = []
    
    while r < rmax:
        ntheta = int(2 * np.pi * (r + np.sqrt(optreg) / 2) / 
                     np.sqrt(optreg)) * 2 + 1
        dtheta = 2 * np.pi / ntheta + 1e-10
        if r > 200:
            deltar = dr0 + 20 / np.pi * np.arctan(r / 200. - 1)
        
        for itheta in range(ntheta):
            theta1.append(itheta * dtheta)
            theta2.append((itheta + 1) * dtheta)
            r1.append(r)
            r2.append(r + deltar)

        totzones += ntheta
        r += deltar

    return [totzones, np.asarray(r1), np.asarray(r2),
            np.asarray(theta1), np.asarray(theta2)]

