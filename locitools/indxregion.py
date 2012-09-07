#!/usr/bin/env python
#
# Original filename: annulus_indices.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: August 2011
# 
# Summary:  Compuate the indices within each region bounded by
# input values in r and theta.
# 
import numpy as np

def indxregion(rmin, rmax, r1, r2, theta1, theta2, optreg, 
               rindex, r, theta, fwhm=6, nfwhm_ex=0, innerfrac=0.15):

    """
    Function indxregion computes the indices of the subtraction and
    optimization regions within the input radial and azimuthal
    boundaries.  Most of the arguments are passed directly from
    the inputs to annulus_indices
    
    1.  rmin:  float, the minimum allowable radius
    2.  rmax:  float, maximum radius in r (argument 8)
    2.  r1:  1D array, the minimum radius of each subtraction region
    3.  r2:  1D array, the maximum radius of each subtraction region
    4.  theta1:  1D array, the minimum angle of each subtraction
        region.  In the range [0, 2\pi).
    5.  theta2:  1D array, the maximum angle of each subtraction
        region.  In the range [0, 2\pi).
    6.  optreg:  integer, the size (in pixels) of the optimization regions
    7.  rindex:  1D array, the indices of the radius array sorted by
        radius.  Output from sorted_arrays.
    8.  r:  1D array, the radii at rindex.  Output from sorted_arrays.
    9.  theta:  1D array, the angles at rindex.  Output from
        sorted_arrays.

    Optional arguments:
    10. fwhm:  float, PSF FWHM in pixels.  Only used if optimization
        regions exclude an area around the PSF.  Default 6.
    11. nfwhm_ex:  float, buffer region around subtraction zone to be
        ignored when calculating LOCI coefficients
    12. innerfrac:  float, amount that optimization regions extend
        interior to subtraction regions.  Default 15%; i.e., 15% of the
        radial extent of the optimization regions is inward.

    Returns a list of four variables: [optindx, subindx, nopt, nsub]
    iopt, isub are 2D arrays.  First dimension runs over the zone ID,
    while the second dimension runs over the pixel indices.
    nopt, nsub are 1D arrays with the number of elements in each
    optimization and subtraction region, respectively.

    """

    ##################################################################
    # Optimization region will be larger than the subtraction
    # region in both r and theta, and may exclude one fwhm about
    # region.  First, set the outer boundary in r and theta.
    ##################################################################
    
    dtheta_opt = min(np.sqrt(optreg) / r2, np.pi / 2)
    dtheta_opt = max(dtheta_opt, theta2 - theta1)
    dr_opt = np.sqrt(2 * optreg / dtheta_opt + r1**2) - r1

    r_opt1 = max(rmin, r1 - innerfrac * dr_opt)
    r_opt2 = np.sqrt(2 * optreg / dtheta_opt + r_opt1**2)


    ##################################################################
    # Compute the indices of the subtraction regions.
    ##################################################################

    elems_sub = np.all([theta >= theta1, theta < theta2,
                        r >= r1, r < r2], axis=0)       
    subindx = np.extract(elems_sub, rindex)

    thetaopt1 = 0.5 * (theta1 + theta2 - dtheta_opt)
    thetaopt2 = thetaopt1 + dtheta_opt
    
    ##################################################################
    # Optimization regions are larger in theta--need to be careful
    # of the case when adding or subtracting dtheta pushes the value
    # out of the range [0, 2pi]
    ##################################################################

    if thetaopt1 < 0:
        theta_ok = np.any([theta - 2 * np.pi >= thetaopt1,
                           theta < thetaopt2], axis=0)
    elif thetaopt2 > 2 * np.pi:
        theta_ok = np.any([theta >= thetaopt1,
                           theta + 2 * np.pi < thetaopt2],
                          axis=0)
    else:
        theta_ok = np.all([theta >= thetaopt1,
                           theta < thetaopt2], axis=0)
        
    r_ok = np.all([r >= r_opt1, r < r_opt2], axis=0)
        
    ##################################################################
    # Now we need to exclude nfwhm_ex fwhm if nfwhm_ex > 0.
    # Use phi and rho for angle and radius.
    ##################################################################

    if nfwhm_ex > 0:
        dphi = (theta2 - theta1 + nfwhm_ex * 2 * fwhm) / r1
        phi1 = 0.5 * (theta1 + theta2 - dphi)
        phi2 = phi1 + dphi
    
        if phi1 < 0:
            phi_bad = np.any([theta - 2 * np.pi >= phi1,
                              theta < phi2], axis=0)
        elif phi1 > 2 * np.pi:
            phi_bad = np.any([theta >= phi1,
                              theta + 2 * np.pi < phi2], axis=0)
        else:
            phi_bad = np.all([theta >= phi1,
                              theta < phi2], axis=0)

        phi_ok = np.logical_not(phi_bad)
        rho_ok = np.any([r < r1 - nfwhm_ex * fwhm,
                         r > r2 + nfwhm_ex * fwhm], axis=0)
        farenough = np.any([phi_ok, rho_ok], axis=0)
    
    ##################################################################
    # Finally, extract the indices of the optimization regions.
    ##################################################################

    if nfwhm_ex > 0:
         elems_opt = np.all([theta_ok, r_ok, farenough], axis=0)

    else:
        elems_opt = np.all([theta_ok, r_ok], axis=0)
    
    optindx = np.extract(elems_opt, rindex)
    nopt = optindx.size
    nsub = subindx.size

    return [optindx[:], subindx[:], nopt, nsub]
                

