# Original filename: dewarp.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: March 2011
# 
# Summary:  Dewarp, recenter, and rotate an image.  
# 

import numpy
import scipy
import scipy.ndimage
import scipy.interpolate
import pyfits
import math

def distortion_map(dimen=2048, mjd=55000):

    """
    Function distortion_map takes two arguments:
    1.  A 2048 x 2048 array of flux values
    2.  The MJD of observations, used to decide which distortion map to use

    Optional argument:
    2.  A 2 element list for the new center [y, x], default [1024, 1024]
    3.  Angle in radians through which to rotate clockwise, default 0
    4.  Minimum proportion of data from good pixels
            (if at least this proportion of the interpolated value
            is from bad pixels, set pixel to NaN), default 0.5
    5.  Integer dimension of output array (default 2897)

    dewarp applies the distortion correction and recenters the
    resulting image in a 2897x2897 array.  The coordinates of the
    central point are defined in the input 2048x2048 array.  The
    function dewarp returns an HDU with the dewarped, rotated, and
    recentered array.  
    """
    
    dimy = dimen
    dimx = dimen

    #################################################################
    # Coefficients determined from the MCMC
    #################################################################

    # This is the map prior to May 2011. MJD 55682=May 1st, 2011
    if mjd < 55680:
        a = [-5.914548e-01, -3.835090e-03, -4.091949e-05, 
             1.056099e+02, -2.330969e-02, -2.250246e-03, 
             -1.624182e-02, 2.437204e-04, -3.810423e-03]
        b = [1.022675e+02, -1.490726e-02, -2.589800e-03, 
             5.355218e-01, 2.314851e-03, 5.392667e-04, 
             2.279097e-02, -8.767285e-04, 1.290849e-03]
    # Plate scale differs by 8% with the optical secondary, used
    # in September 2012 (September 10 = MJD 56180)
    elif np.abs(mjd - 56180) < 10:
        a = [-4.793903e-01, -8.112123e-04, -1.529206e-03, 
             9.756364e+01, -4.754364e-03, -2.093333e-03, 
             3.598333e-03, -1.665588e-04, -1.997566e-03]
        b = [9.448591e+01, -2.187263e-03, -2.283869e-03, 
             6.725795e-01, 1.142461e-02, -4.920316e-04, 
             1.253295e-02, -2.193661e-03, -2.881631e-05]
        
    # This is the map after May 2011 except for the case above.
    else:
        a = [-6.394663e-01, 7.587927e-04, -5.303889e-05, 
             1.034162e+02, -5.243642e-03, -2.053183e-03, 
             -1.668922e-03, -5.230824e-05, -2.206181e-03]
        b = [1.001763e+02, 1.363092e-02, -1.966445e-03, 
             6.793916e-01, -4.658817e-04, -1.361540e-04, 
             1.789453e-02, -2.063529e-03, 5.683193e-07]
        
    #################################################################
    # Define a larger rectilinear grid for interpolating
    # Offset the center to line up with the input coordinates after
    # the dewarping--this accounts for the linear correction, and
    # is very accurate for sources near the center of the FOV
    #################################################################

    x = numpy.linspace(0, dimen - 1., dimen) - dimen // 2
    y = numpy.linspace(0, dimen - 1., dimen) - dimen // 2
    x, y = numpy.meshgrid(x, y)
    
    #################################################################
    # Convert units to arcseconds, pixel scale will be 9.5 mas
    #################################################################
    
    x *= 9.5e-3
    y *= 9.5e-3

    #################################################################
    # Correct for the distortion.  I group terms to save a few flops.
    #################################################################
    
    x_new = (a[0] + (a[1] + a[2] * y) * y) * y \
            + (a[3] + (a[4] + a[5] * y) * y) * x \
            + (a[6] + a[7] * y + a[8] * x) * x**2

    y_new = (b[0] + (b[1] + b[2] * y) * y) * y \
            + (b[3] + (b[4] + b[5] * y) * y) * x \
            + (b[6] + b[7] * y + b[8] * x) * x**2

    x_new += dimx // 2
    y_new += dimy // 2

    #################################################################
    # Return the map, which we can use to interpolate each frame.
    #################################################################

    return [y_new, x_new]

