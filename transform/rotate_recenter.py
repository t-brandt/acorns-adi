# Original filename: rotate_recenter.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: March 2011
# 
# Summary:  Recenter and rotate an image.  
# 

import numpy as np
import scipy.ndimage
import pyfits as pyf
import re
import warnings

def rotate_recenter(frame, flux, center=None, theta=0, newdimen=None,
                    writefiles=False, output_dir="."):
    
    """
    Function rotate_recenter takes one argument:
    1.  An array of flux values

    Optional arguments:
    2.  A 2 element list for the new center [y, x], default [dimy//2, dimx//2]
    3.  Angle in radians through which to rotate clockwise, default 0
    4.  Boolian 

    rotate_recenter rotates the image about the center given and
    recenters it in an output array contained within an HDU.
    """ 

    assert len(flux.shape) == 2, "Input array must be two-dimensional."
    dimy, dimx = flux.shape
    if newdimen is None:
        newdimen = max(dimy, dimx)
    if center is None:
        center = [dimy // 2, dimx // 2]
   
    #################################################################
    # Define a larger rectilinear grid for interpolating
    # Offset the center to line up with the input coordinates after
    # the dewarping--this accounts for the linear correction, and
    # is very accurate for sources near the center of the FOV
    #################################################################

    x = np.linspace(0, newdimen - 1., newdimen) + center[1] - newdimen // 2
    y = np.linspace(0, newdimen - 1., newdimen) + center[0] - newdimen // 2
    x, y = np.meshgrid(x, y)

    #################################################################
    # Rotate by theta about the given center.
    #################################################################

    if theta != 0:
        x_new = np.ndarray(x.shape)
        y_new = np.ndarray(y.shape)
        
        x_c = x[newdimen // 2, newdimen // 2]
        y_c = y[newdimen // 2, newdimen // 2]
        x -= x_c
        y -= y_c
        
        x_new[:, :] = x * np.cos(theta) - y * np.sin(theta)
        y_new[:, :] = x * np.sin(theta) + y * np.cos(theta)

        x = x_new + x_c
        y = y_new + y_c
    
    #################################################################
    # Interpolate, write files with an "_r.fits" extension if
    # requested. 
    #################################################################
    
    flux = scipy.ndimage.map_coordinates(flux, [y, x], order=1)

    if writefiles:
        fluxout = pyf.HDUList()
        header = pyf.open(frame)[0].header
        flux_hdu = pyf.PrimaryHDU(flux, header)
        fluxout.append(flux_hdu)
        outname = re.sub(".*/", output_dir + "/", frame)
        outname = re.sub(".fits", "_r.fits", outname)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fluxout.writeto(outname, clobber = True)
                fluxout.close()
        except IOError, err:
            print err
    
    return flux

