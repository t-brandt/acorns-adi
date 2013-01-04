# Original filename: dewarp.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: March 2011
# 
# Summary:  Dewarp, recenter, and rotate an image.  
# 

import re
import pyfits as pyf
import scipy.ndimage
import time
import warnings

def distortion_interp_flux(flux, y, x):

    flux[:, :] = scipy.ndimage.map_coordinates(flux, [y, x], order=1)

    return flux

def distortion_interp_frame(frame, y, x, storeall=True, output_dir="."):

    frame_ds = re.sub(".fits", "_ds.fits", frame)
    flux = pyf.open(frame_ds)[-1].data
    flux_dw = scipy.ndimage.map_coordinates(flux, [y, x], order=1)
        
    header = pyf.open(frame_ds)[0].header
    flux_hdu = pyf.PrimaryHDU(flux_dw, header)
    fluxout = pyf.HDUList()
    fluxout.append(flux_hdu)
    
    try:
        outname = re.sub(".fits", "_dw.fits", frame)
        outname = re.sub(".*/", output_dir + "/", outname)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            fluxout.writeto(outname, clobber=True)
            fluxout.close()
    except IOError, err:
        print err
        sys.exit(1)

    if storeall:
        return flux_dw
    else:
        return None
    

