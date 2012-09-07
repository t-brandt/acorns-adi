# Original filename: dewarp.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: March 2011
# 
# Summary:  Dewarp, recenter, and rotate an image.  
# 

import numpy
import re
import pyfits
import scipy.ndimage
import time

def distortion_interp_flux(flux, y, x):

    flux[:, :] = scipy.ndimage.map_coordinates(flux, [y, x], order=1)

    return flux

def distortion_interp_frame(frame, y, x, storeall=True, output_dir="."):

    frame_ds = re.sub(".fits", "_ds.fits", frame)
    flux = pyfits.open(frame_ds)[0].data
    flux_dw = scipy.ndimage.map_coordinates(flux, [y, x], order=1)
        
    if storeall:
        return flux_dw
    else:
        header = pyfits.open(frame_ds)[0].header
        flux_hdu = pyfits.PrimaryHDU(flux_dw, header)
        fluxout = pyfits.HDUList()
        fluxout.append(flux_hdu)
    
        try:
            outname = re.sub(".fits", "_dw.fits", frame)
            outname = re.sub(".*HICA", output_dir + "/" + "HICA", outname)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fluxout.writeto(outname, clobber=True)
                fluxout.close()
        except IOError, err:
            print err
            sys.exit(1)

    return
    

