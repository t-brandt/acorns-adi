# Original filename: get_pa.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: March 2011
#
#
# Summary: Fetch the position angle from a FITS header.  This routine
# uses the local sidereal time, RA, Dec, and latitude to get the
# parallactic angle of an object.  All quantities except the latitude
# of the telescope (Subaru's by default) are read from the FITS header.
#

import numpy as np
import pyfits as pyf
import math
import ephem
import re
import string

def format_time(t):
    h, m, s = string.split(t, ':')
    h, m, s = [float(h), float(m), float(s)]
    if h < 12:
        h += 24
    return h + m / 60 + s / 3600

def get_pa(frame, lat=19.825556):

    """
    Function get_pa takes one argument:
    1.  The FITS HDU with basic position variables entered

    Optional argument:
    2.  Telescope altitude in degrees (default 19.825556 for Subaru)
    
    get_pa gets the position angle for an alt-az telescope with the
    image rotator off.  Note, this routine currently assumes J2000
    (FK5) coordinates.  It returns a single floating point number,
    the parallactic angle in degrees.  
    """

    #################################################################
    # Fetch the equatorial coordinates in J2000, precess to epoch
    # of observation (keyword DATE-OBS)
    #################################################################

    fullframe = re.sub("-C.*fits", ".fits", frame)
    if fullframe != frame:
        coadd_id = re.sub(".*HICA[0-9]*-C", "", frame)
        coadd_id = re.sub(".fits", "", coadd_id)
        coadd_id = int(coadd_id)
    else:
        coadd_id = 1

    fits = pyf.open(fullframe, "readonly")[0]
    ra = fits.header['RA']
    dec = fits.header['DEC']
    newepoch = fits.header['DATE-OBS']
    oldepoch = str(fits.header['EQUINOX'])
    
    try:    
        delay_hr = 0.5 * (format_time(pyf.open(frame)[0].header['P_HSTEND']) +
                          format_time(pyf.open(frame)[0].header['P_HSTSTR']))
        delay_hr -= format_time(fits.header['P_HSTSTR'])
    except:
        delay_hr = (coadd_id - 0.5) * (fits.header['EXP1TIME'] + 1.2) / 3600

    #delay_hr = coadd_id * fits.header['EXP1TIME'] / 3600
    # Convert from solar to sidereal time
    delay_hr *= 1.0027

    pos = ephem.Equatorial(str(ra), str(dec), epoch = oldepoch)
    newepoch = re.sub("-", "/", newepoch)
    newpos = ephem.Equatorial(pos, epoch = newepoch)

    #################################################################
    # Convert ra and dec from sexagesimal to degrees
    #################################################################    
    
    ralist = str(newpos.ra).split(':')
    declist = str(newpos.dec).split(':')

    ra_new = (float(ralist[0]) + float(ralist[1]) / 60
              + float(ralist[2]) / 3600) * 15
    dec_new = float(declist[0]) + float(declist[1]) / 60 \
              + float(declist[2]) / 3600
    #print ra, dec, newpos.ra, newpos.dec, ra_new, dec_new

    #################################################################
    # Fetch the local sidereal time of the observation, calculate
    # hour angle in degrees
    #################################################################    

    lst_list = fits.header['LST'].split(':')
    lst = float(lst_list[0]) + float(lst_list[1]) / 60 \
           + float(lst_list[2]) / 3600
    ha = (lst + delay_hr) * 15 - ra_new
    
    d2r = math.pi / 180
    r2d = 180 / math.pi

    #################################################################
    # Use the four-parts formula from spherical trigonometry to
    # determine the parallactic angle.  This is copied from
    # parangle.pro by Tim Robishaw.
    #################################################################    
    
    pa = -r2d * np.arctan2(-np.sin(d2r * ha),
                           np.cos(d2r * dec_new) * np.tan(d2r * lat)
                           - np.sin(d2r * dec_new) * np.cos(d2r * ha))
    if dec_new > lat:
        pa = (pa + 360) % 360
    
    return pa + 180

