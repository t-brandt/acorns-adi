#!/usr/bin/env python
#
# Original filename: makeheader.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: May 2012
#
# Summary: Create a header for the final image files containing the
# reduction parameters
#

import pyfits as pyf
import re

def makeheader(flux, header, filesetup, adipar, locipar):

    newhead = pyf.PrimaryHDU(flux).header

    newhead.update('Version', 1.0, 'Version of ACORNS-ADI')
    newhead.add_comment('*', before='Version')
    newhead.add_comment('********** Basic Information **********', before='Version')
    newhead.add_comment('*', before='Version')
    newhead.update('Version', 1.0, 'Version of ACORNS-ADI')
    newhead.update('Object', header['OBJECT'], 'Target Object Name')
    newhead.update('RA', header['RA'], 'Right Ascension of Pointing (J2000)')
    newhead.update('Dec', header['DEC'], 'Declination of Pointing (J2000)')
    newhead.update('Date', header['DATE'], 'Observation Date')
    newhead.update('NFrames', len(filesetup.framelist), 'Number of Science Frames Used')

    t_tot = 0
    for frame in filesetup.framelist:
        fullframe = re.sub("-C.*fits", ".fits", frame)
        t_tot += pyf.open(fullframe)[0].header['EXPTIME']
        
    newhead.update('T_Int', t_tot / 60, 'Total Integration Time (minutes)')

    
    newhead.update('Destripe', 'Median', 'Destriping Method')
    if adipar.bias_only:
        newhead.update('SciPix', 'False', 'Use Science Pixels for Destriping?')
    else:
        newhead.update('SciPix', 'True', 'Use Science Pixels for Destriping?')
        
    newhead.add_comment('*', before='Destripe')
    newhead.add_comment('********** Destriping Parameters **********', before='Destripe')
    newhead.add_comment('*', before='Destripe')

    newhead.update('Cen_Mask', adipar.r_ex, 'Radius of Central Masked Region (pixels)')
    newhead.update('Reduc', adipar.adi, 'ADI Reduction Method')

    if adipar.adi == 'LOCI':

        newhead.update('FWHM', locipar.fwhm, 'PSF Full Width at Half Maximum (pixels)')
        newhead.add_comment('*', before='FWHM')
        newhead.add_comment('********** LOCI Reduction Parameters **********', before='FWHM')
        newhead.add_comment('*', before='FWHM')
        newhead.update('nFWHM', locipar.nfwhm, 'Angular Protection Zone in Units of PSF FWHM')
        newhead.update('nPSF', locipar.npsf, 'PSF Footprints per Optimization Region')
        newhead.update('dr_in', locipar.innerfrac, 'Fraction of dr_sub Interior to r_sub')
        newhead.update('dr_min', locipar.dr0, 'Minimum Radial Increment')
        newhead.update('smooth', locipar.smooth, 'Size of 2D Median Smoothing Filter for LOCI')
        newhead.update('ngroup', 1 + int((len(filesetup.framelist) - 1) / locipar.max_n), 'Number of Groups of LOCI frames')
        newhead.update('Feedback', locipar.feedback, 'Feedback coefficient for LOCI refinement')


    newhead.update('CRPIX1', flux.shape[1] // 2 + 1, 'Object Centroid in X (indexed from 1)')
    newhead.add_comment('*', before='CRPIX1')
    newhead.add_comment('********** Coordinates **********', before='CRPIX1')
    newhead.add_comment('*', before='CRPIX1')

    newhead.update('CRPIX2', flux.shape[0] // 2 + 1, 'Object Centroid in Y (indexed from 1)')
    newhead.update('CRVAL1', header['CRVAL1'], 'Origin in X')
    newhead.update('CRVAL2', header['CRVAL2'], 'Origin in Y')
    newhead.update('CDELT1', 2.638E-06, 'Pixel scale in X (degrees)')
    newhead.update('CDELT2', 2.638E-06, 'Pixel scale in Y (degrees)')
    
    newhead.update('CTYPE1', 'RA---TAN', 'Pixel coordinate system')
    newhead.update('CTYPE2', 'DEC---TAN', 'Pixel coordinate system')
    newhead.update('CUNIT1', 'DEGREE')
    newhead.update('CUNIT2', 'DEGREE')
    
    newhead.update('CD1_1', -2.638E-06, 'Pixel scale in X (degrees)')
    newhead.update('CD1_2', 0, '')
    newhead.update('CD2_1', 0, '')
    newhead.update('CD2_2', 2.638E-06, 'Pixel scale in Y (degrees)')

    


    
    return newhead

