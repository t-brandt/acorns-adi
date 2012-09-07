#!/usr/bin/env python
#
# Original filename: finecenter.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: May 2012
# 
# Summary:  Interactively refine the centroid of an image sequence
# 

from speckle_centroid import speckle_centroid
from easygui import *
from pylab import *
import pyfits as pyf


def finecenter(flux, objname, output_dir):

    """

    function finecenter interactively refines the centroid of an image
    sequence.  The routine averages a sequence of frames, plots the
    central 80x80 pixels (0.8x0.8 arcseconds for HiCIAO) and labels
    the center.  The user then inputs the offset, and, when satisfied,
    clicks 'OK'.  The final image is saved in the output directory.
    The function takes three inputs:
    
    1.  A 3D flux array; the first index should run over frame number
    2.  The object name, for naming the output file.  This should be
        a string.
    3.  The output directory.  This should be a string, and the directory
        should exist and have write permissions.

    The function returns the user-determined offset in the centroid as
    [y_offset, x_offset].

    """

    ################################################################
    # Extract, plot the central portion of the flux array.  Try to
    # centroid between pairs of speckles (local maxima in intensity)
    # to form a guess as to the absolute centroid.  This doesn't
    # always improve things; the user will soon check interactively.
    ################################################################

    dimy, dimx = flux.shape
    y, x = speckle_centroid('', flux, center=[dimy // 2, dimx // 2])

    dimy, dimx = flux.shape
    di = min(40, dimy // 2 - 1, dimx // 2 - 1)
    subarr = flux[dimy // 2 - di:dimy // 2 + di + 1,
                  dimx // 2 - di:dimx // 2 + di + 1]
    
    grid = np.arange(di * 2 + 1)
    grid_y, grid_x = np.meshgrid(grid, grid)

    yc, xc = [0., 0.]
    
    while 1:
        
        ################################################################
        # Plot the central part of the image with bullseye-type
        # annotations until the user indicates he/she is satisfied.
        ################################################################

        r = np.sqrt((grid_y - di - y + 100 - yc)**2 +
                    (grid_x - di - x + 100 - xc)**2)
        figure(figsize=(8,8))
        imshow(np.sqrt(subarr), interpolation='bilinear', origin='lower')
        contour(grid_x, grid_y, r, [4, 4, 15, 25, 35],
                linewidths=(4, 1, 3, 3, 3,),
                colors=('k', 'm', 'k', 'k', 'k'),
                linestyles=('solid', 'solid', 'dashed', 'dashed', 'dashed'))
        plot([di + xc + x - 100], [di + y + yc - 100],
             color='y', marker='+', markersize=8, mew=2)
        axis('off')
        savefig(output_dir + '/' + objname + '_center_verify.png')
        show(block=False)
        clf()
        
        shift = enterbox(msg='Check/refine the pipeline absolute center.\n' +
                         'Enter an offset to apply in the format dx, dy.\n' +
                         'Graphic has dimensions ' + str(subarr.shape[0]) +
                         ' x ' + str(subarr.shape[1]) + '.\n' +
                         'Click OK to keep the current center.',
                         title='Fine Centroiding',
                         default=str(xc) + ',' + str(yc))
        try:
            yc_tmp = float(shift.split(',')[1])
            xc_tmp = float(shift.split(',')[0])
            if yc_tmp == yc and xc_tmp == xc:
                if ynbox(msg='Keep the center shown?',
                         title='Fine Centroiding'): 
                    break
            else:
                yc, xc = [yc_tmp, xc_tmp]
        except:
            msgbox(msg='Invalid format, please try again.',
                   title='Fine Centroiding') 
            
    ################################################################
    # Return the user-determined offset.
    ################################################################
            
    return [yc + y - 100, xc + x - 100]
            

