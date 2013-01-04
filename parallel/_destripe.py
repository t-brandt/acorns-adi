#!/usr/bin/env python
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: June 2011
# 
# Summary:  Parallel wrapper for function destripe in module destripe
# 

from progressbar import ProgressBar
import sys
from destripe import destripe
import multiprocessing
import pyfits as pyf
import numpy as np
import re
from parallel import *

def _destripe(filesetup, flat, hotpix, memory, adipar, 
              write_files=False, clean=True, storeall=True,
              full_destripe=True, do_horiz=True, extraclean=True,
              PDI=False):
    """
    Function _destripe is essentially a parallelized wrapper for the
    function destripe in the destripe module.  _destripe takes five
    required and four optional arguments:
    1.  A list of raw frames (as filenames)
    2.  The file name of the flatfield
    3.  The file name of a hot pixel mask (0 for bad pixels, 1 for good)
    4.  The total system RAM (in bytes)
    5.  The ADI parameters, as output by adiparam

    Optional arguments:
    6.  Write destriped files to new fits files (default false)
    7.  Directory to write destriped files (default '.', unused
               unless (4) is true)
    8.  Use only the reference pixels to estimate the bias.  This is
               useful for crowded fields.  (default false)
    9.  Set masked pixels to the median of their neighbors (default true)
    10. Perform full analysis apppropriate to Hawaii-IIRG detector
               (defualt true)
    11. Calibrate the zero point in different readout channels (default true)
    12. Mask outlier pixels as identified with a large median filter (default
               true)
    13. Use two readout channels, as in HiCIAO's PDI mode (default false)

    destripe estimates the bias at each pixel, subtracts it, and produces
    a flatfielded image.  32 bit integer data is converted to 32 bit
    floating point.
    """

    if not (storeall or write_files):
        print "Error: attempting to run destripe without saving files to either disk or memory"
        sys.exit(1)
        
    print '\nDestriping Frames in Parallel.  Selected Options:'
    if write_files:
        print '    Writing Destriped Files to ' + filesetup.reduce_dir
    else:
        print '    Storing Results in Memory, Not Writing Destriped Files'
    if full_destripe:
        print '    Performing Full Hawaii-IIRG Bias Correction'
        if adipar.bias_only:
            print '    Using Reference Pixels Only to Remove Vertical Stripes'
        else:
            print '    Using Peripheral Science Pixels to Remove Vertical Stripes'
            if adipar.r_ex > 0:
                print '    Not using pixels out to ' + str(int(adipar.r_ex)) + \
                      " in radius."
    else:
        print '    Flat-fielding Only'
    if clean:
        print '    Replacing Masked Pixel Values with Median of Neighbors'
    else:
        print '    Replacing Masked Pixel Values with NaN'

    ######################################################################
    # Create a big array to store all HiCIAO frames.
    # Total required memory ~17MB/frame.
    ######################################################################
    
    nframes = len(filesetup.framelist)
    fits = pyf.open(filesetup.framelist[0], "readonly")
    if storeall:
        flux = np.ndarray((nframes, fits[-1].data.shape[0],
                              fits[-1].data.shape[1]), np.float32)
    flatdata = np.ndarray(flat[-1].data.shape, np.float32)
    flatdata[:] = flat[-1].data

    ######################################################################
    # Mask pixels with low gain to exclude vignetted regions from the fit
    ######################################################################

    if hotpix is not None:
        np.putmask(hotpix[-1].data, flat[-1].data < 0.2, 0)
        np.putmask(hotpix[-1].data, np.isnan(flat[-1].data), 0)
        dimy = fits[-1].data.shape[0]
        dimx = fits[-1].data.shape[1]
        hotpix[-1].data[0:4, :] = 1
        hotpix[-1].data[dimy - 4:dimy, :] = 1
        hotpix[-1].data[:, 0:4] = 1
        hotpix[-1].data[:, dimx - 4:dimx] = 1
    
    fullframe = re.sub("-C.*fits", ".fits", filesetup.framelist[0])
    try:
        x1 = pyf.open(fullframe)[0].header['PRD-MIN1'] - 1
        y1 = pyf.open(fullframe)[0].header['PRD-MIN2'] - 1
    except:
        x1 = y1 = 0

    #hotpix[-1].data[y1 + 96, x1 + 107] = 0
    if hotpix is not None:
        pixmask = np.where(hotpix[-1].data[y1:y1 + dimy, x1:x1 + dimx] == 0)

    ######################################################################
    # Set up ncpus workers, one for each thread.
    ######################################################################

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    
    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for i in range(ncpus) ]
    for w in consumers:
        w.start()

    ######################################################################
    # Submit the jobs with a poison pill for each thread.
    ######################################################################

    for i in range(nframes):
        tasks.put(Task(i, destripe,
                       (filesetup.framelist[i], flat[-1].data, pixmask, 
                        write_files, filesetup.reduce_dir, 
                        adipar.bias_only, clean, storeall, adipar.r_ex,
                        extraclean, full_destripe, do_horiz, PDI)))
    for i in range(ncpus):
         tasks.put(None)

    ######################################################################
    # Jobs are hanging; this is a messy fix.
    ######################################################################
    
    success = np.zeros(nframes, bool)

    for i in range(nframes):
        
        # Anyone alive?
        someone_alive = False
        for w in consumers:
            if w.is_alive():
                someone_alive = True
        if someone_alive == False:
            #print "All processes are dead."
            continue

        p.render((i + 1) * 100 / nframes,
                 'Destriping Frame {0}'.format(i + 1))
        index, result = results.get()
            
        if storeall:
            flux[index] = result
            success[index] = True

    # Did we lose some frames?  Finish them up.

    ibad = np.extract(np.logical_not(success), np.arange(nframes))
    for i in range(len(ibad)):
        p.render((i + 1 + nframes - len(ibad)) * 100 / nframes,
                 'Destriping Frame {0}'.format(i + 1 + nframes - len(ibad)))
        result = destripe(filesetup.framelist[ibad[i]], flat[-1].data, pixmask, 
                          write_files, filesetup.reduce_dir, 
                          adipar.bias_only, clean, storeall, adipar.r_ex,
                          extraclean, full_destripe, do_horiz, PDI)
        if storeall:
            flux[ibad[i]] = result



    ######################################################################
    # Return all of the destriped data; first index is frame ID.
    ######################################################################

    if storeall:
        return flux
    else:
        return None

