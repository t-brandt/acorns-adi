#!/usr/bin/env python
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: April 2012
# 
# Summary:  Parallel wrapper for function centroid in module cc_centroid
# 

from progressbar import ProgressBar
import sys
import numpy as np
import pyfits as pyf
from centroid import *
import multiprocessing
from parallel import *
import utils

def _cc_centroid(framelist, flux=None, ref_psf=None, psf_dir='psfref',
                 usemask=True, side=None):

    """
    Function _centroid reads in the reference PSF components and then
    calls the function cc_centroid on each frame.
    _centroid takes two arguments:
    1.  A list of HiCIAO frames (as filenames). Frames should be
        destriped, flat-fielded, and dewarped.
    2.  An array with all of the destriped HiCIAO data (optional)

    cc_centroid estimates the centroid by fitting the PSF components
    at many offsets, making a map of chi2 as a function of offset, and
    centroiding the chi2 map.
    """

    nframes = len(framelist)

    ####################################################################
    # Read in reference data
    ####################################################################

    im = pyf.open(psf_dir + '/pcacomp_0.fits')[0].data
    if ref_psf is None:
        refimage = np.ndarray((3, im.shape[0], im.shape[1]))
    else:
        refimage = np.ndarray((4, im.shape[0], im.shape[1]))
        refimage[-1] = utils.arr_resize(ref_psf, newdim=im.shape[0], padval=0)
               
    refimage[0] = im
    refimage[1] = pyf.open(psf_dir + '/pcacomp_1.fits')[0].data
    refimage[2] = pyf.open(psf_dir + '/pcacomp_2.fits')[0].data
    
    ######################################################################
    # Set up ncpus workers, one for each thread.
    ######################################################################

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
        if flux is not None:
            tasks.put(Task(i, cc_centroid, (refimage, flux[i], framelist[i],
                                            usemask, side)))
        else:
            tasks.put(Task(i, cc_centroid, (refimage, flux, framelist[i],
                                            usemask, side)))
    for i in range(ncpus):
        tasks.put(None)

    ######################################################################
    # Return centroids; first index is frame ID.
    ######################################################################

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    success = np.ndarray(nframes, bool)
    centers = np.ndarray((nframes, 2))
    dr_rms = np.ndarray(nframes)

    for i in range(nframes):        
        p.render((i + 1) * 100 / nframes,
                 'Centroid {0} of {1}'.format(i+1, nframes))
        index, result = results.get()
        if result is not None:
            centers[index, 0], centers[index, 1], dr_rms[index] = result
            success[index] = True
        else:
            centers[index] = [0, 0]
            success[index] = False

    print 'Successfully centroided {0} of {1} frames'.format(np.sum(success),
                                                             len(success))
    
    if np.sum(1 - success) > 0:
        print 'Algorithm failed on the following {0} frames:'.format(np.sum(1 - success))
        for i in range(len(success)):
            if not success[i]:
                print framelist[i]        

    return [success, centers, np.median(dr_rms)]

