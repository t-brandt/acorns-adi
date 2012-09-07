#!/usr/bin/env python
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: June 2011
# 
# Summary:  Parallel wrapper for function centeroflight in module centroid
# 

from progressbar import ProgressBar
import sys
import numpy
from centroid import *
from parallel import *
import multiprocessing

def _speckle_centroid(framelist, flux=None, centerguess=None, storeall=True):

    """
    Function _centeroflight is essentially a parallelized wrapper for the
    function centeroflight.  _centeroflight takes two arguments:
    1.  A list of HiCIAO frames (as filenames)
    2.  An array with all of the destriped HiCIAO data

    centeroflight estimates the centroid of the image by iteratively
    computing the center of light.  
    """

    nframes = len(framelist)
    x = numpy.ndarray(nframes, numpy.float32)
    y = numpy.ndarray(nframes, numpy.float32)

    if centerguess is None:
        centerguess = numpy.ones(nframes, 2) * flux.shape[1] / 2.
    

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
        tasks.put(Task(i, speckle_centroid, (framelist[i], flux[i, :, :], 
                                             centerguess[i, :])))
    for i in range(ncpus):
        tasks.put(None)

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')

    for i in range(nframes):
       p.render((i + 1) * 100 / nframes,
                'Speckle Centroid {0} of {1}'.format(i + 1, nframes))
       index, result = results.get()
       y[index], x[index] = result

    return [y, x]

