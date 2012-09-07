#!/usr/bin/env python
#

import sys
import numpy as np
from combine import *
from parallel import *
from progressbar import ProgressBar

def _radialsub(filesetup, flux, mode='median', center=None, rmax=None, 
               smoothwidth=3):

    ######################################################################
    # Print parameters.  Are we asked to recenter, rotate, or both?  
    ######################################################################

    nframes = len(filesetup.framelist)
    print '\nSubtracting a Radial Profile from %d Frames' % nframes
         
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

    profile, rsort, rindex, r_ref = radprof(flux[0], mode=mode)
    flux[0] -= profile
    for i in range(1, nframes):
        tasks.put(Task(i, radprof, (flux[i], mode, rsort, rindex, r_ref,
                                    center, rmax, smoothwidth, None)))
    for i in range(ncpus):
        tasks.put(None)
                       
    ######################################################################
    # Fetch and return all of the transformed data; first index is frame ID.
    ######################################################################

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    
    for i in range(1, nframes):
        p.render((i + 1) * 100 / nframes,
                 'Transforming Frame {0}'.format(i + 1))
        index, result = results.get()
        flux[index] -= result[0]

    #np.putmask(flux, np.logical_not(np.isfinite(flux)), 0)

    return flux
    

