#!/usr/bin/env python
#

import sys
import numpy as np
from transform import *
from parallel import *
from progressbar import ProgressBar

def _rotate_recenter(filesetup, flux, storeall=True, centers=None, 
                     theta=None, newdimen=None, write_files=False):

    ######################################################################
    # Print parameters.  Are we asked to recenter, rotate, or both?  
    ######################################################################

    if filesetup is not None:
        nframes = len(filesetup.framelist)
    else:
        nframes = flux.shape[0]
    rotate = True
    recenter = True
    
    if theta is None:
        rotate = False
        theta = np.zeros(nframes)
    if centers is None:
        recenter = False
        centers = np.ones((nframes, 2))
        centers[:, 0] = flux.shape[1] // 2
        centers[:, 1] = flux.shape[2] // 2
    if newdimen is None:
        #newdimen = max(flux.shape[1], flux.shape[2])
        newdimen = flux.shape[1]

    if not (rotate or recenter):
        print "Called _rotate_recenter with no new centers or rotation angles."
        return flux
    elif rotate and recenter:
        print '\nRotating and Recentering %d Frames' % nframes
    elif rotate:
        print '\nRotating %d Frames' % nframes
    else:
        print '\nRecentering %d Frames' % nframes
         
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
        if filesetup is not None:
            tasks.put(Task(i, rotate_recenter, (filesetup.framelist[i],
                                                flux[i], centers[i], theta[i],
                                                newdimen, write_files,
                                                filesetup.reduce_dir)))
        else:
            tasks.put(Task(i, rotate_recenter, ('', flux[i], centers[i],
                                                theta[i], newdimen,
                                                write_files, '')))
            
    for i in range(ncpus):
        tasks.put(None)
                       
    ######################################################################
    # Fetch and return all of the transformed data; first index is frame ID.
    ######################################################################

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    
    if newdimen != flux.shape[1]:
        flux2 = np.zeros((flux.shape[0], newdimen, newdimen), np.float32)

    for i in range(nframes):
        p.render((i + 1) * 100 / nframes,
                 'Transforming Frame {0}'.format(i + 1))
        index, result = results.get()
        if newdimen != flux.shape[1]:
            flux2[index] = result
        else:
            flux[index] = result

    if newdimen != flux.shape[1]:
        flux = flux2

    #np.putmask(flux, np.logical_not(np.isfinite(flux)), 0)

    return flux
    

