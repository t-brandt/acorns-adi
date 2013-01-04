#!/usr/bin/env python
#

from progressbar import ProgressBar
from transform import *
import re
import pyfits as pyf
from parallel import *

def _dewarp(filesetup, memory, flux=None, storeall=True, write_files=False):

    nframes = len(filesetup.framelist)
    print '\nDewarping %d Frames' % nframes

    fullframe = re.sub("-C.*fits", ".fits", filesetup.framelist[0])
    if flux is not None:
        dimen = flux.shape[1]
    else:
        dimen = pyf.open(fullframe)[0].header['NAXIS1']
        write_files = True

    mjd = pyf.open(fullframe)[0].header['MJD']
    
    ######################################################################
    # Calculate the pixel index map for the distortion correction.
    ######################################################################

    y, x = distortion_map(dimen, mjd=mjd)

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
        if write_files:
            frame = re.sub('.*/', filesetup.reduce_dir + '/',
                           filesetup.framelist[i])
            tasks.put(Task(i, distortion_interp_frame, (frame, y, x, storeall,
                                                        filesetup.reduce_dir)))
        else:
            tasks.put(Task(i, distortion_interp_flux, (flux[i], y, x)))
    for i in range(ncpus):
        tasks.put(None)
                       
    ######################################################################
    # Fetch and return all of the dewarped data; first index is frame ID.
    ######################################################################

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    for i in range(nframes):
        p.render((i + 1) * 100 / nframes,
                 'Dewarping Frame {0}'.format(i + 1))
        index, result = results.get()
        if storeall:
            flux[index] = result

    if storeall:
        return flux
    else:
        return None
        

