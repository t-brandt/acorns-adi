from progressbar import ProgressBar
import numpy as np
import scipy as sp
import time, multiprocessing
from parallel import *
import combine

def absdiff(coef, refs, flux):
    return np.sum(np.abs(flux - sp.dot(coef, refs)))

def pca_sub(refs, flux, usesq=False):
    
    coef = sp.linalg.lstsq(refs.T, flux)[0]
    if usesq:
        flux -= sp.dot(coef, refs)
        return flux
    else:
        coef_new = sp.optimize.fmin_powell(absdiff, coef, (refs, flux),
                                           disp=False)
        flux -= sp.dot(coef_new, refs)
    return flux

def pca(flux, ncomp=30):

    nframes = flux.shape[0]
    oldshape = flux.shape
    
    noise = np.std(flux, axis=0) + 0.1
    dimy, dimx = oldshape[1], oldshape[2]
    #noise[dimy // 2 - 15:dimy // 2 + 16, dimx // 2 - 15:dimx // 2 + 16] += 1e10
    noise = np.reshape(noise, -1)

    flux = np.reshape(flux, (nframes, -1))

    # First step: subtract the mean

    meanflux = np.mean(flux, axis=0)
    #noise = combine.radial_std(np.reshape(meanflux, (oldshape[1], oldshape[2])))
    #noise = np.reshape(noise, -1)

    for i in range(nframes):
        flux[i, :] -= meanflux
        flux[i, :] /= noise
    np.putmask(flux, np.isnan(flux), 0)
    np.putmask(meanflux, np.isnan(meanflux), 0)

    # Next, do a singular value decomposition

    u, s, v = sp.linalg.svd(flux.T, full_matrices=False)

    #for i in range(nframes):
    #    flux[i, :] += meanflux / noise

    #print s[0:50]

    u = u.T
    
    refs = np.ndarray((ncomp + 2, flux.shape[1]))
    refs[0, :] = meanflux[:] / noise[:]
    refs[1, :] = 1. / noise[:]
    refs[2:, :] = u[0:ncomp, :]
    
    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for j in range(ncpus) ]
    for w in consumers:
        w.start()

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    print "\nDoing PCA subtraction on {0} frames with {1} references.".format(nframes, ncomp)

    for i in range(nframes):
        tasks.put(Task(i, pca_sub, (refs, flux[i, :])))
    for i in range(ncpus):
        tasks.put(None)

    for i in range(nframes):
        p.render((i + 1) * 100 / nframes,
                 'Processing Frame {0}'.format(i + 1))
        #coef = sp.linalg.solve(sp.dot(refs, refs.T),
        #                       sp.dot(refs, flux[i, :].T))
        #for j in range(coef.size):
        #    flux[i, :] -= coef[j] * refs[j, :]
        index, result = results.get()
        flux[index, :] = result
    
    for i in range(u.shape[0]):
        u[i, :] *= noise[:]

    for i in range(nframes):
        flux[i, :] *= noise[:]
        
    for i in range(nframes):
        flux[i, :] += meanflux

    flux = np.reshape(flux, oldshape)

    return flux, np.reshape(u, oldshape)

