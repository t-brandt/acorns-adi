from progressbar import ProgressBar
import numpy as np
import scipy as sp
import pyfits as pyf
import time, multiprocessing
from parallel import *
import combine
from svd import stochastic_svd
import os.path

def absdiff(coef, refs, flux):
    return np.sum(np.abs(flux - sp.dot(coef, refs)))

def pca_sub(refs, flux, usesq=True):
    
    coef = sp.linalg.lstsq(refs.T, flux)[0]
    if usesq:
        flux -= sp.dot(coef, refs)
        return flux
    else:
        coef_new = sp.optimize.fmin_powell(absdiff, coef, (refs, flux),
                                           disp=False)
        flux -= sp.dot(coef_new, refs)
            
    return flux

def pca(flux, ncomp=5, nread=0, dosub=True, pcadir='.'):

    np.seterr(all='ignore')              
    nframes = flux.shape[0]
    
    oldshape = flux.shape
    
    noise = np.std(flux, axis=0) + 0.1
    dimy, dimx = oldshape[1], oldshape[2]
    x, y = np.meshgrid(np.arange(dimy) - dimy // 2, np.arange(dimx) - dimx // 2)
    #noise += 1e5 * (x**2 + y**2 < 15**2)
    noise = np.reshape(noise, -1)
    
    flux = np.reshape(flux, (nframes, -1))
    
    # First step: subtract the mean
    
    meanflux = np.mean(flux, axis=0)
    #np.putmask(noise, meanflux > 3e4, 3e4)
    for i in range(nframes):
        flux[i] -= meanflux
        flux[i] /= noise
        np.putmask(flux[i], np.isnan(flux[i]), 0)
    np.putmask(meanflux, np.isnan(meanflux), 0)

    # Next, do a singular value decomposition
    #u, s, v = sp.linalg.svd(flux.T, full_matrices=False)
    #u = u.T[:ncomp]
    if ncomp > 0:
        u = stochastic_svd(flux, rank=ncomp, power_iter=3).T[:ncomp]

    if not dosub:
        for i in range(nframes):
            flux[i] *= noise
            flux[i] += meanflux
        if ncomp > 0:
            for i in range(ncomp):
                u[i] *= noise
            
    pca_arr = np.zeros((ncomp + nread + 2, oldshape[1], oldshape[2]),
                         np.float32)
    if ncomp > 0:
        pca_arr[:ncomp] = np.reshape(u, (ncomp, oldshape[1], oldshape[2]))
        
    pca_arr[ncomp] = 1 
    pca_arr[ncomp + 1] = np.reshape(meanflux, (oldshape[1], oldshape[2]))

    if nread > 0:
        meanflux2 = pyf.open(pcadir + '/meanflux.fits')[0].data
        pcashape = meanflux2.shape
        if pca_arr.shape[1] > pcashape[0]:
            n = (pca_arr.shape[1] - pcashape[0]) // 2
            pca_arr[ncomp + 1, n:-n, n:-n] = meanflux2
        else:
            n = (pcashape[0] - pca_arr.shape[1]) // 2
            pca_arr[ncomp + 1] = meanflux2[n:-n, n:-n]
            

    for i in range(nread):
        if pca_arr.shape[1] > pcashape[0]:
            pca_arr[ncomp + i + 2, n:-n, n:-n] = pyf.open(pcadir + '/pca_comp_alt' +
                                                          str(i) + '.fits')[0].data
        else:
            pca_arr[ncomp + i + 2] = pyf.open(pcadir + '/pca_comp_alt' +
                                              str(i) + '.fits')[0].data[n:-n, n:-n]

    if not dosub:
        return np.reshape(flux, oldshape), pca_arr

    pca_arr = np.reshape(pca_arr, (pca_arr.shape[0], -1))
    for i in range(ncomp, pca_arr.shape[0]):
        pca_arr[i] /= noise
        
    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for j in range(ncpus) ]
    for w in consumers:
        w.start()

    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    print "\nDoing PCA subtraction on {0} frames with {1} references.".format(nframes, ncomp + nread + 2)

    for i in range(nframes):
        tasks.put(Task(i, pca_sub, (pca_arr, flux[i])))
    for i in range(ncpus):
        tasks.put(None)

    for i in range(nframes):
        p.render((i + 1) * 100 / nframes,
                 'Processing Frame {0}'.format(i + 1))
        index, result = results.get()

        flux[index] = result

    for i in range(pca_arr.shape[0]):
        pca_arr[i] *= noise

    for i in range(nframes):
        flux[i] *= noise
        
    flux = np.reshape(flux, oldshape)
    
    return flux, np.reshape(pca_arr, (pca_arr.shape[0], oldshape[1], oldshape[2]))

