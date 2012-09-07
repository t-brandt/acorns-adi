from progressbar import ProgressBar
import numpy as np
import scipy as sp
import pyfits as pyf
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

def pca(flux, ncomp=5, nread=0, dosub=True):

    np.seterr(all='ignore')              
    nframes = flux.shape[0]

    pcadir = '/scr/endor1/tbrandt/SEEDS/psfref'
    if nread > 0:
        meanflux = pyf.open(pcadir + '/meanflux.fits')[0].data
        pcashape = meanflux.shape
        if not dosub:
            pca_arr = np.ndarray(
            
        n = (flux.shape[1] - pcashape[0]) // 2
        
        fluxsub = flux[:, n:-n, n:-n]
        noise = np.std(fluxsub, axis=0) + 0.1
        x = np.linspace(0, pcashape[1] - 1, pcashape[1]) - pcashape[1] // 2
        y = np.linspace(0, pcashape[0] - 1, pcashape[0]) - pcashape[0] // 2
        x, y = np.meshgrid(x, y)
        r = np.sqrt(x**2 + y**2)
        
        #np.putmask(noise, r > 20, 1e5)
        
        fluxsub = np.reshape(fluxsub, (nframes, -1))
        noise = np.reshape(noise, -1)
        for i in range(nframes):
            fluxsub[i] /= noise

        refs = np.ndarray((ncomp + 2, pcashape[0] * pcashape[1]))
        
        refs[0] = np.reshape(meanflux, -1) / noise
        refs[1] = np.ones(pcashape[0] * pcashape[1]) / noise

        for i in range(ncomp):
            pcacomp = pyf.open(pcadir + '/pca_comp' + str(i) + '.fits')[0].data
            refs[i + 2] = np.reshape(pcacomp, -1) / noise

        #print np.sum(np.isnan(fluxsub)), np.sum(np.isnan(refs))
        #np.putmask(fluxsub, np.isnan(fluxsub), 0)
        #np.putmask(refs, np.isnan(refs), 0)        

    else:
        oldshape = flux.shape
        
        noise = np.std(flux, axis=0) + 0.1
        dimy, dimx = oldshape[1], oldshape[2]
        noise = np.reshape(noise, -1)

        flux = np.reshape(flux, (nframes, -1))
        
        # First step: subtract the mean
        
        meanflux = np.mean(flux, axis=0)
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

        if not dosub:
            for i in range(nframes):
                flux[i] += meanflux
                flux[i] *= noise
            for i in range(ncomp):
                u[i] *= noise
            return np.reshape(flux, oldshape), \
                   np.reshape(u[:ncomp], (ncomp, oldshape[1], oldshape[2])
    
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
        if readrefs:
            tasks.put(Task(i, pca_sub, (refs, fluxsub[i, :])))
        else:
            tasks.put(Task(i, pca_sub, (refs, flux[i, :])))
    for i in range(ncpus):
        tasks.put(None)

    for i in range(nframes):
        p.render((i + 1) * 100 / nframes,
                 'Processing Frame {0}'.format(i + 1))
        index, result = results.get()

        if readrefs:
            fluxsub[index, :] = result
        else:
            flux[index, :] = result

    if readrefs:
        for i in range(nframes):
            fluxsub[i] *= noise
        fluxsub = np.reshape(fluxsub, (nframes, pcashape[0], pcashape[1]))
        flux[:, n:-n, n:-n] = fluxsub
        return flux
    
    else:
        for i in range(u.shape[0]):
            u[i, :] *= noise[:]

        for i in range(nframes):
            flux[i, :] *= noise[:]
        
        #for i in range(nframes):
        #    flux[i, :] += meanflux

        flux = np.reshape(flux, oldshape)

        return flux, np.reshape(u, oldshape)

