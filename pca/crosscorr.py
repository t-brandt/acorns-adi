from parallel import *
import locitools
import numpy as np
import math
import scipy
import scipy.stats


def get_indices(n, corr_arr):
    nframes = corr_arr.shape[0]
    indxval = np.zeros((nframes, int(n)), np.int)
    for i in range(nframes):
        mincorr = scipy.stats.scoreatpercentile(corr_arr[i],
                                                per=100 - 100. * n / nframes)
        indxval[i] = np.extract(corr_arr[i] >= mincorr,
                                np.arange(nframes))[0:int(n)]
    return indxval
        

def corr(foo, flux):
    #print flux.shape
    nframes = flux.shape[0]
    corr_arr = np.ndarray((nframes, nframes), np.float)

    for i in range(nframes):
        for j in range(i, nframes):
            corr_arr[i, j] = np.sum(flux[i] * flux[j])
            corr_arr[j, i] = corr_arr[i, j]

    return corr_arr


def allcorr(rad, flux, n=None, locipar=None):

    rindex, rsort, thetasort = locitools.sorted_arrays(flux.shape)
    nframes = flux.shape[0]
    if n is None:
        n = nframes
    flux2 = np.ndarray(flux.shape, np.float64)
    flux2[:, :, :] = flux[:, :, :]
    flux2 = np.reshape(flux2, (nframes, -1))

    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for j in range(ncpus) ]
    for w in consumers:
        w.start()

    allrad = range(0, int(rad[-1] + 10))
    
    for i in allrad:
        # Circle in entirely contained within domain
        if i < np.sqrt(flux2.shape[1]) / 2:
            imin = min(math.pi * i**2, flux2.shape[1] - 10)
            imax = min(math.pi * (i + 1)**2, flux2.shape[1])
        # Circle is outside domain
        elif i > np.sqrt(flux2.shape[1]) / np.sqrt(2):
            imin = min(math.pi * i**2, flux2.shape[1] - 10)
            imax = min(math.pi * (i + 1)**2, flux2.shape[1])
        else:
            R = np.sqrt(flux2.shape[1]) / 2
            imin = 4 * R * np.sqrt(i**2 - R**2)
            imin += i**2 * (4 * np.arcsin(R / i) - math.pi)
            imax = 4 * R * np.sqrt((i + 1)**2 - R**2)
            imax += (i + 1)**2 * (4 * np.arcsin(R / (i + 1)) - math.pi)
            imax = min(imax, flux2.shape[1])
        
        #print flux2.shape, type(flux2)
        #print i, len((flux2[:, rindex[imin:imax]]))
        tasks.put(Task(i, corr, (None, flux2[:, rindex[imin:imax]])))
    for i in range(ncpus):
        tasks.put(None)

    corr_arr = np.ndarray((len(allrad), nframes, nframes), np.float)
      
    for i in allrad:        
        index, result = results.get()
        corr_arr[index] = result

    ret_arr = np.ndarray((len(rad), nframes, nframes), np.float)

    for i in range(len(rad)):
        if locipar is not None:
            optreg = math.pi * (locipar.fwhm / 2)**2 * locipar.npsf
            inval = int(np.sqrt(optreg) * locipar.innerfrac * 1.5)
            outval = int(np.sqrt(optreg) * (1 - locipar.innerfrac) * 1.5)
        else:
            inval = 15
            outval = 30

        imin = max(0, int(rad[i] - inval))
        ret_arr[i] = np.sum(corr_arr[imin:int(rad[i] + outval)], axis=0)
        diag = np.sqrt(ret_arr[i, range(nframes), range(nframes)])
        for j in range(nframes):
            ret_arr[i, j, :] /= diag[j]
            ret_arr[i, :, j] /= diag[j]
            
    ret_arr = np.abs(ret_arr)

    indices = np.ndarray((len(rad), nframes, int(n)), np.int)
    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for j in range(ncpus) ]
    for w in consumers:
        w.start()

    for i in range(len(rad)):
        tasks.put(Task(i, get_indices, (n, ret_arr[i])))
    for i in range(ncpus):
        tasks.put(None)

    for i in range(len(rad)):        
        index, result = results.get()
        indices[index] = result

    #print indices.shape
    
    return indices



