import scipy as sp
import numpy as np
from parallel import *
from progressbar import ProgressBar

def doubmatmult(A, b):
    tmp = np.dot(A, b)
    return np.dot(A.T, tmp)

def randmatmult(A, n):
    o = np.random.normal(0., 1., (A.shape[0], n))
    return np.dot(A.T, o)

def gramschmidt(A):
    V = A.T
    # Gram-Schmidt
    # Normalize each column, then subtract its components from other columns
    for i in range(V.shape[0]):
        V[i] /= np.sqrt(np.sum(V[i]**2))
        for j in range(i + 1, V.shape[0]):
            V[j] -= V[i] * np.sum(V[i] * V[j])
    return V.T
    
# Algorithm is presented in Halko, Martinsson, & Tropp, arXiv:0909.4061

def stochastic_svd(A, rank=20, extra_dims=None, power_iter=2,
                   chunksize=50):

    assert len(A.shape) == 2, "Stochastic SVD requires a 2-dimensional array as input"
    nframes, npts = A.shape
    A = np.asarray(A)
    print "\nStarting Stochastic SVD"
    
    if extra_dims is None:
        samples = max(10, 2 * rank) # use more samples than requested factors, to improve accuracy
    else:
        samples = rank + int(extra_dims)    

    # Multiply by a random matrix
    
    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    y = np.zeros((npts, samples))
    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for i in range(ncpus) ]
    for w in consumers:
        w.start()
        
    for i1 in range(0, nframes, chunksize):
        i2 = min(i1 + chunksize, nframes)
        tasks.put(Task(i1, randmatmult, (A[i1:i2], samples)))
            
    for i1 in range(0, nframes, chunksize):
        i2 = min(i1 + chunksize, nframes)
        p.render(i2 * 100 / nframes, 'Multiplying by Random Matrix')
        index, result = results.get()
        y += result

    y = gramschmidt(y)
  
    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    
    for iter in range(power_iter):
        yold = y.copy()
        y[:] = 0.
        
        for i1 in range(0, nframes, chunksize):
            i2 = min(i1 + chunksize, nframes)
            tasks.put(Task(i1, doubmatmult, (A[i1:i2], yold)))    

        for i1 in range(0, nframes, chunksize):
            i2 = min(i1 + chunksize, nframes)
            p.render((i2 + iter * nframes) * 100 / nframes / power_iter,
                     'Performing Power Iteration {0} of {1}'.format(iter + 1, power_iter))
            index, result = results.get()
            y += result
        del yold
        y = gramschmidt(y)

    
    y = y.T

    b = np.zeros((y.shape[0], nframes))
        
    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')
    
    for i1 in range(0, nframes, chunksize):
        i2 = min(i1 + chunksize, nframes)
        tasks.put(Task(i1, np.dot, (y, A[i1:i2].T)))

    for j in range(ncpus):
        tasks.put(None)
    
    for i1 in range(0, nframes, chunksize):
        ii = min(i1 + chunksize, nframes)
        p.render(ii * 100 / nframes, 'Multiplying by Original Matrix')
        index, result = results.get()
        i2 = min(index + chunksize, nframes)
        b[:, index:i2] = result

    print 'Taking SVD of Small Matrix'
    
    u, v, s = sp.linalg.svd(b)
    u = np.dot(y.T, u)

    return u

