#!/usr/bin/env python
#
# Original filename: calc_phot.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: Dec 2012
# 
# Summary:  Calculate the photometric scaling factor to convert the
# sensitivity maps into contrast maps
#

import re
import numpy as np
import pyfits as pyf
import multiprocessing
from scipy import signal, optimize
from scipy import linalg
from parallel import _destripe, _rotate_recenter
from parallel import *
import utils

def errfunc(p, x, y, image, mask):
    x2 = (x - p[1])**2 / p[3]
    y2 = (y - p[0])**2 / p[2]
    xy = (x - p[1]) * (y - p[0]) / p[4]

    #print (image - p[5] * np.exp(-x2 - y2 - xy)) * mask
    return np.reshape((image - p[5] * np.exp(-x2 - y2 - xy)) * mask, -1)

def calc_phot(filesetup, adipar, flat, hotpix, mem, window):
    
    print 'Calibrating Photometry of Sensitivity Map'

    ##############################################################
    # First calibrate everything using the standard routines
    ##############################################################
    
    science_frames = [frame for frame in filesetup.framelist]
    filesetup.framelist = filesetup.photlist
    #ref_flux = utils.read_files(filesetup, ext='')
    ref_flux = _destripe(filesetup, flat, hotpix, mem, adipar,
                         write_files=False, storeall=True,
                         full_destripe=adipar.full_destripe, extraclean=False,
                         do_horiz=adipar.full_destripe)
    filesetup.framelist = science_frames
    
    ##############################################################
    # Now convolve (in parallel) with the same aperture used on
    # the science frames (given as the argument 'window', and take
    # the maximum in each photometric reference frame.
    ##############################################################
    
    tasks = multiprocessing.Queue()
    results = multiprocessing.Queue()
    ncpus = multiprocessing.cpu_count()
    consumers = [ Consumer(tasks, results)
                  for j in range(ncpus) ]
    for w in consumers:
        w.start()

    smoothflux = ref_flux.copy()

    for i in range(ref_flux.shape[0]):
        tasks.put(Task(i, signal.convolve2d, (ref_flux[i], window, 'same')))
    for i in range(ncpus):
        tasks.put(None)

    for i in range(ref_flux.shape[0]):
        index, result = results.get()
        smoothflux[index] = result
    
    y = np.arange(ref_flux.shape[1])
    x = np.arange(ref_flux.shape[2])
    x, y = np.meshgrid(x, y)
    centers = np.zeros((ref_flux.shape[0], 2))

    for i in range(ref_flux.shape[0]):
        yc = y[np.where(smoothflux[i] == smoothflux[i].max())][0]
        xc = x[np.where(smoothflux[i] == smoothflux[i].max())][0]
        imcen = ref_flux[i, yc - 5:yc + 6, xc - 5:xc + 6]
        ycen = y[yc - 5:yc + 6, xc - 5:xc + 6] - yc
        xcen = x[yc - 5:yc + 6, xc - 5:xc + 6] - xc
        mask = np.sqrt(xcen**2 + ycen**2) < 4.01

        p0 = [0., 0., 50., 50., 100., imcen.max()]
        p1, success = optimize.leastsq(errfunc, p0[:],
                                       args=(xcen, ycen, imcen, mask))
        yc += p1[0]
        xc += p1[1]

        centers[i] = [yc, xc]

    #print centers[:]
    imcen = _rotate_recenter(None, ref_flux, storeall=True, centers=centers,
                             newdimen=201, write_files=False)

    ref_psf = np.median(imcen, axis=0)
    
    #outim = pyf.HDUList(pyf.PrimaryHDU(ref_psf))
    #outim.writeto('test_ref.fits', clobber=True)
            
    smoothflux = np.reshape(smoothflux, (smoothflux.shape[0], -1))
    smoothflux = smoothflux.max(axis=1)

    try:
        t_arr = []
        mjd_arr = []
        names = []
        ndfilt = {'ND10': 1.023e-1, 'ND1': 9.138e-3,
                  'ND0.1': 6.904e-4, 'ND0.01': 1.756e-4, 'OPEN': 1.0}

        for frame in filesetup.photlist:
            fullframe = re.sub('-C[0-9]*.fits', '.fits', frame)
            if fullframe == frame:
                t_norm = float(pyf.open(frame)[0].header['EXPTIME'])
            else:
                t_norm = float(pyf.open(fullframe)[0].header['EXP1TIME'])

            t_norm *= ndfilt[pyf.open(fullframe)[0].header['FILTER02']]
            t_arr.append(t_norm)
            mjd_arr.append(float(pyf.open(fullframe)[0].header['MJD']))
            names.append(fullframe)

        t_arr = np.asarray(t_arr)
        mjd_arr = np.asarray(mjd_arr)
        
        ##########################################################
        # Values for ND transmission in HiCIAO in the H band.
        # ***BEWARE*** of trusting the resulting values from
        # another instrument, or from HiCIAO in the K band.
        # Unfortunately, the header keywords are not as
        # standardized as they should be.
        ##########################################################

        n = len(science_frames)
        fullframe = re.sub('-C[0-9]*.fits', '.fits', science_frames[n // 2])

        ##########################################################
        # Are we using HiCIAO coadded frames or not?  If so, we
        # want the exp1time keyword.
        ##########################################################

        try:
            if fullframe == science_frames[n // 2]:
                exptime = float(pyf.open(fullframe)[0].header['EXPTIME'])
            else:
                exptime = float(pyf.open(fullframe)[0].header['EXP1TIME'])
        except:
            exptime = 1
            print "Warning: Unable to read keyword EXPTIME to scale photometry"

        smoothflux /= t_arr / exptime
    except:
        print 'Warning: unable to read keyword EXPTIME to scale photometry'

    #print smoothflux[:]
    #mjd_arr -= mjd_arr[0]
    #mjd_arr *= 86000
    
    #A = np.ones((smoothflux.shape[0], 3))
    #A[:, 1] = mjd_arr
    #A[:, 2] = mjd_arr**2

    #coef1 = linalg.lstsq(A, centers[:, 0])[0]
    #print coef1
    #c1 = np.dot(A, coef1)
    #coef2 = linalg.lstsq(A, centers[:, 1])[0]
    #c2 = np.dot(A, coef2)
    #print coef2

    #for i in range(smoothflux.shape[0]):
    #    print '%7.1f %8.3e %7.2f %7.2f %7.2f %7.2f' % (mjd_arr[i], smoothflux[i], centers[i, 0], centers[i, 1], c1[i], c2[i])
    #    
    #print np.mean(smoothflux), np.median(smoothflux), np.std(smoothflux)
    
    return [np.median(smoothflux), ref_psf]
