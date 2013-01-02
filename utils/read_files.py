#!/usr/bin/env python
#
# Original filename: config.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: August 2011
# 
# Summary: Read all files from framelist with a given extension.
# Return a structure with all data.
#

import re
import os.path
import pyfits as pyf
import numpy
from progressbar import ProgressBar
import sys

def read_files(filesetup, ext="", newdimen=None):

    ###################################################################
    # Check for the existence of each file and that it is the same
    # shape as the first file in filesetup.framelist.  Look in both the data
    # directory and in the output directory.  
    ###################################################################

    nframes = len(filesetup.framelist)
    p = ProgressBar('green', width=30, block='=', lastblock='>', empty=' ')

    for i in range(len(filesetup.framelist)):
        frame = filesetup.framelist[i]
        p.render((i + 1) * 100 / nframes, 'Reading File {0}'.format(i + 1))
        
        newfile1 = re.sub(".fits", ext + ".fits", frame)     
        newfile2 = re.sub(".*HICA", filesetup.reduce_dir + "/HICA", newfile1)

        if os.path.isfile(newfile2):
            framedata = pyf.open(newfile2)[-1].data
        elif os.path.isfile(newfile1):
            framedata = pyf.open(newfile1)[-1].data
        else:
            print "Error: failed to read data from " + frame + "."
            sys.exit(1)

    ###################################################################
    # Size array appropriately, place data into newly created array.
    ###################################################################

        refshape = framedata.shape
        if i == 0 and newdimen == None:
            data = numpy.zeros((nframes, refshape[0], refshape[1]),
                               numpy.float32)
        elif i == 0:
            data = numpy.zeros((nframes, newdimen, newdimen), numpy.float32)
        
        if newdimen is None or newdimen == refshape[0]:
            data[i] = framedata
        else:
            y1 = refshape[0] // 2 - newdimen // 2
            x1 = refshape[1] // 2 - newdimen // 2
            if y1 > 0 and x1 > 0:
                data[i] = framedata[y1:y1 + newdimen, x1:x1 + newdimen]
            elif y1 < 0 and x1 < 0:
                y1 *= -1
                x1 *= -1
                data[i, y1:y1 + refshape[0], x1:x1 + refshape[1]] = framedata
            else:
                print "read_files assumes a square array when resizing."
                sys.exit()
                
    return data

    


