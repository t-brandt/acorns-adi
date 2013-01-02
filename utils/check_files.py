#!/usr/bin/env python
#
# Original filename: config.py
#
# Author: Tim Brandt
# Email: tbrandt@astro.princeton.edu
# Date: August 2011
# 
# Summary: Set configuration parameters to sensible values.
#

import re
import os.path
import pyfits
import numpy

def check_files(filesetup, ext=""):

    ###################################################################
    # Check for the existence of each file.  Look in both the data
    # directory and in the output directory.  
    ###################################################################

    filethere = numpy.ndarray(len(filesetup.framelist), numpy.bool)
    filethere[:] = False
    refshape = pyfits.open(filesetup.framelist[0])[-1].data.shape

    for i in range(len(filesetup.framelist)):
        frame = filesetup.framelist[i]
        
        newfile = re.sub(".fits", ext + ".fits", frame)
        if os.path.isfile(newfile):
            filethere[i] = True

        newfile = re.sub(".*HICA", filesetup.reduce_dir + "/HICA", newfile) 
        if os.path.isfile(newfile):
            filethere[i] = True

    return filethere

    


