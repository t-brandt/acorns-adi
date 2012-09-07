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
from subprocess import *
import multiprocessing
import numpy as np

def config(nframes, framesize):

    ###################################################################
    # Fetch the total amount of physical system memory in bytes.
    # This is the second entry on the second line of the standard
    # output of the 'free' command.
    ###################################################################

    print "\nGetting system parameters, setting pipeline execution parameters..."
    osver = Popen(["uname", "-a"], stdout=PIPE).stdout.read()
    if osver.startswith("Linux"):
        print "You are running Linux."
    elif osver.startswith("Darwin"):
        print "You are running Mac OS-X."
    else:
        print "Your operating system is not recognized."

    if osver.startswith("Linux"):
        mem = Popen(["free", "-b"], stdout=PIPE).stdout.read()
        mem = int(mem.split('\n')[1].split()[1])
    elif osver.startswith("Darwin"):
        mem = Popen(["vm_stat"], stdout=PIPE).stdout.read().split('\n')
        blocksize = re.search('.*size of ([0-9]+) bytes.*', mem[0]).group(1)
        totmem = 0.
        for line in mem:
            if np.any(["Pages free:" in line, "Pages active:" in line,
                       "Pages inactive:" in line, "Pages speculative:" in line,
                       "Pages wired down:" in line]):
                totmem += float(line.split(':')[1]) * float(blocksize)
        mem = int(totmem)
        
    ncpus = multiprocessing.cpu_count()
    hostname = Popen("hostname", stdout=PIPE).stdout.read().split()[0]

    print "\n    You are running on " + hostname + "."
    print "    You have " + str(mem / 2**20) + " megabytes of memory and " + \
          str(ncpus) + " threads available."

    datasize = framesize * nframes * 4
    print "    The dataset consists of " + str(nframes) + " frames, " + \
          str(datasize * 100 / mem) + "% of your physical RAM."

    storeall = False
    
    if datasize * 100 / mem < 60:
        storeall = True
        print "    --> You have enough RAM to store all data."
        print "        The pipeline will not need to write all intermediate files."
    else:
        print "    --> You do not have enough RAM to store all data."
        print "        The pipeline will need to write all intermediate files"
        print "        and do the reduction in pieces."

    return mem, ncpus, storeall



