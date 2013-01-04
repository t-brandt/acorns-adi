#!/usr/bin/env python

from easygui import *
import glob
import os
import re
import sys
import pyfits

def pick_dir(msg=None, title="Choose Directory", default=None):
    if msg != None:
        msgbox(msg)
    while 1:
        dirname = diropenbox(title, default=default)
        if dirname == None or not os.path.isdir(dirname):
            if not ccbox(msg="Invalid choice.  Choose again, or exit?",
                         title="Invalid Directory", 
                         choices=["Continue", "Exit"]):
                sys.exit(1)
        else:
            break
    return dirname

def choosefiles(msg="Choose the Raw Files to Use from the Following List.", 
                title="Raw File Selection", data_dir="", filelist=[""]):

    shortlist = [re.sub(data_dir + "/", "", frame) for frame in filelist
                 if os.path.isfile(re.sub("-C[0-9]*", "", frame))]
    nocoadd = [re.sub("-C[0-9]*", "", frame) for frame in filelist 
               if os.path.isfile(frame)]
    n = len(nocoadd)
    if n == 0:
        msgbox("Error:  full files (" + prefix + "*.fits) not found.\n" + 
               "Full files contain header information necessary for\n" + 
               "the ADI data reduction.  Please use a directory with\n" +
               "full files.\n")
        sys.exit(1)

    if 'HICA' in filelist[0]:   
        headerlist = [pyfits.open(frame)[0].header for frame in nocoadd]
        objlist = [header['OBJECT'] for header in headerlist]
        maskid = [header['P_FLID'] for header in headerlist]
        dateid = [header['DATE-OBS'] for header in headerlist]
        exp1time = [header['EXP1TIME'] for header in headerlist]
        exptime = [header['EXPTIME'] for header in headerlist]
        coadd = [header['COADD'] for header in headerlist]
        filtid = [header['FILTER01'] for header in headerlist]
        ndfilt = [header['FILTER02'] for header in headerlist]
        modeid = [header['P_FMID'] for header in headerlist]
        dimen = [''.join([str(header['NAXIS' + str(i + 1)]) + 'x'
                          for i in range(header['NAXIS'])])
                 for header in headerlist]
        
        header = "\n\n\n" + "-" * 250 + "\nFrame ID" + " " * 27 + \
        "Object          Date                   Mode    " + \
        "Mask              Filter    ND Filter    Exp Time      " + \
        "N Coadd    T/Coadd        Dimensions\n" + "-" * 250
        fulllist = [shortlist[i] + '          ' + 
                    objlist[i] + '          ' + 
                    dateid[i] + '         ' + 
                    modeid[i] + '         ' + 
                    maskid[i] + '         ' + 
                    filtid[i] + '         ' + 
                    ndfilt[i] + '         ' + 
                    ('%.2f s        ' % exptime[i]) +
                    ('%d              ' % coadd[i]) + 
                    ('%.2f s        ' % exp1time[i]) + 
                    dimen[i][:-1]
                    for i in range(n)]
    else:
        try:
            headerlist = [pyfits.open(frame)[0].header for frame in nocoadd]
            header1list = [pyfits.open(frame)[-1].header for frame in nocoadd]
            objlist = [header['OBJECT'] for header in headerlist]
            dateid = [header['DATE-OBS'] for header in headerlist]
            exptime = [header['EXPTIME'] for header in headerlist]
            coadd = [header['COADDS'] for header in headerlist]
            filtid = [header['FILTER1'] for header in headerlist]
            ndfilt = [header['FILTER2'] for header in headerlist]
            dimen = [''.join([str(header['NAXIS' + str(i + 1)]) + 'x'
                              for i in range(header['NAXIS'])])
                     for header in header1list]
            
            header = "\n\n\n" + "-" * 250 + "\nFrame ID" + " " * 27 + \
                     "Object            Date                   " + \
                     "Filter 1                  Filter 2                    " + \
                     "Exp Time      N Coadd    Dimensions\n" + "-" * 250
            fulllist = [shortlist[i] + '          ' + 
                        objlist[i] + '          ' + 
                        dateid[i] + '         ' + 
                        filtid[i] + '         ' + 
                        ndfilt[i] + '         ' + 
                        ('%.2f s        ' % exptime[i]) +
                        ('%d              ' % coadd[i]) +
                        dimen[i][:-1]
                        for i in range(n)]
        except:
            try:
                headerlist = [pyfits.open(frame)[0].header for frame in nocoadd]
                header1list = [pyfits.open(frame)[-1].header for frame in nocoadd]
                objlist = [header['OBJECT'] for header in headerlist]
                dateid = [header['DATE-OBS'] for header in headerlist]
                exptime = [header['EXPTIME'] for header in headerlist]
                dimen = [''.join([str(header['NAXIS' + str(i + 1)]) + 'x'
                                  for i in range(header['NAXIS'])])
                         for header in header1list]
                
                header = "\n\n\n" + "-" * 250 + "\nFrame ID" + " " * 27 + \
                         "Object            Date                   " + \
                         "Exp Time      Dimensions\n" + "-" * 250
                fulllist = [shortlist[i] + '          ' + 
                            objlist[i] + '          ' + 
                            dateid[i] + '         ' + 
                            ('%.2f s        ' % exptime[i]) +
                            dimen[i][:-1]
                            for i in range(n)]
            except:
                header1list = [pyfits.open(frame)[-1].header for frame in nocoadd]
                header = "\n\n\n" + "-" * 250 + "\nFrame ID" + " " * 27 + \
                         "Dimensions\n" + "-" * 250
                fulllist = [shortlist[i] + '          ' +
                            dimen[i][:-1]
                            for i in range(n)]
        
    
    while 1:
        usefiles = multchoicebox(msg=msg + header, title=title, 
                                 choices=fulllist)
        if usefiles == None:
            return None
        if choicebox(msg="Use the following " + str(len(usefiles)) + 
                     " files?", title="Confirm", choices=usefiles) != None:
            break
        else:
            if not boolbox(msg="Reselect files or exit?", title="Exit?",
                           choices=["Reselect", "Exit"]):
                sys.exit(1)
    return [data_dir + "/" + frame.split()[0] for frame in usefiles]
            

def pick_file(msg=None, title="Choose File", default=".", filetypes=["*.fits"]):
    if msg != None:
        msgbox(msg)
    while 1:
        filename = fileopenbox(title=title, default=default, 
                               filetypes=filetypes)
        if filename == None or not os.path.isfile(filename):
            if not ccbox(msg="Invalid choice.  Choose again, or exit?",
                         title="Invalid File", 
                         choices=["Continue", "Exit"]):
                sys.exit(1)
        else:
            break
    return filename

class FileSetup(object):

    """

    """
    
    def __init__(self, oldval=None, prefix="HICA"): 

        default = None
        if 'data_dir' in dir(oldval):
            default = oldval.data_dir
        self.data_dir = pick_dir(msg="Select the directory with the raw data.",
                                 title="Raw Data Directory", default=default)
        
        ##################################################################
        # Are there coadd components in data_dir?  If not, don't
        # bother to ask whether to use them.  
        ##################################################################

        files_ok = False
        while not files_ok:
            self.coadd = False
            if len(glob.glob(self.data_dir + "/" + prefix + "*-C[0-9].fits")) > 0:
                self.coadd = ynbox(msg="Will you be using individual coadd " +
                                  "components?\nIf unsure, answer 'No'.")
            if self.coadd:
                # We just tested for the coadded files' existence.
                self.framelist = glob.glob(self.data_dir + 
                                            "/" + prefix + "*-C[0-9].fits")
                for i in range(2, 5):
                    self.framelist += glob.glob(self.data_dir + "/" + prefix +
                                                "*-C" + "[0-9]" * i + ".fits")
                usefiles = choosefiles(data_dir=self.data_dir, 
                                       filelist=self.framelist)
                if usefiles != None:
                    self.framelist = usefiles
                    files_ok = True
            else:
                self.framelist = glob.glob(self.data_dir + "/" + prefix + 
                                           "*[0-9].fits")
                if len(self.framelist) > 0:
                    usefiles = choosefiles(data_dir=self.data_dir, 
                                           filelist=self.framelist)
                    if usefiles != None:
                        self.framelist = usefiles
                        files_ok = True
                                          
                else:
                    if not ccbox(msg="I couldn't find any raw files.",
                              title="Error", choices=["Continue", "Exit"]):
                        sys.exit()
                    
            if not files_ok:
                default = None
                if 'data_dir' in dir(oldval):
                    default = oldval.data_dir
                self.data_dir = pick_dir(msg="Re-select the directory " +
                                         "with the raw data.", 
                                         title="Raw Data Directory",
                                         default=default)
          
        default = self.data_dir
        self.reduce_dir = pick_dir(msg="Next, select a directory to hold " +
                                   "intermediate files.",
                                   title="Intermediate File Directory",
                                   default=default)
        default = self.reduce_dir
        self.output_dir = pick_dir(msg="Finally, select a directory for " +
                                   "the final data products.",
                                   title="Output File Directory",
                                   default=default)

        default = None
        if 'flat' in dir(oldval):
            default = oldval.flat
        self.flat = pick_file(msg="Select the flatfield file.\n" +
                              "As of October 2011, various flats are " + 
                              "available at\n" + 
                              "www.astro.princeton.edu/~tbrandt/seeds/flats/",
                              title="Open Flat", default=default)
        
        default = None
        if 'pixmask' in dir(oldval):
            default = oldval.pixmask
        else:
            default = self.flat
        self.pixmask = None
        while self.pixmask == None:
            if ccbox("Select the hot pixel mask file.\n" +
                     "Optional, but strongly recommended.\n" +
                     "As of October 2011, mask files are available at\n" + 
                     "www.astro.princeton.edu/~tbrandt/seeds/pixmask/",
                     title="Hot Pixel Mask", 
                     choices=["Continue", "Do Not Use"]):
                self.pixmask = pick_file(title="Open Hot Pixel Mask",
                                         default=default)
            else:
                if boolbox(msg="Are you sure you don't want to mask hot pixels?",
                         title="Confirm Decision"):
                    break

        self.scale_phot = ccbox("Would you like to scale the final contrast.\n" +
                     "curves to the central star's flux using unsaturated.\n" +
                     "photometric reference frames?\n",
                     title="Scale Contrast Curves?", 
                     choices=["Yes", "No"])
        if self.scale_phot:
            files_ok = False
            while not files_ok:
                self.phot_dir = pick_dir(msg="Select the directory " +
                                         "with photometric reference frames.",
                                         title="Photometric Reference Frames",
                                         default=self.data_dir)
                self.photlist = glob.glob(self.phot_dir + "/" + prefix +
                                          "*[0-9].fits")
                if len(self.photlist) > 0:
                    usefiles = choosefiles(data_dir=self.phot_dir, 
                                           filelist=self.photlist)
                    if usefiles != None:
                        self.photlist = usefiles
                        files_ok = True
                    else:
                        if not ccbox(msg="I couldn't find any matching files.",
                                     title="Error",
                                     choices=["Retry", "Do Not Scale Contrast"]):
                            self.scale_phot = False
                            files_ok = True
        
        msgbox("Finished setting up files and paths.")
            
    def display(self):
        
        text = "Input files read from " + self.data_dir 
        text += "\nIntermediate files written to " + self.reduce_dir
        text += "\nFinal data products written to " + self.output_dir
        text += "\nUsing flat file " + self.flat
        if self.pixmask == None:
            text += "\nWarning: not masking hot pixels."
        else:
            text += "\nMasking hot pixels flagged in " + self.pixmask
        
        text += "\nFull list of raw files to be used:\n\n"
        for frame in self.framelist:
            text += re.sub(self.data_dir + "/", "     ", frame) + "\n"
        
        textbox(msg="Current File and Directory Setup:", 
                title="Display File Setup",
                text=text)
        

