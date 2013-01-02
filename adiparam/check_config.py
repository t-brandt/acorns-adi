from easygui import *
from os.path import isfile
import sys
import os
from lociparam import LOCIparam
from adiparam import ADIparam
from filesetup import FileSetup
import pickle

def GetConfig(dirfile='./dirinfo', adifile='./adipar', locifile='./locipar',
              verify=True, prefix='HICA'):
    
    if not isfile(dirfile) and isfile(adifile) and isfile(locifile):
        if not ccbox(msg="The ADI Data Reduction Pipeline is not " + 
                     "fully configured.\n" + "Would you like to start " + 
                     "the configuration wizard?", title="Start Wizard",
                     choices=["Start Wizard", "Exit"]):
            sys.exit(1)
    if isfile(dirfile):
        filesetup = pickle.load(open(dirfile, 'r'))
    else:
        filesetup = FileSetup(prefix=prefix)
        pickle.dump(filesetup, open(dirfile, 'w'))
    if isfile(adifile):
        adipar = pickle.load(open(adifile, 'r'))
    else:
        adipar = ADIparam()
        pickle.dump(adipar, open(adifile, 'w'))
    if isfile(locifile):
        locipar = pickle.load(open(locifile, 'r'))
    else:
        locipar = LOCIparam()
        pickle.dump(locipar, open(locifile, 'w'))

    if not verify:
        return filesetup, adipar, locipar

    choice = ""
    while choice != "Run Program":
        choice = buttonbox(msg="The ADI Data Reduction Pipeline is fully " +
                           "configured.\nPlease select a category to " +
                           "review or change its configuration.\n" +
                           "Click 'Run Program' to begin the data reduction.",
                           title="Confirm Configuration",
                           choices=["Run Program", "File Setup", 
                                    "ADI Parameters", "LOCI Parameters", 
                                    "Exit"])
        if choice == "File Setup":
            filesetup.display()
            if not boolbox(msg="Would you like to keep the current " +
                           "configuration?", title="Keep Configuration?",
                           choices=["Yes", "No, Reconfigure"]):
                filesetup = FileSetup(filesetup, prefix=prefix)
                pickle.dump(filesetup, open(dirfile, 'w'))
        elif choice == "ADI Parameters":
            adipar.display()
            if not boolbox(msg="Would you like to keep the current " +
                           "configuration?", title="Keep Configuration?",
                           choices=["Yes", "No, Reconfigure"]):
                adipar = ADIparam()
                pickle.dump(adipar, open(adifile, 'w'))
            
        elif choice == "LOCI Parameters":
            locipar.display()
            if not boolbox(msg="Would you like to keep the current " +
                           "configuration?", title="Keep Configuration?",
                           choices=["Yes", "No, Reconfigure"]):
                locipar = LOCIparam()
                pickle.dump(locipar, open(locifile, 'w'))
        elif choice == "Exit":
            sys.exit(1)
                
    return filesetup, adipar, locipar

