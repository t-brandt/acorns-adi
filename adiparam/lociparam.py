from easygui import *

def getval(msg="", title="", default="", minval=0, maxval=1):
    val_ok = False
    while not val_ok:
        val = enterbox(msg=msg, title=title, default=default)
        if val == None:
            if ynbox(msg="Really Exit?"):
                sys.exit(1)
        else: 
            try:
                val = float(val)
                if not (val >= minval and val <= maxval):
                    msgbox("That is not a valid value.")
                else:
                    val_ok = True
            except:
                msgbox("That is not a valid value.")
    return float(val)

class LOCIparam(object):

    """

    """
    
    def __init__(self):
        if not ccbox(msg="Now you will select the LOCI parameters.",
                     title="LOCI Parameters Wizard", 
                     choices=["Continue", "Exit"]):
            sys.exit(1)
            
        self.rmax = getval(msg="Enter maximum LOCI radius in pixels.\n" +
                           "Recommended values: 300-500 pixels, or 3-5''.\n" +
                           "Maximum radius must be between 50 and 1000 pixels.",
                           title="Max LOCI Radius", 
                           default=1000, minval=50, maxval=1000)
        self.fwhm = getval(msg="Enter PSF fwhm as measured from an " + 
                           "unsaturated frame.\n" +
                           "If you cannot measure the fwhm, assuming " +
                           "a value of 6 is recommended.\n" +
                            "The fwhm must be between 2 and 12 pixels.",
                           title="PSF FWHM", 
                           default=6, minval=2, maxval=12)

        self.nfwhm = getval(msg="Enter the angular protection zone in units " +
                            "of the PSF fwhm.\n" +
                            "Recommended values for point source detection: " +
                            "between 0.5 and 0.75.\n" + 
                            "Value must be between 0.25 and 10.",
                            title="Angular Protection Zone", 
                            default=0.7, minval=0.25, maxval=10)

        self.npsf = getval(msg="Enter the number of PSF footprints "+ 
                           "(=pi*fwhm^2/4)\n" +
                           "to be used in the optimization regions.\n" +
                           "Recommended values: 150-400.  " +
                           "Value must be between 50 and 10000.",
                           title="Angular Protection Zone",
                           default=200, minval=50, maxval=10000)

        #self.smooth = int(buttonbox(msg="Enter smoothing of reference PSFs " +
        #                            "for calculation of LOCI coefficients.\n" +
        #                            "Subtraction will be done with the " +
        #                            "original images.\n" +
        #                            "1=no smoothing, 3=3x3 median filter, " + 
        #                            "5=5x5 median filter.  Recommended:  1",
        #                            title="Reference Smoothing",
        #                            choices=["1", "3", "5"]))
        self.smooth = 1
        
        self.innerfrac = getval(msg="Enter the fraction of the optimization " +
                                "region interior to the subtraction region.\n" +
                                "Recommended values: 0.05 to 0.2, value " +
                                "must be between 0 and 0.5.\n",
                                title="LOCI Inner Fraction", 
                                default=0.15, minval=0., maxval=0.5)
        
        self.dr0 = getval(msg="Enter the increment in r0, the subtraction " +
                          "radius.  Constrained to be constant for now.\n" + 
                          "Recommended values: 5-25.  " +
                          "Must be between 1 and 50.",
                          title="Subtraction Radius Increment",
                          default=5, minval=1, maxval=50)

        self.feedback = getval(msg="Enter the feedback coefficient to by " +
                               "which to overcorrect in LOCI.\n" +
                               "Recommended, required values: 0-2.\n" +
                               "Setting this to zero is equivalent to " +
                               "normal LOCI.", title="Overcorrection Factor",
                               default=0, minval=0, maxval=3)
        
        self.max_n = getval(msg="Enter the maximum number of comparison " +
                            "frames for each LOCI subtraction.\n" +
                            "Recommended values: 100-200, no more than " +
                            "the number of PSF footprints in\n" +
                            "each optimization region.\n" +
                            "Must be between 50 and twice Npsf.",
                            title="Maximum Number of LOCI Comparison Frames",
                            default=self.npsf, minval=50, maxval=self.npsf)
        
    def display(self):
            
        text = "Performing LOCI out to a radius of " + \
            str(int(self.rmax)) + " pixels."
        text += "\nAssuming a PSF full width at half maximum of " + \
            ('%.2f' % self.fwhm) + " pixels."
        text += "\nUsing an angular protection zone of " + \
            ('%.2f' % self.nfwhm) + " PSF FWHM."
        text += "\nUsing optimization regions of " + \
            str(int(self.npsf)) + " PSF footprints."
        #text += "\nSmoothing reference images with " + \
        #    str(int(self.smooth)) + "x" + str(int(self.smooth)) + \
        #    "median filters\n   to calcuate for LOCI coefficients"
        text += "\nUsing optimization regions extending " + \
            str('%d' % int(self.innerfrac * 100)) + "% inward in radius"
        text += "\nUsing a radial increment of " + \
            str(int(self.dr0)) + " pixels."
        text += "\nOvercorrecting LOCI by a factor " + ('%.2f' % self.feedback)
        text += "\nUsing a maximum of " + str(int(self.max_n)) + \
                "comparison frames for each LOCI subtraction."
            
        textbox(msg="Current LOCI Configuration:",
                title="Display LOCI Configuration",
                text=text)

