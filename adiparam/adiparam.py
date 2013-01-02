from easygui import *
import sys

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

class ADIparam(object):

    """

    """
    
    def __init__(self):
        if not ccbox(msg="Now you will set up the ADI parameters.",
                     title="ADI Parameters Wizard", 
                                  choices=["Continue", "Exit"]):
            sys.exit(1)

        self.full_destripe = boolbox(msg="Perform full intensity calibration" +
                                     "for a Hawaii-IIRG detector,\n" +
                                     "or simply flat-field?\n", 
                                     title="Hawaii-IIRG Corrections?",
                                     choices=["Full Hawaii-IIRG Analysis",
                                              "Flatfield Only"])
        self.bias_only = True
        self.r_ex = 0
        
        if self.full_destripe:
            if boolbox(msg="Use only reference pixels to destripe,\n" +
                       "or use peripheral science pixels as well?\n\n" +
                       "Reference pixels strongly recommended for " +
                       "crowded fields,\n" +
                       "peripheral science pixels recommended otherwise.", 
                       title="Destriping Parameters",
                       choices=["Reference Only", "Reference + Science"]):
                self.bias_only = True
                self.r_ex = 0
            else:
                self.bias_only = False

                self.r_ex = getval(msg="Enter radius to exclude from " +
                                   "the calculation\nof the vertical stripes.\n" +
                                   "Zero = use full array.\n" +
                                   "Recommended: ~300 if weakly saturated,\n" +
                                   "  ~600 if strongly saturated.\n" +
                                   "Minimum 0, maximum 800.",
                                   title="Exclusion Radius",
                                   default=0, minval=0, maxval=800)
                
        self.dewarp = boolbox(msg="Apply the HiCIAO distortion correction?\n" +
                              "If not, no correction will be applied.\n" +
                              "(as of Dec 2012, only the HiCIAO distortion " +
                              "correction is available).\n",
                              title="Apply the HiCIAO Distortion Correction?",
                              choices=["Yes", "No"])
        
        self.do_centroid = boolbox(msg="Will you need to register and " +
                                   "centroid the frames?\n" +
                                   "If not, ACORNS will assume that the " +
                                   "the centroids are the central pixels.\n",
                                   title="Register Frames?",
                                   choices=["Yes", "No"])

        if self.do_centroid:
            msgbox("Next, choose how to centroid the frames.  " + 
                   "This is important, and tricky.  Recommendations:\n" + 
                   "     Cross-correlation if there is no mask,\n" + 
                   "     Moffat profile fit to halo if there is a mask,\n" +
                   "     Center of light if the frames are weakly saturated.\n")

            centroid = buttonbox(msg="How would you like to centroid/" +
                                 "register the frames?",
                                 title="Registration Method",
                                 choices=["Cross-Correlation", 
                                          "Moffat Fit",
                                          "Center of Light"])
            if centroid == "Cross-Correlation":
                self.center = 'crosscorr'
            elif centroid == "Moffat Fit":
                self.center = 'moffat'
            else:
                self.center = 'centeroflight' 
                
        adi = buttonbox(msg="How would you like to process the dataset?\n" +
                        "You may use LOCI, or subtract the median PSF.\n" +
                        "LOCI is far more effective at finding point sources.",
                        title="ADI Method",
                        choices=["LOCI", "median PSF"])
        if adi.startswith("median"):
            self.adi = 'median'
        else:
            self.adi = 'LOCI'


    def display(self):
        
        text = "Intensity Calibration Configuration:  "
        if not self.full_destripe:
            text += "Flatfield only"
        elif self.bias_only:
            text += "Full Hawaii-IIRG analysis using only reference pixels"
        else:
            text += "Full Hawaii-IIRG analysis using peripheral science pixels"

            text += "\n    Not using pixels out to a radius of " + \
                    str(int(self.r_ex))

        if self.dewarp:
            text += "\nApplying the HiCIAO distortion correction."
        else:
            text += "\nNot correcting for field distortion."
        
        text += "\nCentering method:  "
        if not self.do_centroid:
            text += "None Requested"
        elif self.center == 'centeroflight':
            text += "center of light"
        elif self.center == 'moffat':
            text += "Moffat profile fit to halo"
        elif self.center == 'crosscorr':
            text += "cross-correlation to reference frame"
        
        text += "\nADI Processing:  "
        if self.adi == 'median':
            text += "median PSF subtraction"
        elif self.adi == 'LOCI':
            text += "LOCI"

        textbox(msg="Current ADI Configuration:",
                title="Display ADI Configuration",
                text=text)

