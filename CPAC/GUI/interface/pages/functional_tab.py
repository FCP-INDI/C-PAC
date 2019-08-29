import wx
import wx.html
import os
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class FunctionalPreProcessing(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/func.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/nuisance.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/nuisance.html')
#            else:
#                self.LoadFile('html/functional.html')
#        except:
#            self.LoadFile('html/functional.html')

    def get_counter(self):
        return self.counter


class TimeSeriesOptions(wx.ScrolledWindow):
    
    def __init__(self, parent, counter =0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.page = GenericClass(self, "Time Series Options")
        self.counter = counter 
                
        self.page.add(label="Perform Slice Time Correction ",
                      control=control.CHOICE_BOX,
                      name='slice_timing_correction',
                      type=dtype.LSTR,
                      comment="Interpolate voxel time courses so they are "
                              "sampled at the same time points.",
                      values=["On", "Off", "On/Off"],
                      wkf_switch=True)

        self.page.add(label="TR (in seconds) ",
                      control=control.TEXT_BOX,
                      name='TR',
                      type=dtype.NUM,
                      values="None",
                      validator=CharValidator("no-alpha"),
                      comment="Specify the TR (in seconds) at which images "
                              "were acquired."
                              "\n\nDefault is None- TR information is then "
                              "read from scan parameters in the data "
                              "configuration file, or the image file header "
                              "if there is no scan information in the data "
                              "configuration.\n\nNote: the selection chosen "
                              "here applies to all scans of all "
                              "participants.")

        self.page.add(label="Slice Acquisition Pattern ",
                      control=control.CHOICE_BOX,
                      name='slice_timing_pattern',
                      type=dtype.STR,
                      comment="Acquisition strategy for acquiring image "
                              "slices.\n\nSlice acquisition information is "
                              "read from scan parameters in the data "
                              "configuration file- if this is not provided, "
                              "then this option will apply.\n\nNote: the "
                              "selection here applies to all scans of all "
                              "participants.",
                      values=["Use NIFTI Header", "alt+z", "alt+z2",
                              "alt-z", "alt-z2", "seq+z", "seq-z"],
                      wkf_switch=True)

        self.page.add(label="First Timepoint ",
                      control=control.INT_CTRL, 
                      name='startIdx', 
                      type=dtype.NUM, 
                      comment="First timepoint to include in analysis.\n\n"
                              "Default is 0 (beginning of timeseries).\n\n"
                              "First timepoint selection in the scan "
                              "parameters in the data configuration file, if "
                              "present, will over-ride this selection.\n\n"
                              "Note: the selection here applies to all scans "
                              "of all participants.",
                      values=0)
        
        self.page.add(label="Last Timepoint ",
                      control=control.TEXT_BOX, 
                      name='stopIdx', 
                      type=dtype.NUM, 
                      values="None",
                      validator=CharValidator("no-alpha"),
                      comment="Last timepoint to include in analysis.\n\n"
                              "Default is None or End (end of timeseries).\n"
                              "\nLast timepoint selection in the scan "
                              "parameters in the data configuration file, if "
                              "present, will over-ride this selection.\n\n"
                              "Note: the selection here applies to all scans "
                              "of all participants.")

        self.page.add(label="Volume Registration Two-pass", 
                      control=control.CHOICE_BOX,
                      name='functional_volreg_twopass',
                      type=dtype.BOOL,
                      comment="This options is useful when aligning high-resolution datasets that may need more alignment than a few voxels.",
                      values=["On", "Off"])

        self.page.set_sizer() 
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter


class EPI_DistCorr(wx.ScrolledWindow):
    
    def __init__(self, parent, counter =0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.page = GenericClass(self, "Field Map Distortion Correction "
                                       "Options")
        self.counter = counter 
        fsl = os.environ.get('FSLDIR')
        if fsl == None:
            fsl = "$FSLDIR"
                
        self.page.add(label="Perform Field Map Distortion Correction ",
                      control=control.CHOICE_BOX,
                      name='runEPI_DistCorr',
                      type=dtype.LSTR,
                      comment="Perform field map correction using a single "
                              "phase difference image, a subtraction of the "
                              "two phase images from each echo. Default "
                              "scanner for this method is SIEMENS.",
                      values=["Off", "On", "On/Off"],
                      wkf_switch=True)
                      
        self.page.add(label="Skull-strip the magnitude file with: ",
                      control=control.CHOICE_BOX,
                      name='fmap_distcorr_skullstrip',
                      type=dtype.LSTR,
                      comment="Since the quality of the distortion heavily "
                              "relies on the skull-stripping step, we "
                              "provide a choice of method (AFNI 3dSkullStrip "
                              "or FSL BET).",
                      values=["BET", "3dSkullStrip"],
                      wkf_switch=True)

        self.page.add(label="BET threshold",
                      control=control.TEXT_BOX,
                      name='fmap_distcorr_frac',
                      type=dtype.LNUM,
                      comment="Set the threshold value for the skull-"
                              "stripping of the magnitude file. Depending "
                              "on the data, a tighter extraction may be "
                              "necessary in order to prevent noisy voxels "
                              "from interfering with preparing the field map."
                              "\n\nThe default value is 0.5.",
                      validator=CharValidator("no-alpha"),
                      values="0.5")

        self.page.add(label="AFNI threshold",
                      control=control.TEXT_BOX,
                      name='fmap_distcorr_threshold',
                      type=dtype.NUM,
                      comment="Set the threshold value for the skull-"
                              "stripping of the magnitude file. Depending "
                              "on the data, a tighter extraction may be "
                              "necessary in order to prevent noisy voxels "
                              "from interfering with preparing the field map."
                              "\n\nThe default value is 0.5.",
                      validator=CharValidator("no-alpha"),
                      values="0.6")
                        
        self.page.add(label="DeltaTE, in ms ",
                      control=control.TEXT_BOX,
                      name='fmap_distcorr_deltaTE',
                      type=dtype.LNUM,
                      comment="Set the Delta-TE value, used for preparing "
                              "field map, time delay between the first and "
                              "second echo images. Default value is 2.46 ms.",
                      validator = CharValidator("no-alpha"),
                      values="2.46")
                            
        self.page.add(label="Dwell Time, in s ",
                      control=control.TEXT_BOX,
                      name='fmap_distcorr_dwell_time',
                      type=dtype.LNUM,
                      comment="Set the Dwell Time for the fugue input. "
                              "This is the time between scans, default "
                              "value is 0.0005s.",
                      validator=CharValidator("no-alpha"),
                      values="0.0005")
                                
        self.page.add(label="Dwell to asymmetric ratio ",
                      control=control.TEXT_BOX,
                      name='fmap_distcorr_dwell_asym_ratio',
                      type=dtype.LNUM,
                      comment="Set the asymmetric ratio value for FSL Fugue "
                              "input.",
                      validator=CharValidator("no-alpha"),
                      values="0.93902439")

        self.page.add(label="Phase-encoding direction ",
                      control=control.CHOICE_BOX,
                      name='fmap_distcorr_pedir',
                      type=dtype.STR,
                      comment="Set the phase-encoding direction. The options "
                              "are: x, y, z, -x, -y, -z.",
                      values=["x", "y", "z", "-x", "-y", "-z"])

        self.page.set_sizer() 
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter


class AnatToFuncRegistration(wx.ScrolledWindow):
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "Functional to Anatomical Registration")
        
        fsl = os.environ.get('FSLDIR')
        if fsl == None:
            fsl = "$FSLDIR"

        self.page.add(label="Run Functional to Anatomical Registration ", 
                      control=control.CHOICE_BOX,
                      name='runRegisterFuncToAnat',
                      type=dtype.LSTR,
                      comment="Run Functional to Anatomical Registration",
                      values=["On", "Off"],
                      wkf_switch = True)

        self.page.add(label="Using BB Register ", 
                      control=control.CHOICE_BOX,
                      name='runBBReg',
                      type=dtype.LSTR,
                      comment="Run Functional to Anatomical Registration with BB Register",
                      values=["On", "Off", "On/Off"],
                      wkf_switch = True)
       
        self.page.add(label="Boundary Based Registration Scheduler ", 
                      control=control.COMBO_BOX,
                      name='boundaryBasedRegistrationSchedule',
                      type=dtype.STR,
                      values = str(os.path.join(fsl,"etc/flirtsch/bbr.sch")),
                      comment="Standard FSL 5.0 Scheduler used for Boundary Based Registration.\n\nIt is not necessary to change this path unless you intend to use non-standard MNI registration.")

        self.page.add(label="Use as Functional-to-Anatomical Registration Input ", 
                      control=control.CHOICE_BOX, 
                      name='func_reg_input', 
                      type=dtype.LSTR, 
                      values =["Mean Functional","Selected Functional Volume"],
                      comment="Choose whether to use the mean of the functional/EPI as the input to functional-to-anatomical registration or one of the volumes from the functional 4D timeseries that you choose.")
     
        self.page.add(label="Functional Volume to Use as Input (Selected Functional Volume only) ", 
                      control=control.INT_CTRL, 
                      name='func_reg_input_volume', 
                      type=dtype.NUM, 
                      values = 0,
                      comment="Only for when 'Use as Functional-to-Anatomical Registration Input' is set to 'Selected Functional Volume'. Input the index of which volume from the functional 4D timeseries input file you wish to use as the input for functional-to-anatomical registration.")
 
        self.page.add(label="Functional Masking ",
                      control=control.CHOICE_BOX,
                      name='functionalMasking',
                      type=dtype.LSTR,
                      comment="Choose which tool to be used in functional masking - AFNI 3dAutoMask or FSL BET.",
                      values=["3dAutoMask", "BET", "3dAutoMask & BET"])

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter


class FuncToMNIRegistration(wx.ScrolledWindow):
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "Functional to MNI Registration")
        
        fsl = os.environ.get('FSLDIR')
        if fsl == None:
            fsl = "$FSLDIR"
        
        self.page.add(label="Run Functional to MNI Registration ", 
                      control=control.CHOICE_BOX,
                      name='runRegisterFuncToMNI',
                      type=dtype.LSTR,
                      comment="Register functional images to a standard "
                              "MNI152 template.\n\nThis option must be "
                              "enabled if you wish to calculate any "
                              "derivatives. If set to On [1], only the "
                              "template-space files will be output. If set "
                              "to On/Off [1,0], both template-space and "
                              "native-space files will be output.",
                      values=["On", "On/Off", "Off"],
                      wkf_switch=True)

        self.page.add(label="Functional-to-Template Resolution ", 
                      control=control.CHOICE_BOX,
                      name='resolution_for_func_preproc',
                      type=dtype.STR,
                      values = ["4mm", "3mm", "2mm", "1mm"],
                      comment="The resolution (in mm) to which the " \
                             "preprocessed, registered functional " \
                             "timeseries outputs are written into. Note " \
                             "that selecting a 1 mm or 2 mm resolution " \
                             "might substantially increase your RAM needs- " \
                             "these resolutions should be selected with " \
                             "caution. For most cases, 3 mm or 4 mm " \
                             "resolutions are suggested.")

        self.page.add(label="Functional Derivatives Resolution ",
                      control = control.CHOICE_BOX,
                      name = "resolution_for_func_derivative",
                      type = dtype.STR,
                      values = ["4mm", "3mm", "2mm", "1mm"],
                      comment = "The resolution (in mm) to which the " \
                                "registered derivative outputs are written " \
                                "into.")

        self.page.add(label="Standard Brain only Template (functional resolution) ", 
                      control=control.COMBO_BOX, 
                      name='template_brain_only_for_func', 
                      type=dtype.STR, 
                      values = str(os.path.join(fsl,"data/standard/MNI152_T1_${resolution_for_func_preproc}_brain.nii.gz")),
                      comment="Standard FSL Skull Stripped Template. Used as a reference image for functional registration")
        
        self.page.add(label="Standard Template with Skull (functional resolution) ", 
                      control=control.COMBO_BOX, 
                      name='template_skull_for_func', 
                      type=dtype.STR, 
                      values =  str(os.path.join(fsl,"data/standard/MNI152_T1_${resolution_for_func_preproc}.nii.gz")),
                      comment="Standard FSL Anatomical Brain Image with Skull")
        
        self.page.add(label="Standard Identity Matrix ", 
                      control=control.COMBO_BOX,
                      name='identityMatrix',
                      type=dtype.STR,
                      values = str(os.path.join(fsl,"etc/flirtsch/ident.mat")),
                      comment="Matrix containing all 1's. Used as an identity matrix during registration.\n\nIt is not necessary to change this path unless you intend to use non-standard MNI registration.")
     
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter

