import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class Nuisance(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/nuisance.html'))         
            
    def get_counter(self):
        return self.counter
            

class NuisanceRegression(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        import os
        
        self.counter = counter

        fsl = os.environ.get('FSLDIR')
        if not fsl:
            fsl = "$FSLDIR"
        
        self.page = GenericClass(self, "Nuisance Signal Regression Options")
        
        self.page.add(label="Run Nuisance Signal Regression ", 
                 control=control.CHOICE_BOX, 
                 name='runNuisance', 
                 type=dtype.LSTR, 
                 comment="Run Nuisance Signal Regression", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        
        self.page.add(label="Lateral Ventricles Mask (Standard Space) ", 
                     control=control.COMBO_BOX, 
                     name='lateral_ventricles_mask', 
                     type=dtype.STR, 
                     values = os.path.join(fsl, "data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz"),
                     comment="Standard Lateral Ventricles Binary Mask")

        self.page.add(label = "Select Regressors:",
                      control = control.LISTBOX_COMBO,
                      name = "Regressors",
                      type = dtype.LDICT,
                      values = ['compcor', 'wm','csf','global','pc1','motion','linear','quadratic', 'gm'],
                      comment = "Select which nuisance signal corrections to apply:\n"\
                                "compcor = CompCor\n"\
                                 "wm = White Matter\n"\
                                 "csf = CSF\n"\
                                 "gm = Gray Matter\n"\
                                 "global = Global Mean Signal\n"\
                                 "pc1 = First Principle Component\n"\
                                 "motion = Motion\n"\
                                 "linear = Linear Trend\n"\
                                 "quadratic = Quadratic Trend",
                     size = (300,120),
                     combo_type =1)
                    
        self.page.add(label= "CompCor Components ",
                      control = control.TEXT_BOX,
                      name = "nComponents",
                      type = dtype.LNUM,
                      values = "5",
                      validator = CharValidator("no-alpha"),
                      comment = "Number of Principle Components to calculate when running CompCor. We recommend 5 or 6.")

        self.page.add(label="Use Friston's 24 (Motion Regression) ",
                      control=control.CHOICE_BOX,
                      name='runFristonModel',
                      type=dtype.LSTR,
                      comment="Use the Friston 24-Parameter Model during volume realignment.\n\nIf this option is turned off, only 6 parameters will be used.\n\nThese parameters will also be output as a spreadsheet.",
                      values=["On", "Off", "On/Off"])


        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        

class MedianAngleCorrection(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Median Angle Correction Options")
        
        self.page.add(label="Run Median Angle Correction ", 
                 control=control.CHOICE_BOX, 
                 name='runMedianAngleCorrection', 
                 type=dtype.LSTR, 
                 comment="Correct for the global signal using Median Angle Correction.", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        
        self.page.add(label= "Target Angle (degrees) ",
                      control = control.TEXT_BOX,
                      name = "targetAngleDeg",
                      type = dtype.LNUM,
                      values = "90",
                      validator = CharValidator("no-alpha"),
                      comment = "Target angle used during Median Angle Correction.")
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter


class FilteringSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Temporal Filtering Options")
        
        self.page.add(label="Run Temporal Filtering ", 
                      control=control.CHOICE_BOX,
                      name='runFrequencyFiltering',
                      type=dtype.LSTR,
                      comment="Apply a temporal band-pass filter to functional data.",
                      values=["Off","On","On/Off"],
                      wkf_switch = True)
        
        self.page.add(label = "Band-Pass Filters ",
                      control = control.LISTBOX_COMBO,
                      name = "nuisanceBandpassFreq",
                      type = dtype.LOFL,
                      values = [0.01, 0.1],
                      comment = "Define one or more band-pass filters by clicking the + button.",
                      size = (200,100),
                      combo_type = 2)

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter


class Scrubbing(wx.ScrolledWindow):

    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter

        self.page = GenericClass(self, "Scrubbing Options")

        self.page.add(label="Run Scrubbing ",
                      control=control.CHOICE_BOX,
                      name='runScrubbing',
                      type=dtype.LSTR,
                      comment="Remove volumes exhibiting excessive motion.",
                      values=["Off", "On", "On/Off"],
                      wkf_switch=True)

        self.page.add(label="Framewise Displacement (FD) Calculation ",
                      control=control.CHOICE_BOX,
                      name='fdCalc',
                      type=dtype.LSTR,
                      comment="Choose which Framewise Displacement (FD) "
                              "calculation to apply the threshold to during "
                              "scrubbing.",
                      values=["Jenkinson", "Power"])

        self.page.add(label="Framewise Displacement (FD) Threshold (mm) ",
                      control=control.TEXT_BOX,
                      name='scrubbingThreshold',
                      type=dtype.LNUM,
                      values="0.2",
                      validator=CharValidator("no-alpha"),
                      comment="Specify the maximum acceptable Framewise "
                              "Displacement (FD) in millimeters.\n\nAny "
                              "volume exhibiting FD greater than this value "
                              "will be removed.",
                      size=(100, -1))

        self.page.add(label="Number of Preceeding Volumes to Remove ",
                      control=control.INT_CTRL,
                      name='numRemovePrecedingFrames',
                      type=dtype.NUM,
                      comment="Number of volumes to remove preceeding a "
                              "volume with excessive FD.",
                      values=1)

        self.page.add(label="Number of Subsequent Volumes to Remove ",
                      control=control.INT_CTRL,
                      name='numRemoveSubsequentFrames',
                      type=dtype.NUM,
                      comment="Number of volumes to remove subsequent to "
                              "a volume with excessive FD.",
                      values=2)

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter
