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
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/functional.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/motion.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/motion.html')
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
                
                
        self.page.add(label= "First Timepoint ",
                 control=control.INT_CTRL, 
                 name='startIdx', 
                 type=dtype.NUM, 
                 comment="First timepoint to include in analysis.\n\nDefault is 0 (beginning of timeseries).", 
                 values=0)
        
        self.page.add(label= "Last Timepoint ",
                 control=control.TEXT_BOX, 
                 name='stopIdx', 
                 type=dtype.NUM, 
                 values= "End",
                 validator = CharValidator("no-alpha"),
                 comment="Last timepoint to include in analysis.\n\nDefault is None or End (end of timeseries).")
        
        self.page.add(label= "TR ",
                 control=control.TEXT_BOX, 
                 name='TR', 
                 type=dtype.NUM, 
                 values= "None",
                 validator = CharValidator("no-alpha"),
                 comment="Specify the TR at which images were acquired.\n\nDefault is None (TR information is read from image file header)")
        
        
    
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
        
        
        self.page.add(label="Run Functional to Anatomical Registration:", 
                     control=control.CHOICE_BOX, 
                     name='runRegisterFuncToAnat', 
                     type=dtype.LSTR, 
                     comment="Run Functional to Anatomical Registration", 
                     values=["On","Off"],
                     wkf_switch = True)

        self.page.add(label="Using BB Register:", 
                     control=control.CHOICE_BOX,
                     name='runBBReg', 
                     type=dtype.LSTR, 
                     comment="Run Functional to Anatomical Registration with BB Register", 
                     values=["On","Off"],
                     wkf_switch = True)
        
        self.page.add(label="Functional Standard Resolution:", 
                     control=control.CHOICE_BOX, 
                     name='standardResolution', 
                     type=dtype.STR, 
                     values = ["3mm", "2mm", "1mm"],
                     comment="The resolution (in mm) to which functional images are transformed during registration")
        
        self.page.add(label="Standard Brain only Template (functional resolution):", 
                      control=control.COMBO_BOX, 
                      name='standardResolutionBrain', 
                      type=dtype.STR, 
                      values = str(os.path.join(fsl,"data/standard/MNI152_T1_${standardResolution}_brain.nii.gz")),
                      comment="Standard FSL Skull Stripped Template. Used as a reference image for functional registration")
        
        self.page.add(label="Standard Template with Skull (functional resolution):", 
                      control=control.COMBO_BOX, 
                      name='standard', 
                      type=dtype.STR, 
                      values =  str(os.path.join(fsl,"data/standard/MNI152_T1_$standardResolution.nii.gz")),
                      comment="Standard FSL Anatomical Brain Image with Skull")
        
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
                     comment="Register functional images to a standard MNI152 template.\n\nThis option must be enabled if you wish to calculate any derivatives.", 
                     values=["On","Off"],
                     wkf_switch = True)
        
        self.page.add(label="Standard Identity Matrix ", 
                     control=control.COMBO_BOX, 
                     name='identityMatrix', 
                     type=dtype.STR, 
                     values = str(os.path.join(fsl,"etc/flirtsch/ident.mat")),
                    comment="Matrix containing all 1's. Used as an identity matrix during registration.\n\nIt is not necessary to change this path unless you intend to use non-standard MNI registration.")
                    
        self.page.add(label="Boundary Based Registration Scheduler ", 
                     control=control.COMBO_BOX, 
                     name='boundaryBasedRegistrationSchedule', 
                     type=dtype.STR, 
                     values = str(os.path.join(fsl,"etc/flirtsch/bbr.sch")),
                     comment="Standard FSL 5.0 Scheduler used for Boundary Based Registration.\n\nIt is not necessary to change this path unless you intend to use non-standard MNI registration.")
     
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
