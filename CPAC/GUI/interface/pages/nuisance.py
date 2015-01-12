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
            
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/nuisance.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/nuisance.html')
#            else:
#                self.LoadFile('html/nuisance.html')
#        except:
#            self.LoadFile('html/nuisance.html')
            
            
    def get_counter(self):
        return self.counter
            
class NuisanceCorrection(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        import os
        
        self.counter = counter

        fsl = os.environ.get('FSLDIR')
        if not fsl:
            fsl = "$FSLDIR"
        
        self.page = GenericClass(self, "Nuisance Signal Correction Options")
        
        self.page.add(label="Run Nuisance Signal Correction ", 
                 control=control.CHOICE_BOX, 
                 name='runNuisance', 
                 type=dtype.LSTR, 
                 comment="Run Nuisance Signal Correction", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        
        self.page.add(label="Lateral Ventricles Mask (Standard Space) ", 
                     control=control.COMBO_BOX, 
                     name='lateral_ventricles_mask', 
                     type=dtype.STR, 
                     values = os.path.join(fsl, "data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz"),
                     comment="Standard Lateral Ventricles Binary Mask")

        self.page.add(label = "Corrections:",
                      #control = control.CHECKLISTBOX_COMBO,
                      control = control.LISTBOX_COMBO,
                      name = "Corrections",
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
