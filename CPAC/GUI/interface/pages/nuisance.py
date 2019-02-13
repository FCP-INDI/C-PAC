import wx
import wx.html
from ..utils.generic_class import GenericClass, Control
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
            

class FilteringSettings(wx.ScrolledWindow):
    pass


class Scrubbing(wx.ScrolledWindow):
    pass


class NuisanceRegressionRegressors(Control):
    
    def __init__(self, parent):

        self.ctrl = wx.TextCtrl(parent, id=wx.ID_ANY, 
                                value="")
        

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
                 values=["Off", "On", "On/Off"],
                 wkf_switch = True)
        
        self.page.add(label="Lateral Ventricles Mask (Standard Space) ", 
                     control=control.COMBO_BOX, 
                     name='lateral_ventricles_mask', 
                     type=dtype.STR, 
                     values = os.path.join(fsl, "data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz"),
                     comment="Standard Lateral Ventricles Binary Mask")

        self.page.add(label = "Select Regressors and Censors:",
                      control = NuisanceRegressionRegressors(self),
                      name = "Regressors",
                      type = dtype.LDICT)

        self.page.set_sizer()

        if hasattr(parent, 'get_page_list'):
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


def test_nuisance():

    app = wx.App(False)
    frame = wx.Frame(None, wx.ID_ANY, "Nuisance Regression")
    NuisanceRegression(frame)
    frame.Show(True)
    app.MainLoop()