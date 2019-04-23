import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
import pkg_resources as p

class DualRegression(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/dual_sca.html'))
            
    def get_counter(self):
        return self.counter

class DualRegressionOptions(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Dual Regression Options")
        
        self.page.add(label="Run Dual Regression ", 
                 control=control.CHOICE_BOX, 
                 name='runDualReg', 
                 type=dtype.LSTR, 
                 comment="Run Dual Regression.\n\nRequires that Spatial Regression be enabled under Time Series Extraction.", 
                 values=["Off","On"],
                 wkf_switch = True)
        
        self.page.add(label="Normalize Time Series ", 
                     control=control.CHOICE_BOX, 
                     name='drNorm', 
                     type=dtype.BOOL, 
                     values = ["On", "Off"],
                     comment="Normalize time series before running Dual Regression.")
                

        self.page.set_sizer()
        parent.get_page_list().append(self)
        

    def get_counter(self):
        return self.counter
