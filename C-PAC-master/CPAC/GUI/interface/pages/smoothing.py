import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p

class Smoothing(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/smoothing.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/smoothing.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/smoothing.html')
#            else:
#                self.LoadFile('html/smoothing.html')
#        except:
#            self.LoadFile('html/smoothing.html')
            
            
    def get_counter(self):
        return self.counter
            
class SmoothingSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Spatial Smoothing Options")
        
        self.page.add(label= "Kernel FWHM (in mm) ",
                 control=control.TEXT_BOX, 
                 name='fwhm', 
                 type=dtype.LNUM, 
                 values= "4",
                 validator = CharValidator("no-alpha"),
                 comment="Full Width at Half Maximum of the Gaussian kernel used during spatial smoothing.\n\nCan be a single value or multiple values separated by commas.\n\nNote that spatial smoothing is run as the last step in the individual-level analysis pipeline, such that all derivatives are output both smoothed and unsmoothed.")
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter