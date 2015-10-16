import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p

class Filtering(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/temporal.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/temporal.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/temporal.html')
#            else:
#                self.LoadFile('html/temporal.html')
#        except:
#            self.LoadFile('html/temporal.html')
            
            
    def get_counter(self):
        return self.counter
            

