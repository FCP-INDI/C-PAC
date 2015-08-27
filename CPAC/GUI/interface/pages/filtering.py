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
