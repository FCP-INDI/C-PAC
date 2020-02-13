import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class AROMA_ICA(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib.request import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/aroma.html'))            
            
    def get_counter(self):
        return self.counter


class AromaSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "ICA-AROMA De-Noising Options")
        
        self.page.add(label="Run ICA-AROMA ",
                 control=control.CHOICE_BOX, 
                 name='runICA', 
                 type=dtype.LSTR, 
                 comment="Run ICA-AROMA de-noising.",
                 values=["Off", "On"],
                 wkf_switch = True)

        self.page.add(label="De-noising Type",
                 control=control.CHOICE_BOX, 
                 name='aroma_denoise_type', 
                 type=dtype.STR, 
                 values= ["nonaggr", "aggr"],
                 validator = CharValidator("no-alpha"),
                 comment="Types of denoising strategy:i)nonaggr-patial "
                         "component regression,ii)aggr-aggressive denoising")
                 
        #self.page.add(label=".mat file path",
        #			  control=control.COMBO_BOX,
        #			  name='aroma_mat_file',
        #			  type=dtype.STR,
        #			  comment="Specify the path to the matfile describing the affine registration of the functional data to structural data.")
        
        #self.page.add(label="MELODIC dir",
        #			  control= control.DIR_COMBO_BOX,
        #              name ='melodic_dir',
        #			  type=dtype.STR,
        #			  comment="Specify the path to the melodic directory if you do not specify any other functional file")        
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
            
