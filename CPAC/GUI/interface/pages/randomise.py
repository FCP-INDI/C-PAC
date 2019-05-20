import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class Randomise(wx.html.HtmlWindow):
    def __init__(self, parent, counter=0):
        wx.html.HtmlWindow.__init__(self, parent,
                                    style=wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        self.counter = counter
        self.LoadFile(p.resource_filename('CPAC',
                                          'GUI/resources/html/group_analysis.html'))

    def get_counter(self):
        return self.counter

class RandomiseSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, " FSL-Randomise Settings")
        
        self.page.add(label="Run Randomise ", 
                      control=control.CHOICE_BOX, 
                      name='run_randomise',
                      type=dtype.LSTR, 
                      comment="Run Randomise", 
                      values=["Off", "On"],
                      wkf_switch=True)

        self.page.add(label="Permutations",
                      control=control.INT_CTRL,
                      name='randomise_permutation',
                      type=dtype.NUM,
                      values=500,
                      comment="Number of permutations you would like to use "
                              "when building up the null distribution to "
                              "test against.")

        self.page.add(label="Threshold ",
                      control=control.INT_CTRL,
                      name='randomise_thresh',
                      type=dtype.NUM,
                      comment="Cluster-based thresholding corrected for "
                              "multiple comparisons by using the null "
                              "distribution of the max (across the image) "
                              "cluster mask.",
                      values=5)

        self.page.add(label="Demean ", 
                     control=control.CHOICE_BOX, 
                     name='randomise_demean', 
                     type=dtype.BOOL, 
                     values = ["On", "Off"],
                     comment="Demean data temporally before model fitting.")

        self.page.add(label="Threshold-Free Cluster Enhancement ", 
                     control=control.CHOICE_BOX, 
                     name='randomise_tfce', 
                     type=dtype.BOOL, 
                     values = ["On", "Off"],
                     comment="From the FMRIB FSL-Randomise user guide: TFCE "
                             "(Threshold-Free Cluster Enhancement) is a new "
                             "method for finding 'clusters' in your data with"
                             "out having to define clusters in a binary way. "
                             "Cluster-like structures are enhanced but the "
                             "image remains fundamentally voxelwise.")

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
