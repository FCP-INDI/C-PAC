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
        
        self.page = GenericClass(self, "Randomise")
        
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
                      comment="Number of permutations you would like randomise to run the data.")

        self.page.add(label="Threshold ",
                      control=control.INT_CTRL,
                      name='randomise_thresh',
                      type=dtype.NUM,
                      comment="Cluster-based thresholding corrected for multiple comparisons by using the null distribution of the max (across the image) cluster mask.",
                      values=5)

        self.page.add(label="path to design matrix file",
                      control=control.COMBO_BOX,
                      name='randomise_dmat',
                      type=dtype.STR,
                      comment="Full path to a design matrix file to be used when "
                              "running Randomise",
                      values="None")

        self.page.add(label="path to contrast file",
                      control=control.COMBO_BOX,
                      name='randomise_contrast',
                      type=dtype.STR,
                      comment="Full path to the contrast file to be used when "
                              "running Randomise",
                      values="None")

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter