import wx
import wx.html
import os
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class SCA(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/sca.html'))         
            
    def get_counter(self):
        return self.counter
    


class SCASettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):

        import os

        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Seed-based Correlation Analysis (SCA) Options")
        
        self.page.add(label="Run Seed-based Correlation Analysis ", 
                     control=control.CHOICE_BOX, 
                     name='runSCA', 
                     type=dtype.LSTR, 
                     comment="For each extracted ROI Average time series, " \
                             "CPAC will generate a whole-brain correlation " \
                             "map.\n\nIt should be noted that for a given " \
                             "seed/ROI, SCA maps for ROI Average time " \
                             "series will be the same.", 
                     values=["Off","On"],
                     wkf_switch = True)

        self.page.add(label = "SCA ROI Paths ",
                      control = control.CHECKBOX_GRID,
                      name = "sca_roi_paths",
                      type = 9,
                      values = '',
                      selections = ["Avg","DualReg","MultReg"],
                      comment="Enter paths to region-of-interest (ROI) " \
                              "NIFTI files (.nii or .nii.gz) to be used for "\
                              "time-series extraction, and then select " \
                              "which types of analyses to run.",
                      size = (450, -1))

        self.page.add(label="Normalize Time Series (Dual Regression) ", 
                     control=control.CHOICE_BOX, 
                     name='mrsNorm', 
                     type=dtype.BOOL, 
                     values = ["On", "Off"],
                     comment="Normalize each time series before running " \
                             "Dual Regression SCA.")

        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        


'''
class MultipleRegressionSCA(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Multiple Regression SCA Options")
        
        self.page.add(label="Run Multiple Regression SCA ", 
                 control=control.CHOICE_BOX, 
                 name='runMultRegSCA', 
                 type=dtype.LSTR, 
                 comment="CPAC will enter all extracted time series from ROI Average TSE, ROI Voxelwise TSE, and Spatial Regression into a single multiple regression model and output a single correlation map.", 
                 values=["Off","On"],
                 wkf_switch = True)
        
        self.page.add(label="Demean Time Series ", 
                     control=control.CHOICE_BOX, 
                     name='mrsDemean', 
                     type=dtype.BOOL, 
                     values = ["On", "Off"],
                     comment="Demean each time series before running Multiple Regression SCA.")
                
        self.page.add(label="Normalize Time Series ", 
                     control=control.CHOICE_BOX, 
                     name='mrsNorm', 
                     type=dtype.BOOL, 
                     values = ["On", "Off"],
                     comment="Normalize each time series before running Multiple Regression SCA.")
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
'''

