import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import os
import pkg_resources as p
    
class MDMRSettings(wx.ScrolledWindow):
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "MDMR")

        self.page.add(label="Run MDMR?",
                      control=control.CHOICE_BOX,
                      name="mdmr",
                      type=dtype.LSTR,
                      comment="",
                      values=["Off", "On"])

        self.page.add(label="Mask ROI", 
                     control=control.COMBO_BOX,
                     name='mdmr_roi', 
                     type=dtype.STR, 
                     values=str(""),
                     comment="")
               
        self.page.add(label="Regressors file", 
                     control=control.COMBO_BOX,
                     name="mdmr_regressors", 
                     type=dtype.STR, 
                     values="",
                     comment="")

        self.page.add(label="Regressor Columns", 
                     control=control.TEXT_BOX, 
                     name="mdmr_regressors_cols", 
                     type=dtype.STR, 
                     values="",
                     comment="")

        self.page.add(label="Parallel nodes",
                      control=control.INT_CTRL,
                      name="mdmr_parallel",
                      type=dtype.NUM,
                      comment="",
                      values=1)

        self.page.add(label="Iterations",
                      control=control.INT_CTRL,
                      name="mdmr_iterations",
                      type=dtype.NUM,
                      comment="",
                      values=1)

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
                
