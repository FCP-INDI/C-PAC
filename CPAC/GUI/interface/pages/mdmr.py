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
                
        self.page = GenericClass(self, "Multivariate Distance Matrix Regression")

        self.page.add(label="Run MDMR?",
                      control=control.CHOICE_BOX,
                      name="mdmr",
                      type=dtype.LSTR,
                      comment="Run Multivariate Distance Matrix Regression.",
                      values=["Off", "On"])

        self.page.add(label="Mask ROI", 
                     control=control.COMBO_BOX,
                     name='mdmr_roi_file', 
                     type=dtype.STR, 
                     values=str(""),
                     comment="Path to a mask file in NIFTI format.")
               
        self.page.add(label="Regressors file", 
                     control=control.COMBO_BOX,
                     name="mdmr_regressor_file", 
                     type=dtype.STR, 
                     values="",
                     comment="CSV file containing the variables to regress. "
                             "It may have one line per subject.")
        
        self.page.add(label="Regressors Participant Column Name",
                      control=control.TEXT_BOX,
                      name="mdmr_regressor_participant_column",
                      type=dtype.STR,
                      comment="Name of the participants column in your Regressors "
                              "file.",
                      values="")

        self.page.add(label="Regressors columns", 
                     control=control.TEXT_BOX, 
                     name="mdmr_regressor_columns", 
                     type=dtype.STR,
                     values="",
                     comment="Columns from the CSV file indicating factor "
                             "variables. Other columns will be handled as covariates. "
                             "Separated by commas.")

        self.page.add(label="Permutations",
                      control=control.INT_CTRL,
                      name="mdmr_permutations",
                      type=dtype.NUM,
                      comment="Number of permutations to compute p-values.",
                      values=500)

        self.page.add(label="Parallel nodes",
                      control=control.INT_CTRL,
                      name="mdmr_parallel_nodes",
                      type=dtype.NUM,
                      comment="Number of parallel MDMR to run simultaneously.",
                      values=1)

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
