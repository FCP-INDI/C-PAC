import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
import os
import pkg_resources as p
    
class CWASSettings(wx.ScrolledWindow):
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "Connectome-wide Association Studies")

        self.page.add(label="Run CWAS?",
                      control=control.CHOICE_BOX,
                      name="runCWAS",
                      type=dtype.LSTR,
                      comment="Used to determine if Connectome-wide Association Studies "
                              "will be added to the pipeline or not.",
                      values=["Off", "On"],
                      wkf_switch = True)

        self.page.add(label="Mask ROI File", 
                     control=control.COMBO_BOX,
                     name='cwas_roi_file', 
                     type=dtype.STR, 
                     values=str(""),
                     validation_req=False,
                     comment="Path to a mask file. Voxels outside of the mask will "
                             "be excluded from CWAS.")
               
        self.page.add(label="Regressor file", 
                     control=control.COMBO_BOX,
                     name="cwas_regressor_file", 
                     type=dtype.STR, 
                     values="",
                     comment="Path to a CSV file containing the phenotypic "
                             "regressor.")
        
        self.page.add(label="Regressor Participant Column Name",
                      control=control.TEXT_BOX,
                      name="cwas_regressor_participant_column",
                      type=dtype.STR,
                      comment="Name of the participants column in your "
                              "regressor file.",
                      values="")

        self.page.add(label="Regressor of Interest columns", 
                     control=control.TEXT_BOX, 
                     name="cwas_regressor_columns", 
                     type=dtype.STR,
                     values="",
                     comment="Columns from the CSV file indicating factor "
                             "variables. Other columns will be handled as covariates. "
                             "Separated by commas.")

        self.page.add(label="Permutations",
                      control=control.INT_CTRL,
                      name="cwas_permutations",
                      type=dtype.NUM,
                      comment="Number of permutation tests to run on the "
                              "Pseudo-F statistics.",
                      values=500)

        self.page.add(label="Parallel nodes",
                      control=control.INT_CTRL,
                      name="cwas_parallel_nodes",
                      type=dtype.NUM,
                      comment="Number of Nipype nodes created while "
                              "computing CWAS.  Dependent upon computing resources.",
                      values=1)

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
