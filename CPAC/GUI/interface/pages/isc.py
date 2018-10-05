import os

import pkg_resources as p
import wx
import wx.html

from CPAC.GUI.interface.utils.constants import control, dtype
from CPAC.GUI.interface.utils.generic_class import GenericClass


class ISCSettings(wx.ScrolledWindow):

    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "Inter-subject Correlation")

        self.page.add(label="Run ISC?",
                      control=control.CHOICE_BOX,
                      name="runISC",
                      type=dtype.LSTR,
                      comment="Used to determine if Inter-subject Correlation (ISC) "
                              "will be added to the pipeline or not.",
                      values=["Off", "On"],
                      wkf_switch = True)

        self.page.add(label="Run ISFC?",
                      control=control.CHOICE_BOX,
                      name="runISFC",
                      type=dtype.LSTR,
                      comment="Used to determine if Inter-subject Functional Correlation (ISFC) "
                              "will be added to the pipeline or not.",
                      values=["Off", "On"])

        self.page.add(label="ROI analysis",
                      control=control.CHOICE_BOX,
                      name="isc_level_roi",
                      type=dtype.LSTR,
                      comment="Used to determine if the ISC and ISFC will run in the ROI level.",
                      values=["Off", "On"])

        self.page.add(label="Voxel-wise analysis",
                      control=control.CHOICE_BOX,
                      name="isc_level_voxel",
                      type=dtype.LSTR,
                      comment="Used to determine if the ISC and ISFC will run in the voxel level. "
                              "Depending on the image resolution, it may take several hours "
                              "and consume a great amount of available memory.",
                      values=["Off", "On"])

        self.page.add(label="Voxel-wise analysis: Standard Deviation filter", 
                      control=control.FLOAT_CTRL, 
                      name='isc_level_voxel_std_filter', 
                      type=dtype.NUM, 
                      comment="Filter out voxels that, in the correlation distribution, is greater "
                              "then the informed standard deviation. Zero value will disable the filter.", 
                      values=0.0)

        self.page.add(label="Permutations",
                      control=control.INT_CTRL,
                      name="isc_permutations",
                      type=dtype.NUM,
                      comment="Number of permutation tests to compute the statistics.",
                      values=1000)

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
