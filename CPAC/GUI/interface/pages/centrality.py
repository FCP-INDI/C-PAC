import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class Centrality(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/centrality.html'))            
            
    def get_counter(self):
        return self.counter
            


class CentralitySettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Network Centrality Options")
        
        self.page.add(label="Calculate Network Centrality Measures ", 
                 control=control.CHOICE_BOX, 
                 name='runNetworkCentrality', 
                 type=dtype.LSTR, 
                 comment="Calculate Degree, Eigenvector Centrality, or Functional Connectivity Density.", 
                 values=["Off","On"],
                 wkf_switch = True)

        self.page.add(label="Mask Specification File ", 
                     control=control.COMBO_BOX, 
                     name='templateSpecificationFile', 
                     type=dtype.STR, 
                     values = "",
                     comment="Full path to a NIFTI file describing the " \
                             "mask. Centrality will be calculated for all " \
                             "voxels within the mask.")
        
        # --- Add degree centrality ---
        self.page.add(label = "Degree Centrality Weight Options",
                      control = control.CHECKLIST_BOX,
                      name = "degWeightOptions",
                      type = dtype.LBOOL,
                      values = ['Binarized','Weighted'],
                      comment = "Enable/Disable degree centrality by selecting the connectivity weights")

        self.page.add(label="Degree Centrality Threshold Type", 
                     control=control.CHOICE_BOX, 
                     name='degCorrelationThresholdOption', 
                     type = dtype.LSTR, 
                     comment="Select the type of threshold used when creating the degree centrality adjacency matrix.", 
                     values=["Sparsity threshold", "Significance threshold", "Correlation threshold"])
        
        self.page.add(label="Degree Centrality Threshold Value", 
                     control=control.FLOAT_CTRL, 
                     name='degCorrelationThreshold', 
                     type=dtype.NUM, 
                     comment="Based on the Threshold Type selected above, enter a Threshold Value.\n\nP-value for Significance Threshold\nSparsity value for Sparsity Threshold\nPearson's r value for Correlation Threshold", 
                     values=0.001)
        
        # --- Add eigenvector centrality ---
        self.page.add(label = "Eigenvector Centrality Weight Options",
                      control = control.CHECKLIST_BOX,
                      name = "eigWeightOptions",
                      type = dtype.LBOOL,
                      values = ['Binarized','Weighted'],
                      comment = "Enable/Disable eigenvector centrality by selecting the connectivity weights"
                      )
        self.page.add(label="Eigenvector Centrality Threshold Type", 
                     control=control.CHOICE_BOX, 
                     name='eigCorrelationThresholdOption', 
                     type = dtype.LSTR, 
                     comment="Select the type of threshold used when creating the eigenvector centrality adjacency matrix.", 
                     values=["Sparsity threshold", "Significance threshold", "Correlation threshold"])
        
        self.page.add(label="Eigenvector Centrality Threshold Value", 
                     control=control.FLOAT_CTRL, 
                     name='eigCorrelationThreshold', 
                     type=dtype.NUM, 
                     comment="Based on the Threshold Type selected above, enter a Threshold Value.\n\nP-value for Significance Threshold\nSparsity value for Sparsity Threshold\nPearson's r value for Correlation Threshold", 
                     values=0.001)
        
        # --- Add lFCD ---
        self.page.add(label = "Local Functional Connectivity Density Weight Options",
                      control = control.CHECKLIST_BOX,
                      name = "lfcdWeightOptions",
                      type = dtype.LBOOL,
                      values = ['Binarized','Weighted'],
                      comment = "Enable/Disable lFCD by selecting the connectivity weights"
                      )
        self.page.add(label="Local Functional Connectivity Density Threshold Type ", 
                     control=control.CHOICE_BOX, 
                     name='lfcdCorrelationThresholdOption', 
                     type = dtype.LSTR, 
                     comment="Select the type of threshold used when creating the lFCD adjacency matrix.", 
                     values=["Significance threshold", "Correlation threshold"])
        
        self.page.add(label="Local Functional Connectivity Density Threshold Value", 
                     control=control.FLOAT_CTRL, 
                     name='lfcdCorrelationThreshold', 
                     type=dtype.NUM, 
                     comment="Based on the Threshold Type selected above, enter a Threshold Value.\n\nP-value for Significance Threshold\nSparsity value for Sparsity Threshold\nPearson's r value for Correlation Threshold", 
                     values=0.001)
               
        self.page.add(label="Maximum RAM Use (GB) ", 
                     control=control.FLOAT_CTRL, 
                     name='memoryAllocatedForDegreeCentrality', 
                     type=dtype.NUM, 
                     comment="Maximum amount of RAM (in GB) to be used when calculating Degree Centrality.\n\nCalculating Eigenvector Centrality will require additional memory based on the size of the mask or number of ROI nodes.", 
                     values=2)
        
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
