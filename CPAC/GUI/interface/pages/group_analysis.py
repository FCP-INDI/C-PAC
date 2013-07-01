import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p

class GroupAnalysis(wx.html.HtmlWindow):
    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/group_analysis.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/fsl_ga.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/fsl_ga.html')
#            else:
#                self.LoadFile('html/group_analysis.html')
#        except:
#            self.LoadFile('html/group_analysis.html')
#            
            
    def get_counter(self):
        return self.counter
    
    
class GPASettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "FSL/FEAT Group Analysis Options")
        
        self.page.add(label="Run Group Analysis ", 
                      control=control.CHOICE_BOX, 
                      name='runGroupAnalysis', 
                      type=dtype.LSTR, 
                      comment="Run group analysis using FSL/FEAT.", 
                      values=["Off","On"],
                      wkf_switch = True)
        
        self.page.add(label = "Select Derivatives ",
                    control = control.CHECKLIST_BOX,
                    name = "derivativeList",
                    type = dtype.LSTR,
                    values = ['ROI Average SCA', 
                              'Voxelwise SCA',
                              'Multiple Regression SCA',
                              'Dual Regression',
                              'VMHC',
                              'ALFF', 
                              'f/ALFF', 
                              'ReHo',
                              'Network Centrality'],
                    comment = "Select which derivatives you would like to include when running group analysis.\n\nWhen including Dual Regression, make sure to correct your P-value for the number of maps you are comparing.\n\nWhen including Multiple Regression SCA, you must have more degrees of freedom (subjects) than there were time series.",
                    size = (220,160))
 
        self.page.add(label = "Models to Run ",
                      control = control.LISTBOX_COMBO,
                      name = 'modelFile',
                      type = dtype.LSTR,
                      values = "",
                      comment="Use the + to add FSL models to be run.",
                      size = (400,100),
                      combo_type = 3)

        self.page.add(label="Models Contain F-tests ", 
                 control=control.CHOICE_BOX, 
                 name='fTest', 
                 type=dtype.BOOL, 
                 comment = "Set this option to True if any of the models specified above contain F-tests.", 
                 values=["False","True"])
        
        self.page.add(label="Z threshold ", 
                     control=control.FLOAT_CTRL, 
                     name='zThreshold', 
                     type=dtype.NUM, 
                     comment="Only voxels with a Z-score higher than this value will be considered significant.", 
                     values=2.3)

        self.page.add(label="Cluster Significance Threshold ", 
                     control=control.FLOAT_CTRL, 
                     name='pThreshold', 
                     type=dtype.NUM, 
                     comment="Significance threshold (P-value) to use when doing cluster correction for multiple comparisons.", 
                     values=0.05)

                
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        
        
class BASC(wx.html.HtmlWindow):
    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/basc.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/basc.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/basc.html')
#            else:
#                self.LoadFile('html/basc.html')
#        except:
#            self.LoadFile('html/basc.html')
            
            
    def get_counter(self):
        return self.counter
    
    
    
class BASCSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Bootstrap Analysis of Stable Clusters (BASC)")
        
        self.page.add(label="Run BASC ", 
                      control=control.CHOICE_BOX, 
                      name='runBASC', 
                      type=dtype.LSTR, 
                      comment="Run Bootstrap Analysis of Stable Clusters", 
                      values=["Off","On"],
                      wkf_switch = True)
        
        self.page.add(label="Mask File ", 
                     control=control.COMBO_BOX, 
                     name='bascROIFile', 
                     type=dtype.STR, 
                     values = "None",
                     comment="Full path to a mask file to be used when running BASC. Voxels outside this mask will be excluded from analysis.\n\nIf you do not wish to use a mask, set this field to None.\n\nNote: BASC is very computationally intensive, we strongly recommend you limit your analysis to specific brain areas of interest.")

        self.page.add(label= "Number of Time Series Bootstraps ",
                 control=control.INT_CTRL, 
                 name='bascTimeseriesBootstraps', 
                 type=dtype.NUM, 
                 comment="Number of bootstraps to apply to individual time series.", 
                 values=100)
            
        self.page.add(label= "Number of Dataset Bootstraps ",
                 control=control.INT_CTRL, 
                 name='bascDatasetBootstraps', 
                 type=dtype.NUM, 
                 comment="Number of bootstraps to apply to the original dataset.", 
                 values=100)
        
        self.page.add(label="Correlation Threshold File ", 
                     control=control.COMBO_BOX, 
                     name='bascAffinityThresholdFile', 
                     type=dtype.STR, 
                     values = "",
                     comment="Path to a text file containing correlation threshold for each subject. These thresholds will be applied to the correlation matrix before clustering.\n\nThis file should contain one value per line, with each line corresponding to the subject on the same line in the group analysis subject list file.\n\nIn most cases, the same threshold can be used for all subjects. Different thresholds are useful when subjects have time series of different lengths.")
        
        self.page.add(label= "Number of Clusters ",
                 control=control.INT_CTRL, 
                 name='bascClusters', 
                 type=dtype.NUM, 
                 comment="Number of clusters to create during clustering at both the individual and group levels.", 
                 values=6)
        
        
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        
class CWAS(wx.html.HtmlWindow):
    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/cwas.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/cwas.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/cwas.html')
#            else:
#                self.LoadFile('html/cwas.html')
#        except:
#            self.LoadFile('html/cwas.html')
            
            
    def get_counter(self):
        return self.counter
    
    
    
class CWASSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Connectome-wide Association Studies (CWAS)")
        
        self.page.add(label="Run CWAS:", 
                      control=control.CHOICE_BOX, 
                      name='runCWAS', 
                      type=dtype.LSTR, 
                      comment="Run CWAS", 
                      values=["Off","On"],
                      wkf_switch = True)
        
        self.page.add(label="CWAS ROI File:", 
                     control=control.COMBO_BOX, 
                     name='cwasROIFile', 
                     type=dtype.STR, 
                     values = "/path/to/cwas_mask_file",
                     comment="Path to a mask file. Voxels outside this mask will be excluded from CWAS.")
        
        self.page.add(label="CWAS Regressor File:", 
                     control=control.COMBO_BOX, 
                     name='cwasRegressorFile', 
                     type=dtype.STR, 
                     values= '/path/to/cwas_regressor_file',
                     comment = "Path to a text file containing phenotypic regressor.")
        
        self.page.add(label= "CWAS FSamples",
                 control=control.INT_CTRL, 
                 name='cwasFSamples', 
                 type=dtype.NUM, 
                 comment="Number of permutation tests to run on the Psuedo-F statistic.", 
                 values=5000)
            
        self.page.add(label= "CWAS Parallel Nodes",
                 control=control.INT_CTRL, 
                 name='cwasParallelNodes', 
                 type=dtype.NUM, 
                 comment="Number of NiPype nodes to be created while computing CWAS.\n"\
                         "This number depends on computing resources.", 
                 values=10)
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter