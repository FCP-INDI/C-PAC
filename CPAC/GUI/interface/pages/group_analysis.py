import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class GroupAnalysis(wx.html.HtmlWindow):
    def __init__(self, parent, counter=0):
        wx.html.HtmlWindow.__init__(self, parent,
                                    style=wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        self.counter = counter
        self.LoadFile(p.resource_filename('CPAC',
                                          'GUI/resources/html/group_analysis.html'))

    def get_counter(self):
        return self.counter
    
    
class GPASettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, " FSL/FEAT Group Analysis Options")

        self.page.add(label="Run FSL FEAT ",
                      control=control.CHOICE_BOX,
                      name='run_fsl_feat',
                      type=dtype.LSTR,
                      comment="Run FSL FEAT group-level analysis.",
                      values=["Off", "On"],
                      wkf_switch=True)

        self.page.add(label="Number of Models to Run Simultaneously ",
                      control=control.INT_CTRL,
                      name='numGPAModelsAtOnce',
                      type=dtype.NUM,
                      comment="This number depends on computing resources.",
                      values=1)

        self.page.add(label="Models to Run ",
                      control=control.LISTBOX_COMBO,
                      name='modelConfigs',
                      type=dtype.LSTR,
                      values="",
                      comment="Use the + to add FSL model configuration to "
                              "be run.",
                      size=(400,100),
                      combo_type=3)
                        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        
        
class BASC(wx.html.HtmlWindow):
    def __init__(self, parent, counter=0):

        wx.html.HtmlWindow.__init__(self, parent,
                                    style=wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC',
                                          'GUI/resources/html/basc.html'))

    def get_counter(self):
        return self.counter
    

class BASCSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, " PyBASC - Bootstrapped Analysis of "
                                       "Stable Clusters (BASC)")
        
        self.page.add(label="Run BASC ", 
                      control=control.CHOICE_BOX, 
                      name='run_basc',
                      type=dtype.LSTR, 
                      comment="Run Bootstrap Analysis of Stable Clusters", 
                      values=["Off", "On"],
                      wkf_switch=True)

        self.page.add(label="Output File Resolution ",
                      control=control.CHOICE_BOX,
                      name='basc_resolution',
                      type=dtype.STR,
                      values=["4mm", "3mm", "2mm", "1mm"],
                      comment="")

        self.page.add(label="Maximum Processor Use ",
                      control=control.INT_CTRL,
                      name='basc_proc',
                      type=dtype.NUM,
                      comment="Maximum amount of processors to use while "
                              "performing BASC.",
                      values=2)

        self.page.add(label="Maximum RAM Use (GB) ",
                      control=control.INT_CTRL,
                      name='basc_memory',
                      type=dtype.NUM,
                      comment="Maximum amount of RAM (in GB) to be used when "
                              "running BASC.",
                      values=4)

        self.page.add(label="ROI File ",
                      control=control.COMBO_BOX,
                      name='basc_roi_file',
                      type=dtype.STR,
                      values="None",
                      comment="Full path to a mask file to be used when "
                              "running BASC. Voxels outside this mask will "
                              "be excluded from analysis.\n\nIf you do not "
                              "wish to use a mask, set this field to None."
                              "\n\nNote: BASC is very computationally "
                              "intensive, we strongly recommend you limit "
                              "your analysis to specific brain areas of "
                              "interest.")

        self.page.add(label="Second ROI File ",
                      control=control.COMBO_BOX,
                      name='basc_roi_file_two',
                      type=dtype.STR,
                      values="None",
                      comment="")

        self.page.add(label="Number of Time Series Bootstraps ",
                      control=control.INT_CTRL,
                      name='basc_timeseries_bootstraps',
                      type=dtype.NUM,
                      comment="Number of bootstraps to apply to individual "
                              "time series.",
                      values=100)
            
        self.page.add(label="Number of Dataset Bootstraps ",
                      control=control.INT_CTRL,
                      name='basc_dataset_bootstraps',
                      type=dtype.NUM,
                      comment="Number of bootstraps to apply to the original "
                              "dataset.",
                      values=100)
        
        self.page.add(label="Number of Clusters ",
                      control=control.INT_CTRL,
                      name='basc_clusters',
                      type=dtype.NUM,
                      comment="Number of clusters to create during "
                              "clustering at both the individual and group "
                              "levels.",
                      values=10)

        self.page.add(label="Affinity Threshold ",
                      control=control.TEXT_BOX,
                      name='basc_affinity_threshold',
                      type=dtype.LNUM,
                      values="0.0",
                      validator=CharValidator("no-alpha"),
                      comment="",
                      size=(100, -1))

        self.page.add(label="Output Size ",
                      control=control.INT_CTRL,
                      name='basc_output_size',
                      type=dtype.NUM,
                      comment="",
                      values=800)

        self.page.add(label="Run Cross-Clustering ",
                      control=control.CHOICE_BOX,
                      name='basc_cross_clustering',
                      type=dtype.BOOL,
                      values=["True", "False"],
                      comment="")

        self.page.add(label="Participant Inclusion (Optional) ",
                      control=control.COMBO_BOX,
                      name='basc_inclusion',
                      type=dtype.STR,
                      values="None",
                      comment="Full path to a text file listing which "
                              "participant IDs you want included in the "
                              "analysis, with one ID on each line.\n\nTip: "
                              "A sample group-level participant inclusion "
                              "text file is generated when you first create "
                              "your data configuration.")

        self.page.add(label="Pipeline Inclusion (Optional) ",
                      control=control.COMBO_BOX,
                      name='basc_pipeline',
                      type=dtype.STR,
                      values="None",
                      comment="If there are multiple pipeline output "
                              "directories, and you only want to run BASC on "
                              "one or some of them, you can list them here - "
                              "pipeline names separated by commas (check "
                              "the output directory of your individual-level "
                              "analysis run to see which pipeline "
                              "directories are available).\n\nIf nothing is "
                              "listed, all available pipelines will be run.")

        self.page.add(label="Series/Scan Inclusion (Optional) ",
                      control=control.COMBO_BOX,
                      name='basc_scan_inclusion',
                      type=dtype.STR,
                      values="None",
                      comment="If there are multiple series or scans in any "
                              "of the pipeline outputs for which PyBASC is "
                              "being run, and you only want to run for some "
                              "of them, you can list them here - scan labels "
                              "separated by commas (ex. 'rest_run-1, "
                              "rest_run-3').\n\nIf nothing is listed, all "
                              "available pipelines will be run.")

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
