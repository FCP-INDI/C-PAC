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
    
    
class GeneralGA(wx.ScrolledWindow):

    def __init__(self, parent, counter=0, gpa_settings=None):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        self.page = GenericClass(self, " General Group Analysis Options")

        default_gpa_settings = {}
        default_gpa_settings['pipeline_dir'] = ''
        default_gpa_settings['participant_list'] = 'None'
        default_gpa_settings['output_dir'] = ''
        default_gpa_settings['work_dir'] = ''
        default_gpa_settings['log_dir'] = ''
        default_gpa_settings['FSLDIR'] = 'FSLDIR'

        if not gpa_settings:
            self.gpa_settings = default_gpa_settings
        else:
            for key in default_gpa_settings.keys():
                if key not in gpa_settings.keys():
                    gpa_settings[key] = default_gpa_settings[key]
            self.gpa_settings = gpa_settings

        self.page.add(label="Pipeline Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="pipeline_dir",
                      type=dtype.STR,
                      comment="Individual-level analysis pipeline output "
                              "directory.",
                      values=self.gpa_settings['pipeline_dir'])

        self.page.add(label="Participant Inclusion (Optional) ",
                      control=control.COMBO_BOX,
                      name="participant_list",
                      type=dtype.STR,
                      comment="(Optional)\n\nFull path to a list of "
                              "participants to be included in the model. You "
                              "can use this to easily prune participants "
                              "from your model. In group-level analyses "
                              "involving phenotype files, this allows you to "
                              "prune participants without removing them from "
                              "the phenotype CSV/TSV file. This should be a "
                              "text file with one subject per line. An easy "
                              "way to manually create this file is to copy "
                              "the subjects column from your phenotype file.",
                      values=self.gpa_settings['participant_list'])

        self.page.add(label="Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="output_dir",
                      type=dtype.STR,
                      comment="Output directory for the results of FSL FEAT.",
                      values=self.gpa_settings['output_dir'])

        self.page.add(label="Working Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="work_dir",
                      type=dtype.STR,
                      comment="Much like the working directory for "
                              "individual-level analysis, this is where the "
                              "intermediate and working files will be stored "
                              "during your run.\n\nThis directory can be "
                              "deleted later on. However, saving this "
                              "directory allows the group analysis run to "
                              "skip steps that have been already completed, "
                              "in the case of re-runs.",
                      values=self.gpa_settings['work_dir'])

        self.page.add(label="Log Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="log_dir",
                      type=dtype.STR,
                      comment="Directory to write log information.",
                      values=self.gpa_settings['log_dir'])

        self.page.add(label="FSL Directory ",
                      control=control.COMBO_BOX,
                      name="FSLDIR",
                      type=dtype.STR,
                      comment="Directory of your FSL install. Can be kept as "
                              "'FSLDIR' unless you want to use a custom FSL "
                              "directory.",
                      values=self.gpa_settings['FSLDIR'])

        if not gpa_settings:
            self.gpa_settings = default_gpa_settings
        else:
            for key in default_gpa_settings.keys():
                if key not in gpa_settings.keys():
                    gpa_settings[key] = default_gpa_settings[key]
            self.gpa_settings = gpa_settings

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
            return self.counter

class GPASettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter=0, gpa_settings=None):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, " FSL-FEAT Group Analysis Options")

        default_gpa_settings = {}
        default_gpa_settings['run_fsl_feat'] = 1
        default_gpa_settings['num_models_at_once'] = 1
        default_gpa_settings['model_name'] = ''
        default_gpa_settings['pheno_file'] = ''
        default_gpa_settings['participant_id_label'] = ''
        default_gpa_settings['ev_selections'] = {'categorical': [''],
                                                 'demean': ['']}
        default_gpa_settings['design_formula'] = ''
        default_gpa_settings['mean_mask'] = ''
        default_gpa_settings['custom_roi_mask'] = 'None'
        default_gpa_settings['coding_scheme'] = ''
        default_gpa_settings['derivative_list'] = ''
        default_gpa_settings['sessions_list'] = []
        default_gpa_settings['series_list'] = []
        default_gpa_settings['group_sep'] = ''
        default_gpa_settings['grouping_var'] = 'None'
        default_gpa_settings['z_threshold'] = ''
        default_gpa_settings['p_threshold'] = ''

        if not gpa_settings:
            self.gpa_settings = default_gpa_settings
        else:
            for key in default_gpa_settings.keys():
                if key not in gpa_settings.keys():
                    gpa_settings[key] = default_gpa_settings[key]
            self.gpa_settings = gpa_settings

        self.page.add(label="Run FSL-FEAT ",
                      control=control.CHOICE_BOX,
                      name='run_fsl_feat',
                      type=dtype.LSTR,
                      comment="Run the FSL-FEAT pipeline.",
                      values=["Off", "On"],
                      wkf_switch=True)

        self.page.add(label="Number of Models to Run Simultaneously ",
                      control=control.INT_CTRL,
                      name='num_models_at_once',
                      type=dtype.NUM,
                      comment="This number depends on computing resources.",
                      values=self.gpa_settings['num_models_at_once'])

        self.page.add(label="Model Name ",
                      control=control.TEXT_BOX,
                      name="model_name",
                      type=dtype.STR,
                      comment="Specify a name for the new model. Output and working directories for group analysis, as well as the FLAME model files (.mat, .con, .grp, etc.) will be labeled with this name.",
                      values=self.gpa_settings['model_name'],
                      size=(200, -1))

        self.page.add(label="Phenotype/EV File ",
                      control=control.COMBO_BOX,
                      name="pheno_file",
                      type=dtype.STR,
                      comment="Full path to a .csv or .tsv file containing "
                              "EV information for each subject.",
                      values=self.gpa_settings['pheno_file'])

        self.page.add(label="Participant Column Name ",
                      control=control.TEXT_BOX,
                      name="participant_id_label",
                      type=dtype.STR,
                      comment="Name of the participants column in your EV "
                              "file.",
                      values=self.gpa_settings['participant_id_label'],
                      style=wx.EXPAND | wx.ALL,
                      size=(160, -1))
        
        load_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        load_pheno_btn = wx.Button(self, 2, 'Load Phenotype File',
                                   (220,10), wx.DefaultSize, 0)
        load_panel_sizer.Add(load_pheno_btn)

        self.Bind(wx.EVT_BUTTON, self.populateEVs, id=2)

        self.page.add_pheno_load_panel(load_panel_sizer)

        self.page.add(label="Model Setup ",
                      control=control.GPA_CHECKBOX_GRID,
                      name="model_setup",
                      type=10,
                      values='',
                      comment="A list of EVs from your phenotype file will "
                              "populate in this window. From here, you can "
                              "select whether the EVs should be treated as "
                              "categorical or if they should be demeaned "
                              "(continuous/non-categorical EVs only). "
                              "'MeanFD', 'MeanFD_Jenkinson', 'Measure Mean', "
                              "and 'Custom_ROI_Mean' will also appear in "
                              "this window automatically as options to be "
                              "used as regressors that can be included in "
                              "your model design. Note that the MeanFD and "
                              "mean of measure values are automatically "
                              "calculated and supplied by C-PAC via "
                              "individual-level analysis.",
                      size=(450, -1))

        self.page.add(label="Design Matrix Formula ",
                      control=control.TEXT_BOX,
                      name="design_formula",
                      type=dtype.STR,
                      comment="Specify the formula to describe your model design. Essentially, including EVs in this formula inserts them into the model. The most basic format to include each EV you select would be 'EV + EV + EV + ..', etc. You can also select to include MeanFD, MeanFD_Jenkinson, Measure_Mean, and Custom_ROI_Mean here. See the C-PAC User Guide for more detailed information regarding formatting your design formula.",
                      values= self.gpa_settings['design_formula'],
                      size=(450, -1))

        self.page.add(label="Custom ROI Mean Mask ",
                      control=control.COMBO_BOX,
                      name="custom_roi_mask",
                      type=dtype.STR,
                      comment="Optional: Full path to a NIFTI file containing one or more ROI masks. The means of the masked regions will then be computed for each subject's output and will be included in the model as regressors (one for each ROI in the mask file) if you include 'Custom_ROI_Mean' in the Design Matrix Formula.",
                      values=self.gpa_settings['custom_roi_mask'])

        self.page.add(label = "Select Derivatives ",
                    control = control.CHECKLIST_BOX,
                    name = "derivative_list",
                    type = dtype.LSTR,
                    values = ['ALFF',
                              'f/ALFF',
                              'ReHo',
                              'ROI Average SCA',
                              'Dual Regression',
                              'Multiple Regression SCA',
                              'Network Centrality',
                              'VMHC'],
                    comment = "Select which derivatives you would like to "
                              "include when running group analysis.\n\nWhen "
                              "including Dual Regression, make sure to "
                              "correct your P-value for the number of maps "
                              "you are comparing.\n\nWhen including Multiple "
                              "Regression SCA, you must have more degrees of "
                              "freedom (subjects) than there were time "
                              "series.",
                    size = (350,180))

        self.page.add(label="Coding Scheme ", 
                     control=control.CHOICE_BOX, 
                     name="coding_scheme", 
                     type=dtype.LSTR, 
                     comment="Choose the coding scheme to use when generating your model. 'Treatment' encoding is generally considered the typical scheme. Consult the User Guide for more information.", 
                     values=["Treatment", "Sum"])

        self.page.add(label="Mask for Means Calculation ", 
                 control=control.CHOICE_BOX, 
                 name='mean_mask', 
                 type=dtype.LSTR, 
                 comment = "Choose whether to use a group mask or individual-specific mask when calculating the output means to be used as a regressor.\n\nThis only takes effect if you include the 'Measure_Mean' or 'Custom_ROI_Mean' regressors in your Design Matrix Formula.", 
                 values=["Group Mask","Individual Mask"])

        self.page.add(label="Z threshold ", 
                     control=control.FLOAT_CTRL, 
                     name='z_threshold', 
                     type=dtype.NUM, 
                     comment="Only voxels with a Z-score higher than this value will be considered significant.", 
                     values=2.3)

        self.page.add(label="Cluster Significance Threshold ", 
                     control=control.FLOAT_CTRL, 
                     name='p_threshold', 
                     type=dtype.NUM, 
                     comment="Significance threshold (P-value) to use when doing cluster correction for multiple comparisons.", 
                     values=0.05)

        self.page.add(label="Model Group Variances Separately ",
                      control=control.CHOICE_BOX,
                      name='group_sep',
                      type=dtype.NUM,
                      comment="Specify whether FSL should model the variance for each group separately.\n\nIf this option is enabled, you must specify a grouping variable below.",
                      values=['Off', 'On'])

        self.page.add(label="Grouping Variable ",
                      control=control.TEXT_BOX,
                      name="grouping_var",
                      type=dtype.STR,
                      comment="The name of the EV that should be used to group subjects when modeling variances.\n\nIf you do not wish to model group variances separately, set this value to None.",
                      values=self.gpa_settings['grouping_var'],
                      size=(160, -1))

        self.page.add(label = 'Sessions (Repeated Measures Only) ',
                      control = control.LISTBOX_COMBO,
                      name = 'sessions_list',
                      type = dtype.LSTR,
                      values = self.gpa_settings['sessions_list'],
                      comment = 'Enter the session names in your dataset ' \
                                'that you wish to include within the same ' \
                                'model (this is for repeated measures/' \
                                'within-subject designs).\n\nTip: These ' \
                                'will be the names listed as "unique_id" in '\
                                'the original individual-level participant ' \
                                'list, or the labels in the original data ' \
                                'directories you marked as {session} while ' \
                                'creating the CPAC participant list.',
                      size = (200,100),
                      combo_type = 7)

        self.page.add(label = 'Series/Scans (Repeated Measures Only) ',
                      control = control.LISTBOX_COMBO,
                      name = 'series_list',
                      type = dtype.LSTR,
                      values = self.gpa_settings['series_list'],
                      comment = 'Enter the series names in your dataset ' \
                                'that you wish to include within the same ' \
                                'model (this is for repeated measures/' \
                                'within-subject designs).\n\nTip: These ' \
                                'will be the labels listed under "func:"" '\
                                'in the original individual-level ' \
                                'participant list, or the labels in the ' \
                                'original data directories you marked as ' \
                                '{series} while creating the CPAC ' \
                                'participant list.',
                      size = (200,100),
                      combo_type = 8)

        if 'group_sep' in self.gpa_settings.keys():
            for ctrl in self.page.get_ctrl_list():
                name = ctrl.get_name()
                if name == 'group_sep':
                    if self.gpa_settings['group_sep'] == True:
                        ctrl.set_value('On')
                    elif self.gpa_settings['group_sep'] == False:
                        ctrl.set_value('Off')
                        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_pheno_header(self, pheno_file_obj):
        phenoHeaderString = pheno_file_obj.readline().rstrip('\r\n')
        phenoHeaderString = phenoHeaderString.replace(" ", "_")
        phenoHeaderString = phenoHeaderString.replace("/","_")
        
        if ',' in phenoHeaderString:
            self.phenoHeaderItems = phenoHeaderString.split(',')
            
        
        elif '\t' in phenoHeaderString:
            self.phenoHeaderItems = phenoHeaderString.split('\t')
            
        
        else:
            self.phenoHeaderItems = [phenoHeaderString]
        
        if self.gpa_settings['participant_id_label'] in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(self.gpa_settings['participant_id_label'])
        else:
            print('Header labels found:\n{0}'.format(self.phenoHeaderItems))
            err = 'Please enter the name of the participant ID column as ' \
                  'it is labeled in the phenotype file.'
            print(err)
            errSubID = wx.MessageDialog(self, err,
                                        'Blank/Incorrect Subject Header Input',
                                        wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
            raise Exception

    ''' button: LOAD PHENOTYPE FILE '''
    def populateEVs(self, event):
        # this runs when the user clicks 'Load Phenotype File'

        import os

        if self.gpa_settings is None:
            self.gpa_settings = {}

        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()
            self.gpa_settings[name] = str(ctrl.get_selection())

        ### CHECK PHENOFILE if can open etc.
        # function for file path checking
        def testFile(filepath, paramName):
            try:
                fileTest = open(filepath)
                fileTest.close()
            except:
                errDlgFileTest = wx.MessageDialog(
                    self, 'Error reading file - either it does not exist or '
                          'you do not have read access. \n\n' \
                          'Parameter: %s' % paramName,
                    'File Access Error',
                    wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()
                raise Exception

        '''
        # get participant inclusion list
        self.subs = []
        if '/' in self.gpa_settings['participant_list'] and \
                '.' in self.gpa_settings['participant_list']:
            testFile(self.gpa_settings['participant_list'], 'Participant List')
            subFile = open(
                os.path.abspath(self.gpa_settings['participant_list']))
            sub_IDs = subFile.readlines()
            for sub in sub_IDs:
                self.subs.append(sub.rstrip("\n"))
        '''

        # deal with phenotype file
        testFile(self.gpa_settings['pheno_file'], 'Phenotype/EV File')

        phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']),"rU")
        self.get_pheno_header(phenoFile)

        for ctrl in self.page.get_ctrl_list():

            # update the 'Model Setup' box and populate it with the EVs and
            # their associated checkboxes for categorical and demean
            if ctrl.get_name() == 'model_setup':

                evs_for_checkbox = []
                for EV in self.phenoHeaderItems:
                    evs_for_checkbox.append(EV)

                if "session" in evs_for_checkbox:
                    evs_for_checkbox.remove("session")

                if "series" in evs_for_checkbox:
                    evs_for_checkbox.remove("series")

                ctrl.set_value(evs_for_checkbox)

            # populate the design formula text box with a formula which
            # includes all of the EVs, and two of the measures (MeanFD and
            # the measure/derivative mean) - the user can edit this if they
            # need to, obviously
            if ctrl.get_name() == 'design_formula':

                formula_string = ''

                for EV in self.phenoHeaderItems:
                    if formula_string == '':
                        formula_string = EV
                    else:
                        formula_string = formula_string + ' + ' + EV

                formula_string = formula_string + ' + MeanFD_Jenkinson'

                ctrl.set_value(formula_string)

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

        import os

        self.counter = counter
        
        self.page = GenericClass(self, " PyBASC - Bootstrapped Analysis of "
                                       "Stable Clusters (BASC)")

        fsl = os.environ.get('FSLDIR')
        if fsl == None:
            fsl = "$FSLDIR"

        self.page.add(label="Run BASC ", 
                      control=control.CHOICE_BOX, 
                      name='run_basc',
                      type=dtype.LSTR, 
                      comment="Run Bootstrap Analysis of Stable Clusters", 
                      values=["Off", "On"],
                      wkf_switch=True)

        self.page.add(label="Series/Scan Inclusion (Optional) ",
                      control=control.COMBO_BOX,
                      name='basc_scan_inclusion',
                      type=dtype.STR,
                      values="None",
                      comment="If there are multiple series/scans in any "
                              "of the pipeline outputs for which PyBASC is "
                              "being run, and you only want to run for some "
                              "of them, you can list them here - scan labels "
                              "separated by commas (ex. 'rest_run-1, "
                              "rest_run-3').\n\nIf nothing is listed, PyBASC "
                              "will run once for each scan, for all "
                              "available scans.")

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

        self.page.add(label="Standard Brain only Template (functional resolution) ",
                      control=control.COMBO_BOX,
                      name='template_brain_only_for_func',
                      type=dtype.STR,
                      values = str(os.path.join(fsl,"data/standard/MNI152_T1_${basc_resolution}_brain.nii.gz")),
                      comment="Standard FSL Skull Stripped Template.")

        self.page.add(label="ROI Mask File 1",
                      control=control.COMBO_BOX,
                      name='basc_roi_mask_file',
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

        self.page.add(label="ROI Mask File for Cross-Clustering",
                      control=control.COMBO_BOX,
                      name='basc_cross_cluster_mask_file',
                      type=dtype.STR,
                      values="None",
                      comment="")

        self.page.add(label="Similarity Metric ",
                     control=control.CHOICE_BOX,
                     name='basc_similarity_metric_list',
                     type = dtype.LSTR,
                     comment="",
                     values=['correlation', 'euclidean', 'cityblock',
                             'cosine'])

        self.page.add(label="Number of Time Series Bootstraps ",
                      control=control.INT_CTRL,
                      name='basc_timeseries_bootstrap_list',
                      type=dtype.NUM,
                      comment="Number of bootstraps to apply to individual "
                              "time series.",
                      values=100)
            
        self.page.add(label="Number of Dataset Bootstraps ",
                      control=control.INT_CTRL,
                      name='basc_dataset_bootstrap_list',
                      type=dtype.NUM,
                      comment="Number of bootstraps to apply to the original "
                              "dataset.",
                      values=30)

        self.page.add(label="Number of Clusters ",
                      control=control.INT_CTRL,
                      name='basc_n_clusters_list',
                      type=dtype.NUM,
                      comment="Number of clusters to create during "
                              "clustering at both the individual and group "
                              "levels.",
                      values=2)

        self.page.add(label="Affinity Threshold ",
                      control=control.TEXT_BOX,
                      name='basc_affinity_thresh',
                      type=dtype.LNUM,
                      values="0.0",
                      validator=CharValidator("no-alpha"),
                      comment="",
                      size=(100, -1))

        self.page.add(label="Output Sizes ",
                      control=control.INT_CTRL,
                      name='basc_output_sizes',
                      type=dtype.NUM,
                      comment="",
                      values=800)

        self.page.add(label="Run Cross-Clustering ",
                      control=control.CHOICE_BOX,
                      name='basc_cross_cluster',
                      type=dtype.BOOL,
                      values=["True", "False"],
                      comment="")

        self.page.add(label="Blocklength List ",
                      control=control.INT_CTRL,
                      name='basc_blocklength_list',
                      type=dtype.NUM,
                      comment="",
                      values=1)

        self.page.add(label="Group Dimension Reduce ",
                      control=control.CHOICE_BOX,
                      name='basc_group_dim_reduce',
                      type=dtype.BOOL,
                      values=["False", "True"],
                      comment="")

        self.page.set_sizer()
        parent.get_page_list().append(self)

