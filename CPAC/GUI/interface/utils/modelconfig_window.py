import wx
from . import generic_class
from .constants import control, dtype, substitution_map
import os
import yaml

ID_RUN = 11


class ModelConfig(wx.Frame):

    # this creates the wx.Frame mentioned above in the class declaration
    def __init__(self, parent, gpa_settings=None):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - Create New FSL FEAT Model",
            size=(900, 650))

        default_gpa_settings = {}
        default_gpa_settings['pipeline_dir'] = ''
        default_gpa_settings['participant_list'] = 'None'
        default_gpa_settings['output_dir'] = ''
        default_gpa_settings['work_dir'] = ''
        default_gpa_settings['log_dir'] = ''
        default_gpa_settings['FSLDIR'] = 'FSLDIR'
        default_gpa_settings['num_models_at_once'] = 1
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

        self.parent = parent

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        vertSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)

        self.window = wx.ScrolledWindow(self.panel, size=(-1,300))

        self.page = generic_class.GenericClass(self.window,
                                               " FSL FEAT Group-Level "
                                               "Analysis - Model Builder")

        self.page.add(label="Pipeline Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="pipeline_dir",
                      type=dtype.STR,
                      comment="Individual-level analysis pipeline output "
                              "directory.",
                      values=self.gpa_settings['pipeline_dir'])

        self.page.add(label="Participant List ",
                      control=control.COMBO_BOX,
                      name="participant_list",
                      type=dtype.STR,
                      comment="Full path to a list of participants to be "
                              "included in the model.\n\nThis should be a "
                              "text file with each participant ID on its "
                              "own line.",
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
                      comment="Working directory for the intermediates of the "
                              "FSL FEAT pipeline. Can be deleted afterwards.",
                      values=self.gpa_settings['work_dir'])

        self.page.add(label="Log Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="pipeline_dir",
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

        self.page.add(label="Number of Models to Run Simultaneously ",
                      control=control.INT_CTRL,
                      name='num_models_at_once',
                      type=dtype.NUM,
                      comment="This number depends on computing resources.",
                      values=self.gpa_settings['num_models_at_once'])

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
        load_pheno_btn = wx.Button(self.window, 2, 'Load Phenotype File',
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

        self.page.set_sizer()

        if 'group_sep' in self.gpa_settings.keys():
            for ctrl in self.page.get_ctrl_list():
                name = ctrl.get_name()
                if name == 'group_sep':
                    if self.gpa_settings['group_sep'] == True:
                        ctrl.set_value('On')
                    elif self.gpa_settings['group_sep'] == False:
                        ctrl.set_value('Off')

        mainSizer.Add(self.window, 1, wx.EXPAND)

        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        buffer = wx.StaticText(btnPanel, label="\t\t\t\t\t\t")
        hbox.Add(buffer)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel", (
            220, 10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)

        load = wx.Button(btnPanel, wx.ID_ADD, "Load Config", (
            200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.load, id=wx.ID_ADD)
        hbox.Add(load, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        next = wx.Button(btnPanel, 3, "Build Models", (200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.build_models, id=3)
        hbox.Add(next, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        # reminder: functions bound to buttons require arguments
        #           (self, event)
        btnPanel.SetSizer(hbox)

        mainSizer.Add(
            btnPanel, 0.5,  flag=wx.ALIGN_RIGHT | wx.RIGHT, border=20)

        self.panel.SetSizer(mainSizer)

        self.Show()

        # this fires only if we're coming BACK to this page from the second
        # page, and these parameters are already pre-loaded. this is to
        # automatically repopulate the 'Model Setup' checkbox grid and other
        # settings under it
        if self.gpa_settings['pheno_file'] != '':

            phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']))

            self.get_pheno_header(phenoFile)

            # update the 'Model Setup' box and populate it with the EVs and
            # their associated checkboxes for categorical and demean
            for ctrl in self.page.get_ctrl_list():

                name = ctrl.get_name()

                if name == 'model_setup':
                    ctrl.set_value(self.phenoHeaderItems)
                    ctrl.set_selection(self.gpa_settings['ev_selections'])

                if name == 'coding_scheme':
                    ctrl.set_value(self.gpa_settings['coding_scheme'])

                if name == 'mean_mask':
                    ctrl.set_value(self.gpa_settings['mean_mask'])

                if name == 'z_threshold':
                    try:
                        ctrl.set_value(self.gpa_settings['z_threshold'][0])
                    except IndexError:
                        ctrl.set_value(self.gpa_settings['z_threshold'])

                if name == 'p_threshold':
                    ctrl.set_value(self.gpa_settings['p_threshold'])

                if name == 'group_sep':
                    ctrl.set_value(self.gpa_settings['group_sep'])

                if name == 'grouping_var':
                    grouping_var = self.gpa_settings['grouping_var']

                    if isinstance(grouping_var, list) or "[" in grouping_var:
                        new_grouping_var = ""
                        for cov in grouping_var:
                            new_grouping_var += "{0},".format(cov)
                        new_grouping_var = new_grouping_var.rstrip(",")

                        ctrl.set_value(new_grouping_var)

                if ("list" in name) and (name != "participant_list"):
                    value = self.gpa_settings[name]

                    if isinstance(value, str):
                        value = value.replace("[","").replace("]","")
                        if "\"" in value:
                            value = value.replace("\"","")
                        if "'" in value:
                            value = value.replace("'","")
                        values = value.split(",")
                    else:
                        # instead, is a list- most likely when clicking
                        # "Back" on the modelDesign_window
                        values = value

                    ctrl.set_value(values)

    def cancel(self, event):
        self.Close()

    def display(self, win, msg):
        wx.MessageBox(msg, "Error")
        win.SetBackgroundColour("pink")
        win.SetFocus()
        win.Refresh()
        raise ValueError

    ''' button: LOAD SETTINGS '''
    def load(self, event):
        # when the user clicks 'Load Settings', which loads the
        # self.gpa_settings dictionary - it populates the values for both
        # windows, so when they hit Next, the next window is also populated

        dlg = wx.FileDialog(
            self, message="Choose the config fsl yaml file",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
            style=wx.OPEN | wx.CHANGE_DIR)

        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            config_map = yaml.safe_load(open(path, 'r'))
            s_map = dict((v, k) for k, v in substitution_map.items())

            # load the group analysis .yml config file (in dictionary form)
            # into the self.gpa_settings dictionary which holds all settings
            self.gpa_settings = config_map

            if self.gpa_settings is None:
                errDlgFileTest = wx.MessageDialog(
                    self, "Error reading file - group analysis " \
                          "configuration file appears to be blank.",
                    "File Read Error",
                    wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()
                raise Exception

            # repopulate the model setup checkbox grid, since this has to be
            # done specially
            if 'pheno_file' in self.gpa_settings.keys():

                phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']))
                self.get_pheno_header(phenoFile)

                # update the 'Model Setup' box and populate it with the EVs
                # and their associated checkboxes for categorical and demean
                for ctrl in self.page.get_ctrl_list():

                    if ctrl.get_name() == 'model_setup':
                        ctrl.set_value(self.phenoHeaderItems)
                        ctrl.set_selection(self.gpa_settings['ev_selections'])

            # populate the rest of the controls
            for ctrl in self.page.get_ctrl_list():

                name = ctrl.get_name()
                value = config_map.get(name)
                dtype = ctrl.get_datatype()

                # the model setup checkbox grid is the only one that doesn't
                # get repopulated the standard way. instead it is repopulated
                # by the code directly above

                if ("list" in name) and (name != "participant_list"):

                    mapped_values = [s_map.get(item)
                                 for item in value if s_map.get(item) != None]
                    if not mapped_values:
                        value = [str(item) for item in value]
                    else:
                        value = mapped_values

                    new_derlist = []

                    for val in value:
                        new_derlist.append(val)

                    if len(new_derlist) > 0:
                        ctrl.set_value(new_derlist)
                    else:
                        ctrl.set_value(None)

                elif name == 'z_threshold' or name == 'p_threshold':
                    try:
                        value = value[0]
                        ctrl.set_value(value)
                    except TypeError:
                        # if the user has put it in as a float and not a list
                        ctrl.set_value(str(value))

                elif name == 'group_sep':
                    value = s_map.get(value)
                    ctrl.set_value(value)

                elif name == 'grouping_var':
                    if isinstance(value, list) or "[" in value:
                        grouping_var = ""
                        for cov in value:
                            grouping_var += "{0},".format(cov)
                        grouping_var = grouping_var.rstrip(",")
                    else:
                        grouping_var = value

                    ctrl.set_value(grouping_var)

                elif name != 'model_setup' and name != 'derivative_list':
                    ctrl.set_value(value)

            dlg.Destroy()

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

    def build_models(self, event):
        from CPAC.pipeline.cpac_group_runner import build_feat_models
        from CPAC.utils.create_fsl_flame_preset import write_config_dct_to_yaml

        dialog_msg = 'Building your FSL-FEAT models. Check the terminal ' \
                     'window for details and progress.'
        dialog_title = 'Building models..'
        bld_dialog = wx.MessageDialog(self, dialog_msg, dialog_title,
                                      wx.OK | wx.ICON_INFORMATION)
        bld_dialog.ShowModal()
        bld_dialog.Destroy()

        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()
            val = str(ctrl.get_selection())
            if val in substitution_map.keys():
                val = substitution_map[val]
            if isinstance(val, list):
                new_val = []
                for v in val:
                    if v in substitution_map.keys():
                        v = substitution_map[v]
                    new_val.append(v)
                val = new_val
            elif isinstance(val, str):
                if 'derivative' in name:
                    val = val.replace('[', '').replace(']', '').replace(' ', '').replace("'", "")
                    val = val.split(',')
                    new_val = []
                    for v in val:
                        if v in substitution_map.keys():
                            v = substitution_map[v]
                        new_val.append(v)
                    val = new_val

            self.gpa_settings[name] = val

        group_config_path = os.path.join(self.gpa_settings['output_dir'],
                                         'group_config_{0}.yml'.format(self.gpa_settings['model_name']))
        write_config_dct_to_yaml(self.gpa_settings, group_config_path)
        retval = build_feat_models(group_config_path)

        if retval == 0:
            self.Close()
