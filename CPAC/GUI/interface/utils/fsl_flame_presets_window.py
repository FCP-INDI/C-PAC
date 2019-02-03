import wx
import generic_class
from .constants import control, dtype, substitution_map
import os
import yaml

from ....utils import create_fsl_flame_preset

ID_RUN = 11


class FlamePresetsOne(wx.Frame):

    # this creates the wx.Frame mentioned above in the class declaration
    def __init__(self, parent, gpa_settings=None):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - FSL FEAT Presets",
            size=(900, 550))

        self.parent = parent

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        vertSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)
        self.window = wx.ScrolledWindow(self.panel, size=(-1, 255))
        self.page = generic_class.GenericClass(self.window,
                                               " FSL FEAT Group-Level "
                                               "Analysis - Model Presets")

        if not gpa_settings:
            # if this window is being opened for the first time
            self.gpa_settings = {}
            self.gpa_settings["flame_preset"] = ""
            self.gpa_settings["pipeline_dir"] = ""
            self.gpa_settings["participant_list"] = "None"
            self.gpa_settings['derivative_list'] = ""
            self.gpa_settings['z_threshold'] = 2.3
            self.gpa_settings['p_threshold'] = 0.05
            self.gpa_settings["model_name"] = ""
            self.gpa_settings["output_dir"] = ""
        else:
            # if we're coming back from the "next" window
            self.gpa_settings = gpa_settings

        self.page.add(label="Choose Preset: ",
                      control=control.CHOICE_BOX,
                      name="flame_preset",
                      type=dtype.LSTR,
                      comment="Select the type of preset you'd like to "
                              "generate. The preset generator will prompt "
                              "you for more information relevant to the "
                              "type of preset you selected on the next "
                              "window (except for the Single Group Average "
                              "(One-Sample T-Test), which needs no "
                              "additional information).",
                      values=["Single Group Average (One-Sample T-Test)",
                              "Single Group Average with Additional Covariate",
                              "Unpaired Two-Group Difference (Two-Sample Unpaired T-Test)",
                              "Paired Two-Group Difference (Two-Sample Paired T-Test)",
                              "Tripled Two-Group Difference ('Tripled' T-Test)"])

        self.page.add(label="Pipeline Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="pipeline_dir",
                      type=dtype.STR,
                      comment="Full path to the individual-level pipeline "
                              "output directory you wish to run FSL-FEAT for. "
                              "\n\nThis will be the path to a directory titled "
                              "'pipeline_{name}'.",
                      values=self.gpa_settings['pipeline_dir'])

        self.page.add(label="[Optional] Participant Inclusion List ",
                      control=control.COMBO_BOX,
                      name="participant_list",
                      type=dtype.STR,
                      comment="[Optional] Full path to the group-level "
                              "analysis participant list text file. Use this "
                              "to quickly prune participants from your "
                              "analysis. This should be a text "
                              "file with each participant ID you "
                              "want included in the model, on each line.",
                      values=self.gpa_settings['participant_list'])

        self.page.add(label="Select Derivatives ",
                      control=control.CHECKLIST_BOX,
                      name="derivative_list",
                      type=dtype.LSTR,
                      values=['ALFF',
                              'f/ALFF',
                              'ReHo',
                              'ROI Average SCA',
                              'Dual Regression',
                              'Multiple Regression SCA',
                              'Network Centrality',
                              'VMHC'],
                      comment="Select which derivatives you would like to "
                              "include when running group analysis.\n\nWhen "
                              "including Dual Regression, make sure to "
                              "correct your P-value for the number of maps "
                              "you are comparing.\n\nWhen including "
                              "Multiple Regression SCA, you must have more "
                              "degrees of freedom (subjects) than there were "
                              "time series.",
                      size=(350,180))

        self.page.add(label="Z threshold ",
                      control=control.FLOAT_CTRL,
                      name='z_threshold',
                      type=dtype.NUM,
                      comment="Only voxels with a Z-score higher than this "
                              "value will be considered significant.",
                      values=self.gpa_settings['z_threshold'])

        self.page.add(label="Cluster Significance Threshold (P-value) ",
                      control=control.FLOAT_CTRL,
                      name='p_threshold',
                      type=dtype.NUM,
                      comment="Significance threshold (P-value) to use when "
                              "doing cluster correction for multiple "
                              "comparisons.",
                      values=self.gpa_settings['p_threshold'])

        self.page.add(label="Model Name ",
                      control=control.TEXT_BOX,
                      name="model_name",
                      type=dtype.STR,
                      comment="Specify a name for the new model. Output and "
                              "working directories for group analysis, as "
                              "well as the FLAME model files (.mat, .con, "
                              ".grp, etc.) will be labeled with this name.",
                      values=self.gpa_settings['model_name'],
                      size=(200, -1))

        self.page.add(label="Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="output_dir",
                      type=dtype.STR,
                      comment="Full path to the directory where CPAC should "
                              "place the model files (.mat, .con, .grp) and "
                              "the outputs of group analysis.",
                      values=self.gpa_settings['output_dir'])

        if gpa_settings:
            # manually re-set the preset, participant list, and derivatives
            for ctrl in self.page.get_ctrl_list():

                name = ctrl.get_name()
                if name == 'z_threshold':
                    ctrl.set_value(self.gpa_settings['z_threshold'])

                elif name == 'p_threshold':
                    ctrl.set_value(self.gpa_settings['p_threshold'])

                elif ("list" in name) and (name != "participant_list"):

                    value = self.gpa_settings[name]
                    if isinstance(value, str):
                        value = value.replace("[", "").replace("]", "")
                        if "\"" in value:
                            value = value.replace("\"", "")
                        if "'" in value:
                            value = value.replace("'", "")
                        values = value.split(",")
                    else:
                        # instead, is a list- most likely when clicking
                        # "Back" on the modelDesign_window
                        values = value
                    new_derlist = []

                    for val in values:
                        new_derlist.append(val)
                    ctrl.set_value(new_derlist)
                    ctrl.set_selection(new_derlist)

                else:
                    ctrl.set_value(self.gpa_settings[name])

        self.page.set_sizer()
        mainSizer.Add(self.window, 1, wx.EXPAND)

        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        buffer = wx.StaticText(btnPanel, label="\t\t\t\t\t\t")
        hbox.Add(buffer)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel", (220, 10),
                           wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)

        next = wx.Button(btnPanel, 3, "OK", (200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.load_next_stage, id=3)
        hbox.Add(next, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        # reminder: functions bound to buttons require arguments
        #           (self, event)
        btnPanel.SetSizer(hbox)

        mainSizer.Add(
            btnPanel, 0.5,  flag=wx.ALIGN_RIGHT | wx.RIGHT, border=20)

        self.panel.SetSizer(mainSizer)

        self.Show()

    def gather_form_data(self):
        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()

            if ("list" in name) and (name != "participant_list"):
                # for the options that are actual lists
                #   ex. derivative_list, sessions_list, etc.
                self.gpa_settings[name] = []
                for option in list(ctrl.get_selection()):
                    self.gpa_settings[name].append(option)
            else:
                self.gpa_settings[name] = str(ctrl.get_selection())

    def cancel(self, event):
        self.Close()

    def display(self, win, msg):
        wx.MessageBox(msg, "Error")
        win.SetBackgroundColour("pink")
        win.SetFocus()
        win.Refresh()
        raise ValueError

    ''' button: OK '''
    def load_next_stage(self, event):

        self.gather_form_data()

        if "Covariate" in self.gpa_settings["flame_preset"] or \
                "Unpaired" in self.gpa_settings["flame_preset"]:
            # open the next window!
            FlamePresetsTwoPheno(self.parent, self.gpa_settings)
        elif "Two-Sample Paired" in self.gpa_settings["flame_preset"] or \
                "Tripled" in self.gpa_settings["flame_preset"]:
            # open the next window!
            FlamePresetsTwoConditions(self.parent, self.gpa_settings)
        elif "One-Sample" in self.gpa_settings["flame_preset"]:
            # no additional info is needed, and we can run the preset
            # generation straight away
            print "Generating FSL FEAT/FLAME model configuration...\n"
            create_fsl_flame_preset.run(self.gpa_settings["pipeline_dir"], 
                                        self.gpa_settings["derivative_list"],
                                        self.gpa_settings["z_threshold"],
                                        self.gpa_settings["p_threshold"],
                                        "single_grp",
                                        self.gpa_settings["participant_list"],
                                        output_dir=self.gpa_settings["output_dir"],
                                        model_name=self.gpa_settings["model_name"])

            yaml_path = os.path.join(self.gpa_settings["output_dir"],
                                     self.gpa_settings["model_name"],
                                     "group_config_{0}.yml"
                                     "".format(
                                         self.gpa_settings["model_name"]))

            dialog_msg = 'Generated your FSL-FEAT preset. Check the terminal ' \
                         'window for details.\n\nGroup config file created:\n' \
                         '{0}\n\nYou can load this group configuration file into ' \
                         'the Pipelines box and either run group-level analysis ' \
                         'or edit the model (under General Settings and FSL-FEAT ' \
                         'Settings).'.format(yaml_path)
            dialog_title = 'FSL-FEAT Preset Generated'
            bld_dialog = wx.MessageDialog(self, dialog_msg, dialog_title,
                                      wx.OK | wx.ICON_INFORMATION)
            bld_dialog.ShowModal()
            bld_dialog.Destroy()

        self.Close()


class FlamePresetsTwoPheno(wx.Frame):
    def __init__(self, parent, gpa_settings):

        wx.Frame.__init__(self, parent=parent,
                          title="CPAC - FSL FEAT Presets",
                          size=(700, 275))

        self.parent = parent
        self.gpa_settings = gpa_settings

        if "pheno_file" not in self.gpa_settings.keys():
            self.gpa_settings["pheno_file"] = ""
        if "participant_id_label" not in self.gpa_settings.keys():
            self.gpa_settings["participant_id_label"] = ""
        if "covariate" not in self.gpa_settings.keys():
            self.gpa_settings["covariate"] = ""

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = wx.Panel(self)
        self.window = wx.ScrolledWindow(self.panel)

        self.page = generic_class.GenericClass(self.window,
                                               " FSL FEAT Group-Level "
                                               "Analysis - Model Presets")

        self.page.add(label="Phenotype/EV File ",
                      control=control.COMBO_BOX,
                      name="pheno_file",
                      type=dtype.STR,
                      comment="Full path to a .csv file containing EV "
                              "information for each subject.",
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

        if "Covariate" in self.gpa_settings["flame_preset"]:
            self.page.add(label="Phenotype covariate to include: ",
                          control=control.TEXT_BOX,
                          name='covariate',
                          type=dtype.STR,
                          values=self.gpa_settings['covariate'],
                          comment="For the additional covariate for the "
                                  "single group average, enter the column "
                                  "name of the covariate in the phenotype "
                                  "file provided.",
                          size=(160, -1))
        elif "Unpaired" in self.gpa_settings["flame_preset"]:
            self.page.add(label="Two groups from pheno to compare: ",
                          control=control.TEXT_BOX,
                          name='covariate',
                          type=dtype.STR,
                          values=self.gpa_settings['covariate'],
                          comment="Enter the two column names of the group "
                                  "variables in the phenotype provided, "
                                  "separated by a comma. If the two groups "
                                  "are encoded in a single column, enter the "
                                  "name of that one column.",
                          size=(320, -1))

        self.page.set_sizer()

        mainSizer.Add(self.window, 1, wx.EXPAND)

        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        buffer = wx.StaticText(btnPanel, label="\t\t\t\t\t\t")
        hbox.Add(buffer)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel", (220, 10),
                           wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "< Back", (220, 10),
                           wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.go_back, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)

        next = wx.Button(btnPanel, 3, "Generate Model", (200, -1),
                         wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.click_OK, id=3)
        hbox.Add(next, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        # reminder: functions bound to buttons require arguments
        #           (self, event)
        btnPanel.SetSizer(hbox)

        mainSizer.Add(
            btnPanel, 0.5,  flag=wx.ALIGN_RIGHT | wx.RIGHT, border=20)

        self.panel.SetSizer(mainSizer)

        # TODO
        # text blurb depending on the specific preset

        self.Show()

    def gather_form_data(self):
        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()
            self.gpa_settings[name] = str(ctrl.get_selection())
            if "covariate" in name:
                # control for any spaces in the string if there are two
                # covariate/column names listed, i.e. "male, female"
                self.gpa_settings["covariate"] = \
                    self.gpa_settings["covariate"].replace(", ", ",")

    def cancel(self, event):
        self.Close()

    def go_back(self, event):
        self.gather_form_data()
        FlamePresetsOne(self.parent, self.gpa_settings)
        self.Close()

    def testFile(self, filepath, paramName):
        try:
            fileTest = open(filepath)
            fileTest.close()
        except:
            errDlgFileTest = wx.MessageDialog(
                self, 'Error reading file - either it does not exist ' \
                      'or you do not have read access. \n\n' \
                      'Parameter: %s' % paramName,
                'File Access Error',
                wx.OK | wx.ICON_ERROR)
            errDlgFileTest.ShowModal()
            errDlgFileTest.Destroy()
            raise Exception

    def test_pheno_contents(self):
        with open(os.path.abspath(self.gpa_settings['pheno_file']), "rU") as phenoFile:
            phenoHeaderString = phenoFile.readline().rstrip('\r\n')
            if '.csv' in self.gpa_settings['pheno_file'] or \
                    '.CSV' in self.gpa_settings['pheno_file']:
                self.phenoHeaderItems = phenoHeaderString.split(',')
            elif '.tsv' in self.gpa_settings['pheno_file'] or \
                    '.TSV' in self.gpa_settings['pheno_file']:
                self.phenoHeaderItems = phenoHeaderString.split('\t')
            else:
                errSubID = wx.MessageDialog(
                    self, 'This does not seem to be a valid phenotype file.',
                    'Invalid Phenotype File',
                    wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()
                raise Exception

        if self.gpa_settings['participant_id_label'] in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(self.gpa_settings['participant_id_label'])
        else:
            errSubID = wx.MessageDialog(
                self, 'Please enter the name of the participant ID column'
                ' as it is labeled in the phenotype file.',
                'Blank/Incorrect Participant Header Input',
                wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
            raise Exception

        for col_name in self.gpa_settings["covariate"].split(","):
            if col_name not in self.phenoHeaderItems and \
                    col_name.replace(" ", "") not in self.phenoHeaderItems:
                err = wx.MessageDialog(self, "The covariate name entered "
                                             "does not exist in the pheno"
                                             "type file provided.\n\nName: "
                                             "{0}".format(col_name),
                                       'Incorrect Covariate Input',
                                       wx.OK | wx.ICON_ERROR)
                err.ShowModal()
                err.Destroy()
                raise Exception

    def substitute_derivative_names(self):
        # change the human-friendly strings of the derivative names to the
        # CPAC output directory derivative names using constants.py
        new_deriv_list = []
        for deriv_string in self.gpa_settings["derivative_list"]:
            new_deriv_list.append(substitution_map.get(deriv_string))
        self.gpa_settings["derivative_list"] = new_deriv_list

    def click_OK(self, event):

        # gather data
        self.gather_form_data()
        self.substitute_derivative_names()

        # check pheno file
        self.testFile(self.gpa_settings['pheno_file'], 'Phenotype/EV File')
        self.test_pheno_contents()

        # which preset?
        if "Covariate" in self.gpa_settings["flame_preset"]:
            preset = "single_grp_cov"
        elif "Unpaired" in self.gpa_settings["flame_preset"]:
            preset = "unpaired_two"

        # generate the preset files
        print "Generating FSL FEAT/FLAME model configuration...\n"
        create_fsl_flame_preset.run(self.gpa_settings["pipeline_dir"],
                                    self.gpa_settings["derivative_list"],
                                    self.gpa_settings["z_threshold"],
                                    self.gpa_settings["p_threshold"],
                                    preset,
                                    self.gpa_settings["participant_list"],
                                    self.gpa_settings["pheno_file"],
                                    self.gpa_settings["participant_id_label"],
                                    output_dir=self.gpa_settings[
                                        "output_dir"],
                                    model_name=self.gpa_settings[
                                        "model_name"],
                                    covariate=self.gpa_settings["covariate"])

        yaml_path = os.path.join(self.gpa_settings["output_dir"],
                                 self.gpa_settings["model_name"],
                                 "group_config_{0}.yml"
                                 "".format(self.gpa_settings["model_name"]))

        dialog_msg = 'Generated your FSL-FEAT preset. Check the terminal ' \
                     'window for details.\n\nGroup config file created:\n' \
                     '{0}\n\nYou can load this group configuration file into ' \
                     'the Pipelines box and either run group-level analysis ' \
                     'or edit the model (under General Settings and FSL-FEAT ' \
                     'Settings).'.format(yaml_path)
        dialog_title = 'FSL-FEAT Preset Generated'
        bld_dialog = wx.MessageDialog(self, dialog_msg, dialog_title,
                                      wx.OK | wx.ICON_INFORMATION)
        bld_dialog.ShowModal()
        bld_dialog.Destroy()

        self.Close()


class FlamePresetsTwoConditions(wx.Frame):
    def __init__(self, parent, gpa_settings):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - FSL FEAT Presets",
            size=(550, 325))

        self.parent = parent
        self.gpa_settings = gpa_settings

        if "condition_type" not in self.gpa_settings.keys():
            self.gpa_settings["condition_type"] = "Sessions"
        if "conditions_list" not in self.gpa_settings.keys():
            self.gpa_settings["conditions_list"] = []

        if "Two-Sample Paired" in self.gpa_settings["flame_preset"]:
            num_condition = "two"
        elif "Tripled" in self.gpa_settings["flame_preset"]:
            num_condition = "three"

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = wx.Panel(self)
        self.window = wx.ScrolledWindow(self.panel)

        self.page = generic_class.GenericClass(self.window,
                                               " Group-Level Analysis - FSL "
                                               "FEAT Presets")

        # TODO
        # text blurb depending on the specific preset

        self.page.add(label="Conditions: Sessions or Series? ",
                      control=control.CHOICE_BOX,
                      name="condition_type",
                      type=dtype.LSTR,
                      comment="Choose whether the {0} paired conditions you "
                              "wish to compare using the {1} are Sessions or "
                              "Series/Scans.\n\nFor example, if each "
                              "participant has {0} conditions each, are they "
                              "separated by session, or by multiple "
                              "functional series or scans within each a "
                              "single session each?".format(num_condition,
                                                            self.gpa_settings["flame_preset"]),
                      values=["Sessions", "Series"])

        self.page.add(label='Session or Series/Scan IDs: ',
                      control=control.LISTBOX_COMBO,
                      name='conditions_list',
                      type=dtype.LSTR,
                      values=self.gpa_settings['conditions_list'],
                      comment="Enter the {0} session or series/scan ID names "
                              "you wish to compare.\n\nIf they are sessions, "
                              "they can be found under 'unique_id' in the "
                              "individual-level data configuration YAML file "
                              "for the data in question, or in the "
                              "'participant_session' IDs that you would see "
                              "in the individual-level analysis output dir"
                              "ectory (ex. '3005_1' would be participant "
                              "3005, session 1).\n\nIf they are series/"
                              "scans, they will be the functional scan "
                              "names, which can be found in the data config"
                              "uration YAML file nested under "
                              "'func'.".format(num_condition),
                      size=(200, 100),
                      combo_type=7)

        self.page.set_sizer()

        mainSizer.Add(self.window, 1, wx.EXPAND)

        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        buffer = wx.StaticText(btnPanel, label="\t\t\t\t\t\t")
        hbox.Add(buffer)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel", (220, 10),
                           wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "< Back", (220, 10),
                           wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.go_back, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)

        next = wx.Button(btnPanel, 3, "Generate Model", (200, -1),
                         wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.click_OK, id=3)
        hbox.Add(next, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        # reminder: functions bound to buttons require arguments
        #           (self, event)
        btnPanel.SetSizer(hbox)

        mainSizer.Add(
            btnPanel, 0.5,  flag=wx.ALIGN_RIGHT | wx.RIGHT, border=20)

        self.panel.SetSizer(mainSizer)

        # TODO
        # text blurb depending on the specific preset

        for ctrl in self.page.get_ctrl_list():

            name = ctrl.get_name()
            value = self.gpa_settings[name]

            if "conditions_list" in name:
                if isinstance(value, str):
                    value = value.replace("[", "").replace("]", "")
                    if "\"" in value:
                        value = value.replace("\"", "")
                    if "'" in value:
                        value = value.replace("'", "")
                    values = value.split(",")
                else:
                    # instead, is a list- most likely when clicking
                    # "Back" on the modelDesign_window
                    values = value

                new_derlist = []

                for val in values:
                    new_derlist.append(val)

                ctrl.set_value(new_derlist)
            elif "condition_type" in name:
                ctrl.set_value(value)

        self.Show()

    def gather_form_data(self):
        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()
            self.gpa_settings[name] = ctrl.get_selection()

    def cancel(self, event):
        self.Close()

    def go_back(self, event):
        self.gather_form_data()
        FlamePresetsOne(self.parent, self.gpa_settings)
        self.Close()

    def substitute_derivative_names(self):
        # change the human-friendly strings of the derivative names to the
        # CPAC output directory derivative names using constants.py
        new_deriv_list = []
        for deriv_string in self.gpa_settings["derivative_list"]:
            new_deriv_list.append(substitution_map.get(deriv_string))
        self.gpa_settings["derivative_list"] = new_deriv_list

    def click_OK(self, event):
        # gather data
        self.gather_form_data()
        self.substitute_derivative_names()

        # which preset?
        if "Two-Sample Paired" in self.gpa_settings["flame_preset"]:
            preset = "paired_two"
        elif "Tripled" in self.gpa_settings["flame_preset"]:
            preset = "tripled_two"

        # generate the preset files
        print "Generating FSL FEAT/FLAME model configuration...\n"
        create_fsl_flame_preset.run(self.gpa_settings["pipeline_dir"],
                                    self.gpa_settings["derivative_list"],
                                    self.gpa_settings["z_threshold"],
                                    self.gpa_settings["p_threshold"],
                                    preset,
                                    self.gpa_settings["participant_list"],
                                    output_dir=self.gpa_settings[
                                        "output_dir"],
                                    model_name=self.gpa_settings[
                                        "model_name"],
                                    covariate=self.gpa_settings[
                                        "conditions_list"],
                                    condition_type=self.gpa_settings[
                                        "condition_type"])

        yaml_path = os.path.join(self.gpa_settings["output_dir"],
                                 self.gpa_settings["model_name"],
                                 "group_config_{0}.yml"
                                 "".format(self.gpa_settings["model_name"]))

        dialog_msg = 'Generated your FSL-FEAT preset. Check the terminal ' \
                     'window for details.\n\nGroup config file created:\n' \
                     '{0}\n\nYou can load this group configuration file into ' \
                     'the Pipelines box and either run group-level analysis ' \
                     'or edit the model (under General Settings and FSL-FEAT ' \
                     'Settings).'.format(yaml_path)
        dialog_title = 'FSL-FEAT Preset Generated'
        bld_dialog = wx.MessageDialog(self, dialog_msg, dialog_title,
                                      wx.OK | wx.ICON_INFORMATION)
        bld_dialog.ShowModal()
        bld_dialog.Destroy()

        self.Close()
