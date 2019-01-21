# CPAC/GUI/interface/windows/dataconfig_window.py
#
#

'''
This module starts the data configuration GUI for building a subject list
'''

# Import packages
import wx
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
import os
import yaml
import CPAC
import pkg_resources as p
import sys

# Init variables
ID_RUN_EXT = 11
ID_RUN_MEXT = 12


class DataConfig(wx.Frame):

    # Init method
    def __init__(self, parent):

        wx.Frame.__init__(self, parent, title="CPAC - Data Configuration "
                                              "Setup",
                          size=(1040, 620))
        
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        
        self.panel = wx.Panel(self)
        
        self.window = wx.ScrolledWindow(self.panel)
        
        self.page = GenericClass(self.window, "Data Configuration Setup")

        self.page.add(label="Data format ",
                      control=control.CHOICE_BOX,
                      name='dataFormat',
                      type=dtype.BOOL,
                      comment="Select if data is organized using BIDS "
                              "standard or a custom format.",
                      values=["BIDS", "Custom"],
                      wkf_switch=True)

        self.page.add(label="[BIDS only] BIDS Base Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="bidsBaseDir",
                      type=dtype.STR,
                      comment="Base directory of BIDS-organized data.\n"
                              "BIDS Data Format only.\n\nThis should be the "
                              "path to the overarching directory containing "
                              "the entire dataset.",
                      values="")

        self.page.add(label="[Custom only] Anatomical File Path Template ",
                      control=control.TEXT_BOX,
                      name="anatomicalTemplate",
                      type=dtype.STR,
                      comment="File Path Template for Anatomical Files\n"
                              "Custom Data Format only.\n\n"
                              "Place tags for the appropriate data directory "
                              "levels with the tags {site}, {participant}, "
                              "and {session}. Only {participant} is "
                              "required.\n\nExamples:\n"
                              "/data/{site}/{participant}/{session}/anat"
                              "/mprage.nii.gz\n"
                              "/data/{site}/{participant}/anat.nii.gz\n\n"
                              "See the User Guide for more detailed "
                              "instructions.",
                      values="",
                      style=wx.EXPAND | wx.ALL,
                      size=(532,-1))
        
        self.page.add(label="[Custom only] Functional File Path Template ",
                 control=control.TEXT_BOX,
                 name="functionalTemplate",
                 type=dtype.STR,
                 comment = "File Path Template for Functional Files\n"
                           "Custom Data Format only.\n\n"
                           "Place tags for the appropriate data directory "
                           "levels with the tags {site}, {participant}, "
                           "{session}, and {series}. Only {participant} is "
                           "required.\n\nExamples:\n"
                           "/data/{site}/{participant}/{session}/func/"
                           "{series}_bold.nii.gz\n"
                           "/data/{site}/{participant}/{series}/func.nii.gz"
                           "\n\n"
                           "See the User Guide for more detailed "
                           "instructions.",
                 values ="",
                 style= wx.EXPAND | wx.ALL,
                 size = (532,-1))

        self.page.add(label = "Save Config Files Here: ",
                      control = control.DIR_COMBO_BOX,
                      name = "outputSubjectListLocation",
                      type = dtype.STR,
                      comment="Directory where CPAC should place data "
                              "configuration files.",
                      values="")

        self.page.add(label = "Participant List Name ",
                      control = control.TEXT_BOX,
                      name = "subjectListName",
                      type = dtype.STR,
                      comment = "A label to be appended to the generated " \
                                "participant list files.",
                      values = "",
                      style= wx.EXPAND | wx.ALL,
                      size = (300,-1))

        self.page.add(label="(Optional) Which Anatomical Scan? ",
                      control=control.TEXT_BOX,
                      name="anatomical_scan",
                      type=dtype.STR,
                      comment="Scan/Run ID for the Anatomical Scan\n\n"
                              "Sometimes, there are multiple anatomical "
                              "scans for each participant in a dataset.\n\n"
                              "If this is the case, you can choose which "
                              "anatomical scan to use for this participant "
                              "by entering the identifier that makes the "
                              "scan unique.\n\nExamples:\n\nBIDS dataset\n"
                              "../anat/sub-001_run-1_T1w.nii.gz\n"
                              "../anat/sub-001_run-2_T1w.nii.gz\n"
                              "Pick the second with 'run-2'.\n\n"
                              "Custom dataset\n"
                              "Example use case: let's say most anatomicals "
                              "in your dataset are '../mprage.nii.gz', but "
                              "some participants only have '../anat1.nii.gz' "
                              "and '../anat2.nii.gz'. You want the "
                              "mprage.nii.gz files included, but only the "
                              "anat2.nii.gz in the others.\n\nPlace a "
                              "wildcard (*) in the anatomical filepath "
                              "template above (../*.nii.gz), then enter "
                              "'anat2' in this field to 'break the tie' for "
                              "participants that have the 'anat1' and "
                              "'anat2' scans.",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532,-1))

        # Add AWS credentials path
        self.page.add(label="(Optional) AWS credentials file ",
                 control=control.COMBO_BOX,
                 name='awsCredentialsFile',
                 type=dtype.COMBO,
                 comment='Required if downloading data from a non-public S3 '\
                         'bucket on Amazon Web Services instead of using '\
                         'local files.',
                 values='None')

        self.page.add(label="(Optional) [Custom only] Scan Parameters File ",
                 control=control.COMBO_BOX,
                 name="scanParametersCSV",
                 type=dtype.COMBO,
                 comment="For Slice Timing Correction.\nCustom Data Format "
                         "only.\n\n"
                         "Path to a .csv file (if not using BIDS-format "
                         "JSON files) containing information about scan "
                         "acquisition parameters.\n\nFor instructions on "
                         "how to create this file, see the User Guide.\n\n"
                         "If 'None' is specified, CPAC will look for scan "
                         "parameters information provided in the pipeline "
                         "configuration file.",
                 values="None")

        self.page.add(label="(Optional) Brain Mask File Path Template ",
                      control=control.TEXT_BOX,
                      name="brain_mask_template",
                      type=dtype.STR,
                      comment="File Path Template for brain mask files.\n"
                              "For anatomical skull-stripping.\nNote: can "
                              "use this for BIDS data directories as "
                              "well, as the specification for anatomical "
                              "brain masks is still subject to change.\n\n"
                              "Place tags for the appropriate data "
                              "directory levels with the tags {site}, "
                              "{participant}, and {session}. Only "
                              "{participant} is required.\n\nExamples:\n"
                              "/data/{site}/{participant}/{session}/fmap/"
                              "phase.nii.gz\n/data/{site}/{participant}/"
                              "{session}/{participant}_{session}_"
                              "brain-mask.nii.gz",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532,-1))

        self.page.add(label="(Optional) [Custom only] Field Map Phase File "
                            "Path Template ",
                      control=control.TEXT_BOX,
                      name="fieldMapPhase",
                      type=dtype.STR,
                      comment="File Path Template for Field Map Phase "
                              "files\nFor field-map based distortion "
                              "correction.\nCustom Data Format only.\n\n"
                              "Place tags for the appropriate data "
                              "directory levels with the tags {site}, "
                              "{participant}, and {session}. Only "
                              "{participant} is required.\n\nExamples:\n"
                              "/data/{site}/{participant}/{session}/fmap/"
                              "phase.nii.gz\n/data/{site}/{participant}/"
                              "{session}/{participant}_{session}_"
                              "phase.nii.gz",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532,-1))

        self.page.add(label="(Optional) [Custom only] Field Map Magnitude "
                            "File Path Template ",
                      control=control.TEXT_BOX,
                      name="fieldMapMagnitude",
                      type=dtype.STR,
                      comment="File Path Template for Field Map Magnitude "
                              "files\nFor field-map based distortion "
                              "correction.\nCustom Data Format only.\n\n"
                              "Place tags for the appropriate data "
                              "directory levels with the tags {site}, "
                              "{participant}, and {session}. Only"
                              "{participant} is required.\n\nExamples:\n"
                              "/data/{site}/{participant}/{session}/fmap/"
                              "magnitude.nii.gz\n/data/{site}/{participant}/"
                              "{session}/{participant}_{session}_"
                              "magnitude.nii.gz",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532,-1))

        self.page.add(label="(Optional) Include: Subjects ",
                 control=control.COMBO_BOX, 
                 name = "subjectList", 
                 type =dtype.COMBO,
                 comment = "Include only a sub-set of the participants "
                           "present in the folders defined above.\n\nList "
                           "participants in this box (e.g., sub101, sub102) "
                           "or provide the path to a text file with one "
                           "participant ID on each line.\n\nIf 'None' is "
                           "specified, CPAC will include all participants."
                           "\n\nNote: for BIDS datasets, only list the ID "
                           "value that comes after the label tag. For "
                           "example, list '5001' to include what is labeled "
                           "as 'sub-5001'.",
                 values = "None")
        
        self.page.add(label="(Optional) Exclude: Subjects ",
                 control=control.COMBO_BOX, 
                 name = "exclusionSubjectList", 
                 type = dtype.COMBO, 
                 comment = "Exclude a sub-set of the participants present in "
                           "the folders defined above.\n\nList participants "
                           "in this box (e.g., sub101, sub102) or provide "
                           "the path to a text file with one participant ID "
                           "on each line.\n\nIf 'None' is specified, CPAC "
                           "will not exclude any participants.",
                 values = "None")
        
        self.page.add(label= "(Optional) Include: Sites ",
                 control =control.COMBO_BOX,
                 name = "siteList",
                 type =dtype.COMBO,
                 comment = "Include only a sub-set of the sites present in "
                           "the folders defined above.\n\nList sites in this "
                           "box (e.g., NYU, UCLA) or provide the path to a "
                           "text file with one site name on each line.\n\n"
                           "If 'None' is specified, CPAC will include all "
                           "sites.",
                 values ="None",
                 style= wx.EXPAND | wx.ALL,
                 size = (532,-1))

        self.page.add(label= "(Optional) Exclude: Sites ",
                 control =control.COMBO_BOX,
                 name = "exclusionSiteList",
                 type =dtype.COMBO,
                 comment = "Exclude a sub-set of the sites present in the "
                           "folders defined above.\n\nList sites in this box "
                           "(e.g., NYU, UCLA) or provide the path to a text "
                           "file with one site name on each line.\n\n"
                           "If 'None' is specified, CPAC will include all "
                           "sites.",
                 values ="None",
                 style= wx.EXPAND | wx.ALL,
                 size = (532,-1))

        self.page.add(label="(Optional) Include: Sessions ",
                      control=control.COMBO_BOX,
                      name="sessionList",
                      type=dtype.COMBO,
                      comment="Include only a sub-set of the sessions "
                              "present in the folders defined above.\n\nList "
                              "sessions in this box (e.g., session-1, "
                              "session-2) or provide the path to a text file "
                              "with one session name on each line.\n\nIf "
                              "'None' is specified, CPAC will include all "
                              "sessions.",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532, -1))

        self.page.add(label="(Optional) Exclude: Sessions ",
                      control=control.COMBO_BOX,
                      name="exclusionSessionList",
                      type=dtype.COMBO,
                      comment="Exclude a sub-set of the sessions present in "
                              "the folders defined above.\n\nList sessions "
                              "in this box (e.g., session-1, session-2) or "
                              "provide the path to a text file with one "
                              "session name on each line.\n\nIf 'None' is "
                              "specified, CPAC will include all sessions.",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532, -1))

        self.page.add(label="(Optional) Include: Series ",
                      control=control.COMBO_BOX,
                      name="scanList",
                      type=dtype.COMBO,
                      comment="Include only a sub-set of the series present "
                              "in the folders defined above.\n\nList series "
                              "in this box (e.g., func-1, func-2) or provide "
                              "the path to a text file with one series name "
                              "on each line.\n\nIf 'None' is specified, CPAC "
                              "will include all series.",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532, -1))

        self.page.add(label="(Optional) Exclude: Series ",
                      control=control.COMBO_BOX,
                      name="exclusionScanList",
                      type=dtype.COMBO,
                      comment="Exclude a sub-set of the series present "
                              "in the folders defined above.\n\nList series "
                              "in this box (e.g., func-1, func-2) or provide "
                              "the path to a text file with one series name "
                              "on each line.\n\nIf 'None' is specified, CPAC "
                              "will include all series.",
                      values="None",
                      style=wx.EXPAND | wx.ALL,
                      size=(532, -1))

        self.page.set_sizer()
         
        mainSizer.Add(self.window, 1, wx.EXPAND)
        
        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        buffer2 = wx.StaticText(btnPanel, label = "\t")
        hbox.Add(buffer2)
    
        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel",(220,10), wx.DefaultSize, 0)

        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add( cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)
        
        load = wx.Button(btnPanel, wx.ID_ADD, "Load Preset", (280,10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.load, id=wx.ID_ADD)
        hbox.Add(load, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)
        
        save = wx.Button(btnPanel, wx.ID_SAVE, "Save Preset", (280,10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, lambda event: self.save(event,'save'), id=wx.ID_SAVE)
        hbox.Add(save, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        run_ext = wx.Button(btnPanel, ID_RUN_EXT, "Generate Data Config", (280,10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, lambda event: self.save(event,'run'), id=ID_RUN_EXT)
        hbox.Add( run_ext, 1, flag=wx.LEFT | wx.BOTTOM, border=10)
    
        btnPanel.SetSizer(hbox)
        
        mainSizer.Add(btnPanel, 0.5,  flag=wx.ALIGN_RIGHT|wx.RIGHT, border=20)
        
        self.panel.SetSizer(mainSizer)
        
        self.Show()
        
    def cancel(self, event):
        self.Close()

    # Generate the subject list from config
    def run(self, config):
        
        # Import packages
        import CPAC

        # Try to build subject list from config
        try:
            # Load in configuration file
            config_map = yaml.load(open(config, 'r'))

            # Extract arguments for supplementary files
            sublist_outdir = config_map.get('outputSubjectListLocation')
            sublist_name = config_map.get('subjectListName')

            # Make backwards-compatible
            if isinstance(sublist_name, list):
                sublist_name = sublist_name[0]

            sublist_name = "data_config_{0}.yml".format(sublist_name)

            # Get subject list output path
            out_location = os.path.join(sublist_outdir, sublist_name)

            # Build the subject list from the data config
            CPAC.utils.build_data_config.run(config)

            # check GUI's data config list dialog box for duplicate names
            while True:
                parent = self.Parent
                map = parent.get_sublist_map()
                map[sublist_name] = out_location
                parent.listbox2.Append(sublist_name)
                ret = 1
                break

                # TODO: fix this- should warn user and force a decision if
                # TODO: another data config in the dialog box already exists
                # TODO: with the same name
                #if map.get(sublist_name) == None:
                #    map[sublist_name] = out_location
                #    parent.listbox2.Append(sublist_name)
                #    ret = 1
                #    break
                #else:
                #    dlg3 = wx.MessageDialog(self, 'A data config with '
                #                                  'this name already '
                #                                  'exists.', 'Error!',
                #                            wx.OK | wx.ICON_ERROR)
                #    dlg3.ShowModal()
                #    dlg3.Destroy()

            return ret

        # Import error if CPAC not available
        except ImportError as exc:
            wx.MessageBox("Error importing CPAC. Unable to run extract data "
                          "tool.", "Error")
            print "Error importing CPAC"
            print exc
            return -1
        # Problem reading in data from disk
        except IOError as exc:
            print "Error loading data config file", exc
            return -1

    def save(self, event, flag):
        # Save data config
        
        config_list = []
        config_dict = {}

        def display(win, msg):
            wx.MessageBox(msg, "Error")
            win.SetBackgroundColour("pink")
            win.SetFocus()
            win.Refresh()
            raise ValueError

        key_order = ['dataFormat', 'bidsBaseDir', 'anatomicalTemplate',
                     'functionalTemplate', 'outputSubjectListLocation',
                     'subjectListName', 'anatomical_scan',
                     'awsCredentialsFile', 'scanParametersCSV',
                     'brain_mask_template', 'fieldMapPhase',
                     'fieldMapMagnitude', 'subjectList',
                     'exclusionSubjectList', 'siteList', 'exclusionSiteList',
                     'sessionList', 'exclusionSessionList', 'scanList',
                     'exclusionScanList']

        path_fields = ['scanParametersCSV', 'awsCredentialsFile',
                       'outputSubjectListLocation']

        try:
            for ctrl in self.page.get_ctrl_list():
                win = ctrl.get_ctrl()
                value = str(ctrl.get_selection())
                value = value.strip()
                name = ctrl.get_name()
                dtype = ctrl.get_datatype()
                help = ctrl.get_help()

                if name == 'subjectListName':
                    subject_list_name = value

                if len(value) == 0:
                    if name != "bidsBaseDir" and \
                                    name != "anatomicalTemplate" and \
                                    name != "functionalTemplate":
                        display(win,"%s field must contain some text!"%ctrl.get_name())

                if name in path_fields:
                    if "s3://" not in value and len(value) > 0 and \
                                    value is not None and \
                                    "None" not in value and \
                                    "none" not in value:
                        if not os.path.exists(value):
                            display(win, "%s field contains incorrect path. "
                                         "Please update the path!"
                                    % ctrl.get_name())
         
                config_list.append((name, value, dtype))
                config_dict[name] = (value, dtype, help)

            # some final checks
            if "BIDS" in config_dict["dataFormat"][0]:
                if len(config_dict["anatomicalTemplate"][0]) > 0 or \
                                len(config_dict["functionalTemplate"][0]) > 0:
                    if ("None" in config_dict["anatomicalTemplate"][0] or "none" in config_dict["anatomicalTemplate"][0]) and \
                            ("None" in config_dict["functionalTemplate"][0] or "none" in config_dict["functionalTemplate"][0]):
                        pass
                    else:
                        err = wx.MessageDialog(self, "Custom filepath template "
                                                     "provided, but data format "
                                                     "is set to BIDS instead of "
                                                     "Custom.",
                                                     'Error!',
                                                     wx.OK | wx.ICON_ERROR)
                        err.ShowModal()
                        err.Destroy()
                        return

                elif "s3://" not in config_dict["bidsBaseDir"][0] and \
                        not os.path.exists(config_dict["bidsBaseDir"][0]):
                    err = wx.MessageDialog(self, "Data format is set to "
                                                 "BIDS, but no BIDS base "
                                                 "directory is set, or the "
                                                 "BIDS directory does not "
                                                 "exist.",
                                                 'Error!',
                                                 wx.OK | wx.ICON_ERROR)
                    err.ShowModal()
                    err.Destroy()
                    return

            elif "Custom" in config_dict["dataFormat"][0]:
                if "/" in str(config_dict["bidsBaseDir"]):
                    err = wx.MessageDialog(self, "BIDS base directory "
                                                 "provided, but data format "
                                                 "is set to Custom instead "
                                                 "of BIDS.",
                                                 'Error!',
                                                 wx.OK | wx.ICON_ERROR)
                    err.ShowModal()
                    err.Destroy()
                    return

                if len(config_dict["anatomicalTemplate"][0]) == 0:
                    err = wx.MessageDialog(self, "Custom data format "
                                                 "selected, but no custom "
                                                 "anatomical filepath "
                                                 "template provided.",
                                                 'Error!',
                                                 wx.OK | wx.ICON_ERROR)
                    err.ShowModal()
                    err.Destroy()
                    return
                
        except Exception, e:
            errdlg = wx.MessageDialog(self, "Could not save your "
                                            "participant list information."
                                            "\n\n%s" % e,
                               'Error!',
                           wx.OK | wx.ICON_ERROR)
            errdlg.ShowModal()
            errdlg.Destroy()
            print e
            return
            
        else:
            dlg = wx.FileDialog(
                self, message="Save file as ...", 
                defaultDir=str(config_dict["outputSubjectListLocation"]),
                defaultFile="data_settings_{0}.yaml".format(subject_list_name),
                wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml", 
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                dlg.Destroy()

                with open(path, "wt") as f:
                    f.write(
                        "# CPAC Data Settings File\n# Version {0}"
                        "\n".format(CPAC.__version__))
                    f.write(
                        "#\n# http://fcp-indi.github.io for more info.\n#\n"
                        "# Use this file to generate the data configuration "
                        "(participant list) YAML file by loading it via the "
                        "'Load Preset' button in the Data Configuration "
                        "Builder UI, or via command line by providing it to "
                        "the cpac_data_config_setup.py script with the "
                        "--data_settings_file input flag.\n\n\n")

                    for key in key_order:
                        value = config_dict[key][0]
                        help = config_dict[key][2]

                        # if "/" in val or "%s" in val or 'None' in val or \
                        #         key == 'subjectListName':
                        #     value = val
                        # else:
                        #     value =[item.strip() for item in val.split(',')]

                        help = help.replace("\n", "\n# ")

                        f.write("# {0}\n".format(help))
                        f.write("{0}: {1}\n\n\n".format(key, value))

                print "\nSaving data settings file:\n{0}\n".format(path)
                
                if flag == 'run':
                    if self.run(path) > 0:
                        self.Close()

    # Load in data configuration file
    def load(self, event):

        dlg = wx.FileDialog(
        self, message="Choose the config yaml file",
            defaultDir=os.getcwd(), 
            defaultFile="",
            wildcard= "YAML files(*.yaml, *.yml)|*.yaml;*.yml",
            style=wx.OPEN | wx.CHANGE_DIR)
        # Once user click's OK
        if dlg.ShowModal() == wx.ID_OK:
            # Try and load in the data config file to GUI
            try:
                path = dlg.GetPath()
                # Try and load in file contents
                try:
                    config_map = yaml.load(open(os.path.realpath(path),'r'))
                # Otherwise, report error
                except IOError as exc:
                    err_msg = 'File %s does not exist. Check and try again. '\
                              'Error:\n%s' %(path, exc)
                    raise Exception(err_msg)
                except Exception as exc:
                    err_msg = 'Unable to load in the specified file: %s'\
                              'Error:\n%s' %(path, exc)
                    raise Exception(err_msg)

                # If it's a dictionary, check it has anat template key
                if type(config_map) == dict:
                    if not config_map.has_key('anatomicalTemplate'):
                        err_msg = 'File is not a data settings '\
                                  'file. It might be a pipeline '\
                                  'configuration file.'
                        raise Exception(err_msg)
                # It didn't load in as a dictionary, report error
                else:
                    err_msg = 'File is not a data settings '\
                              'file. It might be a data config/participant ' \
                              'list file.'
                    raise Exception(err_msg)

                # Populate GUI fields
                for ctrl in self.page.get_ctrl_list():
                    name = ctrl.get_name()
                    value = config_map.get(name)
                    dtype = ctrl.get_datatype()

                    if not value:
                        val = "None"
                    elif isinstance(value, list):
                        val = None
                        for v in value:
                            if val:
                                val = val + ',' + str(v)
                            else:
                                val = str(v)
                    else:
                        val = value

                    if "None" in val or "none" in val:
                        ctrl.set_value("None ")
                    else:
                        ctrl.set_value(str(val))

            # There was an error loading parameters, report it
            except Exception as exc:
                err_msg = 'CPAC could not load your participant list ' \
                          'information. Check the formatting of your ' \
                          'data_settings YAML file.'\
                          '\n\nIssue info:\n%s' % exc
                errdlg = wx.MessageDialog(self, err_msg, 'Error!',
                                          wx.OK | wx.ICON_ERROR)
                errdlg.ShowModal()
                errdlg.Destroy()

            # Close dialog
            dlg.Destroy()
