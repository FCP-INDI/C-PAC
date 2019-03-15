import wx
from CPAC.GUI.interface.utils.constants import substitution_map
import pkg_resources as p

from CPAC.GUI.interface.pages import (
    AnatomicalPreprocessing,
    Segmentation, Registration, FunctionalPreProcessing,
    SkullStripProcessing, SkullStripOptions, AFNI_options, BET_options,
    Scrubbing, AnatToFuncRegistration, FuncToMNIRegistration,
    VMHC, VMHCSettings, ReHo, ReHoSettings,
    SCA, SCASettings,
    Settings, ComputerSettings, DirectorySettings,
    Nuisance, NuisanceRegression, MedianAngleCorrection,
    CentralitySettings, Centrality,
    ALFF, ALFFSettings,
    AfterWarping, AfterWarpingOptions,
    FilteringSettings,
    TimeSeries, EPI_DistCorr, ROITimeseries,
    GroupAnalysis, GeneralGA, GPASettings, MDMRSettings, ISCSettings, RandomiseSettings,
    TimeSeriesOptions, BASCSettings,
    AROMA_ICA, AromaSettings
)

from ..utils.constants import dtype as data_types

ID_SUBMIT = 6


def gen_checkboxgrid_config_string(label, value):

    # example inputs
    #     label (string): "tsa_roi_paths" (the name of the input control)
    #     value (dictionary):
    #         {'/path/to/roi.nii.gz': ['Voxel', 'Avg']}

    string = ""

    if len(value.keys()) == 0:
        string = string + label + ": None\n"
        return string
    else:
        string = string + label + ":\n"

    flag = 0                  

    for entry in value.keys():

        # each "entry" is a filepath to an ROI .nii.gz file
        # each "value[entry]" is a list of selected analysis types (strings)

        if flag == 0:
            string = string + "  - "
            flag = 1
        else:
            string = string + "    "

        string = string + entry + ": "

        firstnum = 0

        selection_string = str(value[entry])
        selection_string = selection_string.replace("'", "")
        selection_string = selection_string.replace("[", "")
        selection_string = selection_string.replace("]", "")

        string = string + selection_string

        string = string + "\n"

    return string


class Mybook(wx.Treebook):

    def __init__(self, parent, ind=True):
        wx.Treebook.__init__(self, parent, wx.ID_ANY, style=
                             wx.BK_DEFAULT)
        self.page_list = []

        # create the page windows as children of the notebook

        if ind:
            page1 = Settings(self)
            page2 = ComputerSettings(self)
            page3 = DirectorySettings(self)

            page4 = SkullStripProcessing(self)
            page5 = SkullStripOptions(self)
            page6 = AFNI_options(self)
            page7 = BET_options(self)

            page8 = AnatomicalPreprocessing(self)
            page9 = Registration(self, 1)
            page10 = Segmentation(self, 2)

            page11 = FunctionalPreProcessing(self)
            page12 = TimeSeriesOptions(self)
            page13 = EPI_DistCorr(self)
            page14 = AnatToFuncRegistration(self, 5)
            page15 = FuncToMNIRegistration(self, 6)

            page16= Nuisance(self)
            ica_aroma = AromaSettings(self)
            page17= NuisanceRegression(self, 7)
            page18= MedianAngleCorrection(self, 8)

            page19 = FilteringSettings(self, 9)

            page20 = TimeSeries(self)
            page22 = ROITimeseries(self)

            page27 = SCA(self)
            page28 = SCASettings(self)

            page31 = VMHC(self)
            page32 = VMHCSettings(self)

            page33 = ALFF(self)
            page34 = ALFFSettings(self)

            page35 = ReHo(self)
            page36 = ReHoSettings(self)

            page37 = Centrality(self)
            page38 = CentralitySettings(self)

            page39 = AfterWarping(self)
            page40 = AfterWarpingOptions(self)

            # add the pages to the notebook with the label to show on the tab
            self.AddPage(page1, "Environment Setup", wx.ID_ANY)
            self.AddSubPage(page2, "Computer Settings", wx.ID_ANY)
            self.AddSubPage(page3, "Output Settings", wx.ID_ANY)

            self.AddPage(page4, "Skull-Strip", wx.ID_ANY)
            self.AddSubPage(page5, "Skull-Stripping ", wx.ID_ANY)
            self.AddSubPage(page6, "3dSkullStrip options", wx.ID_ANY)
            self.AddSubPage(page7, "FSL BET options", wx.ID_ANY)

            self.AddPage(page8, "Anatomical Preprocessing", wx.ID_ANY)
            self.AddSubPage(page9, "Anatomical Registration", wx.ID_ANY)
            self.AddSubPage(page10, "Tissue Segmentation", wx.ID_ANY)
            self.AddPage(page11, "Functional Preprocessing", wx.ID_ANY)
            self.AddSubPage(page12, "Time Series Options", wx.ID_ANY)
            self.AddSubPage(page13, "Distortion Correction", wx.ID_ANY)
            self.AddSubPage(page14, "Functional to Anatomical Registration",
                            wx.ID_ANY)
            self.AddSubPage(page15, "Functional to MNI Registration",
                            wx.ID_ANY)

            self.AddPage(page16, "Nuisance", wx.ID_ANY)
            self.AddSubPage(ica_aroma, "ICA-AROMA De-noising", wx.ID_ANY)
            self.AddSubPage(page17, "Nuisance Regression", wx.ID_ANY)
            self.AddSubPage(page18, "Median Angle Correction", wx.ID_ANY)

            self.AddPage(page20, "Time Series Extraction (TSE)", wx.ID_ANY)
            self.AddSubPage(page22, "Region-of-Interest TSE Options",
                            wx.ID_ANY)

            self.AddPage(page27, "Seed-based Correlation Analysis (SCA)",
                         wx.ID_ANY)
            self.AddSubPage(page28, "SCA Options", wx.ID_ANY)

            self.AddPage(page31, "Voxel-mirrored Homotopic Connectivity",
                         wx.ID_ANY)
            self.AddSubPage(page32, "VMHC Settings", wx.ID_ANY)

            self.AddPage(page33, "ALFF and f/ALFF", wx.ID_ANY)
            self.AddSubPage(page34, "ALFF and f/ALFF Options", wx.ID_ANY)

            self.AddPage(page35, "Regional Homogeneity (ReHo)", wx.ID_ANY)
            self.AddSubPage(page36, "ReHo Options", wx.ID_ANY)

            self.AddPage(page37, "Network Centrality", wx.ID_ANY)
            self.AddSubPage(page38, "Network Centrality Options", wx.ID_ANY)

            self.AddPage(page39, "After Warping", wx.ID_ANY)
            self.AddSubPage(page40, "After Warping Options", wx.ID_ANY)

        else:
            page45 = GroupAnalysis(self)
            page46 = GeneralGA(self)
            page47 = GPASettings(self)
            page48 = RandomiseSettings(self)
            page49 = BASCSettings(self)
            page50 = MDMRSettings(self)
            page51 = ISCSettings(self)

            self.AddPage(page45, "Group Analysis Settings", wx.ID_ANY)
            self.AddSubPage(page46, "General Settings", wx.ID_ANY)
            self.AddSubPage(page47, "FSL FEAT Settings", wx.ID_ANY)
            self.AddSubPage(page48, "FSL Randomise Settings", wx.ID_ANY)
            self.AddSubPage(page49, "PyBASC Settings", wx.ID_ANY)
            self.AddSubPage(page50, "MDMR Settings", wx.ID_ANY)
            self.AddSubPage(page51, "ISC Settings", wx.ID_ANY)

        self.Bind(wx.EVT_TREEBOOK_PAGE_CHANGED, self.OnPageChanged)
        self.Bind(wx.EVT_TREEBOOK_PAGE_CHANGING, self.OnPageChanging)

        # This is a workaround for a sizing bug on Mac...
        wx.FutureCall(100, self.AdjustSize)
        
        self.SetSelection(1)

        self.Refresh()

    def OnPageChanged(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        event.Skip()

    def OnPageChanging(self, event):
        old = event.GetOldSelection()
        new = event.GetSelection()
        sel = self.GetSelection()
        event.Skip()

    def AdjustSize(self):
        self.GetTreeCtrl().InvalidateBestSize()
        self.SendSizeEvent()

    def get_page_list(self):
        return self.page_list

    def get_page(self, index):
        return self.page_list[index]


class MainFrame(wx.Frame):

    def __init__(self, parent, option='save', path="", pipeline_id="",
                 ind=True, size=(1200, 520)):
        wx.Frame.__init__(
            self, parent=parent, title="CPAC Pipeline Configuration",
            size=size)

        self.ind = ind

        # Here we create a panel and a notebook on the panel
        self.p = wx.Panel(self)
        self.nb = Mybook(self.p, ind=ind)
        self.path = path
        self.pipeline_id = pipeline_id
        self.option = option
        self.parent = parent

        btnPanel = wx.Panel(self.p, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        submit = wx.Button(btnPanel, wx.ID_SAVE, "Save",
            (280, 10), wx.DefaultSize, 0)
        hbox.Add(submit, 0.6, wx.ALIGN_RIGHT | wx.ALL, 5)
        self.Bind(wx.EVT_BUTTON, self.submit_item, id=wx.ID_SAVE)

        if ind:
            testConfig = wx.Button(btnPanel, wx.ID_PREVIEW,
                                   "Test Configuration", (350, 10),
                                   wx.DefaultSize, 0)
            hbox.Add(testConfig, 0, wx.ALIGN_RIGHT | wx.ALL, 5)
            self.Bind(wx.EVT_BUTTON, self.testConfig, id=wx.ID_PREVIEW)
        
        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel",
            (220, 10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, wx.ALIGN_RIGHT | wx.ALL, 5)
        
        btnPanel.SetSizer(hbox)

        # finally, put the notebook in a sizer for the panel to manage
        # the layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.nb, 1, wx.EXPAND)
        sizer.Add(btnPanel, 0.6, wx.EXPAND | wx.RIGHT, 20)

        self.p.SetSizer(sizer)

        self.Layout()
        self.Show()
        if option == 'edit' or option == 'load':
            self.load()

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
        
        if self.config_map['participant_id_label'] in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(self.config_map['participant_id_label'])
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

    def OpenPage(self, page):

        tree = self.nb.GetTreeCtrl()

        def GetItemByName(tree_ctrl_instance, search_text):
            retval = None
            root_list = [tree_ctrl_instance.GetRootItem()]
            for root_child in root_list:
                item, cookie = tree_ctrl_instance.GetFirstChild(root_child)
                while item.IsOk():
                    if tree_ctrl_instance.GetItemText(item) == search_text:
                        retval = item
                        break
                    if tree_ctrl_instance.ItemHasChildren(item):
                        root_list.append(item)
                    item, cookie = tree_ctrl_instance.GetNextChild(root_child, cookie)
            return retval
            
        page_item = GetItemByName(tree, page)
        tree.SelectItem(page_item)
        return page_item

    def load(self):
        import os
        import yaml

        try:
            config_file_map = yaml.load(open(self.path, 'r'))
        except:
            raise Exception("Error importing file - %s , Make"
                            " sure it is in correct yaml format")
        self.config_map = config_file_map

        # repopulate the model setup checkbox grid, since this has to be
        # done specially
        pheno_exists = False
        if 'pheno_file' in config_file_map.keys():
            if config_file_map['pheno_file']:
                if os.path.isfile(config_file_map['pheno_file']):
                    pheno_exists = True
                    phenoFile = open(
                        os.path.abspath(config_file_map['pheno_file']))
                    self.get_pheno_header(phenoFile)

        for page in self.nb.get_page_list():
            ctrl_list = page.page.get_ctrl_list()

            for ctrl in ctrl_list:
                name = ctrl.get_name()
                val = config_file_map.get(str(name))
                sample_list = ctrl.get_values()
                s_map = dict((v, k)
                            for k, v in substitution_map.iteritems())

                if name == 'Regressors':
                    ctrl.set_value(val)
                    continue

                if name == 'model_setup':
                    # update the 'Model Setup' box and populate it with the 
                    # EVs and their associated checkboxes for categorical 
                    # and demean
                    if pheno_exists:
                        ctrl.set_value(self.phenoHeaderItems)
                        ctrl.set_selection(config_file_map['ev_selections'])

                if val is not None:
                    if ("list" in name) and (name != "participant_list") and ('basc' not in name):
                        try:
                            mapped_vals = [s_map.get(item) for item in val if s_map.get(item) != None]
                        except TypeError:
                            mapped_vals = None

                        if not mapped_vals:
                            val = [str(item) for item in val]
                        else:
                            val = mapped_vals

                        new_derlist = []
                        for val in val:
                            new_derlist.append(val)

                        if len(new_derlist) > 0:
                            ctrl.set_value(new_derlist)
                        else:
                            ctrl.set_value(None)
                    elif name == 'z_threshold' or name == 'p_threshold':
                        try:
                            val = val[0]
                            ctrl.set_value(val)
                        except TypeError:
                            # if the user has put it in as a float and not a list
                            ctrl.set_value(str(val))  

                    elif name == 'group_sep':
                        val = s_map.get(val)
                        ctrl.set_value(val)

                    elif name == 'grouping_var':
                        if isinstance(val, list) or "[" in val:
                            grouping_var = ""
                            for cov in val:
                                grouping_var += "{0},".format(cov)
                            grouping_var = grouping_var.rstrip(",")
                        else:
                            grouping_var = val

                        ctrl.set_value(grouping_var)

                                      # and this line is just sad
                    elif name not in ('model_setup', 'derivative_list',
                                      'Regressors', 'nuisanceBandpassFreq',
                                      'tsa_roi_paths', 'sca_roi_paths'):
                        try:
                            ctrl.set_value(str(val))
                        except ValueError:
                            ctrl.set_value(int(val[0]))

                    if isinstance(val, list):
                        if ctrl.get_datatype() == 8:
                            value = []
                            for item in val:
                                data = ""
                                for k, v in item.iteritems():
                                    if v == 1 and k in sample_list:
                                        if data:
                                            data = data + "," + k
                                        else:
                                            data = k
                                value.append(data)

                        elif ctrl.get_datatype() == 6:
                            value = []
                            for v in val:
                                value.append(str(v))

                        elif ctrl.get_datatype() == 3:
                            value = [sample_list[i]
                                     for i, x in enumerate(val) if x == True]

                        elif ctrl.get_datatype() == 4:
                            if 1 in val and 0 in val:
                                val = [10]
                                
                            if 'ANTS' in val and 'FSL' in val:
                                val = [11]

                            if '3dAutoMask' in val and 'BET' in val:
                                val = [12]
                            if 'AFNI' in val and 'BET' in val:
                                val = [12]
                            
                            value = [s_map.get(item)
                                         for item in val if s_map.get(item) != None]
                            if not value:
                                value = [str(item) for item in val]

                        elif ctrl.get_datatype() == 5 and \
                            ctrl.get_type() == 6:
                                value = [sample_list[v] for v in val]

                        elif ctrl.get_datatype() == 9:
                            value = val[0] # pass the dictionary straight up
                                
                        else:
                            value = None
                            for v in val:
                                if value:
                                    value = value + "," + str(v)
                                else:
                                    value = str(v)
                    else:
                        if ctrl.get_datatype() == 2 and \
                            ctrl.get_type() == 0 and \
                                str(val) not in sample_list:
                                value = sample_list[val]
                        else:
                            value = str(val)
                else:
                    value = ""

                ctrl.set_value(value)

    # Test the subject list
    def test_sublist(self, sublist):
        '''
        Instance method to test a subject list for errors

        Parameters
        ----------
        self : MainFrame (wx.Frame object)
            the method is aware of the instance as self
        sublist : list (dict)
            a C-PAC-formatted subject list (yaml list of dictionaries)

        Returns
        -------
        pass_flg : boolean
            flag which indicates whether the subject list passed testing
        '''

        # Import packages
        import os
        import tempfile
        import nibabel as nb

        from CPAC.utils.datasource import check_for_s3

        # Init variables
        err_str = ''
        err_msg = ''
        not_found_flg = False
        bad_dim_flg = False
        pass_flg = False
        checked_s3 = False
        s3_str = 's3://'

        # Check to ensure the user is providing an actual subject
        # list and not some other kind of file
        try:
            subInfo = sublist[0]
        except:
            msg = 'ERROR: Subject list file not in proper format - ' \
                  'check if you loaded the correct file? \n\n'\
                  'Error name: config_window_0001'
            errDlg4 = wx.MessageDialog(self, msg, 'Subject List Error',
                                       wx.OK | wx.ICON_ERROR)
            errDlg4.ShowModal()
            errDlg4.Destroy()
            # Raise Exception
            raise Exception

        # Another check to ensure the actual subject list was generated
        # properly and that it will work
        if 'subject_id' not in subInfo:
            msg = 'ERROR: Subject list file not in proper format - '\
                  'check if you loaded the correct file? \n\n'\
                  'Error name: config_window_0002'
            errDlg3 = wx.MessageDialog(self, msg , 'Subject List Error',
                                       wx.OK | wx.ICON_ERROR)
            errDlg3.ShowModal()
            errDlg3.Destroy()
            # Raise Exception
            raise Exception

        # Iterate and test each subject's files
        for sub in sublist:
            anat_file = sub['anat']

            func_files = None
            if 'func' in sub:
                func_files = sub['func']
            elif 'rest' in sub:
                func_files = sub['rest']

            checked_anat_s3 = False

            if not anat_file:
                err = "\n\n[!] Could not read in at least one of your anatom"\
                      "ical input files. Please double-check the formatting "\
                      "of your participant list YAML file.\n\n"
                raise Exception(err)

            if func_files and not isinstance(func_files, dict):
                err = "\n\n[!] The functional files in the participant " \
                      "list YAML should be listed with a scan name key and " \
                      "a file path value.\n\nFor example:\nfunc_1: " \
                      "/path/to/func_1.nii.gz\n\n"
                raise Exception(err)

            if anat_file.lower().startswith(s3_str):
                dl_dir = tempfile.mkdtemp()
                try:
                    creds_path = sub['creds_path']
                except KeyError:
                    # if no creds path is provided, it could be that the user
                    # is downloading public data - leave it to downstream to
                    # handle creds issues
                    creds_path = None
                anat_file = check_for_s3(anat_file, creds_path, dl_dir=dl_dir)
                checked_anat_s3 = True
            # Check if anatomical file exists
            if os.path.exists(anat_file):
                try:
                    img = nb.load(anat_file)
                except Exception as e:
                    print(e)
                    continue
                hdr = img.get_header()
                dims = hdr.get_data_shape()
                # Check to make sure it has the proper dimensions
                if len(dims) != 3:
                    bad_dim_flg = True
                    err_str_suffix = 'Anat file not 3-dimensional: %s\n' \
                                     % anat_file
                    err_str = err_str + err_str_suffix
            # Anat file doesnt exist
            else:
                not_found_flg = True
                err_str_suffix = 'File not found: %s\n' % anat_file
                err_str = err_str + err_str_suffix
            # If we're just checking s3 files, remove the temporarily 
            # downloaded
            if checked_anat_s3:
                try:
                    os.remove(anat_file)
                except:
                    pass

            if func_files:
                for func_file in func_files.values():
                    checked_s3 = False
                    if '.nii' not in func_file:
                        # probably a JSON file
                        continue

                    if func_file.lower().startswith(s3_str):
                        dl_dir = tempfile.mkdtemp()
                        try:
                            creds_path = sub['creds_path']
                        except KeyError:
                            # if no creds path is provided, it could be that the 
                            # user is downloading public data - leave it to down-
                            # stream to handle creds issues
                            creds_path = None
                        func_file = check_for_s3(func_file, creds_path,
                                                dl_dir=dl_dir, img_type='func')
                        checked_s3 = True

                    if os.path.exists(func_file):
                        try:
                            img = nb.load(func_file)
                        except Exception as e:
                            print(e)
                            continue
                        hdr = img.get_header()
                        dims = hdr.get_data_shape()
                        # Check to make sure it has the proper dimensions
                        if len(dims) != 4:
                            bad_dim_flg = True
                            err_str_suffix = 'Func file not 4-dimensional: %s\n' \
                                            % func_file
                            err_str = err_str + err_str_suffix
                    # Functional file doesnt exist
                    else:
                        not_found_flg = True
                        err_str_suffix = 'File not found: %s\n' % func_file
                        err_str = err_str + err_str_suffix
                    # If we're just checking s3 files, remove the temporarily 
                    # downloaded
                    if checked_s3:
                        try:
                            os.remove(func_file)
                        except:
                            pass

            # Check flags for error message
            if not_found_flg:
                err_msg = 'One or more of your input files are missing.\n'
            if bad_dim_flg:
                err_msg = ''.join([err_msg, 'One or more of your input '
                                            'images have improper '
                                            'dimensionality\n'])
            # If err_msg was populated, display in window
            if err_msg:
                err_msg = 'ERROR: ' + err_msg + \
                          'See terminal output for more details'
                errDlgFileTest = wx.MessageDialog(self,
                                                  err_msg,
                                                  'Pipeline Not Ready',
                                                  wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()
                raise Exception(err_str)
            else:
                pass_flg = True

        # Return the flag
        return pass_flg

    # Test pipeline config file
    def testConfig(self, event):
        '''
        This function runs when the user clicks the "Test Configuration"
        button in the pipeline configuration window.
        
        It prompts the user for a sample subject list (i.e. one that they will
        be using with the config they are building). Then it builds the
        pipeline but does not run it. It then reports whether or not the
        config will run or not depending on if the pipeline gets built
        successfully.
        '''

        # Import packages
        import os
        import yaml
        from CPAC.utils import Configuration

        from CPAC.pipeline.cpac_pipeline import prep_workflow
        from CPAC.pipeline.cpac_runner import build_strategies

        def display(win, msg, changeBg=True):
            wx.MessageBox(msg, "Error", style=wx.OK | wx.ICON_ERROR)
            if changeBg:
                win.SetBackgroundColour("pink")
            win.SetFocus()
            win.Refresh()

        # Collect a sample subject list and parse it in
        testDlg0 = wx.MessageDialog(
            self, 'This tool will run a quick check on the current pipeline '
                  'configuration. Click OK to provide a subject list you '
                  'will be using with this setup.',
            'Subject List',
            wx.OK | wx.ICON_INFORMATION)
        testDlg0.ShowModal()
        testDlg0.Destroy()
        
        dlg = wx.FileDialog(
            self, message="Choose the CPAC Subject list file",
            defaultDir=os.getcwd(), 
            defaultFile="CPAC_subject_list.yml",
            wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
            style=wx.OPEN | wx.CHANGE_DIR)
        
        if dlg.ShowModal() == wx.ID_OK:
            subListPath = dlg.GetPath()
        
        # Load and test the subject list
        print('Checking data configuration: {0}...'.format(subListPath))
        sublist = yaml.load(open(os.path.realpath(subListPath), 'r'))
        sub_flg = self.test_sublist(sublist)
        if not sub_flg:
            raise Exception
        print('Data configuration looks good!')
        # Following code reads in the parameters and selections from the
        # pipeline configuration window and populate the config_list

        config_list = []
        wf_counter = []

        for page in self.nb.get_page_list():

            switch = page.page.get_switch()

            ctrl_list = page.page.get_ctrl_list()
            validate = False

            if switch:
                switch_val = str(switch.get_selection()).lower()

                if switch_val == 'on' or switch_val == 'true' or \
                    switch_val == '1':

                    validate = True
                    wf_counter.append(page.get_counter())

            for ctrl in ctrl_list:

                # option_name will be the selection name as it is written
                # as the dictionary key of the config.yml dictionary
                option_name = ctrl.get_name()

                #validating
                if (switch == None or validate) and ctrl.get_validation() \
                    and (option_name != 'derivativeList') and \
                        (option_name != 'modelConfigs'):

                    win = ctrl.get_ctrl()
                    
                    if isinstance(ctrl.get_selection(), list):
                        value = ctrl.get_selection()
                        if not value and 'session' not in option_name and 'series' not in option_name:
                            display(
                                win, "\"%s\" field is empty or the items are " \
                                     "not checked!" % ctrl.get_pretty_name(), False)
                            return

                    elif (option_name == "tsa_roi_paths") or \
                             (option_name == "sca_roi_paths"):

                        # fires if the control is the checkbox grid for
                        # multiple paths assigned to multiple options
                        # (i.e. timeseries analysis)

                        config_list.append(ctrl)
                        continue

                    else:
                        value = str(ctrl.get_selection())

                    if len(value) == 0:
                        display(win, "\"%s\" field is empty!" % ctrl.get_pretty_name())
                        return
                        
                    if '/' in value and '$' not in value and not \
                        isinstance(value, list):

                        if not os.path.exists(ctrl.get_selection()) and \
                                        value != 'On/Off':
                            display(
                                win, "\"%s\" field contains incorrect path. " \
                                "Please update the path!" % ctrl.get_pretty_name())
                            return
                    
                config_list.append(ctrl)

        # Write out a pipeline_config file, read it in and then delete it
        # (Will revise the data structure of the config files later so this
        # can just pass the data structure instead of doing it this way)
        try:
            test_cfg_yml = '/tmp/test_config.yml'
            self.write(test_cfg_yml, config_list)
            c = Configuration(yaml.load(open(os.path.realpath(test_cfg_yml), 'r')))
            os.remove(test_cfg_yml)
        except:
            errDlg2 = wx.MessageDialog(
                self, 'A problem occurred with preparing the pipeline test run. \n\n' \
                      'Please ensure you have rights access to the directories you' \
                      ' have chosen for the CPAC working, crash, and output folders.',
                'Test Configuration Error',
                wx.OK | wx.ICON_ERROR)
            errDlg2.ShowModal()
            errDlg2.Destroy()

        if (1 in c.runNuisance) or (c.Regressors != None):
            strategies = sorted(build_strategies(c))
        else:
            strategies = None

        # Run the actual pipeline building prep and see if it works or not
        testDlg1 = wx.MessageDialog(
            self, 'Click OK to run the test. This should take only a few seconds.',
            'Running Test',
            wx.OK | wx.ICON_INFORMATION)
        testDlg1.ShowModal()

        # Check file paths first
        
        # Just getting proper names of config file parameters
        try:
            params_file = open(p.resource_filename('CPAC', 'GUI/resources/config_parameters.txt'), "r")
        except:
            print "Error: Could not open configuration parameter file.", "\n"
            raise Exception            

        paramInfo = params_file.read().split('\n')
        paramList = []

        for param in paramInfo:
            if param != '':
                paramList.append(param.split(','))

        # function for file path checking
        def testFile(filepath, paramName, switch):
            try:
                if (1 in switch) and (filepath != None):
                    fileTest = open(filepath)
                    fileTest.close()
            except:
                testDlg1.Destroy()
                
                for param in paramList:
                    if param[0] == paramName:
                        paramTitle = param[1]
                        paramGroup = param[2]
                        break
                    
                errDlgFileTest = wx.MessageDialog(
                    self, 'Error reading file - either it does not exist or '\
                          'you do not have read access. \n\n' \
                          'Parameter: %s \n' \
                          'In tab: %s \n\n' \
                          'Path: %s' % (paramTitle, paramGroup, filepath),
                    'Pipeline Not Ready',
                    wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()

        # Check S3 output bucket access if writing to S3
        output_dir = c.outputDirectory
        s3_str = 's3://'
        if output_dir.lower().startswith(s3_str):
            output_dir_sp = output_dir.split('/')
            output_dir_sp[0] = output_dir_sp[0].lower()
            output_dir = '/'.join(output_dir_sp)

        if type(output_dir) is str and output_dir.lower().startswith(s3_str):
            from indi_aws import fetch_creds
            creds_path = c.awsOutputBucketCredentials
            bucket_name = output_dir.split(s3_str)[1].split('/')[0]
            try:
                bucket = fetch_creds.return_bucket(creds_path, bucket_name)
                print 'Connection with output bucket "%s" successful!' % bucket_name
            except Exception as exc:
                err_msg = 'Unable to access output S3 bucket: "%s" with '\
                          'credentials in: "%s". Check bucket name '\
                          'and credentials file and try again'\
                          % (bucket_name, creds_path)
                testDlg1.Destroy()

                errDlg1 = wx.MessageDialog(self, err_msg, 'Pipeline Not Ready',
                                           wx.OK | wx.ICON_ERROR)
                errDlg1.ShowModal()
                errDlg1.Destroy()
                return

        testFile(c.template_brain_only_for_anat, \
                     'template_brain_only_for_anat',[1])
        testFile(c.template_skull_for_anat,'template_skull_for_anat',[1])
        testFile(c.PRIORS_WHITE,'PRIORS_WHITE',c.runSegmentationPreprocessing)
        testFile(c.PRIORS_GRAY,'PRIORS_GRAY',c.runSegmentationPreprocessing)
        testFile(c.PRIORS_CSF,'PRIORS_CSF',c.runSegmentationPreprocessing)
        testFile(c.template_brain_only_for_func, \
                     'template_brain_only_for_func',c.runRegisterFuncToMNI)
        testFile(c.template_skull_for_func,'template_skull_for_func', \
                     c.runRegisterFuncToMNI)
        testFile(c.identityMatrix,'identityMatrix',c.runRegisterFuncToMNI)
        testFile(c.boundaryBasedRegistrationSchedule, \
                     'boundaryBasedRegistrationSchedule', \
                     c.runRegisterFuncToAnat)
        testFile(c.lateral_ventricles_mask,'lateral_ventricles_mask', \
                     c.runNuisance)
        testFile(c.template_symmetric_brain_only, \
                     'template_symmetric_brain_only',c.runVMHC)
        testFile(c.template_symmetric_skull,'template_symmetric_skull', \
                     c.runVMHC)
        testFile(c.dilated_symmetric_brain_mask, \
                     'dilated_symmetric_brain_mask',c.runVMHC)
        testFile(c.configFileTwomm,'configFileTwomm',c.runVMHC)
        testFile(c.templateSpecificationFile,'templateSpecificationFile', \
                     c.runNetworkCentrality)

        if c.tsa_roi_paths and type(c.tsa_roi_paths[0]) == dict:
            for roi_path in c.tsa_roi_paths[0].keys():
                testFile(roi_path, "tsa_roi_paths", c.runROITimeseries)
        if c.sca_roi_paths and type(c.sca_roi_paths[0]) == dict:
            for roi_path in c.sca_roi_paths[0].keys():
                testFile(roi_path, "sca_roi_paths", c.runSCA)
        try:
            # Run the pipeline building
            prep_workflow(sublist[0], c, strategies, 0)

        except Exception as xxx:
            print xxx
            print "an exception occurred"
            
            testDlg1.Destroy()
            
            errDlg1 = wx.MessageDialog(
                self, 'There are issues with the current configuration ' \
                      'which need to be resolved - please check to make ' \
                      'sure the options you are running have the proper ' \
                      'pre-requisites selected.\n\nIssue Info:\n%s' \
                      % str(xxx),
                'Pipeline Not Ready',
                wx.OK | wx.ICON_ERROR)
            errDlg1.ShowModal()
            errDlg1.Destroy()
            
        else:
            testDlg1.Destroy()
            
            okDlg1 = wx.MessageDialog(
                self, 'The current configuration will run successfully. You '\
                      'can safely save and run this setup!',
                'Pipeline Ready',
                wx.OK | wx.ICON_INFORMATION)
            okDlg1.ShowModal()
            okDlg1.Destroy()

    def submit_item(self, event):
        import os
        import linecache

        def display(win, msg, changeBg=True):
            wx.MessageBox(msg, "Error", style=wx.OK | wx.ICON_ERROR)
            if changeBg:
                win.SetBackgroundColour("pink")
            win.SetFocus()
            win.Refresh()

        config_list = []
        hash_val = 0
        wf_counter = []

        for page in self.nb.get_page_list():

            switch = page.page.get_switch()

            ctrl_list = page.page.get_ctrl_list()
            validate = False

            if switch:
                switch_val = str(switch.get_selection()).lower()

                if switch_val == 'on' or switch_val == 'true' or \
                    switch_val == '1':

                    validate = True
                    wf_counter.append(page.get_counter())

            for ctrl in ctrl_list:

                # option_name will be the selection name as it is written
                # as the dictionary key of the config.yml dictionary
                option_name = ctrl.get_name()

                # validating
                if (switch == None or validate) and ctrl.get_validation() \
                    and option_name not in ['derivativeList', 'modelConfigs', 'sessions_list', 'series_list']:
                
                    win = ctrl.get_ctrl()
                    
                    if isinstance(ctrl.get_selection(), list):

                        # fires if the control is a ListBoxCombo (type 7 in
                        # generic_class.py)

                        value = ctrl.get_selection()

                        if not value:
                            display(
                                win, "\"%s\" field is empty or the items are " \
                                     "not checked!" % ctrl.get_name(), False)
                            return

                    elif (option_name == "tsa_roi_paths") or \
                             (option_name == "sca_roi_paths"):

                        # fires if the control is the checkbox grid for
                        # multiple paths assigned to multiple options
                        # (i.e. timeseries analysis)

                        config_list.append(ctrl)
                        continue

                    elif option_name == "templateSpecificationFile":

                        # let's make sure this is a NIFTI file
                        if not ctrl.get_selection().endswith(".nii") and \
                                not ctrl.get_selection().endswith(".nii.gz"):
                            display(
                                win, "The Mask Specification File field "
                                     "must contain a NIFTI file (ending in "
                                     ".nii or .nii.gz).", False)
                            return

                    else:

                        value = str(ctrl.get_selection())

                    if len(value) == 0:
                        display(win, "\"%s\" field is empty!" % ctrl.get_pretty_name())
                        return

                    if '/' in value and '$' not in value and \
                        not isinstance(value, list):

                        if not os.path.exists(ctrl.get_selection()) and \
                                        value != 'On/Off':
                            display(
                                win, "\"%s\" field contains incorrect path. "
                                     "Please update the path!"
                                     % ctrl.get_pretty_name())
                            return

                config_list.append(ctrl)

        # Get the user's CPAC pipeline name for use in this script
        # ....or the pipeline directory for group configs
        pipeline_name = None
        pipeline_dir = None
        for config in config_list:
            if config.get_name() == 'pipelineName':
                pipeline_name = config.get_selection()
                if len(pipeline_name) == 0:
                    noNameDlg = wx.MessageDialog(
                        self, 'Please enter a pipeline name.',
                        'Error!',
                        wx.OK | wx.ICON_ERROR)
                    noNameDlg.ShowModal()
                    noNameDlg.Destroy()
                    return
            if config.get_name() == 'pipeline_dir':
                pipeline_dir = config.get_selection()
                if len(pipeline_dir) == 0:
                    noNameDlg = wx.MessageDialog(
                        self, 'Please enter a pipeline directory.',
                        'Error!',
                        wx.OK | wx.ICON_ERROR)
                    noNameDlg.ShowModal()
                    noNameDlg.Destroy()
                    return

        if self.ind:
            dlg = wx.FileDialog(
                self, message="Save CPAC configuration file as ...",
                defaultDir=os.getcwd(),
                defaultFile=("pipeline_config_{0}".format(pipeline_name)),
                wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
                style=wx.SAVE)
        elif not self.ind and pipeline_dir:
            pipeline_name = pipeline_dir.rstrip('/').split('/')[-1]
            dlg = wx.FileDialog(
                self, message="Save CPAC configuration file as ...",
                defaultDir=os.getcwd(),
                defaultFile=("group_config_{0}".format(pipeline_name)),
                wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
                style=wx.SAVE)

        if dlg.ShowModal() == wx.ID_OK:
            self.path = dlg.GetPath()

            # Strips any user-input file extension and enforces .yml as
            # the extension
            self.path = os.path.splitext(self.path)[0] + '.yml'

            self.write(self.path, config_list)
            dlg.Destroy()
            if self.option == 'save':

                # this runs if you hit 'Save' from within the pipeline config
                # editor AND the editor was opened from the main window by
                # clicking 'New' instead of 'Edit'

                ### this is the old code for generating random city names
                ### to name pipeline configs. remove at some point?
                #for counter in wf_counter:
                #    if counter != 0:
                #        hash_val += 2 ** counter
                #print "wf_counter -- ", wf_counter
                #print "hashval --> ", hash_val
                #pipeline_id = linecache.getline(p.resource_filename('CPAC', \
                #       'GUI/resources/pipeline_names.py'), hash_val)

                if os.path.exists(self.path):
                    self.update_listbox(pipeline_name)

            else:

                # this runs if you hit 'Save' from within the pipeline config
                # editor AND the editor was opened from the main window by
                # clicking 'Edit' instead of 'New'

                pipeline_map = self.parent.get_pipeline_map()

                if self.option == 'load':
                    # this runs if your pipeline config is being updated
                    pipeline_map[pipeline_name] = self.path
                    self.Parent.listbox.Append(pipeline_name)

                elif pipeline_map.get(pipeline_name) != None:
                    # this runs if you hit Edit, change your pipeline config
                    # file BUT keep the Pipeline Name the same and save it
                    pipeline_map[pipeline_name] = self.path

                else:
                    # this runs if you hit Edit, change your pipeline config
                    # AND also change the Pipeline Name and save it with the
                    # new path - this adds the new pipeline to the listbox on
                    # the main CPAC window
                    pipeline_map[pipeline_name] = self.path
                    self.Parent.listbox.Append(pipeline_name)
                

            self.SetFocus()
            self.Close()

    def cancel(self, event):
        self.Close()

    def update_listbox(self, value):

        if len(value) > 0:
            self.pipeline_id = value
            pipeline_map = self.parent.get_pipeline_map()
            if pipeline_map.get(self.pipeline_id) == None:
                pipeline_map[self.pipeline_id] = self.path
                self.Parent.listbox.Append(self.pipeline_id)

            else:
                dlg2 = wx.MessageDialog(
                    self, 'Pipeline already exists. Please enter a new name',
                    'Error!',
                    wx.OK | wx.ICON_ERROR)
                dlg2.ShowModal()
                dlg2.Destroy()

    def write(self, path, config_list):

        import ast
        import CPAC

        try:

            f = open(path, 'w')

            print >>f, "# CPAC Pipeline Configuration YAML file"
            print >>f, "# Version %s\n#" % str(CPAC.__version__)
            print >>f, "# http://fcp-indi.github.io for more info.\n#"
            print >>f, "# Tip: This file can be edited manually with a " \
                       "text editor for quick modifications.\n\n"

            for item in config_list:

                label = item.get_name()
                value = item.get_selection()
                dtype = item.get_datatype()
                item_type = item.get_type()

                sample_list = item.get_values()

                comment = item.get_help()

                for line in comment.split("\n"):
                    if line:
                        print>>f, "#", line

                # prints setting names and values (ex. " runAnatomicalProcessing: [1] ") into the
                # pipeline_config file, using a different set of code depending on the data type

                # parameters that are strings (ex. " False " or a path)
                if dtype == 0 or dtype == 1:

                    print >>f, label, ": ", str(value)
                    print >>f,"\n"

                # parameters that are integers
                elif dtype == 2:

                    # Add check for ReHo cluster
                    if label == 'clusterSize':
                        #print 'Using ReHo cluster size of ', value
                        pass
                    elif item_type == 0:
                        value = sample_list.index(value)
                    else:
                        if substitution_map.get(value) != None:
                            value = substitution_map.get(value)
                        elif value != 'None':
                            try:
                                value = ast.literal_eval(str(value))
                            except Exception as e:
                                raise Exception("could not parse value: "
                                                "{0}\n\n{1}"
                                                "\n".format(str(value), e))
                    
                    print >>f, label, ": ", value
                    print >>f,"\n"

                # parameters that are lists (ex. " [False, False] ")
                elif dtype == 3:

                    map = ast.literal_eval(str(value))
                    values = []
                    for x in range(0, len(map.keys())):
                        values.append(False)
                    for k, v in map.iteritems():
                        item, idx = k
                        values[idx] = v

                    print>>f, label, ": ", values
                    print>>f,"\n"

                # parameters that are switches (ex. [0] or [1] )
                elif dtype == 4:

                    values=[]

                    if isinstance(value, list):
                        value = ast.literal_eval(str(value))
                    else:
                        value = str(value).split(",")

                    for val in value:
                        val = val.strip()
                        sval = substitution_map.get(val)
                        if sval != None:
                            values.append(sval)
                        else:
                            values.append(val)
                            
                    if values == [10]:
                        values = [1,0]
                    elif values == [11]:
                        values = ['ANTS','FSL']
                    elif values == [12]:
                        values = ['3dAutoMask','BET']

                    print>>f, label, ": ", values
                    print>>f, "\n"

                # parameters that are bracketed numbers (int or float)
                elif dtype == 5:

                    ### parse user input   ### can't use internal function
                    # type() here???
                    if value.find(',') != -1:
                        lvalue = value.split(',')
                    elif value.find(';') != -1:
                        lvalue = value.split(';')
                    elif value.find(':') != -1:
                        lvalue = value.split(':')
                    else:
                        lvalue = [value]

                    if value.find('.') != -1:
                        lvalue = [float(item) for item in lvalue]
                    elif len(value) > 0:
                        try:
                            new_lvalue = [int(item) for item in lvalue]
                        except ValueError:
                            # this trips only if user inputs a percentage for
                            # de-spiking/scrubbing motion threshold
                            new_lvalue = [str(item) for item in lvalue]
                        lvalue = new_lvalue
                    else:
                        lvalue = 0
                    
                    print>>f, label, ":", lvalue   ###
                    print>>f, "\n"

                # parameters that are ? (bandpass filter specs)
                elif dtype == 6:

                    values = []
 
                    for val in ast.literal_eval(str(value)):
                        try:
                            val=ast.literal_eval(str(val))
                        except Exception as err:
                            print "Exception trying to translate: " \
                                  "%s, %s, %s, %s"%(label,str(value),val,err)
                            print "value type: %s"%(type(val))
                        values.append(ast.literal_eval(str(val)))

                    print>>f, label, ":", values
                    print>>f, "\n"

                # parameters that are whole words
                #     ALSO: the Nuisance Corrections lists                
                elif dtype == 8:

                    if type(value) == list:

                        import yaml
                        yml = yaml.dump({
                            label: value
                        }, default_flow_style=False)

                        print>>f, yml

                    else:

                        print>>f, label,":"

                        value = ast.literal_eval(str(value))

                        for val in value:
                            val = val.split(',')
                            f.write("  - ")
                            flag = 0
                            for sample in sample_list:
                                if flag == 0:
                                    space = ""
                                    flag = 1
                                else:
                                    space = "    "
                                if sample in val:
                                    print>>f, space, sample, ": ", 1
                                else:
                                    print>>f, space, sample, ": ", 0

                        print >>f, "\n"

                elif dtype == 9:

                    # checkbox grid (ROI extraction etc.)
                    string = gen_checkboxgrid_config_string(label, value)

                    print >>f, string
                    print >>f, "\n"

                elif label == 'derivative_list':
                    # this takes the user selection in the derivative
                    # list and matches it with the output directory
                    # folder name for each chosen derivative via the
                    # substitution map in constants.py

                    # go over each string in the list
                    value = []
                    for val in ast.literal_eval(str(item[1])):
                        if substitution_map.get(val) != None:
                            if substitution_map.get(val) not in value:
                                value.append(substitution_map.get(val))
                        elif val != 'None':
                            if ast.literal_eval(val) not in value:
                                value.append(ast.literal_eval(val))
                    print>>f, label, ":", value
                    print>>f, "\n"

                elif label == 'model_setup':
                    # basically, ctrl is checkbox_grid in this case, and
                    # get_selection goes to generic_class.py first, which links
                    # it to the custom GetGridSelection() function in the
                    # checkbox_grid class in custom_control.py
                    print>>f, 'ev_selections:', value
                    print>>f, "\n"

                else:
                    value = ast.literal_eval(str(value))

                    print>>f, label, ":", value
                    print>>f, "\n"

            f.close()

        except Exception, e:
            print e
            print "Error Writing the pipeline configuration file %s" % path
            raise


