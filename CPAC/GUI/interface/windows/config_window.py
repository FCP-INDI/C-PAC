import wx
from CPAC.GUI.interface.utils.constants import substitution_map
import pkg_resources as p
from CPAC.GUI.interface.pages import WorkflowConfig, Motion, AnatomicalPreprocessing, \
    DerivativesConfig, Segmentation,  Registration, FunctionalPreProcessing,\
    MotionOptions, Scrubbing, AnatToFuncRegistration, FuncToMNIRegistration,\
    VMHC, VMHCSettings, ReHo, ReHoSettings, \
    SCA, SCASettings, MultipleRegressionSCA,\
    Settings, ComputerSettings, DirectorySettings, \
    Nuisance, NuisanceCorrection, MedianAngleCorrection,\
    CentralitySettings, Centrality,\
    ALFF, ALFFSettings,\
    Smoothing, SmoothingSettings,\
    Filtering, FilteringSettings,\
    TimeSeries, ROITimeseries, VOXELTimeseries, \
    SpatialRegression, GenerateSeeds, VerticesTimeSeries,\
    GroupAnalysis, GPASettings, BASCSettings,\
    BASC, CWAS, CWASSettings,\
    DualRegression, DualRegressionOptions, TimeSeriesOptions


ID_SUBMIT = 6


class Mybook(wx.Treebook):

    def __init__(self, parent):
        wx.Treebook.__init__(self, parent, wx.ID_ANY, style=
                             wx.BK_DEFAULT)
        self.page_list = []

        # create the page windows as children of the notebook
        page1 = Settings(self)
        page2 = ComputerSettings(self)
        page3 = DirectorySettings(self)
        page4 = WorkflowConfig(self)
        page47 = DerivativesConfig(self)

        page5 = AnatomicalPreprocessing(self)
        page6 = Registration(self, 1)
        page7 = Segmentation(self, 2)

        page8 = FunctionalPreProcessing(self)
        page9 = TimeSeriesOptions(self)
        page10 = AnatToFuncRegistration(self, 5)
        page11 = FuncToMNIRegistration(self, 6)

        page12 = Nuisance(self)
        page13 = NuisanceCorrection(self, 7)
        page14 = MedianAngleCorrection(self, 8)

        page15 = Filtering(self)
        page16 = FilteringSettings(self, 9)

        page17 = Motion(self)
        page18 = MotionOptions(self)
        page19 = Scrubbing(self, 4)

        page20 = TimeSeries(self)
        page21 = GenerateSeeds(self)
        page22 = ROITimeseries(self)
        page23 = VOXELTimeseries(self)
        page24 = VerticesTimeSeries(self)
        page25 = SpatialRegression(self)

        page26 = SCA(self)
        page27 = SCASettings(self)
        page28 = MultipleRegressionSCA(self)

        page29 = DualRegression(self)
        page30 = DualRegressionOptions(self)

        page31 = VMHC(self)
        page32 = VMHCSettings(self)

        page33 = ALFF(self)
        page34 = ALFFSettings(self)

        page35 = ReHo(self)
        page36 = ReHoSettings(self)

        page37 = Centrality(self)
        page38 = CentralitySettings(self)
        
        page39 = Smoothing(self)
        page40 = SmoothingSettings(self)

        page41 = BASC(self)
        page42 = BASCSettings(self)

        page43 = CWAS(self)
        page44 = CWASSettings(self)

        page45 = GroupAnalysis(self)
        page46 = GPASettings(self)

        # add the pages to the notebook with the label to show on the tab
        self.AddPage(page1, "Environment Setup", wx.ID_ANY)
        self.AddSubPage(page2, "Computer Settings", wx.ID_ANY)
        self.AddSubPage(page3, "Output Settings", wx.ID_ANY)
        self.AddSubPage(page4, "Preprocessing Workflow Options", wx.ID_ANY)
        self.AddSubPage(page47, "Derivatives Settings", wx.ID_ANY)

        self.AddPage(page5, "Anatomical Preprocessing", wx.ID_ANY)
        self.AddSubPage(page6, "Anatomical Registration", wx.ID_ANY)
        self.AddSubPage(page7, "Tissue Segmentation", wx.ID_ANY)

        self.AddPage(page8, "Functional Preprocessing", wx.ID_ANY)
        self.AddSubPage(page9, "Time Series Options", wx.ID_ANY)
        self.AddSubPage(page10, "Functional to Anatomical Registration", wx.ID_ANY)
        self.AddSubPage(page11, "Functional to MNI Registration", wx.ID_ANY)

        self.AddPage(page12, "Nuisance", wx.ID_ANY)
        self.AddSubPage(page13, "Nuisance Correction", wx.ID_ANY)
        self.AddSubPage(page14, "Median Angle Correction", wx.ID_ANY)

        self.AddPage(page15, "Temporal Filtering", wx.ID_ANY)
        self.AddSubPage(page16, "Temporal Filtering Options", wx.ID_ANY)

        self.AddPage(page17, "Motion Correction", wx.ID_ANY)
        self.AddSubPage(page18, "Motion Correction Options", wx.ID_ANY)
        self.AddSubPage(page19, "Scrubbing Options", wx.ID_ANY)

        self.AddPage(page20, "Time Series Extraction (TSE)", wx.ID_ANY)
        self.AddSubPage(page21, "Define New Seeds", wx.ID_ANY)
        self.AddSubPage(page22, "ROI Average TSE", wx.ID_ANY)
        self.AddSubPage(page23, "ROI Voxelwise TSE", wx.ID_ANY)
        self.AddSubPage(page24, "Surface Vertices TSE", wx.ID_ANY)
        self.AddSubPage(page25, "Spatial Regression", wx.ID_ANY)

        self.AddPage(page26, "Seed-based Correlation Analysis (SCA)", wx.ID_ANY)
        self.AddSubPage(page27, "SCA Options", wx.ID_ANY)
        self.AddSubPage(page28, "Mutiple Regression SCA Options", wx.ID_ANY)

        self.AddPage(page29, "Dual Regression", wx.ID_ANY)
        self.AddSubPage(page30, "Dual Regression Options", wx.ID_ANY)

        self.AddPage(page31, "Voxel-mirrored Homotopic Connectivity", wx.ID_ANY)
        self.AddSubPage(page32, "VMHC Settings", wx.ID_ANY)

        self.AddPage(page33, "ALFF and f/ALFF", wx.ID_ANY)
        self.AddSubPage(page34, "ALFF and f/ALFF Options", wx.ID_ANY)

        self.AddPage(page35, "Regional Homogeneity (ReHo)", wx.ID_ANY)
        self.AddSubPage(page36, "ReHo Options", wx.ID_ANY)

        self.AddPage(page37, "Network Centrality", wx.ID_ANY)
        self.AddSubPage(page38, "Network Centrality Options", wx.ID_ANY)
        
        self.AddPage(page39, "Spatial Smoothing", wx.ID_ANY)
        self.AddSubPage(page40, "Spatial Smoothing Options", wx.ID_ANY)

        self.AddPage(page41, "Bootstrap Analysis of Stable Clusters", wx.ID_ANY)
        self.AddSubPage(page42, "BASC Settings", wx.ID_ANY)

        self.AddPage(page43, "CWAS", wx.ID_ANY)
        self.AddSubPage(page44, "CWAS Settings", wx.ID_ANY)

        self.AddPage(page45, "Group Analysis", wx.ID_ANY)
        self.AddSubPage(page46, "Group Analysis Settings", wx.ID_ANY)

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

    def __init__(self, parent, option='save', path="", pipeline_id=""):
        wx.Frame.__init__(
            self, parent=parent, title="CPAC Pipeline Configuration", size=(1200, 520))

        # Here we create a panel and a notebook on the panel
        self.p = wx.Panel(self)
        self.nb = Mybook(self.p)
        self.path = path
        self.pipeline_id = pipeline_id
        self.option = option
        self.parent = parent

        btnPanel = wx.Panel(self.p, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        submit = wx.Button(btnPanel, wx.ID_SAVE, "Save", (
            280, 10), wx.DefaultSize, 0)
        hbox.Add(submit, 0.6, wx.ALIGN_RIGHT | wx.ALL, 5)
        self.Bind(wx.EVT_BUTTON, self.submit_item, id=wx.ID_SAVE)
        
        testConfig = wx.Button(btnPanel, wx.ID_PREVIEW, "Test Configuration", (
            350, 10), wx.DefaultSize, 0)
        hbox.Add(testConfig, 0, wx.ALIGN_RIGHT | wx.ALL, 5)
        self.Bind(wx.EVT_BUTTON, self.testConfig, id=wx.ID_PREVIEW)
        
        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel", (
            220, 10), wx.DefaultSize, 0)
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


    def load(self):

        import yaml
        try:
            config_file_map = yaml.load(open(self.path, 'r'))
        except:
            raise Exception("Error importing file - %s , Make"
                      " sure it is in correct yaml format")


        #for config in config_file_map:
        #    print "\n\n config: ", config, " selection: ", config_file_map[config]



        for page in self.nb.get_page_list():

            ctrl_list = page.page.get_ctrl_list()

            for ctrl in ctrl_list:

                name = ctrl.get_name()

                val = config_file_map.get(str(name))

                #print "loading ctrl ->", name, "->", val
                sample_list = ctrl.get_values()
                #print "sample_list -->", sample_list
                s_map = dict((v, k)
                            for k, v in substitution_map.iteritems())
                if val:
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
                            
                            value = [s_map.get(item)
                                         for item in val if s_map.get(item) != None]
                            if not value:
                                value = [ str(item) for item in val]
                            
                                
                        elif ctrl.get_datatype() == 5 and ctrl.get_type() == 6:
                                value = [sample_list[v] for v in val]
                                
                        else:
                            value = None
                            for v in val:
                                if value:
                                    value = value + "," + str(v)
                                else:
                                    value = str(v)
                    else:
                        if ctrl.get_datatype() == 2 and ctrl.get_type() == 0 and\
                        str(val) not in sample_list:
                                value = sample_list[val]
                        else:
                            value = str(val)
                else:
                    value = ""

                #print "setting value in ctrl -->", value
                #print "type -->", type(value)
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
        import nibabel as nb
        import os
        
        # Init variables
        err_str = ''
        err_msg = ''
        not_found_flg = False
        bad_dim_flg = False
        pass_flg = False
        
        # Check to ensure the user is providing an actual subject
        # list and not some other kind of file
        try:
            subInfo = sublist[0]
        except:
            errDlg4 = wx.MessageDialog(
                self, 'ERROR: Subject list file not in proper format - check if you' \
                        ' loaded the correct file? \n\n' \
                        'Error name: config_window_0001',
                'Subject List Error',
                wx.OK | wx.ICON_ERROR)
            errDlg4.ShowModal()
            errDlg4.Destroy()
    
            raise Exception  
            
        # Another check to ensure the actual subject list was generated
        # properly and that it will work
        if 'subject_id' not in subInfo:
            errDlg3 = wx.MessageDialog(
                self, 'ERROR: Subject list file not in proper format - check if you' \
                        ' loaded the correct file? \n\n' \
                        'Error name: config_window_0002',
                'Subject List Error',
                wx.OK | wx.ICON_ERROR)
            errDlg3.ShowModal()
            errDlg3.Destroy()
    
            raise Exception
        
        # Iterate and test each subject's files
        for sub in sublist:
            anat_file = sub['anat']
            func_files = sub['rest']
            # Check if anatomical file exists
            if os.path.exists(anat_file):
                img = nb.load(anat_file)
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
            # For each functional file
            for func_file in func_files.values():
                # Check if functional file exists
                if os.path.exists(func_file):
                    img = nb.load(func_file)
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
            # Check flags for error message
            if not_found_flg:
                err_msg = 'One or more of your input files are missing.\n'
            if bad_dim_flg:
                err_msg = err_msg + 'One or more of your input images have '\
                          'improper dimensionality\n'
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
        
        import os
        import yaml
        from CPAC.utils import Configuration
        
        from CPAC.pipeline.cpac_pipeline import prep_workflow
        from CPAC.pipeline.cpac_runner import build_strategies
        
        def display(win, msg, changeBg=True):
            wx.MessageBox(msg, "Error")
            if changeBg:
                win.SetBackgroundColour("pink")
            win.SetFocus()
            win.Refresh()
        
        # Collect a sample subject list and parse it in
        testDlg0 = wx.MessageDialog(
            self, 'This tool will run a quick check on the current pipeline configuration.' \
                  ' Click OK to provide a subject list you will be using with this setup.',
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
        sublist = yaml.load(open(os.path.realpath(subListPath), 'r'))
        sub_flg = self.test_sublist(sublist)
        if not sub_flg:
            raise Exception
        
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

                if switch_val == 'on' or switch_val == 'true' or switch_val == '1':
                    validate = True
                    wf_counter.append(page.get_counter())

            for ctrl in ctrl_list:

                # option_name will be the selection name as it is written
                # as the dictionary key of the config.yml dictionary
                option_name = ctrl.get_name()
                
                #validating
                if (switch == None or validate) and ctrl.get_validation() \
                    and (option_name != 'derivativeList') and (option_name != 'modelConfigs'):

                    win = ctrl.get_ctrl()
                    
                    if isinstance(ctrl.get_selection(), list):
                        value = ctrl.get_selection()
                        if not value:
                            display(
                                win, "%s field is empty or the items are not checked!" % ctrl.get_name(), False)
                            return
                    else:
                        value = str(ctrl.get_selection())

                    if len(value) == 0:
                        display(win, "%s field is empty!" % ctrl.get_name())
                        return
                        
                    if '/' in value and '$' not in value and not isinstance(value, list):

                        if not os.path.exists(ctrl.get_selection()) and value != 'On/Off':
                            display(
                                win, "%s field contains incorrect path. Please update the path!" % ctrl.get_name())
                            return
                    
                config_list.append(ctrl)
                
        

        # Get the user's CPAC output directory for use in this script
        for config in config_list:

            if config.get_name() == 'outputDirectory':
                outDir = config.get_selection()
        
        
        # Write out a pipeline_config file, read it in and then delete it
        # (Will revise the data structure of the config files later so this
        # can just pass the data structure instead of doing it this way)
        try:
            
            self.write(outDir + 'testConfig.yml', config_list)
            c = Configuration(yaml.load(open(os.path.realpath(outDir + 'testConfig.yml'), 'r')))
        
            os.remove(outDir + 'testConfig.yml')
        
        except:
        
            errDlg2 = wx.MessageDialog(
                self, 'A problem occurred with preparing the pipeline test run. \n\n' \
                      'Please ensure you have rights access to the directories you' \
                      ' have chosen for the CPAC working, crash, and output folders.',
                'Test Configuration Error',
                wx.OK | wx.ICON_ERROR)
            errDlg2.ShowModal()
            errDlg2.Destroy()


        
        if (1 in c.runNuisance) or (c.Corrections != None):
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
                    self, 'Error reading file - either it does not exist or you' \
                          ' do not have read access. \n\n' \
                          'Parameter: %s \n' \
                          'In tab: %s \n\n' \
                          'Path: %s' % (paramTitle, paramGroup, filepath),
                    'Pipeline Not Ready',
                    wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()
        
        
        testFile(c.template_brain_only_for_anat,'template_brain_only_for_anat',c.runRegistrationPreprocessing)
        testFile(c.template_skull_for_anat,'template_skull_for_anat',c.runRegistrationPreprocessing)
        testFile(c.PRIORS_WHITE,'PRIORS_WHITE',c.runSegmentationPreprocessing)
        testFile(c.PRIORS_GRAY,'PRIORS_GRAY',c.runSegmentationPreprocessing)
        testFile(c.PRIORS_CSF,'PRIORS_CSF',c.runSegmentationPreprocessing)
        testFile(c.template_brain_only_for_func,'template_brain_only_for_func',c.runRegisterFuncToMNI)
        testFile(c.template_skull_for_func,'template_skull_for_func',c.runRegisterFuncToMNI)
        testFile(c.identityMatrix,'identityMatrix',c.runRegisterFuncToMNI)
        testFile(c.boundaryBasedRegistrationSchedule,'boundaryBasedRegistrationSchedule',c.runRegisterFuncToAnat)
        testFile(c.lateral_ventricles_mask,'lateral_ventricles_mask',c.runNuisance)
        testFile(c.seedSpecificationFile,'seedSpecificationFile',[1])
        testFile(c.roiSpecificationFile,'roiSpecificationFile',c.runROITimeseries)
        testFile(c.roiSpecificationFileForSCA,'roiSpecificationFileForSCA',c.runROITimeseries)
        testFile(c.maskSpecificationFile,'maskSpecificationFile',c.runVoxelTimeseries)
        testFile(c.maskSpecificationFileForSCA,'maskSpecificationFileForSCA',c.runVoxelTimeseries)
        testFile(c.spatialPatternMaps,'spatialPatternMaps',c.runSpatialRegression)
        testFile(c.template_symmetric_brain_only,'template_symmetric_brain_only',c.runVMHC)
        testFile(c.template_symmetric_skull,'template_symmetric_skull',c.runVMHC)
        testFile(c.dilated_symmetric_brain_mask,'dilated_symmetric_brain_mask',c.runVMHC)
        testFile(c.configFileTwomm,'configFileTwomm',c.runVMHC)
        testFile(c.templateSpecificationFile,'templateSpecificationFile',c.runNetworkCentrality)
        testFile(c.bascAffinityThresholdFile,'bascAffinityThresholdFile',c.runBASC)
        testFile(c.cwasROIFile,'cwasROIFile',c.runCWAS)
        testFile(c.cwasRegressorFile,'cwasRegressorFile',c.runCWAS)
             
            
        try:
            
            # Run the pipeline building           
            prep_workflow(sublist[0], c, strategies, 0)

        except Exception as xxx:

            print xxx
            print "an exception occured"
            
            testDlg1.Destroy()
            
            errDlg1 = wx.MessageDialog(
                self, 'There are issues with the current configuration ' \
                      'which need to be resolved - please check to make ' \
                      'sure the options you are running have the proper ' \
                      'pre-requisites selected.\n\nIssue Info:\n%s' % xxx,
                'Pipeline Not Ready',
                wx.OK | wx.ICON_ERROR)
            errDlg1.ShowModal()
            errDlg1.Destroy()
            
        else:
            
            testDlg1.Destroy()
            
            okDlg1 = wx.MessageDialog(
                self, 'The current configuration will run successfully. You can safely' \
                      ' save and run this setup!',
                'Pipeline Ready',
                wx.OK | wx.ICON_INFORMATION)
            okDlg1.ShowModal()
            okDlg1.Destroy()



    def submit_item(self, event):
        import os
        import linecache

        def display(win, msg, changeBg=True):
            wx.MessageBox(msg, "Error")
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
                #print "switch_val ---->", switch_val
                if switch_val == 'on' or switch_val == 'true' or switch_val == '1':
                    validate = True
                    wf_counter.append(page.get_counter())

            for ctrl in ctrl_list:

                # option_name will be the selection name as it is written
                # as the dictionary key of the config.yml dictionary
                option_name = ctrl.get_name()
                # validating
                if (switch == None or validate) and ctrl.get_validation() \
                    and (option_name != 'derivativeList') and (option_name != 'modelConfigs'):
                
                    win = ctrl.get_ctrl()
                    
                    if isinstance(ctrl.get_selection(), list):

                        # fires if the control is a ListBoxCombo (type 7 in
                        # generic_class.py)

                        value = ctrl.get_selection()

                        if not value:
                            display(
                                win, "%s field is empty or the items are not checked!" % ctrl.get_name(), False)
                            return
                    else:
                        value = str(ctrl.get_selection())

                    if len(value) == 0:
                        display(win, "%s field is empty!" % ctrl.get_name())
                        return

                    if '/' in value and '$' not in value and not isinstance(value, list):

                        if not os.path.exists(ctrl.get_selection()) and value != 'On/Off':
                            display(
                                win, "%s field contains incorrect path. Please update the path!" % ctrl.get_name())
                            return

                config_list.append(ctrl)


        # Get the user's CPAC pipeline name for use in this script
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

        dlg = wx.FileDialog(
            self, message="Save CPAC configuration file as ...", defaultDir=os.getcwd(),
            defaultFile=("pipeline_config_%s" % pipeline_name), wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml", style=wx.SAVE)


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

        try:
            f = open(path, 'w')

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
                        print 'Using ReHo cluster size of ', value
                    elif item_type == 0:
                        value = sample_list.index(value)
                    else:
                        if substitution_map.get(value) != None:
                            value = substitution_map.get(value)
                        elif value != 'None':
                            value = ast.literal_eval(str(value))
                    
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
                    print>>f,"\n"
                

                # parameters that are bracketed numbers (int or float)
                elif dtype == 5:

                    ### parse user input   ### can't use internal function type() here???
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
                        lvalue = [int(item) for item in lvalue]
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
                            print "Exception trying to translate: %s, %s, %s, %s"%(label,str(value),val,err)
                            print "value type: %s"%(type(val))
                        values.append(ast.literal_eval(str(val)))

                    print>>f, label, ":", values
                    print>>f, "\n"


                # parameters that are whole words
                #     ALSO: the Nuisance Corrections lists                
                elif dtype == 8:

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


                else:

                    value = ast.literal_eval(str(value))

                    print>>f, label, ":", value
                    print>>f, "\n"


            f.close()

        except Exception, e:
            print e
            print "Error Writing the pipeline configuration file %s" % path
            raise


