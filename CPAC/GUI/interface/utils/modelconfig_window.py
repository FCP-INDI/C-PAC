import wx
import generic_class
from .constants import control, dtype, substitution_map
import os
import yaml

import modelDesign_window


ID_RUN = 11


class ModelConfig(wx.Frame):

    # this creates the wx.Frame mentioned above in the class declaration
    def __init__(self, parent, gpa_settings=None):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - Create New FSL Model", size=(900, 650))


        if gpa_settings == None:
            self.gpa_settings = {}
            self.gpa_settings['subject_list'] = ''
            self.gpa_settings['pheno_file'] = ''
            self.gpa_settings['subject_id_label'] = ''
            self.gpa_settings['design_formula'] = ''
            self.gpa_settings['mean_mask'] = ''
            self.gpa_settings['custom_roi_mask'] = 'None'
            self.gpa_settings['coding_scheme'] = ''
            self.gpa_settings['use_zscore'] = True
            self.gpa_settings['derivative_list'] = ''
            self.gpa_settings['repeated_measures'] = ''
            self.gpa_settings['group_sep'] = ''
            self.gpa_settings['grouping_var'] = 'None'
            self.gpa_settings['z_threshold'] = ''
            self.gpa_settings['p_threshold'] = ''
        else:
            self.gpa_settings = gpa_settings
        

        self.parent = parent

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        vertSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)

        self.window = wx.ScrolledWindow(self.panel, size=(-1,300))
        


        self.page = generic_class.GenericClass(self.window, " FSL Model Setup")

        self.page.add(label="Subject List ",
                      control=control.COMBO_BOX,
                      name="subject_list",
                      type=dtype.STR,
                      comment="Full path to a list of subjects to be included in the model.\n\nThis should be a text file with one subject per line.\n\nTip 1: A list in this format contaning all subjects run through CPAC was generated along with the main CPAC subject list (see subject_list_group_analysis.txt).\n\nTIp 2: An easy way to manually create this file is to copy the subjects column from your Regressor/EV spreadsheet.",
                      values=self.gpa_settings['subject_list'])

        self.page.add(label="Phenotype/EV File ",
                      control=control.COMBO_BOX,
                      name="pheno_file",
                      type=dtype.STR,
                      comment="Full path to a .csv file containing EV information for each subject.\n\nTip: A file in this format (containing a single column listing all subjects run through CPAC) was generated along with the main CPAC subject list (see template_phenotypic.csv).",
                      values=self.gpa_settings['pheno_file'])

        self.page.add(label="Subjects Column Name ",
                      control=control.TEXT_BOX,
                      name="subject_id_label",
                      type=dtype.STR,
                      comment="Name of the subjects column in your EV file.",
                      values=self.gpa_settings['subject_id_label'],
                      style=wx.EXPAND | wx.ALL,
                      size=(160, -1))

        
        load_panel_sizer = wx.BoxSizer(wx.HORIZONTAL)
        load_pheno_btn = wx.Button(self.window, 2, 'Load Phenotype File', (220,10), wx.DefaultSize, 0)
        load_panel_sizer.Add(load_pheno_btn)

        self.Bind(wx.EVT_BUTTON, self.populateEVs, id=2)


        self.page.add_pheno_load_panel(load_panel_sizer)



        # experimental checkbox row stuff
        self.page.add(label = "Model Setup ",
                      control = control.CHECKBOX_GRID,
                      name = "model_setup",
                      type = 9,#dtype.LBOOL,
                      values = '',
                      comment="A list of EVs from your phenotype file will populate in this window. From here, you can select whether the EVs should be treated as categorical or if they should be demeaned (continuous/non-categorical EVs only). 'MeanFD', 'MeanFD_Jenkinson', 'Measure Mean', and 'Custom_ROI_Mean' will also appear in this window automatically as options to be used as regressors that can be included in your model design. Note that the MeanFD and mean of measure values are automatically calculated and supplied by C-PAC via individual-level analysis.",
                      size = (450, -1))

        self.page.add(label="Design Matrix Formula ",
                      control=control.TEXT_BOX,
                      name="design_formula",
                      type=dtype.STR,
                      comment="Specify the formula to describe your model design. Essentially, including EVs in this formula inserts them into the model. The most basic format to include each EV you select would be 'EV + EV + EV + ..', etc. You can also select to include MeanFD, MeanFD_Jenkinson, Measure_Mean, and Custom_ROI_Mean here. See the C-PAC User Guide for more detailed information regarding formatting your design formula.",
                      values= self.gpa_settings['design_formula'],
                      size=(450, -1))

        self.page.add(label="Measure Mean Generation ", 
                 control=control.CHOICE_BOX, 
                 name='mean_mask', 
                 type=dtype.LSTR, 
                 comment = "Choose whether to use a group mask or individual-specific mask when calculating the output means to be used as a regressor.\n\nThis only takes effect if you include the 'Measure_Mean' regressor in your Design Matrix Formula.", 
                 values=["Group Mask","Individual Mask"])

        self.page.add(label="Custom ROI Mean Mask ",
                      control=control.COMBO_BOX,
                      name="custom_roi_mask",
                      type=dtype.STR,
                      comment="Optional: Full path to a NIFTI file containing one or more ROI masks. The means of the masked regions will then be computed for each subject's output and will be included in the model as regressors (one for each ROI in the mask file) if you include 'Custom_ROI_Mean' in the Design Matrix Formula.",
                      values=self.gpa_settings['custom_roi_mask'])

        self.page.add(label="Use z-score Standardized Derivatives ", 
                     control=control.CHOICE_BOX, 
                     name='use_zscore', 
                     type=dtype.BOOL, 
                     comment="Run the group analysis model on the z-score " \
                             "standardized version of the derivatives you " \
                             "choose in the list below.",
                     values=["True","False"])

        self.page.add(label = "Select Derivatives ",
                    control = control.CHECKLIST_BOX,
                    name = "derivative_list",
                    type = dtype.LSTR,
                    values = ['ALFF',
                              'ALFF (smoothed)',
                              'f/ALFF',
                              'f/ALFF (smoothed)',
                              'ReHo',
                              'ReHo (smoothed)',
                              'ROI Average SCA',
                              'ROI Average SCA (smoothed)',
                              'Voxelwise SCA',
                              'Voxelwise SCA (smoothed)',
                              'Dual Regression',
                              'Dual Regression (smoothed)',
                              'Multiple Regression SCA',
                              'Multiple Regression SCA (smoothed)',
                              'Network Centrality',
                              'Network Centrality (smoothed)',
                              'VMHC (z-score std only)',
                              'VMHC z-stat (z-score std only)'],
                    comment = "Select which derivatives you would like to include when running group analysis.\n\nWhen including Dual Regression, make sure to correct your P-value for the number of maps you are comparing.\n\nWhen including Multiple Regression SCA, you must have more degrees of freedom (subjects) than there were time series.",
                    size = (350,160))

        self.page.add(label="Coding Scheme ", 
                     control=control.CHOICE_BOX, 
                     name="coding_scheme", 
                     type=dtype.LSTR, 
                     comment="Choose the coding scheme to use when generating your model. 'Treatment' encoding is generally considered the typical scheme. Consult the User Guide for more information.", 
                     values=["Treatment", "Sum"])
                     
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

        self.page.add(label="Run Repeated Measures ", 
                     control=control.CHOICE_BOX, 
                     name='repeated_measures', 
                     type=dtype.BOOL, 
                     comment="Run repeated measures to compare different " \
                             "scans (must use the group analysis subject " \
                             "list and phenotypic file formatted for " \
                             "repeated measures.", 
                     values=["False","True"])
        
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
        
        load = wx.Button(btnPanel, wx.ID_ADD, "Load Settings", (
            200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.load, id=wx.ID_ADD)
        hbox.Add(load, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)
        
        next = wx.Button(btnPanel, 3, "Next >", (200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.load_next_stage, id=3)
        hbox.Add(next, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        # reminder: functions bound to buttons require arguments
        #           (self, event)
        

        btnPanel.SetSizer(hbox)



        #text_sizer = wx.BoxSizer(wx.HORIZONTAL)
        #measure_text = wx.StaticText(self.window, label='Note: Regressor options \'MeanFD\' and \'Measure_Mean\' are automatically demeaned prior to being inserted into the model.')
        #text_sizer.Add(measure_text)

        #mainSizer.Add(text_sizer)


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

            phenoHeaderString = phenoFile.readline().rstrip('\r\n')
            phenoHeaderItems = phenoHeaderString.split(',')
            phenoHeaderItems.remove(self.gpa_settings['subject_id_label'])

            # update the 'Model Setup' box and populate it with the EVs and
            # their associated checkboxes for categorical and demean
            for ctrl in self.page.get_ctrl_list():

                name = ctrl.get_name()

                if name == 'model_setup':
                    ctrl.set_value(phenoHeaderItems)
                    ctrl.set_selection(self.gpa_settings['ev_selections'])

                if name == 'coding_scheme':
                    ctrl.set_value(self.gpa_settings['coding_scheme'])

                if name == 'mean_mask':
                    ctrl.set_value(self.gpa_settings['mean_mask'])

                if name == 'repeated_measures':
                    ctrl.set_value(self.gpa_settings['repeated_measures'])

                if name == 'z_threshold':
                    ctrl.set_value(self.gpa_settings['z_threshold'][0])

                if name == 'p_threshold':
                    ctrl.set_value(self.gpa_settings['p_threshold'])

                if name == 'use_zscore':
                    ctrl.set_value(self.gpa_settings['use_zscore'])
                    
                if name == 'group_sep':
                    ctrl.set_value(self.gpa_settings['group_sep'])

                if name == 'grouping_var':
                    ctrl.set_value(self.gpa_settings['grouping_var'])

                if name == 'derivative_list':

                    value = self.gpa_settings['derivative_list']

                    if isinstance(value, str):
                        value = value.replace("['","").replace("']","").split("', '")

                    new_derlist = []

                    # remove the _z if they are there, just so it can
                    # repopulate the listbox through the substitution map
                    for val in value:
                        if "_z" in val:
                            val = val.replace("_z","")
                            new_derlist.append(val)
                        else:
                            new_derlist.append(val)                           

                    ctrl.set_value(new_derlist)



    def cancel(self, event):
        self.Close()

    def display(self, win, msg):
        wx.MessageBox(msg, "Error")
        win.SetBackgroundColour("pink")
        win.SetFocus()
        win.Refresh()
        raise ValueError


    def load_pheno(self,event):
        pass

    

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

            config_map = yaml.load(open(path, 'r'))
            s_map = dict((v, k) for k, v in substitution_map.iteritems())

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

                phenoHeaderString = phenoFile.readline().rstrip('\r\n')
                phenoHeaderItems = phenoHeaderString.split(',')
                phenoHeaderItems.remove(self.gpa_settings['subject_id_label'])

                # update the 'Model Setup' box and populate it with the EVs and
                # their associated checkboxes for categorical and demean
                for ctrl in self.page.get_ctrl_list():

                    if ctrl.get_name() == 'model_setup':
                        ctrl.set_value(phenoHeaderItems)
                        ctrl.set_selection(self.gpa_settings['ev_selections'])


            # populate the rest of the controls
            for ctrl in self.page.get_ctrl_list():

                name = ctrl.get_name()
                value = config_map.get(name)
                dtype = ctrl.get_datatype()

                # the model setup checkbox grid is the only one that doesn't
                # get repopulated the standard way. instead it is repopulated
                # by the code directly above

                if name == 'derivative_list':
                    value = [s_map.get(item)
                                 for item in value if s_map.get(item) != None]
                    if not value:
                        value = [str(item) for item in value]
                    
                    new_derlist = []

                    for val in value:
                        if "_z" in val:
                            val = val.replace("_z","")
                            new_derlist.append(val)
                        else:
                            new_derlist.append(val)

                    ctrl.set_value(new_derlist)

                elif name == 'repeated_measures' or name == 'use_zscore':
                    ctrl.set_value(str(value))

                elif name == 'z_threshold' or name == 'p_threshold':
                    value = value[0]
                    ctrl.set_value(value)
                    
                elif name == 'group_sep':
                    value = s_map.get(value)
                    ctrl.set_value(value)                

                elif name != 'model_setup' and name != 'derivative_list':
                    ctrl.set_value(value)
                

            dlg.Destroy()




    def read_phenotypic(self, pheno_file, ev_selections):

        import csv

        ph = pheno_file

        # Read in the phenotypic CSV file into a dictionary named pheno_dict
        # while preserving the header fields as they correspond to the data
        p_reader = csv.DictReader(open(os.path.abspath(ph), 'rU'), skipinitialspace=True)

        #pheno_dict_list = []
        
        # dictionary to store the data in a format Patsy can use
        # i.e. a dictionary where each header is a key, and the value is a
        # list of all of that header's values
        pheno_data_dict = {}

        for line in p_reader:

            for key in line.keys():

                if key not in pheno_data_dict.keys():
                    pheno_data_dict[key] = []

                # create a list within one of the dictionary values for that
                # EV if it is categorical; formats this list into a form
                # Patsy can understand regarding categoricals:
                #     example: { ADHD: ['adhd1', 'adhd1', 'adhd2', 'adhd1'] }
                #                instead of just [1, 1, 2, 1], etc.
                if 'categorical' in ev_selections.keys():
                    if key in ev_selections['categorical']:
                        pheno_data_dict[key].append(key + str(line[key]))

                    else:
                        pheno_data_dict[key].append(line[key])

                else:
                    pheno_data_dict[key].append(line[key])

   
            #pheno_dict_list.append(line)
        
            # pheno_dict_list is a list of dictionaries of phenotype header items
            # matched to their values, which also includes subject IDs
            
            # i.e. [{'header1': 'value', 'header2': 'value'}, {'header1': 'value', 'header2': 'value'}, ..]
            
            # these dictionaries are UNORDERED, i.e. header items ARE NOT ORDERED


        return pheno_data_dict


          
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
                    self, 'Error reading file - either it does not exist or you' \
                          ' do not have read access. \n\n' \
                          'Parameter: %s' % paramName,
                    'File Access Error',
                    wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()
                raise Exception
                
        
        testFile(self.gpa_settings['subject_list'], 'Subject List')
        testFile(self.gpa_settings['pheno_file'], 'Phenotype/EV File')

        subFile = open(os.path.abspath(self.gpa_settings['subject_list']))
        phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']),"rU")

        phenoHeaderString = phenoFile.readline().rstrip('\r\n')
        self.phenoHeaderItems = phenoHeaderString.split(',')

        if self.gpa_settings['subject_id_label'] in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(self.gpa_settings['subject_id_label'])
        else:
            errSubID = wx.MessageDialog(
                self, 'Please enter the name of the subject ID column' \
                ' as it is labeled in the phenotype file.',
                'Blank/Incorrect Subject Header Input',
                wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
            raise Exception
            
            
        # some more checks
        sub_IDs = subFile.readlines()
        self.subs = []
        
        for sub in sub_IDs:
            self.subs.append(sub.rstrip("\n"))        
        
        pheno_rows = phenoFile.readlines()
        
        for row in pheno_rows:
        
            # check if the pheno file produces any rows such as ",,,,," due
            # to odd file formatting issues. if so, ignore this row. if there
            # are values present in the row, continue as normal
            if ",," not in row:
                    
                # if it finds a sub from the subject list in the current row
                # taken from the pheno, move on. if it goes through the entire
                # subject list and never finds a match, kick off the "else"
                # clause below containing the error message
                for sub in self.subs:
                
                    # for repeated measures-formatted files
                    if "," in sub:
                    
                        # make the comma separator an underscore to match the
                        # repeated measures-formatted pheno file
                        if sub.replace(",","_") in row:
                            break
                            
                    # for normal
                    else:
                
                        if sub in row:
                            break
                            
                else:
                    errSubID = wx.MessageDialog(
                        self, "Your phenotype file contains a subject ID " \
                        "that is not present in your group analysis " \
                        "subject list.\n\nPhenotype file row with subject " \
                        "ID not in subject list:\n%s" \
                        % row,
                        "Subject Not In List",
                        wx.OK | wx.ICON_ERROR)
                    errSubID.ShowModal()
                    errSubID.Destroy()
                    raise Exception


        for ctrl in self.page.get_ctrl_list():

            # update the 'Model Setup' box and populate it with the EVs and
            # their associated checkboxes for categorical and demean
            if ctrl.get_name() == 'model_setup':
                ctrl.set_value(self.phenoHeaderItems)

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



    ''' button: NEXT '''
    def load_next_stage(self, event):

        import patsy

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
                    self, 'Error reading file - either it does not exist ' \
                          'or you do not have read access. \n\n' \
                          'Parameter: %s' % paramName,
                    'File Access Error',
                    wx.OK | wx.ICON_ERROR)
                errDlgFileTest.ShowModal()
                errDlgFileTest.Destroy()
                raise Exception
                
        
        testFile(self.gpa_settings['subject_list'], 'Subject List')
        testFile(self.gpa_settings['pheno_file'], 'Phenotype/EV File')

     
        phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']),"rU")

        phenoHeaderString = phenoFile.readline().rstrip('\r\n')
        self.phenoHeaderItems = phenoHeaderString.split(',')

        if self.gpa_settings['subject_id_label'] in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(self.gpa_settings['subject_id_label'])
        else:
            errSubID = wx.MessageDialog(
                self, 'Please enter the name of the subject ID column' \
                ' as it is labeled in the phenotype file.',
                'Blank/Incorrect Subject Header Input',
                wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
            raise Exception



        for ctrl in self.page.get_ctrl_list():
            
            name = ctrl.get_name()

            # get the design matrix formula
            if name == 'design_formula':

                self.gpa_settings['design_formula'] = str(ctrl.get_selection())

            # get the EV categorical + demean grid selections
            elif name == 'model_setup':

                # basically, ctrl is checkbox_grid in this case, and
                # get_selection goes to generic_class.py first, which links
                # it to the custom GetGridSelection() function in the
                # checkbox_grid class in custom_control.py
                self.gpa_settings['ev_selections'] = ctrl.get_selection()
                
            elif name == 'group_sep':

                self.gpa_settings['group_sep'] = ctrl.get_selection()


            elif name == 'grouping_var':

                self.gpa_settings['grouping_var'] = ctrl.get_selection()



            if name == 'derivative_list':

                # grab this for below
                derlist_ctrl = ctrl

            else:
                self.gpa_settings[name] = str(ctrl.get_selection())



        self.gpa_settings['derivative_list'] = []

        for derivative in list(derlist_ctrl.get_selection()):
            if self.gpa_settings['use_zscore'] == "True":
                
                self.gpa_settings['derivative_list'].append(derivative + "_z")
                
            else:
            
                self.gpa_settings['derivative_list'].append(derivative)





        self.pheno_data_dict = self.read_phenotypic(self.gpa_settings['pheno_file'], self.gpa_settings['ev_selections'])



        try:
            phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']))
        except:
            print '\n\n[!] CPAC says: The phenotype file path provided ' \
                    'couldn\'t be opened - either it does not exist or ' \
                    'there are access restrictions.\n'
            print 'Phenotype file provided: '
            print self.gpa_settings['pheno_file'], '\n\n'
            raise IOError



        # validate design formula and build Available Contrasts list
        var_list_for_contrasts = []
        EVs_to_test = []
        EVs_to_include = []

        # take the user-provided design formula and break down the included
        # terms into a list, and use this to create the list of available
        # contrasts
                

        formula = self.gpa_settings['design_formula']
            

        # need to cycle through the EVs inside parentheses just to make
        # sure they are valid

        # THEN you have to treat the entire parentheses thing as one EV when
        # it comes to including it in the list for contrasts

        formula_strip = formula.replace('+',' ')
        formula_strip = formula_strip.replace('-',' ')
        formula_strip = formula_strip.replace('**(', '**')
        formula_strip = formula_strip.replace(')**', '**')
        formula_strip = formula_strip.replace('(',' ')
        formula_strip = formula_strip.replace(')',' ')
        EVs_to_test = formula_strip.split()



        # ensure the design formula only has valid EVs in it
        for EV in EVs_to_test:

            # ensure ** interactions have a valid EV on one side and a number
            # on the other
            if '**' in EV:

                both_sides = EV.split('**')

                int_check = 0

                for side in both_sides:

                    if side.isdigit():
                        int_check = 1
                    else:
                        if (side not in self.pheno_data_dict.keys()) and \
                            side != 'MeanFD' and side != 'MeanFD_Jenkinson' \
                            and side != 'Measure_Mean' and \
                            side != 'Custom_ROI_Mean':

                            errmsg = 'CPAC says: The regressor \'%s\' you ' \
                                     'entered within the design formula as ' \
                                     'part of the interaction \'%s\' is not ' \
                                     'a valid EV option.\n\nPlease enter ' \
                                     'only the EVs in your phenotype file ' \
                                     'or the MeanFD, MeanFD_Jenkinson, ' \
                                     'Custom_ROI_Mean, or Measure_Mean ' \
                                     'options.' \
                                     % (side,EV)

                            errSubID = wx.MessageDialog(self, errmsg,
                                'Invalid EV', wx.OK | wx.ICON_ERROR)
                            errSubID.ShowModal()
                            errSubID.Destroy()
                
                            raise Exception


                if int_check != 1:

                    errmsg = 'CPAC says: The interaction \'%s\' you ' \
                             'entered within the design formula requires ' \
                             'a number on one side.\n\nExample: ' \
                             '(EV1 + EV2 + EV3)**3\n\nNote: This would be ' \
                             'equivalent to\n(EV1 + EV2 + EV3) * ' \
                             '(EV1 + EV2 + EV3) * (EV1 + EV2 + EV3)' % EV

                    errSubID = wx.MessageDialog(self, errmsg,
                        'Invalid EV', wx.OK | wx.ICON_ERROR)
                    errSubID.ShowModal()
                    errSubID.Destroy()
                
                    raise Exception


                    
            # ensure these interactions are input correctly
            elif (':' in EV) or ('/' in EV) or ('*' in EV):

                if ':' in EV:
                    both_EVs_in_interaction = EV.split(':')

                if '/' in EV:
                    both_EVs_in_interaction = EV.split('/')

                if '*' in EV:
                    both_EVs_in_interaction = EV.split('*')


                for interaction_EV in both_EVs_in_interaction:

                    if (interaction_EV not in self.pheno_data_dict.keys()) and \
                        interaction_EV != 'MeanFD' and \
                        interaction_EV != 'MeanFD_Jenkinson' and \
                        interaction_EV != 'Measure_Mean' and \
                        interaction_EV != 'Custom_ROI_Mean':

                        errmsg = 'CPAC says: The regressor \'%s\' you ' \
                                 'entered within the design formula as ' \
                                 'part of the interaction \'%s\' is not a ' \
                                 'valid EV option.\n\nPlease enter only ' \
                                 'the EVs in your phenotype file or the ' \
                                 'MeanFD, MeanFD_Jenkinson, Custom_ROI_' \
                                 'Mean, or Measure_Mean options.' \
                                 % (interaction_EV,EV)

                        errSubID = wx.MessageDialog(self, errmsg,
                            'Invalid EV', wx.OK | wx.ICON_ERROR)
                        errSubID.ShowModal()
                        errSubID.Destroy()
                
                        raise Exception    

            else:

                if (EV not in self.pheno_data_dict.keys()) and EV != 'MeanFD' \
                    and EV != 'MeanFD_Jenkinson' and EV != 'Measure_Mean' \
                    and EV != 'Custom_ROI_Mean':

                    errmsg = 'CPAC says: The regressor \'%s\' you ' \
                             'entered within the design formula is not ' \
                             'a valid EV option.' \
                             '\n\nPlease enter only the EVs in your phenotype ' \
                             'file or the MeanFD, MeanFD_Jenkinson, ' \
                             'Custom_ROI_Mean, or Measure_Mean options.' \
                             % EV

                    errSubID = wx.MessageDialog(self, errmsg,
                        'Invalid EV', wx.OK | wx.ICON_ERROR)
                    errSubID.ShowModal()
                    errSubID.Destroy()
                
                    raise Exception





        ''' design formula/input parameters checks '''

        if "Custom_ROI_Mean" in formula and \
            (self.gpa_settings['custom_roi_mask'] == None or \
            self.gpa_settings['custom_roi_mask'] == ""):

            err_string = "You included 'Custom_ROI_Mean' as a regressor " \
                         "in your Design Matrix Formula, but you did not " \
                         "specify a Custom ROI Mean Mask file.\n\nPlease " \
                         "either specify a mask file, or remove " \
                         "'Custom_ROI_Mean' from your model."

            errSubID = wx.MessageDialog(self, err_string,
                'No Custom ROI Mean Mask File', wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()

            raise Exception



        if "Custom_ROI_Mean" not in formula and \
            (self.gpa_settings['custom_roi_mask'] != None and \
            self.gpa_settings['custom_roi_mask'] != "" and \
            self.gpa_settings['custom_roi_mask'] != "None" and \
            self.gpa_settings['custom_roi_mask'] != "none"):

            warn_string = "Note: You specified a Custom ROI Mean Mask file, " \
                          "but you did not include 'Custom_ROI_Mean' as a " \
                          "regressor in your Design Matrix Formula.\n\nThe " \
                          "means of the ROIs specified in the file will not " \
                          "be included as regressors unless you include " \
                          "'Custom_ROI_Mean' in your model."

            errSubID = wx.MessageDialog(self, warn_string,
                'No Custom_ROI_Mean Regressor', wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()

            raise Exception
            
            
        if str(self.gpa_settings["use_zscore"]) == "True":
        
            if "Measure_Mean" in formula:
            
                warn_string = "Note: You have included Measure_Mean as a " \
                    "regressor in your model, but you have selected to run " \
                    "the group-level analysis with the z-score standardized "\
                    "version of the outputs.\n\nThe mean of any z-score " \
                    "standardized output will always be zero."

                errSubID = wx.MessageDialog(self, warn_string,
                    'Measure_Mean Included With z-scored Outputs', wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                raise Exception
        
        else:
        
            for deriv in self.gpa_settings["derivative_list"]:
            
                if "VMHC" in deriv:
                
                    warn_string = "Note: You have selected to run group-" \
                        "level analysis using raw outputs (non-z-score " \
                        "standardized), but you have also included VMHC " \
                        "as one of the outputs to include in your model."

                    errSubID = wx.MessageDialog(self, warn_string,
                        'VMHC Cannot Be Included As Raw Output', wx.OK | wx.ICON_ERROR)
                    errSubID.ShowModal()
                    errSubID.Destroy()

                    raise Exception



        # if there is a custom ROI mean mask file provided, and the user
        # includes it as a regressor in their design matrix formula, calculate
        # the number of ROIs in the file and generate the column names so that
        # they can be passed as possible contrast labels

        if "Custom_ROI_Mean" in formula and \
            (self.gpa_settings['custom_roi_mask'] != None and \
            self.gpa_settings['custom_roi_mask'] != "" and \
            self.gpa_settings['custom_roi_mask'] != "None" and \
            self.gpa_settings['custom_roi_mask'] != "none"):

            import commands

            try:
                ROIstats_output = commands.getoutput("3dROIstats -mask %s %s" \
                                  % (self.gpa_settings['custom_roi_mask'], \
                                  self.gpa_settings['custom_roi_mask']))
            except Exception as e:
                print "[!] CPAC says: AFNI 3dROIstats failed for custom ROI" \
                      "Mean Mask file validation. Please ensure you either " \
                      "have AFNI installed and that you created the mask " \
                      "file properly. Consult the User Guide for more " \
                      "information.\n\n"
                print "Error details: %s\n\n" % e
                raise

            ROIstats_list = ROIstats_output.split("\t")

            # calculate the number of ROIs - 3dROIstats output can be split
            # into a list, and the actual ROI means begin at a certain point
            num_rois = (len(ROIstats_list)-3)/2


            custom_roi_labels = []

            for num in range(0,num_rois):
                custom_roi_labels.append("Custom_ROI_Mean_%d" % int(num+1))
                
                
        
        if str(self.gpa_settings["group_sep"]) == "On":
        
            if (self.gpa_settings["grouping_var"] == "None") or \
                (self.gpa_settings["grouping_var"] is None) or \
                (self.gpa_settings["grouping_var"] == "none"):
                
                warn_string = "Note: You have selected to model group " \
                    "variances separately, but you have not specified a " \
                    "grouping variable."

                errSubID = wx.MessageDialog(self, warn_string,
                    'No Grouping Variable Specified', wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                raise Exception
        
            if self.gpa_settings["grouping_var"] not in formula:           
                
                warn_string = "Note: You have specified '%s' as your " \
                    "grouping variable for modeling the group variances " \
                    "separately, but you have not included this variable " \
                    "in your design formula.\n\nPlease include this " \
                    "variable in your design, or choose a different " \
                    "grouping variable." % self.gpa_settings["grouping_var"]

                errSubID = wx.MessageDialog(self, warn_string,
                    'Grouping Variable not in Design', wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                raise Exception




        def read_phenotypic(pheno_file, ev_selections, subject_id_label):

            import csv
            import numpy as np

            ph = pheno_file

            # Read in the phenotypic CSV file into a dictionary named pheno_dict
            # while preserving the header fields as they correspond to the data
            p_reader = csv.DictReader(open(os.path.abspath(ph), 'rU'), skipinitialspace=True)

            # dictionary to store the data in a format Patsy can use
            # i.e. a dictionary where each header is a key, and the value is a
            # list of all of that header's values
            pheno_data_dict = {}

            for line in p_reader:

                # here, each instance of 'line' is really a dictionary where the
                # keys are the pheno headers, and their values are the values of
                # each EV for that one subject - each iteration of this loop is
                # one subject

                for key in line.keys():

                    if key not in pheno_data_dict.keys():
                        pheno_data_dict[key] = []

                    # create a list within one of the dictionary values for that
                    # EV if it is categorical; formats this list into a form
                    # Patsy can understand regarding categoricals:
                    #     example: { ADHD: ['adhd1', 'adhd1', 'adhd0', 'adhd1'] }
                    #                instead of just [1, 1, 0, 1], etc.
                    if 'categorical' in ev_selections.keys():
                        if key in ev_selections['categorical']:
                            pheno_data_dict[key].append(key + str(line[key]))

                        elif key == subject_id_label:
                            pheno_data_dict[key].append(line[key])

                        else:
                            pheno_data_dict[key].append(float(line[key]))

                    elif key == subject_id_label:
                        pheno_data_dict[key].append(line[key])

                    else:
                        pheno_data_dict[key].append(float(line[key]))



            # this needs to run after each list in each key has been fully
            # populated above
            for key in pheno_data_dict.keys():

                # demean the EVs marked for demeaning
                if 'demean' in ev_selections.keys():
                    if key in ev_selections['demean']:

                        new_demeaned_evs = []

                        mean_evs = 0.0

                        # populate a dictionary, a key for each demeanable EV, with
                        # the value being the sum of all the values (which need to be
                        # converted to float first)
                        for val in pheno_data_dict[key]:
                            mean_evs += float(val)

                        # calculate the mean of the current EV in this loop
                        mean_evs = mean_evs / len(pheno_data_dict[key])

                        # remove the EV's mean from each value of this EV
                        # (demean it!)
                        for val in pheno_data_dict[key]:
                            new_demeaned_evs.append(float(val) - mean_evs)

                        # replace
                        pheno_data_dict[key] = new_demeaned_evs


                # converts non-categorical EV lists into NumPy arrays
                # so that Patsy may read them in properly
                if 'categorical' in ev_selections.keys():
                    if key not in ev_selections['categorical']:
            
                        pheno_data_dict[key] = np.array(pheno_data_dict[key])


            return pheno_data_dict


        patsy_formatted_pheno = read_phenotypic(self.gpa_settings['pheno_file'], self.gpa_settings['ev_selections'], self.gpa_settings['subject_id_label'])


        # let's create dummy columns for MeanFD, Measure_Mean, and
        # Custom_ROI_Mask (if included in the Design Matrix Formula) just so we
        # can get an accurate list of EVs Patsy will generate

        def create_regressor_column(regressor):

            # regressor should be a string of the name of the regressor

            import numpy as np

            regressor_list = []

            for key in patsy_formatted_pheno.keys():
               for val in patsy_formatted_pheno[key]:
                   regressor_list.append(0.0)
               break

            regressor_list = np.array(regressor_list)

            patsy_formatted_pheno[regressor] = regressor_list


        if 'MeanFD' in formula:
            create_regressor_column('MeanFD')
        if 'MeanFD_Jenkinson' in formula:
            create_regressor_column('MeanFD_Jenkinson')
        if 'Measure_Mean' in formula:
            create_regressor_column('Measure_Mean')

        if 'Custom_ROI_Mean' in formula:

            add_formula_string = ""

            for col_label in custom_roi_labels:

                create_regressor_column(col_label)

                # create a string of all the new custom ROI regressor column
                # names to be inserted into the design formula, so that Patsy
                # will accept the phenotypic data dictionary that now has these
                # columns
                if add_formula_string == "":
                    add_formula_string = add_formula_string + col_label
                else:
                    add_formula_string = add_formula_string + " + " + col_label
   

            formula = formula.replace("Custom_ROI_Mean",add_formula_string)   



        if 'categorical' in self.gpa_settings['ev_selections']:
            for EV_name in self.gpa_settings['ev_selections']['categorical']:

                if self.gpa_settings['coding_scheme'] == 'Treatment':
                    formula = formula.replace(EV_name, 'C(' + EV_name + ')')
                elif self.gpa_settings['coding_scheme'] == 'Sum':
                    formula = formula.replace(EV_name, 'C(' + EV_name + ', Sum)')



        # create the dmatrix in Patsy just to see what the design matrix
        # columns are going to be 
        try:
            dmatrix = patsy.dmatrix(formula, patsy_formatted_pheno)
        except:
            print '\n\n[!] CPAC says: Design matrix creation wasn\'t ' \
                    'successful - do the terms in your formula correctly ' \
                    'correspond to the EVs listed in your phenotype file?\n'
            print 'Phenotype file provided: '
            print self.gpa_settings['pheno_file'], '\n\n'
            raise Exception


        column_names = dmatrix.design_info.column_names
        
        
        
        subFile = open(os.path.abspath(self.gpa_settings['subject_list']))

        sub_IDs = subFile.readlines()
        self.subs = []
        
        for sub in sub_IDs:
            self.subs.append(sub.rstrip("\n"))        
        

        
        # check to make sure there are more subjects than EVs!!
        if len(column_names) >= len(self.subs):
            err = "There are more (or an equal amount of) EVs currently " \
                  "included in the model than there are subjects in the " \
                  "group analysis subject list. There must be more " \
                  "subjects than EVs in the design.\n\nNumber of subjects: " \
                  "%d\nNumber of EVs: %d\n\nNote: An 'Intercept' " \
                  "column gets added to the design as an EV, so there will " \
                  "be one more EV than you may have specified in your " \
                  "design." % (len(self.subs),len(column_names))
                  
            errSubID = wx.MessageDialog(self, err,
                    "Too Many EVs or Too Few Subjects",
                    wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
                
            raise Exception
        


        raw_column_strings = []
        
        # remove the header formatting Patsy creates for categorical variables
        # because we are going to use var_list_for_contrasts as a label for
        # users to know what contrasts are available to them
        for column in column_names:

            # if using Sum encoding, a column name may look like this:
            #     C(adhd, Sum)[S.adhd0]

            # this loop leaves it with only "adhd0" in this case, for the
            # contrasts list for the next GUI page

            column_string = column

            string_for_removal = ''

            for char in column_string:

                string_for_removal = string_for_removal + char

                if char == '.':
                    column_string = column_string.replace(string_for_removal, '')
                    string_for_removal = ''

            column_string = column_string.replace(']', '')
            
            if ":" in column_string:
                try:
                    column_string = column_string.split("[")[1]
                except:
                    pass
            
            raw_column_strings.append(column_string)
            
            
        
        if str(self.gpa_settings["group_sep"]) == "On":     

            grouping_options = []
            idx = 0
            
            for column_string in raw_column_strings:

                if self.gpa_settings["grouping_var"] in column_string:

                    grouping_variable_info = []

                    grouping_variable_info.append(column_string)
                    grouping_variable_info.append(idx)

                    grouping_options.append(grouping_variable_info)

                    # grouping_var_idx is the column numbers in the design matrix
                    # which holds the grouping variable (and its possible levels)

                idx += 1               


            # all the categorical values/levels of the grouping variable
            grouping_var_levels = []

            for gv_idx in grouping_options:
            
                for subject in dmatrix:
                
                    if self.gpa_settings["grouping_var"] in self.gpa_settings["ev_selections"]["categorical"]:
                        level_num = str(int(subject[gv_idx[1]]))
                    else:
                        level_num = str(subject[gv_idx[1]])

                    level_label = '__' + self.gpa_settings["grouping_var"] + level_num

                    if level_label not in grouping_var_levels:
                        grouping_var_levels.append(level_label)


            # make the new header for the reorganized data
            for column_string in raw_column_strings:
            
                if column_string != "Intercept":
            
                    if self.gpa_settings["grouping_var"] not in column_string:
                        for level in grouping_var_levels:
                            var_list_for_contrasts.append(column_string + level)
                    elif self.gpa_settings["grouping_var"] in column_string:
                        var_list_for_contrasts.append(column_string)


        else:
        
            for column_string in raw_column_strings:

                if column_string != 'Intercept':
                    var_list_for_contrasts.append(column_string)



        # check for repeated measures file formatting!

        group_sublist_file = open(self.gpa_settings['subject_list'], 'r')

        group_sublist_items = group_sublist_file.readlines()

        group_sublist = [line.rstrip('\n') for line in group_sublist_items \
                          if not (line == '\n') and not line.startswith('#')]

        for ga_sub in group_sublist:

            # ga_sub = subject ID taken off the group analysis subject list

            # let's check to make sure the subject list is formatted for
            # repeated measures properly if repeated measures is enabled
            # and vice versa
            if (self.gpa_settings['repeated_measures'] == "True") and \
                (',' not in ga_sub):

                errmsg = "The group analysis subject list is not in the " \
                         "appropriate format for repeated measures. Please " \
                         "use the appropriate format as described in the " \
                         "CPAC User Guide, or turn off Repeated Measures." \
                         "\n\nNote: CPAC generates a properly-formatted " \
                         "group analysis subject list meant for running " \
                         "repeated measures when you create your original " \
                         "subject list. Look for 'subject_list_group_" \
                         "analysis_repeated_measures.txt' in the directory " \
                         "where you created your subject list."

                errSubID = wx.MessageDialog(self, errmsg,
                    'Subject List Format', wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                raise Exception

            elif (self.gpa_settings['repeated_measures'] == "False") and \
                (',' in ga_sub):

                errmsg = "It looks like your group analysis subject list is " \
                         "formatted for running repeated measures, but " \
                         "'Run Repeated Measures' is not enabled."

                errSubID = wx.MessageDialog(self, errmsg,
                    'Subject List Format', wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                raise Exception


        # make sure the sub IDs in the sublist and pheno files match!

        group_pheno_file = open(self.gpa_settings['pheno_file'], 'r')

        group_pheno_lines = group_pheno_file.readlines()

        # gather the subject IDs from the phenotype file
        def get_pheno_subjects(delimiter):

            for item in group_pheno_lines[0].split(delimiter):
                if item == self.gpa_settings['subject_id_label']:
                    index = group_pheno_lines[0].index(item)

            group_pheno_subs = group_pheno_lines[1:len(group_pheno_lines)]

            pheno_subs = []

            for pheno_sub_line in group_pheno_subs:
                pheno_subs.append(pheno_sub_line.split(delimiter)[index])

            return pheno_subs


        pheno_subs = []

        if "," in group_pheno_lines[0]:
            pheno_subs = get_pheno_subjects(",")

        # now make sure the group sublist and pheno subject IDs match, at least
        # for the ones that exist (i.e. may be less sub IDs in the sublist)
        for sublist_subID, pheno_subID in zip(group_sublist, pheno_subs):

            # if group sublist is formatted for repeated measures
            if "," in sublist_subID:
                sublist_subID = sublist_subID.replace(",","_")

            if sublist_subID != pheno_subID:

                if self.gpa_settings['repeated_measures'] == "False":

                    errmsg = "The subject IDs in your group subject list " \
                             "and your phenotype file do not match. Please " \
                             "make sure these have been set up correctly."

                else:

                    errmsg = "The subject IDs in your group subject list " \
                             "and your phenotype file do not match. Please " \
                             "make sure these have been set up correctly." \
                             "\n\nNote: Repeated measures is enabled - does "\
                             "your phenotype file have properly-formatted " \
                             "subject IDs matching your repeated measures " \
                             "group analysis subject list?"

                errSubID = wx.MessageDialog(self, errmsg,
                    'Subject ID Mismatch', wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                raise Exception



        # open the next window!
        modelDesign_window.ModelDesign(self.parent, self.gpa_settings, var_list_for_contrasts)  # !!! may need to pass the actual dmatrix as well


        self.Close()

        
        

