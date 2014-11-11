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
            self, parent=parent, title="CPAC - Create New FSL Model", size=(850, 650))


        if gpa_settings == None:
            self.gpa_settings = {}
            self.gpa_settings['subject_list'] = ''
            self.gpa_settings['pheno_file'] = ''
            self.gpa_settings['subject_id_label'] = ''
            self.gpa_settings['design_formula'] = ''
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
                      comment="A list of EVs from your phenotype file will populate in this window. From here, you can select whether the EVs should be treated as categorical or if they should be demeaned (continuous/non-categorical EVs only). 'MeanFD' and 'Measure Mean' will also appear in this window automatically as options to be used as regressors that can be included in your model design. Note that the MeanFD and mean of measure values are automatically calculated and supplied by C-PAC via individual-level analysis.",
                      size = (450, -1))

        self.page.add(label="Design Matrix Formula ",
                      control=control.TEXT_BOX,
                      name="design_formula",
                      type=dtype.STR,
                      comment="Specify the formula to describe your model design. Essentially, including EVs in this formula inserts them into the model. The most basic format to include each EV you select would be 'EV + EV + EV + ..', etc. You can also select to include MeanFD and Measure_Mean here. See the C-PAC User Guide for more detailed information regarding formatting your design formula.",
                      values= self.gpa_settings['design_formula'],
                      size=(450, -1))

        self.page.add(label="Coding Scheme ", 
                     control=control.CHOICE_BOX, 
                     name="coding_scheme", 
                     type=dtype.LSTR, 
                     comment="Choose the coding scheme to use when generating your model.", 
                     values=["Treatment", "Sum"])



        self.page.set_sizer()

        


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



        text_sizer = wx.BoxSizer(wx.HORIZONTAL)
        measure_text = wx.StaticText(self.window, label='Note: Regressor options \'MeanFD\' and \'Measure_Mean\' are automatically demeaned prior to being inserted into the model.')
        text_sizer.Add(measure_text)


        mainSizer.Add(text_sizer)

        mainSizer.Add(
            btnPanel, 0.5,  flag=wx.ALIGN_RIGHT | wx.RIGHT, border=20)

        self.panel.SetSizer(mainSizer)

        self.Show()


        # this fires only if we're coming BACK to this page from the second
        # page, and these parameters are already pre-loaded. this is to
        # automatically repopulate the 'Model Setup' checkbox grid
        if self.gpa_settings['pheno_file'] != '':

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
                if name != 'model_setup':
                    ctrl.set_value(str(self.gpa_settings[name]))
                

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


          
            
    def populateEVs(self, event):

        # this runs when the user clicks 'Load Phenotype File'
               
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

     
        phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']))

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

                formula_string = formula_string + ' + MeanFD + Measure_Mean'

                ctrl.set_value(formula_string)




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

     
        phenoFile = open(os.path.abspath(self.gpa_settings['pheno_file']))

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
            if name == 'model_setup':

                # basically, ctrl is checkbox_grid in this case, and
                # get_selection goes to generic_class.py first, which links
                # it to the custom GetGridSelection() function in the
                # checkbox_grid class in custom_control.py
                self.gpa_settings['ev_selections'] = ctrl.get_selection()


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
                            side != 'MeanFD' and side != 'Measure_Mean':

                            errmsg = 'CPAC says: The regressor \'%s\' you ' \
                                     'entered within the design formula as ' \
                                     'part of the interaction \'%s\' is not a ' \
                                     'valid EV option.\n\nPlease enter only ' \
                                     'the EVs in your phenotype file or the ' \
                                     'MeanFD or Measure_Mean options.' \
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
                        interaction_EV != 'MeanFD' and interaction_EV != 'Measure_Mean':

                        errmsg = 'CPAC says: The regressor \'%s\' you ' \
                                 'entered within the design formula as ' \
                                 'part of the interaction \'%s\' is not a ' \
                                 'valid EV option.\n\nPlease enter only ' \
                                 'the EVs in your phenotype file or the ' \
                                 'MeanFD or Measure_Mean options.' \
                                 % (interaction_EV,EV)

                        errSubID = wx.MessageDialog(self, errmsg,
                            'Invalid EV', wx.OK | wx.ICON_ERROR)
                        errSubID.ShowModal()
                        errSubID.Destroy()
                
                        raise Exception    

            else:

                if (EV not in self.pheno_data_dict.keys()) and EV != 'MeanFD' \
                    and EV != 'Measure_Mean':

                    errmsg = 'CPAC says: The regressor \'%s\' you ' \
                             'entered within the design formula is not ' \
                             'a valid EV option.' \
                             '\n\nPlease enter only the EVs in your phenotype ' \
                             'file or the MeanFD or Measure_Mean options.' \
                             % EV

                    errSubID = wx.MessageDialog(self, errmsg,
                        'Invalid EV', wx.OK | wx.ICON_ERROR)
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

        # let's create dummy columns for MeanFD and Measure_Mean just so we
        # can get an accurate list of EVs Patsy will generate
        if 'MeanFD' in formula or 'Measure_Mean' in formula:

            import numpy as np

            MeanFD = []
            Measure_Mean = []

            for key in patsy_formatted_pheno.keys():
               for val in patsy_formatted_pheno[key]:
                   MeanFD.append(0.0)
                   Measure_Mean.append(0.0)
               break

            MeanFD = np.array(MeanFD)
            Measure_Mean = np.array(Measure_Mean)

            patsy_formatted_pheno['MeanFD'] = MeanFD
            patsy_formatted_pheno['Measure_Mean'] = Measure_Mean

        print patsy_formatted_pheno


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
        

        # remove the header formatting Patsy creates for categorical variables
        # because we are going to use var_list_for_contrasts as a label for
        # users to know what contrasts are available to them
        for column in column_names:

            column_string = column

            record = 0
            string_for_removal = ''

            for char in column_string:

                if record == 2:
                    string_for_removal = string_for_removal + char
                elif record == 1 and char == '(':
                    string_for_removal = 'C('
                    record = 2
                else:
                    record = 0
                
                if char == 'C' and record == 0:
                    record = 1

                if char == '.':
                    record = 0
                    column_string = column_string.replace(string_for_removal, '')
                    string_for_removal = ''

            column_string = column_string.replace(']', '')

            if column_string != 'Intercept':
                var_list_for_contrasts.append(column_string)


        
        print 'varlist: ', var_list_for_contrasts

        # open the next window!
        modelDesign_window.ModelDesign(self.parent, self.gpa_settings, var_list_for_contrasts)  # !!! may need to pass the actual dmatrix as well


        self.Close()

        
        

