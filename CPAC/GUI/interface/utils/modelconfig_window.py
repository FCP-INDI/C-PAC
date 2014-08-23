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
            self, parent=parent, title="CPAC - Create New FSL Model", size=(850, 600))


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
                      comment="glob",
                      size = (450, -1))

        self.page.add(label="Design Matrix Formula ",
                      control=control.TEXT_BOX,
                      name="design_formula",
                      type=dtype.STR,
                      comment="Specify a descriptor for your model.",
                      values= self.gpa_settings['design_formula'],
                      size=(450, -1))



        self.page.set_sizer()

        


        mainSizer.Add(self.window, 1, wx.EXPAND)


        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

#         run = wx.Button(btnPanel, ID_RUN, "Create Model", (
#             280, -1), wx.DefaultSize, 0)
#         self.Bind(wx.EVT_BUTTON, lambda event: self.save(
#             event, 'run'), id=ID_RUN)
#         hbox.Add(run, 0, flag=wx.LEFT | wx.ALIGN_LEFT, border=10)

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
                


                '''
                if isinstance(value, list):
                    val = None
                    for v in value:
                        if val:
                            val = val + "," + str(v)
                        else:
                            val = str(v)
                    print 'value1: ', val
                else:
                    val = s_map.get(value)
                    if val == None:
                        val = value
                    print 'value2: ', val
                '''


                #ctrl.set_value(str(val))

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
        
        # somehow get in the header from the phenotype
        # also get in the subject column
        
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

        '''
        try:
            dmatrix = patsy.dmatrix(self.gpa_settings['design_formula'], self.pheno_data_dict)
        except:
            print '\n\n[!] CPAC says: Design matrix creation wasn\'t ' \
                    'successful - do the terms in your formula correctly ' \
                    'correspond to the EVs listed in your phenotype file?\n'
            print 'Phenotype file provided: '
            print self.gpa_settings['pheno_file'], '\n\n'
            raise Exception
        '''



        # eat the pheno_data_dict, picking out the header items. remove the
        # subject header, then identify which ones are categorical (by seeing
        # what the user selected), and then providing the names of each
        # categorical option

        var_list_for_contrasts = []

        for header in self.pheno_data_dict.keys():

            if 'categorical' in self.gpa_settings['ev_selections'].keys():
                if header in self.gpa_settings['ev_selections']['categorical']:
                
                    for val in self.pheno_data_dict[header]:
                        if val not in var_list_for_contrasts:
                            var_list_for_contrasts.append(val)

                else:

                    if header != self.gpa_settings['subject_id_label']:
                        var_list_for_contrasts.append(header)

            else:

                if header != self.gpa_settings['subject_id_label']:
                    var_list_for_contrasts.append(header)


        if 'MeanFD' in self.gpa_settings['design_formula']:
            var_list_for_contrasts.append('MeanFD')

        if 'Measure_Mean' in self.gpa_settings['design_formula']:
            var_list_for_contrasts.append('Measure_Mean')

        

        # open the next window!
        modelDesign_window.ModelDesign(self.parent, self.gpa_settings, var_list_for_contrasts)  # !!! may need to pass the actual dmatrix as well


        self.Close()

        
        

