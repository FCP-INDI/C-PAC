import wx
import generic_class
from .constants import control, dtype, substitution_map
import os
import yaml

import modelDesign_window


ID_RUN = 11


class ModelConfig(wx.Frame):

    # this creates the wx.Frame mentioned above in the class declaration
    def __init__(self, parent, sublist='', pheno='', subID='', setup='', design=''):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - Create New FSL Model", size=(800, 600))

        self.parent = parent

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        vertSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)

        self.window = wx.ScrolledWindow(self.panel, size=(-1,300))
        


        self.page = generic_class.GenericClass(self.window, " FSL Model Setup")

        self.page.add(label="Subject List ",
                      control=control.COMBO_BOX,
                      name="subjectListFile",
                      type=dtype.STR,
                      comment="Full path to a list of subjects to be included in the model.\n\nThis should be a text file with one subject per line.\n\nTip 1: A list in this format contaning all subjects run through CPAC was generated along with the main CPAC subject list (see subject_list_group_analysis.txt).\n\nTIp 2: An easy way to manually create this file is to copy the subjects column from your Regressor/EV spreadsheet.",
                      values=sublist)

        self.page.add(label="Phenotype/EV File ",
                      control=control.COMBO_BOX,
                      name="phenotypicFile",
                      type=dtype.STR,
                      comment="Full path to a .csv file containing EV information for each subject.\n\nTip: A file in this format (containing a single column listing all subjects run through CPAC) was generated along with the main CPAC subject list (see template_phenotypic.csv).",
                      values=pheno)

        self.page.add(label="Subjects Column Name ",
                      control=control.TEXT_BOX,
                      name="subjectColumn",
                      type=dtype.STR,
                      comment="Name of the subjects column in your EV file.",
                      values=subID,
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
                      name = "modelSetup",
                      type = 9,#dtype.LBOOL,
                      values = setup,
                      comment="glob",
                      size = (450, -1))

        self.page.add(label="Design Matrix Formula ",
                      control=control.TEXT_BOX,
                      name="designMatrixFormula",
                      type=dtype.STR,
                      comment="Specify a descriptor for your model.",
                      values= design,
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
        
        '''
        populate = wx.Button(btnPanel, wx.ID_ANY, "Save Model",
            (220, 10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.populateEVs, id=wx.ID_ANY)
        hbox.Add(populate, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)
        '''

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

            for ctrl in self.page.get_ctrl_list():
                name = ctrl.get_name()
                value = config_map.get(name)
                dtype = ctrl.get_datatype()
                if isinstance(value, list):
                    val = None
                    for v in value:
                        if val:
                            val = val + "," + str(v)
                        else:
                            val = str(v)
                else:
                    val = s_map.get(value)
                    if val == None:
                        val = value

                # print "setting value in ctrl name, value -->", name, val
                ctrl.set_value(str(val))

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
                if key in ev_selections['categorical']:
                    pheno_data_dict[key].append(key + str(line[key]))

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
            
            if name == 'subjectListFile':
                self.subjectList = str(ctrl.get_selection())
            
            if name == 'phenotypicFile':
                self.phenoFilePath = str(ctrl.get_selection())
                
            if name == 'subjectColumn':
                self.subjectID = str(ctrl.get_selection())
        
                
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
                
        
        testFile(self.subjectList, 'Subject List')
        testFile(self.phenoFilePath, 'Phenotype/EV File')

     
        phenoFile = open(os.path.abspath(self.phenoFilePath))

        phenoHeaderString = phenoFile.readline().rstrip('\r\n')
        self.phenoHeaderItems = phenoHeaderString.split(',')
        
        if self.subjectID in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(self.subjectID)
        else:
            errSubID = wx.MessageDialog(
                self, 'Please enter the name of the subject ID column' \
                ' as it is labeled in the phenotype file.',
                'Blank/Incorrect Subject Header Input',
                wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
            raise Exception
        

        # update the 'Model Setup' box and populate it with the EVs and their
        # associated checkboxes for include, categorical and demean
        for ctrl in self.page.get_ctrl_list():

            if ctrl.get_name() == 'modelSetup':
                ctrl.set_value(self.phenoHeaderItems)




    # by the time you click next, this should store the design formula, which
    # EVs are categorical and which should be demeaned

    # we should also have the subject list filepath and phenotype filepath,
    # as well as the subject ID column name

    # send this to the next stage, stored nicely, and then run patsy and also
    # the categorical breaking-out




    def load_next_stage(self, event):

        import patsy


        for ctrl in self.page.get_ctrl_list():
            
            name = ctrl.get_name()
            
            # get the design matrix formula
            if name == 'designMatrixFormula':

                self.dmatrix_formula = str(ctrl.get_selection())

            # get the EV categorical + demean grid selections
            if name == 'modelSetup':

                # basically, ctrl is checkbox_grid in this case, and
                # get_selection goes to generic_class.py first, which links
                # it to the custom GetGridSelection() function in the
                # checkbox_grid class in custom_control.py
                self.ev_selections = ctrl.get_selection()


        self.pheno_data_dict = self.read_phenotypic(self.phenoFilePath, self.ev_selections)



        try:
            phenoFile = open(os.path.abspath(self.phenoFilePath))
        except:
            print '\n\n[!] CPAC says: The phenotype file path provided ' \
                    'couldn\'t be opened - either it does not exist or ' \
                    'there are access restrictions.\n'
            print 'Phenotype file provided: '
            print self.phenoFilePath, '\n\n'
            raise IOError

        print '\n\nselections: ', self.ev_selections, '\n\n'
        print '\n\nformula: ', self.dmatrix_formula, '\n\n'
        print '\n\ndata: ', self.pheno_data_dict, '\n\n'


        try:
            dmatrix = patsy.dmatrix(self.dmatrix_formula, self.pheno_data_dict)
        except:
            print '\n\n[!] CPAC says: Design matrix creation wasn\'t ' \
                    'successful - do the terms in your formula correctly ' \
                    'correspond to the EVs listed in your phenotype file?\n'
            print 'Phenotype file provided: '
            print self.phenoFilePath, '\n\n'
            raise Exception


        print '\n\ndmatrix: ', dmatrix, '\n\n'
        print dmatrix.design_info



        # eat the pheno_data_dict, picking out the header items. remove the
        # subject header, then identify which ones are categorical (by seeing
        # what the user selected), and then providing the names of each
        # categorical option

        var_list_for_contrasts = []

        for header in self.pheno_data_dict.keys():

            if header in self.ev_selections['categorical']:
                
                for val in self.pheno_data_dict[header]:
                    if val not in var_list_for_contrasts:
                        var_list_for_contrasts.append(val)

            else:

                if header != self.subjectID:
                    var_list_for_contrasts.append(header)

        

        # open the next window!
        modelDesign_window.ModelDesign(self.parent, self.subjectList, self.phenoFilePath, self.subjectID, var_list_for_contrasts, self.ev_selections, self.dmatrix_formula)   # !!! may need to pass the actual dmatrix as well


        self.Close()

        
        

