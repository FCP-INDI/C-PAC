import wx
import generic_class
from .constants import control, dtype, substitution_map
import os
import ast

import modelconfig_window


ID_RUN = 11


class ModelDesign(wx.Frame):

    def __init__(self, parent, gpa_settings, varlist):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - Create New FSL Model", size=(750, 550))

        self.parent = parent

        self.gpa_settings = gpa_settings

        self.contrasts_list = varlist

        if 'contrasts' not in self.gpa_settings.keys():
            self.gpa_settings['contrasts'] = {}

        if 'grouping_var' not in self.gpa_settings.keys():
            self.gpa_settings['grouping_var'] = 'None'

        if 'model_name' not in self.gpa_settings.keys():
            self.gpa_settings['model_name'] = ''

        if 'output_dir' not in self.gpa_settings.keys():
            self.gpa_settings['output_dir'] = ''


        mainSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)

        self.window = wx.ScrolledWindow(self.panel)

        self.page = generic_class.GenericClass(self.window, " FSL Model Design")
        
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
        


        contrasts_text = 'Available contrasts:\n'

        for con in varlist:
            contrasts_text = contrasts_text + '    ' + con

        varlist_sizer = wx.BoxSizer(wx.HORIZONTAL)
        var_list_text = wx.StaticText(self.window, label=str(contrasts_text))
        varlist_sizer.Add(var_list_text)

        self.page.add_pheno_load_panel(varlist_sizer)
    

        self.page.add(label = 'Contrasts ',
                      control = control.LISTBOX_COMBO,
                      name = 'contrastStrings',
                      type = dtype.LSTR,
                      values = self.gpa_settings['contrasts'],
                      comment = '',
                      size = (300,200),
                      combo_type = 4)

        # this sends the list of available contrast names to the 'Add
        # Contrast' dialog box, so that it may do validation immediately when
        # the user enters contrast strings
        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()
            if name == 'contrastStrings':
                ctrl.set_available_contrasts(varlist)


        self.page.add(label="Model Group Variances Seperately ",
                      control=control.CHOICE_BOX,
                      name='modelGroupVariancesSeparately',
                      type=dtype.NUM,
                      comment="Specify whether FSL should model the variance for each group separately.\n\nIf this option is enabled, you must specify a grouping variable below.",
                      values=['Off', 'On'])

        self.page.add(label="Grouping Variable ",
                      control=control.TEXT_BOX,
                      name="groupingVariable",
                      type=dtype.STR,
                      comment="The name of the EV that should be used to group subjects when modeling variances.\n\nIf you do not wish to model group variances separately, set this value to None.",
                      values=self.gpa_settings['grouping_var'],
                      size=(160, -1))
        
        self.page.add(label="Model Name ",
                      control=control.TEXT_BOX,
                      name="modelName",
                      type=dtype.STR,
                      comment="Specify a name for the new model.",
                      values=self.gpa_settings['model_name'],
                      size=(200, -1))

        self.page.add(label="Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="outputModelFilesDirectory",
                      type=dtype.STR,
                      comment="Full path to the directory where CPAC should place model files.",
                      values=self.gpa_settings['output_dir'])


        if 'group_sep' in self.gpa_settings.keys():

            for ctrl in self.page.get_ctrl_list():

                name = ctrl.get_name()

                if name == 'modelGroupVariancesSeparately':

                    if self.gpa_settings['group_sep'] == True:
                        ctrl.set_value('On')
                    elif self.gpa_settings['group_sep'] == False:
                        ctrl.set_value('Off')
            

        self.page.set_sizer()

        mainSizer.Add(self.window, 1, wx.EXPAND)

        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

#         run = wx.Button(btnPanel, ID_RUN, "Create Model", (
#             280, -1), wx.DefaultSize, 0)
#         self.Bind(wx.EVT_BUTTON, lambda event: self.save(
#             event, 'run'), id=ID_RUN)
#         hbox.Add(run, 0, flag=wx.LEFT | wx.ALIGN_LEFT, border=10)

        buffer1 = wx.StaticText(btnPanel, label="\t\t\t\t\t\t")
        hbox.Add(buffer1)

        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel", (
            220, 10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add(cancel, 0, flag=wx.LEFT | wx.BOTTOM, border=5)
        
        back = wx.Button(btnPanel, 1, "< Back", (
            200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, lambda event: self.back(event), id=1)
        hbox.Add(back, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

        save = wx.Button(btnPanel, wx.ID_SAVE, "Save Settings", (
            200, -1), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, lambda event: self.save(
            event, 'save'), id=wx.ID_SAVE)
        hbox.Add(save, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)


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

    def validate(self):

        return True

        '''
        try:

            columns = [v.strip()
                       for v in config_map.get('columnsInModel')[1].split(",")]

            if not columns:
                self.display(config_map.get('columnsInModel')[
                             0], "No columns specified for the model")
                return -1


            for key, val in config_map.iteritems():

                if key != 'groupingVariable' and len(val[1]) == 0:
                    self.display(val[0], "%s field is empty!" % key)

                if '/' in val[1] and val[2]:
                    if not os.path.exists(val[1]):
                        self.display(val[
                                     0], "%s field contains incorrect path. Please enter correct path!" % key)

                if key == 'deMean' or key == 'categoricalVsDirectional':
                    
                    print val
                    
                    value = [int(v) for v in val[1].split(",")]
                    
                    for v in value:
                        if v not in [1, 0]:
                            self.display(val[
                                         0], "Invalid Entry. Only 1 and 0 entry allowed")

                    if len(value) != len(columns):
                        self.display(val[
                                     0], "Number of values in %s do not match specified columns in the model" % key)

                if key == 'groupingVariable':
                    if str(config_map.get('modelGroupVariancesSeparately')[1]) == "On":
                        if len(val[1]) == 0:
                            self.display(val[0], "%s field is empty!" % key)

                        if val[1] not in columns:
                            self.display(val[
                                         0], "Grouping variable/column not a valid column in the model. Please verify the name")

            return 1

        except Exception:
            return -1
        '''



    def parse_contrast(self, contrast_string):

        contrast_string = contrast_string.replace(' ', '')

        if '>' in contrast_string:
            split_contrast = contrast_string.split('>')
        elif '<' in contrast_string:
            split_contrast = contrast_string.split('<')
        elif '+' in contrast_string:
            split_contrast = contrast_string.split('+')
        elif '-' in contrast_string:
            split_contrast = contrast_string.split('-')
        else:
            print '[!] CPAC says: The contrast \' ', contrast_string, ' \' ' \
                  'did not contain any valid operators ( > , < , + , - ).\n'
            raise Exception

        # in the case of the '+' or '-' contrast operators, which result in
        # the split_contrast list containing a blank element ''
        for item in split_contrast:
            if item == '':
                split_contrast.remove(item)

        return split_contrast



    def collect_input(self):

        for ctrl in self.page.get_ctrl_list():

            name = ctrl.get_name()

            if name == 'contrastStrings':

                for option in ctrl.get_listbox_options():

                    # first, make sure the contrasts are valid!
                    contrasts = self.parse_contrast(option)
       
                    # check to make sure the contrast names are contrasts that
                    # are actually valid - this will only really ever happen
                    # if the user hand-edits the config file, the GUI catches
                    # invalid contrasts when entered
                    for contrast in contrasts:
                        if contrast not in self.contrasts_list:

                            errmsg = 'CPAC says: The contrast \'%s\' you ' \
                                'entered within the string \'%s\' is not ' \
                                'one of the available contrast selections.' \
                                '\n\nPlease enter only the contrast labels ' \
                                'listed under \'Available Contrasts\'.' \
                                % (contrast, option)

                            errSubID = wx.MessageDialog(self, errmsg,
                                'Invalid Contrast', wx.OK | wx.ICON_ERROR)
                            errSubID.ShowModal()
                            errSubID.Destroy()
                            raise Exception
                            

                    # then, add them to gpa_settings appropriately
                    if option in ctrl.get_listbox_selections():
                        self.gpa_settings['contrasts'][option] = True
                    else:
                        self.gpa_settings['contrasts'][option] = False


            if name == 'modelGroupVariancesSeparately':

                self.gpa_settings['group_sep'] = ctrl.get_selection()


            if name == 'groupingVariable':

                self.gpa_settings['grouping_var'] = ctrl.get_selection()


            if name == 'modelName':

                self.gpa_settings['model_name'] = ctrl.get_selection()


            if name == 'outputModelFilesDirectory':

                self.gpa_settings['output_dir'] = ctrl.get_selection()




    def back(self, event):

        self.collect_input()

        modelconfig_window.ModelConfig(self.parent, self.gpa_settings)

        self.Close()


        '''
        TO-DO:
        also pass new selections from within this window
        make it so that if modelconfig receives these, it stores them in case
        the user hits next without changing anything

        then it will restore them.
        but if the user changes the pheno or something, you need to warn them
        '''


        '''
        ALSO TO-DO:
        have the contrasts box have a method for restoring the strings
        AND their checked states
        '''



    def save(self, event, flag):

        # runs when user clicks 'Save Settings', saves a YAML (.yml/.yaml)
        # file which contains all of the user input for the group analysis
        # model setup

        self.collect_input()

        print self.gpa_settings

        config_list = []
        config_map = {}

        vals = self.gpa_settings

        config_list.append(('subject_list', vals['subject_list'], 1, \
                                'Full path to a list of subjects to be ' \
                                'included in the model.\n\nThis should be ' \
                                'a text file with one subject per line.\n' \
                                '\nTip 1: A list in this format contaning ' \
                                'all subjects run through CPAC was ' \
                                'generated along with the main CPAC ' \
                                'subject list (see subject_list_group_' \
                                'analysis.txt).\n\nTip 2: An easy way to ' \
                                'manually create this file is to copy the ' \
                                'subjects column from your Regressor/EV ' \
                                'spreadsheet.'))

        config_list.append(('pheno_file', vals['pheno_file'], 1, \
                                'Full path to a .csv file containing EV ' \
                                'information for each subject.\n\nTip: A ' \
                                'file in this format (containing a single ' \
                                'column listing all subjects run through ' \
                                'CPAC) was generated along with the main ' \
                                'CPAC subject list (see template_' \
                                'phenotypic.csv).'))

        config_list.append(('subject_id_label', vals['subject_id_label'], 1, \
                                'Name of the subjects column in your EV ' \
                                'file.'))

        config_list.append(('ev_selections', vals['ev_selections'], 8, \
                                'Specify which EVs from your phenotype ' \
                                'are categorical or numerical. Of those ' \
                                'which are numerical, specify which are to ' \
                                'be demeaned.'))

        config_list.append(('design_formula', vals['design_formula'], 1, \
                                'Formula for the design matrix. The EVs ' \
                                'included in this formula will be included ' \
                                'in the model. <MORE INFO>'))

        config_list.append(('contrasts', vals['contrasts'], 8, \
                                'A dictionary of contrast descriptions, ' \
                                'including which ones to be included in ' \
                                'the model (marked with either True or ' \
                                'False).'))

        config_list.append(('group_sep', vals['group_sep'], 0, \
                                'Specify whether FSL should model the ' \
                                'variance for each group separately.\n\n' \
                                'If this option is enabled, you must ' \
                                'specify a grouping variable below.'))

        config_list.append(('grouping_var', vals['grouping_var'], 1, \
                                'The name of the EV that should be used to ' \
                                'group subjects when modeling variances.\n' \
                                '\nIf you do not wish to model group ' \
                                'variances separately, set this value to ' \
                                'None.'))

        config_list.append(('model_name', vals['model_name'], 1, \
                                'Specify a name for the new model.'))

        config_list.append(('output_dir', vals['output_dir'], 1, \
                                'Full path to the directory where CPAC ' \
                                'should place model files.'))

            


        try:
            if self.validate() == True:
                dlg = wx.FileDialog(self, message="Save file as ...",
                                    defaultDir=os.getcwd(),
                                    defaultFile="gpa_fsl_config_%s.yaml" % self.gpa_settings['model_name'],
                                    wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
                                    style=wx.SAVE)

                if dlg.ShowModal() == wx.ID_OK:

                    path = dlg.GetPath()
                    f = open(path, 'w')
                    dlg.Destroy()

                    for item in config_list:

                        # parse the different data types accordingly

                        # file selector combo
                        if item[2] == 2:
                            value = substitution_map.get(str(item[1]))
                            if value is None:
                                value = ast.literal_eval(item[1])

                        # directory selector combo
                        elif item[2] == 5:
                            value = [v for v in ast.literal_eval(item[1])]

                        # num ctrl
                        elif item[2] == 4:
                            value = [str(v.strip())
                                     for v in item[1].split(',')]

                        # all other data types
                        else:
                            value = str(item[1])

                        # print out 'help' (comments describing values)
                        for lines in item[3].split('\n'):
                            print >> f, '#', lines

                        # print out 'label: value'
                        print >> f, item[0], ': ', value, '\n\n'

                    print '\n\nCPAC says: Saving the group analysis model ' \
                              'configuration file to: ', path, '\n\n'
                    f.close()
                    
                    self.Parent.box2.GetTextCtrl().SetValue(path)
                    self.Close()


        except Exception:
            print '\n\n[!] CPAC says: I couldn\'t save the group analysis ' \
                      'model configuration file! Maybe check if you have ' \
                      'write permissions?\n\nPath you selected: ', path, \
                      '\n\n'
            raise

            
            
