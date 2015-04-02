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
            self, parent=parent, title="CPAC - Create New FSL Model", size=(850, 650))

        self.parent = parent

        self.gpa_settings = gpa_settings

        self.contrasts_list = varlist


        if 'contrasts' not in self.gpa_settings.keys():
            self.gpa_settings['contrasts'] = {}

        if 'custom_contrasts' not in self.gpa_settings.keys():
            self.gpa_settings['custom_contrasts'] = 'None'

        if 'f_tests' not in self.gpa_settings.keys():
            self.gpa_settings['f_tests'] = []

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
        


        # build 'Available contrasts' string
        contrasts_text = 'Available EVs for contrasts:\n'

        con_length = 65

        for con in varlist:
            contrasts_text = contrasts_text + '    ' + con

            if len(contrasts_text) > con_length:
                contrasts_text = contrasts_text + '\n'
                con_length += 50
            

        varlist_sizer = wx.BoxSizer(wx.HORIZONTAL)
        var_list_text = wx.StaticText(self.window, label=str(contrasts_text))
        varlist_sizer.Add(var_list_text)

        self.page.add_pheno_load_panel(varlist_sizer)
    

        self.page.add(label = 'Contrasts ',
                      control = control.LISTBOX_COMBO,
                      name = 'contrasts',
                      type = dtype.LSTR,
                      values = self.gpa_settings['contrasts'],
                      comment = 'Specify your contrasts in this window. For example, if two of your available contrasts are EV1 and EV0, you can enter contrast descriptions such as EV1 > EV0 or EV1+ . Consult the User Guide for more information about describing contrasts. Alternatively, you can provide your own custom-written contrasts matrix in a CSV file in the \'Custom Contrasts Matrix\' field below.',
                      size = (300,200),
                      combo_type = 4)

        # this sends the list of available contrast names to the 'Add
        # Contrast' dialog box, so that it may do validation immediately when
        # the user enters contrast strings
        for ctrl in self.page.get_ctrl_list():
            name = ctrl.get_name()
            if name == 'contrasts':
                ctrl.set_available_contrasts(varlist)


        self.page.add(label = 'f-Tests ',
                      control = control.LISTBOX_COMBO,
                      name = 'f_tests',
                      type = dtype.LSTR,
                      values = self.gpa_settings['f_tests'],
                      comment = 'Optional: Specify your desired f-tests, if any, in this window.',
                      size = (300,120),
                      combo_type = 5)

        self.page.add(label="Custom Contrasts Matrix ",
                      control=control.COMBO_BOX,
                      name="custom_contrasts",
                      type=dtype.STR,
                      comment="Optional: Full path to a CSV file which specifies the contrasts you wish to run in group analysis. This allows you to describe your own custom contrasts matrix if you do not wish to use the contrasts builder above. Consult the User Guide for proper formatting.\n\nIf you wish to use the standard contrast builder, leave this field blank. If you provide a path for this option, CPAC will use your custom contrasts matrix instead, and will use the f-tests described in this custom file only (ignoring those you have input in the f-tests field in the GUI above).\n\nIf you wish to include f-tests, create a new column in your CSV file for each f-test named 'f_test_1', 'f_test_2', .. etc. Then, mark the contrasts you would like to include in each f-test with a 1, and mark the rest 0. Note that you must select at least two contrasts per f-test.",
                      values=str(self.gpa_settings['custom_contrasts']))
        
        self.page.add(label="Model Name ",
                      control=control.TEXT_BOX,
                      name="model_name",
                      type=dtype.STR,
                      comment="Specify a name for the new model. Output and working directories for group analysis, as well as the FLAMEO model files (.mat, .con, .grp, etc.) will be labeled with this name.",
                      values=self.gpa_settings['model_name'],
                      size=(200, -1))

        self.page.add(label="Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="output_dir",
                      type=dtype.STR,
                      comment="Full path to the directory where CPAC should place the model files (.mat, .con, .grp) and the outputs of group analysis.",
                      values=self.gpa_settings['output_dir'])

           

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

            if name == 'contrasts':

                self.gpa_settings['contrasts'] = []

                for option in ctrl.get_selection(): #listbox_options():

                    # first, make sure the contrasts are valid!
                    contrasts = self.parse_contrast(option)
       
                    # check to make sure the contrast names are contrasts that
                    # are actually valid - this will only really ever happen
                    # if the user hand-edits the config file; otherwise the GUI
                    # catches invalid contrasts when entered
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
                            
                    self.gpa_settings['contrasts'].append(option)



            if name == 'f_tests':

                self.gpa_settings['f_tests'] = []

                for option in ctrl.get_selection():

                    cons_in_ftest = []

                    for con in option.split(","):
                        cons_in_ftest.append(con)

                    # check to make sure the contrasts in the f-test lists
                    # are actually valid - this will only really ever happen
                    # if the user hand-edits the config file, the GUI catches
                    # invalid contrasts when entered
                    for contrast in cons_in_ftest:

                        if contrast not in self.window.input_contrasts:

                            errmsg = 'CPAC says: The contrast \'%s\' you ' \
                                'entered within the f-test \'%s\' is not ' \
                                'one of the available contrast selections.' \
                                '\n\nPlease enter only the contrast labels ' \
                                'of the contrasts you have specified.' \
                                % (contrast, option)

                            errSubID = wx.MessageDialog(self, errmsg,
                                'Invalid Contrast', wx.OK | wx.ICON_ERROR)
                            errSubID.ShowModal()
                            errSubID.Destroy()
                            raise Exception
                            

                    # then, add them to gpa_settings appropriately
                    self.gpa_settings['f_tests'].append(option)


            if name == 'custom_contrasts':

                self.gpa_settings['custom_contrasts'] = ctrl.get_selection()

                custom_confile = self.gpa_settings['custom_contrasts']

                if not ((custom_confile == None) or (custom_confile == '') or \
                    ("None" in custom_confile)):

                    if not os.path.exists(custom_confile):
                    
                        errmsg = "You've specified a Custom Contrasts file " \
                                 "for your group model, but this file " \
                                 "cannot be found. Please double-check the " \
                                 "filepath you have entered."

                        errSubID = wx.MessageDialog(self, errmsg,
                                'Invalid Path', wx.OK | wx.ICON_ERROR)
                        errSubID.ShowModal()
                        errSubID.Destroy()
                        raise Exception


            if name == 'model_name':

                self.gpa_settings['model_name'] = ctrl.get_selection()


            if name == 'output_dir':

                self.gpa_settings['output_dir'] = ctrl.get_selection()



    '''
    def get_custom_contrasts(self.window):

        self.collect_input()

        confilepath = self.gpa_settings['custom_contrasts']

        con_names = []

        if (confilepath != None) or (confilepath != '') or \
            ("None" not in confilepath):

            if os.path.exists(confilepath):

                confile = open(confilepath, 'rb')

                for con in confile.readlines():

                    con_names.append(con.split(",")[0])

                # get rid of "Contrasts" header label that came from the first
                # row in the file (if formatted properly)
                del con_names[0]

        return con_names
    '''




    def back(self, event):

        self.collect_input()

        modelconfig_window.ModelConfig(self.parent, self.gpa_settings)

        self.Close()


        '''
        done:
        also pass new selections from within this window (done)
        make it so that if modelconfig receives these, it stores them in case
        the user hits next without changing anything

        then it will restore them.

        to do:
        but if the user changes the pheno or something, you need to warn them
        '''



    def collect_contrasts(self):
        pass




    def save(self, event, flag):

        # runs when user clicks 'Save Settings', saves a YAML (.yml/.yaml)
        # file which contains all of the user input for the group analysis
        # model setup

        self.collect_input()

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
                                'Specify the formula to describe your ' \
                                'model design. Essentially, including EVs ' \
                                'in this formula inserts them into the ' \
                                'model. The most basic format to include ' \
                                'each EV you select would be \'EV + EV + EV ' \
                                '+ ..\', etc. You can also select to ' \
                                'include MeanFD, Measure_Mean, and ' \
                                'Custom_ROI_Mean here. See the C-PAC User ' \
                                'Guide for more detailed information ' \
                                'regarding formatting your design formula.'))

        config_list.append(('mean_mask', vals['mean_mask'], 4, \
                                'Choose whether to use a group mask or ' \
                                'individual-specific mask when calculating ' \
                                'the output means to be used as a ' \
                                'regressor.\n\nThis only takes effect if ' \
                                'you include the \'Measure_Mean\' regressor ' \
                                'in your Design Matrix Formula.'))

        config_list.append(('custom_roi_mask', vals['custom_roi_mask'], 1, \
                                'Full path to a NIFTI file containing one ' \
                                'or more ROI masks. The means of the masked ' \
                                'regions will then be computed for each ' \
                                'subject\'s output and will be included in ' \
                                'the model as regressors (one for each ROI ' \
                                'in the mask file) if you include ' \
                                '\'Custom_ROI_Mean\' in the Design Matrix ' \
                                'Formula.'))

        config_list.append(('use_zscore', vals['use_zscore'], 0, \
                                'Run the group analysis model on the ' \
                                'z-score standardized version of the ' \
                                'derivatives you choose in the list below.'))

        config_list.append(('derivative_list', vals['derivative_list'], 6, \
                                "Choose the derivatives to run the group " \
                                "model on.\n\nThese must be written out " \
                                "as a list, and must be one of the options " \
                                "listed below.\n\nFor z-scored analyses:\n" \
                                "'alff_to_standard_zstd', " \
                                "'alff_to_standard_smooth_zstd', " \
                                "'falff_to_standard_zstd', " \
                                "'falff_to_standard_smooth_zstd', " \
                                "'reho_to_standard_zstd', " \
                                "'reho_to_standard_smooth_zstd', " \
                                "'sca_roi_to_standard_fisher_zstd', " \
                                "'sca_roi_to_standard_smooth_fisher_zstd', " \
                                "'sca_seed_to_standard_fisher_zstd', " \
                                "'sca_seed_to_standard_smooth_fisher_zstd', " \
                                "'vmhc_fisher_zstd', " \
                                "'vmhc_fisher_zstd_zstat_map', " \
                                "'dr_tempreg_maps_zstat_files_to_standard', " \
                                "'dr_tempreg_maps_zstat_files_to_standard_smooth', " \
                                "'sca_tempreg_maps_zstat_files', " \
                                "'sca_tempreg_maps_zstat_files_smooth', " \
                                "'centrality_outputs_zstd', " \
                                "'centrality_outputs_smoothed_zstd'\n\n" \
                                "For raw (non-z-scored) analyses:\n" \
                                "'alff_to_standard', " \
                                "'alff_to_standard_smooth', " \
                                "'falff_to_standard', " \
                                "'falff_to_standard_smooth', " \
                                "'reho_to_standard', " \
                                "'reho_to_standard_smooth', " \
                                "'sca_roi_to_standard', " \
                                "'sca_roi_to_standard_smooth', " \
                                "'sca_seed_to_standard', " \
                                "'sca_seed_to_standard_smooth', " \
                                "'centrality_outputs', " \
                                "'centrality_outputs_smoothed', " \
                                "'dr_tempreg_maps_files_to_standard', " \
                                "'dr_tempreg_maps_files_to_standard_smooth', " \
                                "'sca_tempreg_maps_files', " \
                                "'sca_tempreg_maps_files_smooth'\n\n" \
                                "Example input: derivative_list :  ['alff_to" \
                                "_standard_smooth_zstd', 'sca_roi_to_" \
                                "standard_smooth_fisher_zstd']\n"))


        config_list.append(('coding_scheme', vals['coding_scheme'], 4, \
                                'Choose the coding scheme to use when ' \
                                'generating your model. \'Treatment\' ' \
                                'encoding is generally considered the ' \
                                'typical scheme. Consult the User Guide for ' \
                                'more information.\n\nAvailable options:\n' \
                                '\'Treatment\', \'Sum\'\n'))
                                
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

        config_list.append(('z_threshold', vals['z_threshold'], 4, \
                                'Only voxels with a Z-score higher than ' \
                                'this value will be considered significant.'))

        config_list.append(('p_threshold', vals['p_threshold'], 4, \
                                'Significance threshold (P-value) to use ' \
                                'when doing cluster correction for multiple ' \
                                'comparisons.'))

        config_list.append(('repeated_measures', vals['repeated_measures'], 0, \
                                'Run repeated measures to compare different ' \
                                'scans (must use the group analysis subject ' \
                                'list and phenotypic file formatted for ' \
                                'repeated measures.'))

        config_list.append(('contrasts', vals['contrasts'], 8, \
                                'A list of contrast descriptions.'))

        config_list.append(('f_tests', vals['f_tests'], 8, \
                                'Optional: A list of f-test strings ' \
                                'containing contrasts. If you do not wish ' \
                                'to run f-tests, leave this blank.'))

        config_list.append(('custom_contrasts', vals['custom_contrasts'], 1, \
                                'Optional: Full path to a CSV file which ' \
                                'specifies the contrasts you wish to run in ' \
                                'group analysis. Consult the User Guide for ' \
                                'proper formatting.\nIf you wish to use the ' \
                                'standard contrast builder, leave this ' \
                                'field blank. If you provide a path for ' \
                                'this option, CPAC will use your custom ' \
                                'contrasts matrix instead, and will use the ' \
                                'f-tests described in this custom file only ' \
                                '(ignoring those you have input in the ' \
                                'f-tests field above).\nIf you wish to ' \
                                'include f-tests, create a new column in ' \
                                'your CSV file for each f-test named ' \
                                '\'f_test_1\', \'f_test_2\', .. etc. Then, ' \
                                'mark the contrasts you would like to ' \
                                'include in each f-test with a 1, and mark ' \
                                'the rest 0. Note that you must select at ' \
                                'least two contrasts per f-test.'))

        config_list.append(('model_name', vals['model_name'], 1, \
                                'Specify a name for the new model.'))

        config_list.append(('output_dir', vals['output_dir'], 1, \
                                'Full path to the directory where CPAC ' \
                                'should place model files.'))
            


        try:

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

                        # patchwork code to get this working for now, but
                        # why does coding_scheme print out as a string
                        # everytime, even with the list formatting?
                        # example: ['Treatment'] is fully a string in the
                        # yaml file, including the [' and '], which are
                        # characters in the string
                        if isinstance(item[1], str) and "[" in item[1] and "]" in item[1]:
                            value = item[1].replace("['","")
                            value = value.replace("']","")

                        else:
                            # regular handling
                            value = [str(v.strip())
                                     for v in item[1].split(',')]

                    # all other data types
                    else:
                        value = str(item[1])


                    if item[0] == 'derivative_list':

                        value = []

                        # this takes the user selection in the derivative
                        # list and matches it with the output directory
                        # folder name for each chosen derivative via the
                        # substitution map in constants.py

                        # go over each string in the list
                        for val in ast.literal_eval(str(item[1])):
                            if substitution_map.get(val) != None:
                                value.append(substitution_map.get(val))
                            elif val != 'None':
                                value.append(ast.literal_eval(val))


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

        
        except Exception as e:

            errmsg = '\n\n[!] CPAC says: Couldn\'t save the group analysis ' \
                      'model configuration file! Maybe check if you have ' \
                      'write permissions?\n\nPath you selected: %s\n\n' \
                      'Error details: %s' % (path, e)

            raise Exception(errmsg)
        

            
            
