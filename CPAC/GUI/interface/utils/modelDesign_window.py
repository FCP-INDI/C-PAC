import wx
import generic_class
from .constants import control, dtype, substitution_map
import os
import ast


ID_RUN = 11


class ModelDesign(wx.Frame):

    def __init__(self, parent, subjectList, phenoFilePath, subjectID):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - Create New FSL Model", size=(900, 575))

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)

        self.window = wx.ScrolledWindow(self.panel)

        self.page = generic_class.GenericClass(self.window, " FSL Model Design")
        
        phenoFile = open(os.path.abspath(phenoFilePath))
        phenoHeaderString = phenoFile.readline().rstrip('\r\n')
        self.phenoHeaderItems = phenoHeaderString.split(',')
        
        if subjectID in self.phenoHeaderItems:
            self.phenoHeaderItems.remove(subjectID)
        else:
            errSubID = wx.MessageDialog(
                self, 'Please enter the name of the subject ID column' \
                ' as it is labeled in the phenotype file.',
                'Blank/Incorrect Subject Header Input',
                wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()
            raise Exception
        
        
        
        # experimental checkbox row stuff
        self.page.add(label = "Model Setup ",
                      control = control.CHECKBOX_GRID,
                      name = "modelSetup",
                      type = dtype.LBOOL,
                      values = self.phenoHeaderItems,
                      comment="glob",
                      size = (400, -1))
        
        # end experimental code


        self.page.add(label="Contrast File ",
                      control=control.COMBO_BOX,
                      name="contrastFile",
                      type=dtype.STR,
                      comment="Full path to a .csv file containing contrasts to be applied to this model.\n\nWhen specifying EVs in this file:\n\n- Continuous EVs should appear the same as their corresponding column name in the EV file.\n\n- Categorical EVs must be split into multiple columns (one for each category), with names of the format EVname__N (e.g. diagnosis__1, diagnosis__2, diagnosis__3)\n\nIf you wish to include F-tests in your model, create a column for each desired F-test, with names in the format f_test_1, f_test_2, etc.",
                      values="")

        self.page.add(label="Model Group Variances Seperately ",
                      control=control.CHOICE_BOX,
                      name='modelGroupVariancesSeparately',
                      type=dtype.NUM,
                      comment="Specify whether FSL should model the variance for each group separately.\n\nIf this option is enabled, you must specify a grouping variable below.",
                      values=["Off", "On"])

        self.page.add(label="Grouping Variable ",
                      control=control.TEXT_BOX,
                      name="groupingVariable",
                      type=dtype.STR,
                      comment="The name of the EV that should be used to group subjects when modeling variances.\n\nIf you do not wish to model group variances separately, set this value to None.",
                      values="None",
                      size=(160, -1))
        
        self.page.add(label="Model Name ",
                      control=control.TEXT_BOX,
                      name="modelName",
                      type=dtype.STR,
                      comment="Specify a name for the new model.",
                      values="",
                      size=(200, -1))

        self.page.add(label="Output Directory ",
                      control=control.DIR_COMBO_BOX,
                      name="outputModelFilesDirectory",
                      type=dtype.STR,
                      comment="Full path to the directory where CPAC should place model files.",
                      values="")


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

    def validate(self, config_map):
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



    def save(self, event, flag):

        config_list = []
        config_map = {}
        
        # temporary code to keep backwards-compatibility with old GPA model config
        # files, to be loaded into the new model builder GUI
        for ctrl in self.page.get_ctrl_list():
            
            # this will take the input from the checkbox grid for the "Model Setup"
            # and break it back down into the old "demean", etc. ctrls
            win = ctrl.get_ctrl()
            value = ctrl.get_selection()
            name = ctrl.get_name()
            
            if name == 'modelSetup':
                                
                # in this case, 'value' should be the choiceDict from
                # CheckBoxGrid in custom_control
                
                # turn the lists in 'value' (choiceDict) into strings
                
                columnsInModel = ', '.join(value['include'])
                categoricalVsDirectional = ', '.join(value['categorical'])
                deMean = ', '.join(value['demean'])
                
                print "columns: ", columnsInModel
                print "categorical: ", categoricalVsDirectional
                print "demean: ", deMean
                
                # add these strings into new "ctrls" named accordingly
                config_list.append(('columnsInModel', columnsInModel, 4, 'Specify the names of ', \
                                    'columns in your EV file that you would like to include in ', \
                                    'this model.\n\nColumn names should be separated by commas ', \
                                    'and appear exactly as they do in your EV file.\n\nBy ', \
                                    'clicking the add button on the right, you can also add ', \
                                    'measure generated by CPAC to the list of EVs'))
                config_list.append(('categoricalVsDirectional', categoricalVsDirectional, 5, 'Specify ', \
                                    'whether each of the EVs in this model should be treated as ', \
                                    'categorical or continuous.\n\nTo do this, place a 1 (categorical) ', \
                                    'or 0 (continuous) in the same list position as the corresponding ', \
                                    'EV.\n\nFor example, if the EVs to include were:\nage, sex, ' \
                                    'diagnosis, mean_fd\n\nOne might specify:\n0,1,1,0'))
                config_list.append(('deMean', deMean, 5, 'Specify whether to demean each of the EVs in ', \
                                    'this model.\n\nTo do this, place a 1 (demean) or 0 (don\'t demean) ', \
                                    'in the same list position as the corresponding EV.\n\nFor example, ', \
                                    'if the EVs to include were:\nage, sex, diagnosis, mean_fd\n\nOne ', \
                                    'might specify:\n1,0,0,1\n\nNote that only continuous EV\'s should ', \
                                    'be demeaned.'))
        
        

        for ctrl in self.page.get_ctrl_list():

            #print "validating ctrl-->", ctrl.get_name()
            #print "ctrl.get_selection()", ctrl.get_selection()
            #print "type(ctrl.get_selection())", type(ctrl.get_selection())

            win = ctrl.get_ctrl()
            value = str(ctrl.get_selection())
            name = ctrl.get_name()
            dtype = ctrl.get_datatype()
            validation = ctrl.get_validation()
            help = ctrl.get_help()

            # this is here because modelSetup is added as one of the controls,
            # but its output has to be processed specifically (for loop above
            # this one)
            if name != 'modelSetup':

                config_list.append((name, value, dtype, help))
                config_map[name] = [win, value, validation]



            '''
            if name == 'deMean':
                
                # take in a list of header items from the demean checklist input
                # and convert them into 1's and 0's while retaining their indices
                demeanList = []

                for headerItem in self.phenoHeaderItems:
                    demeanList.append(0)
                        
                for demeanItem in ast.literal_eval(value):
                    demeanList[self.phenoHeaderItems.index(demeanItem)] = 1
                    value = str(demeanList).strip('[]')
                        
                    # 'value' is now a list of 1's and 0's for demean

                # TEMPORARY CODE FOR HARD-CODING OF EVTYPES AS ALL ZEROES
                evTypes = value.replace("1", "0")
                config_list.append((name, evTypes, dtype, help))
            '''

        # TEMPORARY CODE FOR HARD-CODING OF EVTYPES AS ALL ZEROES
        #config_list.append(("categoricalVsDirectional", evTypes, 5, "Placeholder for EV types - should be all zeroes - do not modify."))


        try:
            if self.validate(config_map) > 0:
                dlg = wx.FileDialog(self, message="Save file as ...",
                                    defaultDir=os.getcwd(),
                                    defaultFile="gpa_fsl_config_%s.yaml" % config_map['modelName'](1),
                                    wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
                                    style=wx.SAVE)

                if dlg.ShowModal() == wx.ID_OK:

                    path = dlg.GetPath()
                    f = open(path, 'w')
                    dlg.Destroy()
                    for item in config_list:
                        if item[2] == 2:
                            value = substitution_map.get(str(item[1]))
                            if value is None:
                                value = ast.literal_eval(item[1])
                        elif item[2] == 5:
                            value = [v for v in ast.literal_eval(item[1])]
                        elif item[2] == 4:
                            value = [str(v.strip())
                                     for v in item[1].split(',')]
                        else:
                            value = str(item[1])

                        for lines in item[3].split('\n'):
                            print >> f, "#", lines

                        print >> f, item[0], ": ", value, "\n"

                    print "saving %s" % path
                    f.close()
                    
                    self.Parent.box2.GetTextCtrl().SetValue(path)
                    self.Close()


        except Exception:
            print "error writing temp file "
            raise

            
            