import wx
import generic_class
from .constants import control, dtype, substitution_map
import os
import yaml

import modelDesign_window


ID_RUN = 11


class ModelConfig(wx.Frame):

    def __init__(self, parent):

        wx.Frame.__init__(
            self, parent=parent, title="CPAC - Create New FSL Model", size=(700, 250))

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        self.panel = wx.Panel(self)

        self.window = wx.ScrolledWindow(self.panel)

        self.page = generic_class.GenericClass(self.window, " FSL Model Setup")

        self.page.add(label="Subject List ",
                      control=control.COMBO_BOX,
                      name="subjectListFile",
                      type=dtype.STR,
                      comment="Full path to a list of subjects to be included in the model.\n\nThis should be a text file with one subject per line.\n\nTip 1: A list in this format contaning all subjects run through CPAC was generated along with the main CPAC subject list (see subject_list_group_analysis.txt).\n\nTIp 2: An easy way to manually create this file is to copy the subjects column from your Regressor/EV spreadsheet.",
                      values="")

        self.page.add(label="Phenotype/EV File ",
                      control=control.COMBO_BOX,
                      name="phenotypicFile",
                      type=dtype.STR,
                      comment="Full path to a .csv file containing EV information for each subject.\n\nTip: A file in this format (containing a single column listing all subjects run through CPAC) was generated along with the main CPAC subject list (see template_phenotypic.csv).",
                      values="")

        self.page.add(label="Subjects Column Name ",
                      control=control.TEXT_BOX,
                      name="subjectColumn",
                      type=dtype.STR,
                      comment="Name of the subjects column in your EV file.",
                      values="",
                      style=wx.EXPAND | wx.ALL,
                      size=(160, -1))

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
        
        populate = wx.Button(btnPanel, wx.ID_ANY, "Next >",
            (220, 10), wx.DefaultSize, 0)
        self.Bind(wx.EVT_BUTTON, self.populateEVs, id=wx.ID_ANY)
        hbox.Add(populate, 0.6, flag=wx.LEFT | wx.BOTTOM, border=5)

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

          
            
    def populateEVs(self, event):
        
        # somehow get in the header from the phenotype
        # also get in the subject column
        
        for ctrl in self.page.get_ctrl_list():
            
            name = ctrl.get_name()
            
            if name == 'subjectListFile':
                subjectList = str(ctrl.get_selection())
            
            if name == 'phenotypicFile':
                phenoFilePath = str(ctrl.get_selection())
                
            if name == 'subjectColumn':
                subjectID = str(ctrl.get_selection())
        
                
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
                
        
        testFile(subjectList, 'Subject List')
        testFile(phenoFilePath, 'Phenotype/EV File')
             
        
        modelDesign_window.ModelDesign(self,subjectList,phenoFilePath,subjectID)
        
        

