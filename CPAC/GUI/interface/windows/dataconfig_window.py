import wx
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
import os
import yaml
import pkg_resources as p
import sys

ID_RUN_EXT = 11
ID_RUN_MEXT = 12

class DataConfig(wx.Frame):
        
    
    def __init__(self, parent):

        wx.Frame.__init__(self, parent, title="CPAC - Subject List Setup", size = (820,450))
        
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        
        self.panel = wx.Panel(self)
        
        self.window = wx.ScrolledWindow(self.panel)
        
        self.page = GenericClass(self.window, "Subject List Setup")
        
        self.page.add(label= "Anatomical File Path Template ",
                 control = control.TEXT_BOX,
                 name = "anatomicalTemplate",
                 type = dtype.STR,
                 comment = "File Path Template for Anatomical Files\n\n"
                           "Replace the site- and subject-level directories with %s.\n\n"
                           "See User Guide for more detailed instructions.",
                 values ="",
                 style= wx.EXPAND | wx.ALL,
                 size = (532,-1))
        
        self.page.add(label= "Functional File Path Template ",
                 control = control.TEXT_BOX,
                 name = "functionalTemplate",
                 type = dtype.STR,
                 comment = "File Path Template for Functional Files\n\n"
                           "Replace the site- and subject-level directories with %s.\n\n"
                           "See User Guide for more detailed instructions.",
                 values ="",
                 style= wx.EXPAND | wx.ALL,
                 size = (532,-1))

        self.page.add(label="Subjects to Include (Optional) ", 
                 control=control.COMBO_BOX, 
                 name = "subjectList", 
                 type = dtype.COMBO, 
                 comment = "Include only a sub-set of the subjects present in the folders defined above.\n\n"
                           "List subjects in this box (e.g., sub101, sub102) or provide the path to a\n"
                           "text file with one subject on each line.\n\n"
                           "If 'None' is specified, CPAC will include all subjects.", 
                 values = "None")
        
        self.page.add(label="Subjects to Exclude (Optional) ", 
                 control=control.COMBO_BOX, 
                 name = "exclusionSubjectList", 
                 type = dtype.COMBO, 
                 comment = "Exclude a sub-set of the subjects present in the folders defined above.\n\n"
                           "List subjects in this box (e.g., sub101, sub102) or provide the path to a\n"
                           "text file with one subject on each line.\n\n"
                           "If 'None' is specified, CPAC will not exclude any subjects.", 
                 values = "None")
        
        self.page.add(label= "Sites to Include (Optional) ",
                 control = control.TEXT_BOX,
                 name = "siteList",
                 type = dtype.STR,
                 comment = "Include only a sub-set of the sites present in the folders defined above.\n\n"
                           "List sites in this box (e.g., NYU, UCLA) or provide the path to a text\n"
                           "file with one site on each line.\n\n"
                           "If 'None' is specified, CPAC will include all sites.",
                 values ="None",
                 style= wx.EXPAND | wx.ALL,
                 size = (532,-1))
        
        self.page.add(label="Scan Parameters File (Optional) ", 
                 control=control.COMBO_BOX, 
                 name = "scanParametersCSV", 
                 type = dtype.COMBO, 
                 comment = "Required for Slice Timing Correction.\n\n"
                           "Path to a .csv file containing information about scan acquisition parameters.\n\n"
                           "For instructions on how to create this file, see the User Guide.\n\n"
                           "If 'None' is specified, CPAC will skip Slice Timing Correction.",
                 values = "None")
        
        self.page.add(label = "Output Directory ", 
                      control = control.DIR_COMBO_BOX, 
                      name = "outputSubjectListLocation", 
                      type = dtype.STR, 
                      comment = "Directory where CPAC should place subject list files.",
                      values = "")

        self.page.add(label = "Subject List Name ", 
                      control = control.TEXT_BOX, 
                      name = "subjectListName", 
                      type = dtype.STR, 
                      comment = "A label to be appended to the generated " \
                                "subject list files.",
                      values = "",
                      style= wx.EXPAND | wx.ALL,
                      size = (300,-1))


        self.page.set_sizer()
         
        mainSizer.Add(self.window, 1, wx.EXPAND)
        
        btnPanel = wx.Panel(self.panel, -1)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        
        self.multiscan = wx.CheckBox(btnPanel, -1, label = "Multiscan Data")
        
        if 'linux' in sys.platform:
            hbox.Add(self.multiscan,0, flag=wx.TOP, border=5)
        else:
            hbox.Add(self.multiscan, 0, flag=wx.RIGHT | wx.BOTTOM, border=5)
        
        img = wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/help.png'), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        help = wx.BitmapButton(btnPanel, id=-1, bitmap=img,
                               pos=(10, 20), size = (img.GetWidth()+5, img.GetHeight()+5))
        help.Bind(wx.EVT_BUTTON, self.onHelp)
        
        if 'linux' in sys.platform:
            hbox.Add(help, 0, flag = wx.TOP, border =5)
        else:
            hbox.Add(help, 0, flag=wx.RIGHT | wx.BOTTOM, border=5)
        
        buffer2 = wx.StaticText(btnPanel, label = "\t")
        hbox.Add(buffer2)

        run_ext = wx.Button(btnPanel, ID_RUN_EXT, "Generate Subject Lists", (280,10), wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, lambda event: self.save(event,'run'), id=ID_RUN_EXT)
        hbox.Add( run_ext, 1, flag=wx.LEFT|wx.ALIGN_LEFT, border=10)
        
        buffer = wx.StaticText(btnPanel, label = "\t\t\t\t")
        hbox.Add(buffer)
    
        cancel = wx.Button(btnPanel, wx.ID_CANCEL, "Cancel",(220,10), wx.DefaultSize, 0 )

        self.Bind(wx.EVT_BUTTON, self.cancel, id=wx.ID_CANCEL)
        hbox.Add( cancel, 0, flag=wx.LEFT|wx.BOTTOM, border=5)
        
        load = wx.Button(btnPanel, wx.ID_ADD, "Load Settings", (280,10), wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, self.load, id=wx.ID_ADD)
        hbox.Add(load, 0.6, flag=wx.LEFT|wx.BOTTOM, border=5)
        
        save = wx.Button(btnPanel, wx.ID_SAVE, "Save Settings", (280,10), wx.DefaultSize, 0 )
        self.Bind(wx.EVT_BUTTON, lambda event: self.save(event,'save'), id=wx.ID_SAVE)
        hbox.Add(save, 0.6, flag=wx.LEFT|wx.BOTTOM, border=5)
    
        btnPanel.SetSizer(hbox)
        
        mainSizer.Add(btnPanel, 0.5,  flag=wx.ALIGN_RIGHT|wx.RIGHT, border=20)
        
        self.panel.SetSizer(mainSizer)
        
        self.Show()
        
    def cancel(self, event):
        self.Close()
    
    def onHelp(self, event):
            comment = "Check the box only if the scans have different slice timing infomation."
            wx.TipWindow(self, comment, 500)
        
    def run(self, config):
            
        try:  

            try:
                config_map = yaml.load(open(config, 'r'))
                out_location = os.path.join(\
                               os.path.realpath(config_map.get('outputSubjectListLocation')),\
                               'CPAC_subject_list_%s.yml' % config_map.get('subjectListName')[0])
            except Exception, e:
                print "Error loading data config file", e
                raise 
            
            
            print "executing extract data"
            multiscan = self.multiscan.IsChecked()
            
            import CPAC
            
            if multiscan:
                CPAC.utils.extract_data_multiscan.run(config)
            else:
                CPAC.utils.extract_data.run(config)
            
            while True:
                
                dlg2 = wx.TextEntryDialog(self, 'Please enter a name for the Subject List',
                                                 'Sublist Name', 'CPAC_subject_list_%s.yml' % config_map.get('subjectListName')[0])
                if dlg2.ShowModal() == wx.ID_OK:
                    if len(dlg2.GetValue()) >0:
                        parent = self.Parent
                        map = parent.get_sublist_map()
                        if map.get(dlg2.GetValue()) == None:
                            map[dlg2.GetValue()]= out_location
                            parent.listbox2.Append(dlg2.GetValue())
                            dlg2.Destroy()
                            break
                        else:
                            dlg3 = wx.MessageDialog(self, 'Subject List with this name already exist','Error!',
                                                    wx.OK | wx.ICON_ERROR)
                            dlg3.ShowModal()
                            dlg3.Destroy()
            
            return 1


        
        except ImportError, e:
            wx.MessageBox("Error importing CPAC. Unable to run extract data tool.", "Error") 
            print "Error importing CPAC"
            print e
            return -1
        
        except Exception, e:
            dlg2 = wx.MessageDialog(self, "Error Creating CPAC Subject List.\n%s"%e,
                               'Error!',
                           wx.OK | wx.ICON_ERROR)
            dlg2.ShowModal()
            dlg2.Destroy()
            return -1
         

    def save(self, event, flag):
        
        config_list =[]
        def display(win, msg):
            wx.MessageBox(msg, "Error")
            win.SetBackgroundColour("pink")
            win.SetFocus()
            win.Refresh()
            raise ValueError
        
        try:
            for ctrl in self.page.get_ctrl_list():
                #print "validating ctrl-->", ctrl.get_name()
                win = ctrl.get_ctrl()
                #print "ctrl.get_selection()", ctrl.get_selection()
                #print "type(ctrl.get_selection())", type(ctrl.get_selection())
                        
                value = str(ctrl.get_selection())
                value = value.strip()
                name = ctrl.get_name()
                dtype= ctrl.get_datatype()

                if name == 'subjectListName':
                    subject_list_name = value
                      
                if len(value) == 0:
                    display(win,"%s field must contain some text!"%ctrl.get_name())
                            
                if 'Template' in name:
                    if value.count('%s') != 2:
                        display(win,"Incorrect template, two \'%s\' values are required. One for site and another for"\
                                " subject location in the path. Please refer to example!")
                        
                    if value.startswith('%s'):
                        display(win, "Template cannot start with %s")
                        
                if '/' in value and 'Template' not in name:
                    if not os.path.exists(value):
                        display(win,"%s field contains incorrect path. Please update the path!"%ctrl.get_name())
         
                config_list.append((name, value, dtype))
                
        except Exception, e:

            errdlg = wx.MessageDialog(self, "Could not save your subject " \
                               "list information.\n\n%s" % e,
                               'Error!',
                           wx.OK | wx.ICON_ERROR)
            errdlg.ShowModal()
            errdlg.Destroy()

            print e
            return
            
        else:
        
            dlg = wx.FileDialog(
                self, message="Save file as ...", 
                defaultDir=os.getcwd(), 
                defaultFile="data_config_%s.yaml" % subject_list_name, 
                wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml", 
                style=wx.SAVE)
            
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                dlg.Destroy()
                f = open(path, 'w')
                for ctrl in config_list:
        
                    if "/" in ctrl[1] or "%s" in ctrl[1] or 'None' in ctrl[1]: 
                        value = ctrl[1]
                    else:
                        value =[val.strip() for val in ctrl[1].split(',')]
                    
                    print name, ":", value, "\n"
                    print >>f, ctrl[0], " : ", value, "\n"
                
                f.close()
                print "saving %s"%path
                
                if flag == 'run':
                    if self.run(path) >0:
                        self.Close()
                    
            

    def load(self, event):
            dlg = wx.FileDialog(
            self, message="Choose the config fsl yaml file",
                defaultDir=os.getcwd(), 
                defaultFile="",
                wildcard= "YAML files(*.yaml, *.yml)|*.yaml;*.yml",
                style=wx.OPEN | wx.CHANGE_DIR)
            # Once user click's OK
            if dlg.ShowModal() == wx.ID_OK:
                # Try and load in the data config file to GUI
                try:
                    path = dlg.GetPath()
                    # Check for path existence
                    if os.path.exists(path):
                        path = os.path.realpath(path)
                        # Try and load in file contents
                        try:
                            config_map = yaml.load(open(path, 'r'))
                        except Exception as e:
                            print 'Unable to load in the specified file: %s' \
                                  % path
                            print 'Error:\n%s' % e
                        # If it's a dictionary, check it has anat template key
                        if type(config_map) == dict:
                            if not config_map.has_key('anatomicalTemplate'):
                                err_msg = 'File is not a data configuration '\
                                          'file. It might be a pipeline '\
                                          'configuration file.'
                                raise Exception(err_msg)
                        # It didn't load in as a dictionary, report error
                        else:
                            err_msg = 'File is not a data configuration '\
                                      'file. It might be a subject list file.'
                            raise Exception(err_msg)
                    # Otherwise, report error
                    else:
                        err_msg = 'File %s does not exist. Check and try '\
                                  'again' % path
                        raise Exception(err_msg)
                    # Populate GUI fields
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
                            val = value
                
                        ctrl.set_value(str(val))
                # There was an error loading parameters, report it
                except Exception as e:
                    errdlg = wx.MessageDialog(self, "CPAC could not load " \
                               "your subject list information. Double-" \
                               "check the formatting of your data_config " \
                               "YAML file.\n\nIssue info:\n%s" % e,
                               'Error!',
                           wx.OK | wx.ICON_ERROR)
                    errdlg.ShowModal()
                    errdlg.Destroy()

                # Close dialog
                dlg.Destroy()