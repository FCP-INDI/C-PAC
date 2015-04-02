# Import packages
import wx
from config_window import MainFrame
from dataconfig_window import DataConfig
from ..utils.custom_control import FileSelectorCombo
from ..utils.constants import multiple_value_wfs
import wx.lib.agw.aquabutton as AB
import os
import pkg_resources as p
import sys
from CPAC.utils import Configuration
import yaml

# Init constants
ID_NEW = 1
ID_RENAME = 2
ID_CLEAR = 3
ID_DELETE = 4
ID_CREATE = 5
ID_LOAD = 6
ID_EDIT = 7
ID_ADD = 8
ID_SHOW = 9
ID_DISPLAY = 10
ID_CLEARALL = 11


# ListBox class definition
class ListBox(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(700, 650),  style= wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX)

        # Import packages
        import CPAC
        
        self.CreateStatusBar()
        self.SetStatusText("The Configurable Pipeline for the Analysis of Connectomes (C-PAC) v" + CPAC.__version__)
    
        self.pipeline_map = {}
        self.sublist_map= {}
        
        mainPanel = wx.Panel(self)
        mainPanel.SetBackgroundColour('#E9E3DB')
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        
        outerPanel1 = wx.Panel(mainPanel)
        outerSizer1 = wx.BoxSizer(wx.HORIZONTAL)
        
        outerPanel2 = wx.Panel(mainPanel)
        outerSizer2 = wx.BoxSizer(wx.HORIZONTAL)
              
        innerPanel1 = wx.Panel(outerPanel1)
        innerSizer1 = wx.BoxSizer(wx.HORIZONTAL)
         
        innerPanel2 = wx.Panel(outerPanel1, )
        innerSizer2 = wx.BoxSizer(wx.HORIZONTAL)
                
        lboxPanel1 = wx.Panel(innerPanel1)
        lboxSizer1 = wx.BoxSizer(wx.VERTICAL)
        btnPanel1 = wx.Panel(innerPanel1, -1)
        btnSizer1 = wx.BoxSizer(wx.VERTICAL)
        
        label = wx.StaticText(lboxPanel1, -1, "Pipelines")
        
        if 'linux' in sys.platform:
            label.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
        else:
            label.SetFont(wx.Font(16, wx.SWISS, wx.NORMAL, wx.BOLD))
            
        self.listbox = wx.CheckListBox(lboxPanel1, -1, size = (160,400))
        
        lboxSizer1.Add(label, 0, wx.ALIGN_CENTER)
        lboxSizer1.Add(self.listbox, 1, wx.EXPAND | wx.ALL, 10)
        lboxPanel1.SetSizer(lboxSizer1)
        
        lboxPanel1.SetBackgroundColour('#E9E3DB')
        
        new = wx.Button(btnPanel1, ID_NEW, 'New', size=(90, 30))
        ren = wx.Button(btnPanel1, ID_RENAME, 'Rename', size=(90, 30))
        dlt = wx.Button(btnPanel1, ID_DELETE, 'Delete', size=(90, 30))
        load = wx.Button(btnPanel1, ID_LOAD, 'Load', size=(90,30))
        edit = wx.Button(btnPanel1, ID_EDIT, 'Edit', size=(90,30))
        shw = wx.Button(btnPanel1, ID_DISPLAY, 'View', size=(90,30))
        clr = wx.Button(btnPanel1, ID_CLEAR, 'Clear', size=(90, 30))
    
        self.Bind(wx.EVT_BUTTON, self.NewItem, id=ID_NEW)
        self.Bind(wx.EVT_BUTTON, self.OnRename, id=ID_RENAME)
        self.Bind(wx.EVT_BUTTON, self.OnDelete, id=ID_DELETE)
        self.Bind(wx.EVT_BUTTON, self.AddConfig, id=ID_LOAD)
        self.Bind(wx.EVT_BUTTON, self.OnEdit, id=ID_EDIT)
        self.Bind(wx.EVT_BUTTON, self.OnDisplay, id= ID_DISPLAY)
        self.Bind(wx.EVT_BUTTON, lambda event: self.OnClear(event, 1), id=ID_CLEAR)
        self.Bind(wx.EVT_LISTBOX_DCLICK, self.OnDisplay)        

        if 'linux' in sys.platform:
            btnSizer1.Add((-1,30))
        else:
            btnSizer1.Add((-1, 27))
        
        btnSizer1.Add(new, 0, wx.TOP)
        btnSizer1.Add(load, 0, wx.TOP)
        btnSizer1.Add(edit, 0, wx.TOP)
        btnSizer1.Add(shw, 0, wx.TOP)
        btnSizer1.Add(ren, 0, wx.TOP)
        btnSizer1.Add(dlt, 0, wx.TOP)
        btnSizer1.Add(clr, 0, wx.TOP)
        btnPanel1.SetSizer(btnSizer1)
        
        btnPanel1.SetBackgroundColour('#E9E3DB')
                
        innerSizer1.Add(lboxPanel1, 1, wx.EXPAND | wx.ALL)
        
        if 'linux' in sys.platform:
            innerSizer1.Add(btnPanel1, 1, wx.EXPAND | wx.ALL, 5)
        else:
            innerSizer1.Add(btnPanel1, 1, wx.EXPAND | wx.ALL)
        
        innerPanel1.SetSizer(innerSizer1)
        innerPanel1.SetBackgroundColour('#E9E3DB')
        
        lboxPanel2 = wx.Panel(innerPanel2)
        lboxSizer2 = wx.BoxSizer(wx.VERTICAL)
        btnPanel2 = wx.Panel(innerPanel2, -1)
        btnSizer2 = wx.BoxSizer(wx.VERTICAL)

        label2 = wx.StaticText(lboxPanel2, -1, "Subject Lists")
        
        if 'linux' in sys.platform:
            label2.SetFont(wx.Font(12, wx.SWISS, wx.NORMAL, wx.BOLD))
        else:
            label2.SetFont(wx.Font(16, wx.SWISS, wx.NORMAL, wx.BOLD))
            
        self.listbox2 = wx.CheckListBox(lboxPanel2, -1, size= (160,400)) 
        self.listbox2.Bind(wx.EVT_LISTBOX_DCLICK, self.OnShow)
        lboxSizer2.Add(label2, 0,wx.ALIGN_CENTER)
        lboxSizer2.Add(self.listbox2, 1, wx.EXPAND | wx.ALL, 10)
        lboxPanel2.SetSizer(lboxSizer2)
        
        lboxPanel2.SetBackgroundColour('#E9E3DB')
                
        create = wx.Button(btnPanel2, ID_CREATE, 'New', size=(90, 30))
        add = wx.Button(btnPanel2, ID_ADD, 'Load', size= (90,30))
        show = wx.Button(btnPanel2, ID_SHOW, 'View', size= (90,30))
        clr2 = wx.Button(btnPanel2, ID_CLEARALL, 'Clear', size=(90, 30))
        self.Bind(wx.EVT_BUTTON, self.CreateItem, id=ID_CREATE)
        self.Bind(wx.EVT_BUTTON, self.AddItem, id=ID_ADD)
        self.Bind(wx.EVT_BUTTON, self.OnShow, id= ID_SHOW)
        self.Bind(wx.EVT_BUTTON, lambda event: self.OnClear(event, 2), id=ID_CLEARALL)
        
        if 'linux' in sys.platform:
            btnSizer2.Add((-1,30))
        else:
            btnSizer2.Add((-1, 27))

        # Add buttons to button sizer
        btnSizer2.Add(create, 0, wx.TOP)
        btnSizer2.Add(add, 0, wx.TOP)
        btnSizer2.Add(show, 0, wx.TOP)
        btnSizer2.Add(clr2, 0, wx.TOP)
        btnPanel2.SetSizer(btnSizer2)
        btnPanel2.SetBackgroundColour('#E9E3DB')
        
        innerSizer2.Add(lboxPanel2, 1, wx.EXPAND | wx.ALL)
        
        if 'linux' in sys.platform:
            innerSizer2.Add(btnPanel2, 1, wx.EXPAND | wx.ALL, 5)
        else:
            innerSizer2.Add(btnPanel2, 1, wx.EXPAND | wx.ALL)
        
        innerPanel2.SetSizer(innerSizer2)
        innerPanel2.SetBackgroundColour('#E9E3DB')
        
        outerSizer1.Add(innerPanel2, 1, wx.EXPAND | wx.ALL)
        outerSizer1.Add(innerPanel1, 1, wx.EXPAND | wx.ALL)
        
        outerPanel1.SetSizer(outerSizer1)
        outerPanel1.SetBackgroundColour('#E9E3DB')
        

        self.runCPAC1 = wx.Button(outerPanel2, -1, 'Run Individual Level Analysis')
        self.runCPAC1.Bind(wx.EVT_BUTTON, self.runIndividualAnalysis)
        

        self.runCPAC2 =  wx.Button(outerPanel2, -1, 'Run Group Level Analysis')
        self.runCPAC2.Bind(wx.EVT_BUTTON, self.runGroupLevelAnalysis)
        
#         outerSizer2.Add(self.runCPAC1, 1, wx.EXPAND | wx.RIGHT, 40)
#         outerSizer2.Add(self.runCPAC2, 1, wx.EXPAND | wx.LEFT, 40)

        outerSizer2.Add(self.runCPAC1, 1, wx.RIGHT, 20)
        outerSizer2.Add(self.runCPAC2, 1, wx.LEFT, 20)

        outerPanel2.SetSizer(outerSizer2)
        outerPanel2.SetBackgroundColour('#E9E3DB')
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        text1 = wx.StaticText(mainPanel, -1, "Configure CPAC")
        
        if 'linux' in sys.platform:
            text1.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.BOLD))
        else:
            text1.SetFont(wx.Font(18, wx.SWISS, wx.NORMAL, wx.BOLD))
            
        img = wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/cpac_new_logo.png'), wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        logo = wx.StaticBitmap(mainPanel, -1, img)
        hbox.Add(text1, 1, wx.TOP | wx.EXPAND, 15)
        hbox.Add(logo, 0,wx.ALIGN_RIGHT | wx.RIGHT)
        
        text2 = wx.StaticText(mainPanel, -1, "Run CPAC")
        
        if 'linux' in sys.platform:
            text2.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.BOLD))
        else:
            text2.SetFont(wx.Font(18, wx.SWISS, wx.NORMAL, wx.BOLD))
        
        mainSizer.Add(hbox, 0, wx.EXPAND | wx.ALL,10)
        mainSizer.Add(outerPanel1, 1, wx.EXPAND | wx.ALL, 20)
        mainSizer.Add(wx.StaticLine(mainPanel), 0, wx.EXPAND|wx.TOP|wx.BOTTOM, 10)
        mainSizer.Add(text2, 0, wx.EXPAND| wx.ALL, 5)
        mainSizer.Add(outerPanel2, 0 ,wx.EXPAND | wx.ALL, 20)
        
        mainPanel.SetSizer(mainSizer)

        self.Centre()
        self.Show(True)
        


    def runAnalysis1(self,pipeline, sublist, p):
        
        try:
            
            import CPAC
            CPAC.pipeline.cpac_runner.run(pipeline, sublist, p)
        
        except ImportError, e:
            wx.MessageBox("Error importing CPAC. %s"%e, "Error") 
            print "Error importing CPAC"
            print e



    def runIndividualAnalysis(self, event):

        try:

            if (self.listbox.GetChecked() or self.listbox.GetSelection()!= -1) and \
                (self.listbox2.GetChecked() or self.listbox2.GetSelection()!= -1):
                
                import thread
                import CPAC
                    
                pipelines = self.listbox.GetCheckedStrings()
                sublists = self.listbox2.GetCheckedStrings()
                    
                for s in sublists:
                    sublist = self.sublist_map.get(s)
                    for p in pipelines:
                        pipeline = self.pipeline_map.get(p)
                        print "running for configuration, subject list, pipeline_id -->", \
                              pipeline, sublist, p
                            
                        thread.start_new_thread(self.runAnalysis1, (pipeline, sublist, p))
                    
            else:

                errmsg = 'CPAC says: No pipeline or subject list ' \
                      'selected.'

                errSubID = wx.MessageDialog(self, errmsg, 'Error',
                    wx.OK | wx.ICON_ERROR)
                errSubID.ShowModal()
                errSubID.Destroy()

                print '\n\n[!] ' + errmsg + '\n\n'


                    

        except Exception, e:

            errSubID = wx.MessageDialog(self, str(e), 'Error',
                 wx.OK | wx.ICON_ERROR)
            errSubID.ShowModal()
            errSubID.Destroy()

            print e

                


    def runGroupLevelAnalysis(self, event):

        # Runs group analysis when user clicks "Run Group Level Analysis" in GUI

        print ""
        print "Running CPAC Group Analysis..."
        print ""
        
        if (self.listbox.GetChecked() or self.listbox.GetSelection()!= -1):
            
            pipelines = self.listbox.GetCheckedStrings()
            sublists = self.listbox2.GetCheckedStrings()

            for s in sublists:

                sublist = self.sublist_map.get(s)

                if not os.path.exists(sublist):
                    print '\n\nCPAC says: Please select a subject list ' \
                              'before running group-level analysis. ' \
                              'Thanks!\n\n'
                    raise Exception


                for p in pipelines:

                    pipeline = self.pipeline_map.get(p)

                    if os.path.exists(pipeline):
                        try:
                            import yaml
                            config = yaml.load(open(pipeline, 'r'))
                        except:
                                raise Exception("Error reading config file- %s", config)
                    
                        if config.get('outputDirectory'):
                            derv_path = os.path.join(config.get('outputDirectory'), 'pipeline_%s' % config.get('pipelineName')) #, '*', 'path_files_here' , '*.txt')
                        else:
                            derv_path = ''
                    
                        # Opens the sub-window which prompts the user
                        # for the derivative file paths
                        runGLA(pipeline, sublist, derv_path, p)

                    else:
                        print "pipeline doesn't exist"
                    
                
        else:
            print "No pipeline selected"



    def get_pipeline_map(self):
        return self.pipeline_map
    
    def get_sublist_map(self):
        return self.sublist_map
    
    def get_pipeline_path(self, idx):
        path = self.pipeline_map.get(idx)
        return path

    def NewItem(self, event):
        MainFrame(self, "save")

        
    def OnRename(self, event):
        sel = self.listbox.GetSelection()
        if sel!= -1:
            text = self.listbox.GetString(sel)
            renamed = wx.GetTextFromUser('Rename item', 'Rename dialog', text)
            if renamed != '':
                self.listbox.Delete(sel)
                self.listbox.Insert(renamed, sel)
                self.pipeline_map[renamed]= self.pipeline_map[text]
                del self.pipeline_map[text]
              
                
    def OnDelete(self, event):

        def removekey(d, key):
            r = dict(d)
            del r[key]
            return r

        selections = self.listbox.GetChecked()
        for sel in selections[::-1]:
            if sel != -1:
                text = self.listbox.GetString(sel)
                dlg = wx.MessageDialog(self, 'Do you also want to delete the configuration file from your system?',
                                       text, wx.YES | wx.NO | wx.ICON_INFORMATION)
    
                result = dlg.ShowModal()
                
                if  result == wx.ID_YES:
                    file_path = self.pipeline_map.get(text)
                
                    if file_path and os.path.exists(file_path):
                        os.remove(file_path)
                        
                self.listbox.Delete(sel)
                self.pipeline_map = removekey(self.pipeline_map, text)    
                
                dlg.Destroy()

        
        
    def OnEdit(self, event):
        
        # get 'sel' integer of listbox selection
        sel = self.listbox.GetSelection()
        
        if sel != -1:
            
            # 'text' - name of pipeline config displayed in listbox
            text = str(self.listbox.GetString(sel))
            
            # 'path' - path of pipeline_config.yml
            path = self.get_pipeline_path(text)
            
            if os.path.exists(path):
                # open the pipeline_config editor window
                MainFrame(self, option ="edit", path=path, pipeline_id = text)
            else:
                print "Couldn't find the config file %s "%path
     
     
     
    def OnLoad(self, event):

        # Does this code do anything?? If you are looking for the code that gets called
        # when you click Load, go to 'AddConfig'

        dlg = wx.FileDialog(
            self, message="Choose the config yaml file",
                defaultDir=os.getcwd(), 
                defaultFile="",
                wildcard= "YAML files(*.yaml, *.yml)|*.yaml;*.yml",
                style=wx.OPEN | wx.CHANGE_DIR)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
                
            if len(path)>0 and os.path.exists(path):
                MainFrame(self, option ="load", path=path)


            else:
                dlg = wx.MessageDialog(self, 'Invalid Path',
                                   'Error!',
                                   wx.OK | wx.ICON_ERROR)
                dlg.Destroy()
            
            dlg.Destroy()
                
                
    def OnClear(self, event, flag):
        if flag ==1:
            self.listbox.Clear()
            self.pipeline_map.clear()
        elif flag ==2:
            self.listbox2.Clear()
            self.sublist_map.clear()

    def CreateItem(self, event):
        DataConfig(self)


    def OnShow(self, event):
        import wx.lib.dialogs
        sel = self.listbox2.GetSelection()
        if sel != -1:
            text = self.listbox2.GetString(sel)
            name = self.sublist_map.get(text)
            if name:
                try:
                    f = open(name, "r")
                    msg = f.read()
                    f.close()
                    dlg = wx.lib.dialogs.ScrolledMessageDialog(self, msg, name, size = (800,800))
                    dlg.ShowModal()
                except:
                    print "Cannot open file %s"%(name)
                    
                    
    def OnDisplay(self, event):
        import wx.lib.dialogs
        sel = self.listbox.GetSelection()
        if sel != -1:
            text = self.listbox.GetString(sel)
            name = self.pipeline_map.get(text)
            if name:
                try:
                    f = open(name, "r")
                    msg = f.read()
                    f.close()
                    dlg = wx.lib.dialogs.ScrolledMessageDialog(self, msg, name, size = (800,1000))
                    dlg.ShowModal()
                except:
                    print "Cannot open file %s"%(name)
    
    
    def AddItem(self, event):
        
        dlg = wx.FileDialog(
            self, message="Choose the CPAC Subject list file",
            defaultDir=os.getcwd(), 
            defaultFile="CPAC_subject_list.yml",
            wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
            style=wx.OPEN | wx.CHANGE_DIR)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            # Check tfor path existence
            if os.path.exists(path):
                path = os.path.realpath(path)
                # Try and load in file contents
                try:
                    f_sl = yaml.load(open(path, 'r'))
                except Exception as e:
                    print 'Unable to load in the specified file: %s' % path
                    print 'Error:\n%s' % e
                # If it's not a list, we know it's not a subject list
                if type(f_sl) != list:
                    err_msg = 'File is not a subject list file. It might be a '\
                              'pipeline or data configuration file.'
                    raise Exception(err_msg)
            # Otherwise, report error
            else:
                err_msg = 'File %s does not exist. Check and try again.' \
                          % path
                raise Exception(err_msg)

            while True:
                dlg2 = wx.TextEntryDialog(self, 'Please enter a alias name for the Subject List',
                                     'Sublist Name', os.path.splitext(os.path.basename(path))[0])
            
                if dlg2.ShowModal() == wx.ID_OK:
                    if len(dlg2.GetValue()) >0:
                        if self.sublist_map.get(dlg2.GetValue()) == None:
                            self.sublist_map[dlg2.GetValue()]= path
                            self.listbox2.Append(dlg2.GetValue())
                            dlg2.Destroy()
                            dlg.Destroy()
                            break
                        else:
                            dlg3 = wx.MessageDialog(self, 'Subject List with this name already exist','Error!',
                              wx.OK | wx.ICON_ERROR)
                            dlg3.ShowModal()
                            dlg3.Destroy()
                else:
                    dlg2.Destroy()
                    dlg.Destroy
                    break
                            
                            
                            
    def check_config(self, config):

        ret_val = 1
        
        def display(win, msg, changeBg=True):
            wx.MessageBox(msg, "Error")
            if changeBg:
                win.SetBackgroundColour("pink")
            win.SetFocus()
            win.Refresh()

        try:
            c = yaml.load(open(config, 'r'))
        except:
            dlg = wx.MessageDialog(self, 'Error loading yaml file. Please check the file format',
                                           'Error!',
                                       wx.OK | wx.ICON_ERROR)
            ret_val = -1
            dlg.ShowModal()
            dlg.Destroy()

                
        # the following code checks the loaded pipeline config file for missing parameters (ex. if an old config file is used and new parameters
        # or features have been added) - if missing parameters are detected, it warns the user and informs them of the new defaults
        missingParams = []
        paramList = []

        try:
            params_file = open(p.resource_filename('CPAC', 'GUI/resources/config_parameters.txt'), "r")
        except:
            print "Error: Could not open configuration parameter file.", "\n"
            raise Exception

        paramInfo = params_file.read().split('\n')

        for param in paramInfo:

            if param != '':
                paramList.append(param.split(','))



        for param in paramList:

            try:
                if str(param[0]) not in c:
                    missingParams.append(param)
            except:
                errdlg = wx.MessageDialog(self, "Your pipeline " \
                                          "configuration file could not be " \
                                          "processed properly - please " \
                                          "ensure it is formatted properly." \
                                          "\n\nConfig file: %s" % config,
                                          "Error!",
                                       wx.OK | wx.ICON_ERROR)
                errdlg.ShowModal()
                errdlg.Destroy()
                break
                

        
        if missingParams:

            message = 'The following parameters are missing from your pipeline configuration file:\n\n'

            for param in missingParams:
                message = message + "\"" + str(param[1]) + "\"" + "\n" + "which can be found in tab:" + "\n" +  "\"" + str(param[2]) + "\"\n\n"

            message = message + "Click OK to review the defaults and accept or change the parameters."

            dlg = wx.MessageDialog(self, message, 'Missing Parameters', wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal()
            dlg.Destroy()

            if os.path.exists(config):
                MainFrame(self, option ="load", path=config)
            else:
                print "Couldn't find the config file %s "%config    

            ret_val = -1    

        return ret_val

    def AddConfig(self, event):
        
        # Gets called when you click 'Load' for pipeline config in the GUI
        dlg = wx.FileDialog(
            self, message="Choose the CPAC Configuration file",
            defaultDir=os.getcwd(), 
            defaultFile="",
            wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
            style=wx.OPEN | wx.CHANGE_DIR)
        
        if dlg.ShowModal() == wx.ID_OK:
            # Load config file into memory and verify its not a subject list
            path = dlg.GetPath()
            # Check for path existence
            if os.path.exists(path):
                path = os.path.realpath(path)
                try:
                    f_cfg = yaml.load(open(path, 'r'))
                except Exception as e:
                    print 'Unable to load in the specified file: %s' % path
                    print 'Error:\n%s' % e
                if type(f_cfg) == dict:
                    if not f_cfg.has_key('pipelineName'):
                        err_msg = 'File is not a pipeline configuration '\
                                  'file. It might be a data configuration file.'
                        raise Exception(err_msg)
                else:
                    err_msg = 'File is not a pipeline configuration '\
                              'file. It might be a subject list file.'
                    raise Exception(err_msg)
            # Otherwise, report error
            else:
                err_msg = 'File %s does not exist. Check and try again.' \
                          % path
                raise Exception(err_msg)
            if self.check_config(path) > 0:
                while True:
                    try:
                        c = Configuration(f_cfg)
                    except Exception as e:
                        print '\n\nERROR: Configuration file could not be '\
                              'loaded properly - the file might be '\
                              'access-protected or you might have chosen the '\
                              'wrong file.\n'
                        print 'Error name: main_window_0001\n\n'
                        print 'Exception: %s' % e
                    

                    if c.pipelineName != None:
                            
                            if self.pipeline_map.get(c.pipelineName) == None:

                                # this runs if you click 'Load' on the main
                                # CPAC window, enter a path, and the pipeline
                                # name attribute of the pipeline config file
                                # you are loading does NOT already exist in
                                # the listbox, i.e., the proper condition

                                self.pipeline_map[str(c.pipelineName)] = path
                                self.listbox.Append(str(c.pipelineName))
                                dlg.Destroy()
                                break
                            
                            else:

                                # this runs if you click 'Load' on the main
                                # CPAC window, enter a path, and the pipeline
                                # name attribute of the pipeline config file
                                # you are loading DOES already exist in
                                # the listbox, which is a conflict
                                       
                                dlg3 = wx.MessageDialog(self, 'The \'' \
                                        'Pipeline Name\' attribute of the ' \
                                        'configuration file you are loading' \
                                        ' already exists in one of the' \
                                        ' configuration files listed under' \
                                        ' \'Pipelines\'.\n\nPlease change' \
                                        ' the pipeline name attribute (not' \
                                        ' the filename) from within the' \
                                        ' pipeline editor (under the' \
                                        ' \'Output Settings\' tab in' \
                                        ' \'Environment Setup\'), or load a' \
                                        ' new configuration file.\n\n' \
                                        'Pipeline configuration with' \
                                        ' conflicting name:\n%s' \
                                         % c.pipelineName,
                                               'Conflicting Pipeline Names',
                                           wx.OK | wx.ICON_ERROR)
                                dlg3.ShowModal()
                                dlg3.Destroy()
                                break
                                
                    else:
                        
                        dlg4 = wx.MessageDialog(self, 'Warning: Pipeline name is blank.\n\nPlease edit' \
                                                ' the pipeline_config.yml file in a text editor and' \
                                                ' restore the pipelineName field.',
                                        'Warning',
                                wx.OK | wx.ICON_ERROR)
                        dlg4.ShowModal()
                        dlg4.Destroy()
                        
                        dlg.Destroy
                        break



    '''
    def AddConfigToList(self, config_path):

        path = dlg.GetPath()
        if self.check_config(path) > 0:
            while True:
                    
                try:
                    c = Configuration(yaml.load(open(os.path.realpath(path), 'r')))
                except:
                    print "\n\n" + "ERROR: Configuration file could not be loaded properly - the file " \
                          "might be access-protected or you might have chosen the wrong file." + "\n"
                    print "Error name: main_window_0001" + "\n\n"
                    raise Exception
                    

                if c.pipelineName != None:
                            
                        if self.pipeline_map.get(c.pipelineName) == None:
                            self.pipeline_map[str(c.pipelineName)] = path
                            self.listbox.Append(str(c.pipelineName))
                            dlg.Destroy()
                            break
                            
                        else:        
                                       
                            dlg3 = wx.MessageDialog(self, 'Pipeline already exists. Please enter a new configuration file.',
                                           'Error!',
                                       wx.OK | wx.ICON_ERROR)
                            dlg3.ShowModal()
                            dlg3.Destroy()
                                
                else:
                        
                    dlg4 = wx.MessageDialog(self, 'Warning: Pipeline name is blank.\n\nPlease edit' \
                                            ' the pipeline_config.yml file in a text editor and' \
                                            ' restore the pipelineName field.',
                                    'Warning',
                            wx.OK | wx.ICON_ERROR)
                    dlg4.ShowModal()
                    dlg4.Destroy()
                        
                    dlg.Destroy()
                    break
    '''

                    
                             
class runCPAC(wx.Frame):
 
    def __init__(self, pipeline, sublist, p, pid):
        wx.Frame.__init__(self, None, wx.ID_ANY, "Running CPAC")
                
        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, wx.ID_ANY)
        log = wx.TextCtrl(panel, wx.ID_ANY, size=(300,100),
                          style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
        btn = wx.Button(panel, wx.ID_ANY, 'Kill CPAC!')
        self.Bind(wx.EVT_BUTTON,lambda event: self.onButton(event, pid), btn)
 
        # Add widgets to a sizer        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(log, 1, wx.ALL|wx.EXPAND, 5)
        sizer.Add(btn, 0, wx.ALL|wx.CENTER, 5)
        panel.SetSizer(sizer)
 
        log.AppendText("running for pipeline --> %s \n"%p)
        log.AppendText("pipeline config --> %s \n"%pipeline)
        log.AppendText("subject list --> %s \n"%sublist)
        log.AppendText("process ids ---> %s \n"%pid)

       
    def onButton(self, event, pid):        
        if pid:
            for id in pid:
                print "killing process id -%s"% id
                os.kill(id, 9)
                print "please restart the gui"
        
            self.Close()
        else:
            dlg = wx.MessageDialog(self, 'Unable to retrieve the process id. Please kill the process from command prompt',
                                   'Alert!',
                                   wx.OK | wx.ICON_INFORMATION)
            dlg.Destroy()
       
            
class runGLA(wx.Frame):
    
    # Opens sub window prompting user to input the derivative file path(s).
    # If this is already supplied to the function, the path will show up
    # in the input box already by default.

    # Once the user clicks "Run", group level analysis begins

    def __init__(self, pipeline, sublist, path, name):
        wx.Frame.__init__(self, None, wx.ID_ANY, "Run Group Level Analysis for Pipeline - %s"%name, size = (730,120))
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel = wx.Panel(self)
        
        flexsizer = wx.FlexGridSizer(cols=2, hgap=5, vgap=10)

        img = wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/help.png'), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
       
        label1 = wx.StaticText(panel, -1, label = 'Pipeline Output Directory ')
        self.box1 = FileSelectorCombo(panel, id = wx.ID_ANY,  size = (500, -1))
        self.box1.GetTextCtrl().SetValue(str(path))
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        help1 = wx.BitmapButton(panel, id=-1, bitmap=img,
                                 pos=(10, 20), size = (img.GetWidth()+5, img.GetHeight()+5))
        help1.Bind(wx.EVT_BUTTON, self.OnShowDoc)
        
        hbox1.Add(label1)
        hbox1.Add(help1)
        
        flexsizer.Add(hbox1)
        flexsizer.Add(self.box1, flag = wx.EXPAND | wx.ALL)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        
        button3 = wx.Button(panel, wx.ID_CANCEL, 'Cancel', size =(120,30))
        button3.Bind(wx.EVT_BUTTON, self.onCancel)
        
        button2 = wx.Button(panel, wx.ID_OK, 'Run', size= (120,30))
        button2.Bind(wx.EVT_BUTTON, lambda event: \
                         self.onOK(event, sublist, pipeline) )
        
        hbox.Add(button3, 1, wx.EXPAND, border =5)
        hbox.Add(button2, 1, wx.EXPAND, border =5)
        
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(hbox,0, wx.ALIGN_CENTER, 5)
        panel.SetSizer(sizer)
        
        self.Show()
        
    def onCancel(self, event):
        self.Close()
        
    def runAnalysis(self, pipeline, sublist, path):
        try:
            import CPAC
            CPAC.pipeline.cpac_group_runner.run(pipeline, sublist, path)
        except AttributeError as e:
            print "Exception while running cpac_group_runner"
            #print "%d: %s" % (e.errno, e.strerror)
            print "Error type: ", type(e)
            print "Error args: ", e.args
            print "e: ", e
            print ""
            raise Exception
            
        
    def onOK(self, event, sublist, pipeline):

        # Once the user clicks "Run" in the derivative path file window
        # (from runGLA function), get the filepath and run the
        # "runAnalysis" function
        
        import thread

        if self.box1.GetValue():
            thread.start_new(self.runAnalysis, (pipeline, sublist, self.box1.GetValue()))
            self.Close()
        else:
            wx.MessageBox("Please provide the path to the output directory for the pipeline you want to run group-level analysis for.")
            
    def OnShowDoc(self, event):
        wx.TipWindow(self, "Path to output directory of the pipeline you wish to run group-level analysis for.", 500)

