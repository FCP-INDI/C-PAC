import wx
from config_window import MainFrame
from dataconfig_window import DataConfig
from ..utils.custom_control import FileSelectorCombo
from ..utils.constants import multiple_value_wfs
import wx.lib.agw.aquabutton as AB
import os
import pkg_resources as p
import sys
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

class ListBox(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(700, 650),  style= wx.SYSTEM_MENU | wx.CAPTION | wx.CLOSE_BOX)
        
        self.CreateStatusBar()
        self.SetStatusText("The Configurable Pipeline for the Analysis of Connectomes (C-PAC)")
    
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
                    
                    pipelines = self.listbox.GetCheckedStrings()
                    sublists = self.listbox2.GetCheckedStrings()
                    
                    #self.runCPAC1.SetPulseOnFocus(True)
                    
                    import CPAC
                    
                    for s in sublists:
                        sublist = self.sublist_map.get(s)
                        for p in pipelines:
                            pipeline = self.pipeline_map.get(p)
                            print "running for configuration, subject list, pipeline_id -->", \
                                  pipeline, sublist, p
                            
                            thread.start_new_thread(self.runAnalysis1, (pipeline, sublist, p))

                            
                            import time
                            time.sleep(20)
                            
                            try:
                                import yaml
                                config = yaml.load(open(pipeline, 'r'))
                            except:
                                raise Exception("Error reading config file- %s", config)
                            
                            try:
                                if config.get('outputDirectory'):
                                    pid = [ int(id.strip()) for id in open(os.path.join(config.get('outputDirectory'),\
                                       'pid.txt')).readlines()]
                            except ImportError:
                                print "unable to find the file %s"% os.path.join(config.outputDirectory, 'pid.txt')
                                pid = None
                            except Exception:
                                print "Unable to retrieve process id"
                                pid = None
                            
                            runCPAC(pipeline, sublist, p, pid).Show()

                            #print "Pipeline %s successfully ran for subject list %s"%(p,s)
                    
                else:
                    print "no pipeline and subject list selected"
                    

        except Exception, e:
                print e
                #wx.MessageBox(e, "Error") 
                
    def runGroupLevelAnalysis(self, event):
        print "running Group Analysis"
        
        if (self.listbox.GetChecked() or self.listbox.GetSelection()!= -1):
            pipelines = self.listbox.GetCheckedStrings()
            for p in pipelines:
                pipeline = self.pipeline_map.get(p)
                
                if os.path.exists(pipeline):
                    try:
                        import yaml
                        config = yaml.load(open(pipeline, 'r'))
                    except:
                            raise Exception("Error reading config file- %s", config)
                    
                    if config.get('outputDirectory'):
                        derv_path = os.path.join(config.get('outputDirectory'), 'pipeline_*', '*', 'path_files_here' , '*.txt')
                    else:
                        derv_path = ''
                    
                    runGLA(pipeline, derv_path, p)
                
                else:
                    print "pipeline doesn't exist"
                    
                
        else:
            print "No pipeline selected"

    def get_pipeline_map(self):
        return self.pipeline_map
    
    def get_sublist_map(self):
        return self.sublist_map
    
    def get_pipeline_path(self, id):
        path = self.pipeline_map.get(id)
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
        
        sel = self.listbox.GetSelection()
        if sel != -1:
            text = self.listbox.GetString(sel)
            path = self.get_pipeline_path(text)
            if os.path.exists(path):
                MainFrame(self, option ="edit", path=path, pipeline_id = text)
            else:
                print "Couldn't find the config file %s "%path
     
    def OnLoad(self, event):
        
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
        
        try:
            import yaml
            c = yaml.load(open(config, 'r'))
        except:
            dlg = wx.MessageDialog(self, 'Error loading yaml file. Please check the file format',
                                           'Error!',
                                       wx.OK | wx.ICON_ERROR)
            ret_val = -1
            dlg.ShowModal()
            dlg.Destroy()
        else:
            for wf in multiple_value_wfs:
                if c.get(wf):
                    if len(c.get(wf))>1:
                        dlg = wx.MessageDialog(self, "Configuration with multiple pipeline is not yet accepted by the gui. The multiple"\
                                                      "pipeline is due to workflow - %s"%wf,
                                           'Error!',
                                       wx.OK | wx.ICON_ERROR)
                        dlg.ShowModal()
                        dlg.Destroy()
                        ret_val = -1
            
        return ret_val
    
                            
    def AddConfig(self, event):
        
        dlg = wx.FileDialog(
            self, message="Choose the CPAC Configuration file",
            defaultDir=os.getcwd(), 
            defaultFile="CPAC_subject_list.yml",
            wildcard="YAML files(*.yaml, *.yml)|*.yaml;*.yml",
            style=wx.OPEN | wx.CHANGE_DIR)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if self.check_config(path) > 0:
                while True:
                    dlg2 = wx.TextEntryDialog(self, 'Please enter a unique pipeline id for the configuration',
                                             'Pipeline Id', "")
                    if dlg2.ShowModal() == wx.ID_OK:
                        if len(dlg2.GetValue()) >0:
                            
                                
                            if self.pipeline_map.get(dlg2.GetValue()) == None:
                                self.pipeline_map[dlg2.GetValue()] = path
                                self.listbox.Append(dlg2.GetValue())
                                dlg2.Destroy()
                                dlg.Destroy()
                                break
                            else:        
                                       
                                dlg3 = wx.MessageDialog(self, 'Pipeline already exist. Please enter a new name',
                                               'Error!',
                                           wx.OK | wx.ICON_ERROR)
                                dlg3.ShowModal()
                                dlg3.Destroy()
                    else:
                        dlg2.Destroy()
                        dlg.Destroy
                        break
       
              
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
    
    def __init__(self, pipeline, path, name):
        wx.Frame.__init__(self, None, wx.ID_ANY, "Run Group Level Analysis for Pipeline - %s"%name, size = (680,120))
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel = wx.Panel(self)
        
        flexsizer = wx.FlexGridSizer(cols=2, hgap=5, vgap=10)

        img = wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/help.png'), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
       
        label1 = wx.StaticText(panel, -1, label = 'Derivative Path File ')
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
                         self.onOK(event, pipeline) )
        
        hbox.Add(button3, 1, wx.EXPAND, border =5)
        hbox.Add(button2, 1, wx.EXPAND, border =5)
        
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(hbox,0, wx.ALIGN_CENTER, 5)
        panel.SetSizer(sizer)
        
        self.Show()
        
    def onCancel(self, event):
        self.Close()
        
    def runAnalysis(self, pipeline, path):
        try:
            import CPAC
            CPAC.pipeline.cpac_group_runner.run(pipeline, path)
        except Exception:
            print "Exception while running cpac_group_runner"
            
        
    def onOK(self, event, pipeline):
        
        import thread
        
        if self.box1.GetValue():
            thread.start_new(self.runAnalysis, (pipeline, self.box1.GetValue()))
            self.Close()
        else:
            wx.MessageBox("Please provide the path for the file containing output derivative path for each subject.")
            
    def OnShowDoc(self, event):
        wx.TipWindow(self, "Path to file containing derivative path. \n\nThis should be a text file with one path to derivative per line.", 500)
