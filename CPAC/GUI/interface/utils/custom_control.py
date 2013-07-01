import wx
import wx.combo
import os
from wx.lib.masked import NumCtrl
import modelconfig_window
import wx.lib.agw.balloontip as BT
import pkg_resources as p

class FileSelectorCombo(wx.combo.ComboCtrl):
    def __init__(self, *args, **kw):
        wx.combo.ComboCtrl.__init__(self, *args, **kw)
        bmp = wx.BitmapFromImage(wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/folder3.gif')))
        self.SetButtonBitmaps(bmp, False)
        
    # Overridden from ComboCtrl, called when the combo button is clicked
    def OnButtonClick(self):
        path = ""
        name = ""
        wildcard = "CPAC files (*.gz,*.nii,*.txt,*.mat.*.cnf,*.sch,*.csv)|*gz;*.nii;*.txt;*.cnf;*.sch;*.mat;*.csv"
        if self.GetValue():
            path, name = os.path.split(self.GetValue())
        
        dlg = wx.FileDialog(self, "Choose File", path, name,
                           wildcard= wildcard, style=wx.FD_OPEN|wx.CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.SetValue(dlg.GetPath())
                self.GetTextCtrl().SetValue(dlg.GetPath())
                dlg.Destroy()
        except:
            pass
        
        self.SetFocus()


class DirSelectorCombo(wx.combo.ComboCtrl):
    
    def __init__(self, *args, **kw):
        wx.combo.ComboCtrl.__init__(self, *args, **kw)
        bmp = wx.BitmapFromImage(wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/folder7.gif')))
        self.SetButtonBitmaps(bmp, False)
        
    # Overridden from ComboCtrl, called when the combo button is clicked
    def OnButtonClick(self):
        import os
        
        dlg = wx.DirDialog(self, "Choose a directory:",
                                style= wx.DD_NEW_DIR_BUTTON,
                          defaultPath=os.getcwd())

        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.SetValue(dlg.GetPath())
                self.GetTextCtrl().SetValue(dlg.GetPath())
                dlg.Destroy()
        except:
            pass
        
        
        self.SetFocus()
        
class CheckBox(wx.Frame):
    
    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, title="Select Corrections", size = (200,230))
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self)
        self.ctrl = wx.CheckListBox(panel, id = wx.ID_ANY,
                                    choices = values)
        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(self.ctrl, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)
        
        self.Show()
    
    def onButtonClick(self,event):
        parent = self.Parent
        if self.ctrl.GetCheckedStrings():
            val=""
            for sel in self.ctrl.GetCheckedStrings():
                if val:
                    val = val+ "," + sel
                else:
                    val = sel
            parent.listbox.Append(val)
            self.Close()
        

class TextBoxFrame(wx.Frame):

    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, title="Enter Frequency Cutoffs (in Hz)", size = (300,140))
        
        panel = wx.Panel(self)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15) 
        
        label1 = wx.StaticText(panel, -1, label = 'Low-frequency cutoff')
        self.box1 = NumCtrl(panel, id = wx.ID_ANY, value= values[0],
                            integerWidth=2, fractionWidth = 3, 
                            allowNegative=False, allowNone = True)
        
    
        flexsizer.Add(label1)
        flexsizer.Add(self.box1,0,wx.ALIGN_RIGHT, 5)
        
        label2 = wx.StaticText(panel, -1, label = 'High-frequency cutoff')
        self.box2 = NumCtrl(panel, id = wx.ID_ANY, value= values[1],
                            integerWidth=2, fractionWidth = 3, 
                            allowNegative=False, allowNone = True)
        
        flexsizer.Add(label2, 0, wx.EXPAND, 2)
        flexsizer.Add(self.box2,0, wx.ALIGN_RIGHT, 5)
        
        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)
        
        self.Show()
    
    def onButtonClick(self,event):
        parent = self.Parent
        
        if self.box1.GetValue() and self.box2.GetValue():
            
            if self.box1.GetValue() >= self.box2.GetValue():
                dlg = wx.MessageDialog(self, 'Lower Bound should be less than Upper Bound',
                                       'Error!',
                                   wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
            else:
                val = [self.box1.GetValue() , self.box2.GetValue()]
                parent.listbox.Append(str(val))
                self.Close()
                          
    
class ConfigFslFrame(wx.Frame):
    
    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, title="Specify FSL Model and Subject List", size = (680,200))
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel = wx.Panel(self)
        
        button1 = wx.Button(panel, -1, 'Create New FSL Model', size= (160,50))
        button1.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(button1, 0, wx.ALIGN_CENTER|wx.TOP, border = 15)
        
        flexsizer = wx.FlexGridSizer(cols=2, hgap=5, vgap=10)
        img = wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/help.png'), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
#         hbox1 = wx.BoxSizer(wx.HORIZONTAL)
#         img = wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/help.png'), wx.BITMAP_TYPE_ANY).ConvertToBitmap()
#         help1 = wx.BitmapButton(panel, id=-1, bitmap=img,
#                                  pos=(10, 20), size = (img.GetWidth()+5, img.GetHeight()+5))
#         help1.Bind(wx.EVT_BUTTON, lambda event: \
#                          self.OnShowDoc(event, 1))
        
#         label1 = wx.StaticText(panel, -1, label = 'Model Directory ')
#         self.box1 = DirSelectorCombo(panel, id = wx.ID_ANY, size = (500, -1))
#         hbox1.Add(label1)
#         hbox1.Add(help1)
# 
#         flexsizer.Add(hbox1)
#         flexsizer.Add(self.box1,flag = wx.EXPAND | wx.ALL)
        
        label2 = wx.StaticText(panel, -1, label = 'FSL Model Config')
        self.box2 = FileSelectorCombo(panel, id = wx.ID_ANY,  size = (500, -1))
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        help2 = wx.BitmapButton(panel, id=-1, bitmap=img,
                                 pos=(10, 20), size = (img.GetWidth()+5, img.GetHeight()+5))
        help2.Bind(wx.EVT_BUTTON, lambda event: \
                         self.OnShowDoc(event, 2))
        
        hbox2.Add(label2)
        hbox2.Add(help2)
        flexsizer.Add(hbox2)
        flexsizer.Add(self.box2, flag = wx.EXPAND | wx.ALL)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        button3 = wx.Button(panel, wx.ID_CANCEL, 'Cancel', size =(120,30))
        button3.Bind(wx.EVT_BUTTON, self.onCancel)
        
        button2 = wx.Button(panel, wx.ID_OK, 'OK', size= (120,30))
        button2.Bind(wx.EVT_BUTTON, self.onOK)
        
        hbox.Add(button3, 1, wx.EXPAND, border =5)
        hbox.Add(button2, 1, wx.EXPAND, border =5)
        
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(hbox,0, wx.ALIGN_CENTER, 5)
        panel.SetSizer(sizer)
        
        self.Show()
        
    def onCancel(self, event):
        self.Close()
        
    def onButtonClick(self,event):
        modelconfig_window.ModelConfig(self)

    def onOK(self, event):
        parent = self.Parent
        print parent
        if self.box2.GetValue():
            val = str(self.box2.GetValue())
            parent.listbox.Append(val)
            self.Close()
        else:
             wx.MessageBox("Please provide the path to the fsl model config file.")
        
#         if self.box1.GetValue() and self.box2.GetValue():
#                 val = [str(self.box1.GetValue()) , str(self.box2.GetValue())]
#                 parent.listbox.Append(str(val))
#                 self.Close()
#         else:
#             wx.MessageBox("Please provide the path for the model directory and subject list.")

    def OnShowDoc(self, event, flag):
        if flag == 1:
            wx.TipWindow(self, "Full path to a directory containing files for a single FSL model. All models must include .con, .mat, and .grp files. Models in which an F-Test is specified must also include a .fts file.", 500)
        elif flag == 2:
            #wx.TipWindow(self, "Full path to a subject list to be used with this model.\n\nThis should be a text file with one subject per line.", 500)
            wx.TipWindow(self, "Full path to a CPAC FSL model configuration file to be used.\n\nFor more information, please refer to the user guide.", 500)

class ListBoxCombo(wx.Panel):
    
    def __init__(self, parent, size, validator, style, values, combo_type):
        wx.Panel.__init__(self, parent)
        
        self.ctype = combo_type
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.values = values
        self.listbox = wx.CheckListBox(self, id = wx.ID_ANY, size = size, style = wx.LB_HSCROLL)
        self.listbox.Bind(wx.EVT_LISTBOX_DCLICK, self.onHover)
        bmp = wx.Bitmap(p.resource_filename('CPAC', 'GUI/resources/images/plus12.jpg'), wx.BITMAP_TYPE_ANY)
        self.button = wx.BitmapButton(self, -1, bmp,size = (bmp.GetWidth(), bmp.GetHeight()))# size= (30,30))
        self.button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(self.listbox,wx.EXPAND | wx.ALL, 10)
        sizer.Add(self.button)
        self.SetSizer(sizer)
        
    def onButtonClick(self, event):
        if self.ctype == 3:
            ConfigFslFrame(self, self.values)
        elif self.ctype == 2:    
            TextBoxFrame(self, self.values)
        elif self.ctype == 1:
            CheckBox(self, self.values)
        
    def GetListBoxCtrl(self):
        return self.listbox

    def onHover(self, event):
        sel = self.listbox.GetSelection()
        if sel != -1:
            text = self.listbox.GetString(sel)
            tip = BT.BalloonTip(toptitle= "selected value", message = text,
                                tipstyle= BT.BT_LEAVE)
            tip.SetTarget(self)
            tip.SetBalloonColour(wx.WHITE)
            # Set the font for the balloon title
            tip.SetTitleFont(wx.Font(15, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False))
            # Set the colour for the balloon title
            tip.SetTitleColour(wx.BLACK)
            # Leave the message font as default
            tip.SetMessageFont(wx.Font(9, wx.FONTFAMILY_SWISS, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD, False))
            # Set the message (tip) foreground colour
            tip.SetMessageColour(wx.BLUE)
            
            
            
class ParametersCheckBox(wx.Frame):
    
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Select EV as measures from CPAC", size = (300,130))
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        panel = wx.Panel(self)
        self.ctrl = wx.CheckListBox(panel, id = wx.ID_ANY,
                                    choices = ['MeanFD', 'MeanFD_Jenkinson', 'MeanDVARS'])
        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(self.ctrl, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)
        
        self.Show()
    
    def onButtonClick(self,event):
        parent = self.Parent
        if self.ctrl.GetCheckedStrings():
            val=""
            for sel in self.ctrl.GetCheckedStrings():
                if val:
                    val = val+ "," + sel
                else:
                    val = sel
             
            if parent.GetValue():
                new_val = parent.GetValue() + "," + val
                parent.GetTextCtrl().SetValue(new_val)
            else:
                parent.GetTextCtrl().SetValue(val)
            self.Close()
    
            
# class TextBoxCombo(wx.Panel):
#     def __init__(self, parent, size, validator, style, values):
#         wx.Panel.__init__(self, parent)
#                
#         sizer = wx.FlexGridSizer(cols=2, hgap=5, vgap=10)
#         self.values = values
#         
#         self.ctrl1 = wx.TextCtrl(parent, 
#                                 id = wx.ID_ANY, 
#                                 style= style, 
#                                 validator = validator, 
#                                 size = wx.DefaultSize)
#         
#         bmp = wx.Bitmap(p.resource_filename('CPAC', 'GUI/resources/images/plus12.jpg'), wx.BITMAP_TYPE_ANY)
#         self.button = wx.BitmapButton(self, -1, bmp,size = (bmp.GetWidth(), bmp.GetHeight()))
#         self.button.Bind(wx.EVT_BUTTON, self.onButtonClick)
#         sizer.Add(self.ctrl1)#, wx.EXPAND | wx.ALL, 10)
#         sizer.Add(self.button)
#         self.SetSizer(sizer)
#         
#     def onButtonClick(self, event):
#         print "calling config Fsl Frame"
#         ParametersCheckBox(self, self.values)
# 
#     def GetTextCtrl(self):
#         return self.ctrl1
    
class TextBoxCombo(wx.combo.ComboCtrl):
    def __init__(self, *args, **kw):
        wx.combo.ComboCtrl.__init__(self, *args, **kw)
        bmp = wx.Bitmap(p.resource_filename('CPAC', 'GUI/resources/images/plus12.jpg'), wx.BITMAP_TYPE_ANY)
        self.SetButtonBitmaps(bmp, False)
        
    # Overridden from ComboCtrl, called when the combo button is clicked
    def OnButtonClick(self):
        ParametersCheckBox(self)

