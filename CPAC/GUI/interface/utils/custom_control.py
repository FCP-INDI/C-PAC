import wx
import wx.combo
import os
from wx.lib.masked import NumCtrl
import modelconfig_window, modelDesign_window
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


class FSLModelSelectorCombo(wx.combo.ComboCtrl):
    def __init__(self, *args, **kw):
        wx.combo.ComboCtrl.__init__(self, *args, **kw)
        bmp = wx.BitmapFromImage(wx.Image(p.resource_filename('CPAC', 'GUI/resources/images/folder3.gif')))
        self.SetButtonBitmaps(bmp, False)
        
    # Overridden from ComboCtrl, called when the combo button is clicked
    def OnButtonClick(self):
        path = ""
        name = ""
        wildcard = "YAML files (*.yml,*.yaml)|*.yml;*.yaml"
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
        
        button1 = wx.Button(panel, -1, 'Create or Load FSL Model', size= (210,50))
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
        self.box2 = FSLModelSelectorCombo(panel, id = wx.ID_ANY,  size = (500, -1))
        
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




class ContrastsFrame(wx.Frame):

    def __init__(self, parent, values, avail_cons):

        wx.Frame.__init__(self, parent, title="Add Contrast Description", \
                size = (300,80))
        
        self.avail_cons = avail_cons

        panel = wx.Panel(self)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)
        
        label1 = wx.StaticText(panel, -1, label = 'Contrast')
        self.box1 = wx.TextCtrl(panel, id=wx.ID_ANY, size=(200,-1))
    
        flexsizer.Add(label1)
        flexsizer.Add(self.box1,0,wx.ALIGN_RIGHT, 5)      
        
        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)
        
        self.Show()
    

    def onButtonClick(self, event):

        parent = self.Parent

        add_con = 0
        
        if self.box1.GetValue():
          
            val = self.box1.GetValue()

            # do validation first

            contrasts_in_string = self.parse_contrast(val)

            for contrast in contrasts_in_string:

                if contrast not in self.avail_cons:

                    errmsg = 'CPAC says: The contrast \'%s\' you ' \
                        'entered within the string \'%s\' is not ' \
                        'one of the available contrast selections.' \
                        '\n\nPlease enter only the contrast labels ' \
                        'listed under \'Available Contrasts\'.' \
                        % (contrast, val)

                    errSubID = wx.MessageDialog(self, errmsg,
                        'Invalid Contrast', wx.OK | wx.ICON_ERROR)
                    errSubID.ShowModal()
                    errSubID.Destroy()

                    add_con += 1


            if add_con == 0:

                parent.listbox.Append(str(val))
                parent.options.append(str(val))
                parent.raise_listbox_options()
                self.Close()



    def parse_contrast(self, contrast_string):

        orig_string = contrast_string

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

            errmsg = 'CPAC says: The contrast \'%s\' did not contain any ' \
                     'valid operators.\n\nValid operators: > , < , + , -' \
                     % orig_string

            errCon = wx.MessageDialog(self, errmsg, 'Invalid Operator',
                         wx.OK | wx.ICON_ERROR)
            errCon.ShowModal()
            errCon.Destroy()



        # in the case of the '+' or '-' contrast operators, which result in
        # the split_contrast list containing a blank element ''
        for item in split_contrast:
            if item == '':
                split_contrast.remove(item)

        return split_contrast



class f_test_frame(wx.Frame):
    
    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, title="Select Contrasts for f-Test", size = (280,200))
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

        if len(self.ctrl.GetCheckedStrings()) < 2:
            
            errmsg = "Please select at least two contrasts for your f-test."
            errCon = wx.MessageDialog(self, errmsg, "Not Enough Contrasts",
                     wx.OK | wx.ICON_ERROR)
            errCon.ShowModal()
            errCon.Destroy()

            raise Exception


        if self.ctrl.GetCheckedStrings():
            val=""
            for sel in self.ctrl.GetCheckedStrings():
                if val:
                    val = val + "," + sel
                else:
                    val = sel
            parent.listbox.Append(val)
            parent.options.append(val)
            self.Close()





class ListBoxCombo(wx.Panel):
    
    def __init__(self, parent, size, validator, style, values, combo_type):
        wx.Panel.__init__(self, parent)
        
        self.parent = parent

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

        self.options = []
        #self.listbox_selections = []

        # if it is the Contrasts checklist box in the group analysis model
        # builder GUI
        if self.ctype == 4:

            # if this is a 'load' situation when the user loads their group
            # analysis .yml file and they already have contrasts inserted into
            # the list
            if values:

                selected_contrasts = []

                for val in values:

                    # insert the contrast strings into the GUI's checkbox list
                    self.listbox.Append(str(val))
                    self.options.append(str(val))

                # select the contrast checkboxes that the user checked when
                # they load their gpa config.yml file into the model builder
                #     actually- they should all automatically be checked
                #     because the GUI only saves the ones that are checked
                self.listbox.SetCheckedStrings(values)

                # have to do this to make sure the f-test option list populates
                self.raise_listbox_options()


        if self.ctype == 5:

            # if this is a 'load' situation when the user loads their group
            # analysis .yml file and they already have f-tests inserted into
            # the list
            if values:

                selected_ftests = []

                for val in values:

                    # insert the f-test strings into the GUI's checkbox list
                    self.listbox.Append(str(val))
                    self.options.append(str(val))

                # select the contrast checkboxes that the user checked when
                # they load their gpa config.yml file into the model builder
                self.listbox.SetCheckedStrings(values)

                   
        
    def onButtonClick(self, event):

        if self.ctype == 1:
            CheckBox(self, self.values)
        elif self.ctype == 2:    
            TextBoxFrame(self, self.values)
        elif self.ctype == 3:
            ConfigFslFrame(self, self.values)
        elif self.ctype == 4:
            ContrastsFrame(self, self.values, self.avail_cons)
        elif self.ctype == 5:         

            # here: get the contrasts.csv and populate "self.parent.input_contrasts"
            # if custom_contrasts is a thing:

            input_contrasts = self.parent.input_contrasts

            #custom_cons = self.parent.get_custom_contrasts()

            #if len(custom_cons) > 0:
            #    input_contrasts = custom_cons

            if len(input_contrasts) < 2:

                errmsg = "Please input at least two contrasts before " \
                             "attempting to enter f-test selections."
                errCon = wx.MessageDialog(self, errmsg, "Not Enough Contrasts",
                         wx.OK | wx.ICON_ERROR)
                errCon.ShowModal()
                errCon.Destroy()

            else:

                f_test_frame(self, input_contrasts)

        
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


    def get_listbox_options(self):
        return self.options


    def raise_listbox_options(self):
        self.parent.input_contrasts = self.options


    def set_available_contrasts(self, avail_cons):

        # this is the list of contrast names available to the user to be
        # placed into the contrast strings - this gets passed to
        # ContrastsFrame so it can do string checking immediately
        self.avail_cons = avail_cons


    #def get_listbox_selections(self):
    #    return self.listbox_selections
            
            
            
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
        
        
        
class CheckBoxGrid(wx.Panel):
    
    def __init__(self, parent, idx, values, size):
        wx.Panel.__init__(self, parent, id=idx, size=size)
        #wx.ScrolledWindow.__init__(self, parent, id=idx, size=size, style=wx.VSCROLL)
        
        self.scrollWin = wx.ScrolledWindow(self, pos=(0,25), size=(450,205), style=wx.VSCROLL)  #wx.SUNKEN_BORDER) #wx.SUNKEN_BORDER | wx.VSCROLL)
        self.scrollWin.SetBackgroundColour(wx.WHITE)


        self.values = []
        self.values = values
        
        #wx.StaticText(self, label="Include EV", pos=(250,0))
        cat_label = wx.StaticText(self, label="Categorical", pos=(300,0))
        demean_label = wx.StaticText(self, label="Demean", pos=(400,0))
        

        j = 0
        self.idx = 100
        self.includeCBList = []
        self.categoricalCBList = []
        self.demeanCBList = []
        
        self.includeCBValues = []
        self.categoricalCBValues = []
        self.demeanCBValues = []
        
        self.cbDict = {}
        self.cbValuesDict = {}
        
        # this is for saving checkbox states during operation
        # for returning previous values when boxes come out from
        # being grayed out
        self.tempChoiceDict = {}
        
        # dictionary of 3 lists, each list is a list of chars
        # of either '1' or '0', a list for include, categorical,
        # demean
        self.choiceDict = {}

     
        #set_checkbox_grid_values(self.values)

            
        #self.cbDict['include'] = self.includeCBList
        self.cbDict['categorical'] = self.categoricalCBList
        self.cbDict['demean'] = self.demeanCBList


        '''
        for idNum in range(100,self.maxIDNum):
        
            self.Bind(wx.EVT_CHECKBOX, lambda event: self.onCheck_UpdateValue(event, idNum, wx.FindWindowById(idNum)), wx.FindWindowById(idNum))
        '''



        
    def set_checkbox_grid_values(self, value_list):

        j = 0
        self.idx = 100
        self.includeCBList = []
        self.categoricalCBList = []
        self.demeanCBList = []
        
        self.includeCBValues = []
        self.categoricalCBValues = []
        self.demeanCBValues = []

        self.maxIDNum = (len(value_list)*2)+101

        # clear the checkbox grid panel in case its already populated
        self.scrollWin.DestroyChildren()


        self.grid_sizer = wx.GridBagSizer(wx.VERTICAL)

        self.scrollWin.SetSizer(self.grid_sizer)
        self.scrollWin.SetScrollRate(10,10)
        self.scrollWin.EnableScrolling(True,True)


        row_panel = wx.Panel(self.scrollWin,wx.HORIZONTAL)
        row_panel.SetBackgroundColour(wx.WHITE)

        # iterate over each phenotype header item
        for name in value_list:
           
            # set up the label of each header item
            EV_label = wx.StaticText(row_panel, label=name, pos=(5,j))
                  

            # Categorical checkbox for header item
            self.cb = wx.CheckBox(row_panel, id=self.idx+1, pos=(300,j))
            self.cb.SetValue(False)
            self.categoricalCBList.append(self.cb)
            #self.row_sizer.Add(self.cb, pos=(300,0))
            
            self.cbValuesDict[self.idx+1] = [name, 'categorical', False]
            
            self.cb.Bind(wx.EVT_CHECKBOX, lambda event: self.onCheck_UpdateValue(event))

            
            
            # Demean checkbox for header item
            self.cb = wx.CheckBox(row_panel, id=self.idx+2, pos=(400,j))#, style=wx.CHK_3STATE)
            self.cb.SetValue(False)
            self.demeanCBList.append(self.cb)

            
            self.cbValuesDict[self.idx+2] = [name, 'demean', False]
            
            self.cb.Bind(wx.EVT_CHECKBOX, lambda event: self.onCheck_UpdateValue(event))
                      
                
            # just a nice amount to space the checkboxes out by
            j += 30
            
            # increment IDs
            self.idx += 2



        # automatically include some of the pre-calculated measures from
        # individual-level analysis as labels in the Model Setup checkbox
        # to remind users that they can include these into the design formula

        meanFD_label = wx.StaticText(row_panel, label='MeanFD (demeaned)', pos=(5,j))
        
        meanFDJ_label = wx.StaticText(row_panel, label='MeanFD_Jenkinson (demeaned)', pos=(5,j+30))

        measure_mean_label = wx.StaticText(row_panel, label='Measure_Mean (demeaned)', pos=(5,j+60))

        custom_roi_mean_label = wx.StaticText(row_panel, label='Custom_ROI_Mean (demeaned)', pos=(5,j+90))

        # add the panel that contains all of the rows (labels and checkboxes)
        # to the grid sizer. the grid sizer is necessary for wxPython to know
        # when to provide a scrollbar in the scrollWin object
        self.grid_sizer.Add(row_panel, pos=(0,0))

        w,h = self.grid_sizer.GetMinSize()
        self.scrollWin.SetVirtualSize((w,h))




    def onReload_set_selections(self, ev_selections):

        self.choiceCategoricalList = []
        self.choiceDemeanList = []

        for cb_id in self.cbValuesDict.keys():
            
            cb_name = self.cbValuesDict[cb_id][0]


            if 'categorical' in ev_selections.keys():

                if (cb_name in ev_selections['categorical']) and (cb_id % 2 != 0):

                    cb = wx.FindWindowById(cb_id)
                    cb.SetValue(True)

                    self.choiceCategoricalList.append(cb_name)


            if 'demean' in ev_selections.keys():

                if (cb_name in ev_selections['demean']) and (cb_id % 2 == 0):

                    cb = wx.FindWindowById(cb_id)
                    cb.SetValue(True)

                    self.choiceDemeanList.append(cb_name)



        self.choiceDict['categorical'] = self.choiceCategoricalList                
        self.choiceDict['demean'] = self.choiceDemeanList


        



        
    def onCheck_UpdateValue(self, event):#, idNum):
               
        # somehow take in the self.cbValuesDict[idx] (name, column, value)
        # and then GetValue from all idx, and update value for that idx
        
        # then have another function which returns this entire dict
               
        self.choiceCategoricalList = []
        self.choiceDemeanList = []


        for idNum in range(101,self.maxIDNum):
               
            self.cbValuesDict[idNum][2] = wx.FindWindowById(idNum).GetValue()
        

        for idNum in range(101,self.maxIDNum):           
     
            if self.cbValuesDict[idNum][1] == 'categorical':
                
                if self.cbValuesDict[idNum][2] == True:
                    self.choiceCategoricalList.append(self.cbValuesDict[idNum][0])
                    
                    if ('demean',idNum+1) not in self.tempChoiceDict.keys():
                        self.tempChoiceDict[('demean',idNum+1)] = wx.FindWindowById(idNum+1).GetValue()

                    wx.FindWindowById(idNum+1).Set3StateValue(2)  # set demean to N/A

                else:
                    #self.choiceCategoricalList.append('0')
                    
                    #wx.FindWindowById(idNum+1).Set3StateValue(0)  # undo demean as N/A
                    if ('demean',idNum+1) in self.tempChoiceDict.keys():
                        wx.FindWindowById(idNum+1).SetValue(self.tempChoiceDict[('demean',idNum+1)])
                        del self.tempChoiceDict[('demean',idNum+1)]
                        #wx.FindWindowById(idNum+1).Set3StateValue(0)  # undo demean as N/A
                        #wx.FindWindowById(idNum+1).SetValue(self.tempChoiceDict[('demean',idNum+1)])
                    
            elif self.cbValuesDict[idNum][1] == 'demean':
                
                if self.cbValuesDict[idNum][2] == True:
                    self.choiceDemeanList.append(self.cbValuesDict[idNum][0])
                #else:
                    #self.choiceDemeanList.append('0')
                     
        
                             
        self.choiceDict['categorical'] = self.choiceCategoricalList                
        self.choiceDict['demean'] = self.choiceDemeanList


        
    def GetGridSelection(self):
        return self.choiceDict




