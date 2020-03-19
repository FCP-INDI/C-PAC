import wx
import wx.combo
import os
from wx.lib.masked import NumCtrl
from . import modelconfig_window, modelDesign_window, fsl_flame_presets_window
import wx.lib.agw.balloontip as BT
import pkg_resources as p


class FileSelectorCombo(wx.combo.ComboCtrl):
    def __init__(self, *args, **kw):
        wx.combo.ComboCtrl.__init__(self, *args, **kw)
        bmp = wx.BitmapFromImage(wx.Image(p.resource_filename('CPAC', \
                                     'GUI/resources/images/folder3.gif')))
        self.SetButtonBitmaps(bmp, False)

    # Overridden from ComboCtrl, called when the combo button is clicked
    def OnButtonClick(self):
        path = ""
        name = ""
        wildcard = "CPAC files (*.gz,*.nii,*.txt,*.mat.*.cnf,*.sch,*.csv)|" \
                   "*gz;*.nii;*.txt;*.cnf;*.sch;*.mat;*.csv;*.tsv"
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
        bmp = wx.BitmapFromImage(wx.Image(p.resource_filename('CPAC', \
                                     'GUI/resources/images/folder3.gif')))
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
        bmp = wx.BitmapFromImage(wx.Image(p.resource_filename('CPAC', \
                                     'GUI/resources/images/folder7.gif')))
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
        wx.Frame.__init__(self, parent, title="Select Regressors", \
                              size = (200,230))
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


class StringBoxFrame(wx.Frame):

    def __init__(self, parent, values, title, label):

        wx.Frame.__init__(self, parent, title=title, \
                size = (300,80))

        self.values = values

        panel = wx.Panel(self)

        sizer = wx.BoxSizer(wx.VERTICAL)

        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)

        label1 = wx.StaticText(panel, -1, label=label)
        self.box1 = wx.TextCtrl(panel, id=wx.ID_ANY, size=(200,-1))

        flexsizer.Add(label1)
        flexsizer.Add(self.box1,0,wx.ALIGN_RIGHT, 5)

        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)

        self.Show()

    def onButtonClick(self,event):
        parent = self.Parent
        if self.box1.GetValue():
            parent.listbox.Append(self.box1.GetValue())
            self.Close()


class TextBoxFrame(wx.Frame):

    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, \
                              title="Enter Frequency Cutoffs (in Hz)", \
                              size = (300,140))

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
                dlg = wx.MessageDialog(self, "Lower Bound should be less " \
                                       "than Upper Bound",
                                       'Error!',
                                   wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
            else:
                val = [self.box1.GetValue() , self.box2.GetValue()]
                parent.listbox.Append(str(val))
                self.Close()


class BETCoordinateFrame(wx.Frame):
    def __init__(self, parent, values):
       wx.Frame.__init__(self, parent, \
                             title="Enter Center of gravity coordinates", \
                             size = (450,150))

       panel = wx.Panel(self)

       sizer = wx.BoxSizer(wx.VERTICAL)

       flexsizer = wx.FlexGridSizer(cols=3, hgap=10, vgap=15)

       label1 = wx.StaticText(panel, -1, label = 'x-coordinate')
       self.box1 = NumCtrl(panel, id = wx.ID_ANY, value= values[0],
                           integerWidth=2, fractionWidth = 3,
                           allowNegative=True, allowNone = True)


       flexsizer.Add(label1)
       flexsizer.Add(self.box1,0,wx.ALIGN_RIGHT, 5)

       label2 = wx.StaticText(panel, -1, label = 'y-coordinate')
       self.box2 = NumCtrl(panel, id = wx.ID_ANY, value= values[1],
                           integerWidth=2, fractionWidth = 3,
                           allowNegative=True, allowNone = True)

       flexsizer.Add(label2, 0, wx.EXPAND, 2)
       flexsizer.Add(self.box2,0, wx.ALIGN_LEFT, 5)

       label3 = wx.StaticText(panel, -1, label = 'z-coordinate')
       self.box3 = NumCtrl(panel, id = wx.ID_ANY, value= values[2],
                           integerWidth=2, fractionWidth = 3,
                           allowNegative=True, allowNone = True)

       flexsizer.Add(label3, 0, wx.EXPAND, 3)
       flexsizer.Add(self.box3,0, wx.ALIGN_LEFT, 5)

       button = wx.Button(panel, -1, 'OK', size= (90,30))
       button.Bind(wx.EVT_BUTTON, self.onButtonClick)
       sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
       sizer.Add(button,0, wx.ALIGN_CENTER)
       panel.SetSizer(sizer)

       self.Show()

    def onButtonClick(self,event):
       parent = self.Parent

       if type(self.box1.GetValue() or self.box2.GetValue() or self.box3.GetValue()) is not int:

           dlg = wx.MessageDialog(self, "the values should be int or long" \
                                      'Error!',
                                  wx.OK | wx.ICON_ERROR)
           dlg.ShowModal()
           dlg.Destroy()
       else:
               val = [self.box1.GetValue() , self.box2.GetValue(), self.box3.GetValue()]
               parent.listbox.Append(str(val))
               self.Close()

class ResampleNumBoxFrame(wx.Frame):

    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, \
                              title="Enter Resolution to Resample To", \
                              size = (350,100))

        panel = wx.Panel(self)

        sizer = wx.BoxSizer(wx.VERTICAL)

        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)

        label1 = wx.StaticText(panel, -1, label = 'Resolution (in mm)')
        self.box1 = NumCtrl(panel, id = wx.ID_ANY, value= values[0],
                            integerWidth=2, fractionWidth = 3,
                            allowNegative=False, allowNone = True)


        flexsizer.Add(label1)
        flexsizer.Add(self.box1,0,wx.ALIGN_RIGHT, 5)

        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)

        self.Show()

    def onButtonClick(self,event):
        parent = self.Parent

        if type(self.box1.GetValue()) is not float:
            dlg = wx.MessageDialog(self, "Resolution must be a decimal " \
                                   "value, such as 2.5 or 3.0.",
                                   'Error!',
                               wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            val = self.box1.GetValue()
            parent.listbox.Append(str(val))
            self.Close()


class FilepathBoxFrame(wx.Frame):

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Enter File Path", \
                              size=(520,90))

        panel = wx.Panel(self)

        sizer = wx.BoxSizer(wx.VERTICAL)

        flexsizer = wx.FlexGridSizer(cols=2, hgap=10, vgap=15)

        label1 = wx.StaticText(panel, -1, label = 'File Path')
        self.box1 = FileSelectorCombo(panel, id = wx.ID_ANY, size = (420, -1))


        flexsizer.Add(label1)
        flexsizer.Add(self.box1,0,wx.ALIGN_RIGHT, 5)

        button = wx.Button(panel, -1, 'OK', size= (90,30))
        button.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer.Add(button,0, wx.ALIGN_CENTER)
        panel.SetSizer(sizer)

        self.Show()

    def onButtonClick(self,event):

        import os

        parent = self.Parent

        if self.box1.GetValue():

            val = self.box1.GetValue()

            # check the file just in case. you have to be safe, you just
            # gotta be safe man

            if not os.path.isfile(str(val)):
                errmsg = "That is not a valid filepath."
                errCon = wx.MessageDialog(self, errmsg, "Invalid Filepath",
                         wx.OK | wx.ICON_ERROR)
                errCon.ShowModal()
                errCon.Destroy()

            else:
                parent.add_checkbox_grid_value(val)
                self.Close()


class ConfigFslFrame(wx.Frame):

    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent,
                          title="Specify FSL FEAT Model",
                          size = (680,200))
        sizer_vert = wx.BoxSizer(wx.VERTICAL)
        sizer_horz = wx.BoxSizer(wx.HORIZONTAL)
        panel = wx.Panel(self)

        button0 = wx.Button(panel, -1, 'Choose FEAT Model Preset',
                            size=(210, 50))
        button0.Bind(wx.EVT_BUTTON, self.onPresetClick)
        sizer_horz.Add(button0, 0, wx.ALIGN_CENTER | wx.TOP, border=15)

        button1 = wx.Button(panel, -1, 'FEAT Model Builder/Editor',
                                size= (210, 50))
        button1.Bind(wx.EVT_BUTTON, self.onButtonClick)
        sizer_horz.Add(button1, 1, wx.ALIGN_CENTER | wx.TOP, border=15)

        sizer_vert.Add(sizer_horz, 0, wx.ALIGN_CENTER, border=15)

        flexsizer = wx.FlexGridSizer(cols=2, hgap=5, vgap=10)
        img = wx.Image(p.resource_filename('CPAC',
                                           'GUI/resources/images/help.png'),
                                         wx.BITMAP_TYPE_ANY).ConvertToBitmap()

        label2 = wx.StaticText(panel, -1, label='FSL FEAT Model Config')
        self.box2 = FSLModelSelectorCombo(panel, id = wx.ID_ANY,
                                          size = (500, -1))

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        help2 = wx.BitmapButton(panel, id=-1, bitmap=img,
                                pos=(10, 20), size = (img.GetWidth()+5,
                                img.GetHeight()+5))
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

        sizer_vert.Add(flexsizer, 1, wx.EXPAND | wx.ALL, 10)
        sizer_vert.Add(hbox,0, wx.ALIGN_CENTER, 5)
        panel.SetSizer(sizer_vert)

        self.Show()

    def onCancel(self, event):
        self.Close()

    def onPresetClick(self, event):
        fsl_flame_presets_window.FlamePresetsOne(self)

    def onButtonClick(self, event):
        modelconfig_window.ModelConfig(self)

    def onOK(self, event):
        parent = self.Parent
        if self.box2.GetValue():
            val = str(self.box2.GetValue())
            parent.listbox.Append(val)
            self.Close()
        else:
             wx.MessageBox("Please provide the path to the fsl model "
                           "config file.")

    def OnShowDoc(self, event, flag):
        if flag == 1:
            wx.TipWindow(self, "Full path to a directory containing files "
                               "for a single FSL model. All models must "
                               "include .con, .mat, and .grp files. Models "
                               "in which an F-Test is specified must also "
                               "include a .fts file.", 500)
        elif flag == 2:
            wx.TipWindow(self, "Full path to a CPAC FSL FEAT group analysis "
                               "model configuration YAML file. You can "
                               "either create one using one of the presets, "
                               "or by using the model builder to build one "
                               "from scratch.\n\nFor more information, "
                               "please refer to the user guide.", 500)


class ContrastsFrame(wx.Frame):

    def __init__(self, parent, values, dmatrix_obj):

        wx.Frame.__init__(self, parent, title="Add Contrast Description",
                size = (300,80))

        self.dmatrix_obj = dmatrix_obj
        #self.avail_cons = avail_cons

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

        frame = self
        parent = self.Parent

        add_con = 0

        if self.box1.GetValue():

            from CPAC.GUI.interface.utils.modelDesign_window import check_contrast_equation

            val = self.box1.GetValue()

            ret = check_contrast_equation(frame, self.dmatrix_obj, val)

            if ret == 0:
                parent.listbox.Append(str(val))
                parent.options.append(str(val))
                parent.raise_listbox_options()
                self.Close()



class f_test_frame(wx.Frame):

    def __init__(self, parent, values):
        wx.Frame.__init__(self, parent, title="Select Contrasts for f-Test",
                              size = (280,200))
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
        self.listbox = wx.CheckListBox(self, id = wx.ID_ANY, size = size, \
                                           style = wx.LB_HSCROLL)
        self.listbox.Bind(wx.EVT_LISTBOX_DCLICK, self.onHover)
        bmp = wx.Bitmap(p.resource_filename('CPAC', \
                                          'GUI/resources/images/plus12.jpg'),\
                                          wx.BITMAP_TYPE_ANY)
        self.button = wx.BitmapButton(self, -1, bmp,size = (bmp.GetWidth(), \
                                                            bmp.GetHeight()))
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

                # have to do this to make sure the f-test option list
                # populates
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
            ContrastsFrame(self, self.values, self.dmatrix_obj)
        elif self.ctype == 5:

            # self.parent.input_contrasts will only be populated if the user
            # has input contrasts via ContrastsFrame
            input_contrasts = self.parent.input_contrasts

            if len(input_contrasts) < 2:

                errmsg = "Please input at least two contrasts before " \
                             "attempting to enter f-test selections."
                errCon = wx.MessageDialog(self, errmsg, \
                         "Not Enough Contrasts", wx.OK | wx.ICON_ERROR)
                errCon.ShowModal()
                errCon.Destroy()

            else:
                f_test_frame(self, input_contrasts)

        elif self.ctype == 6:
            ResampleNumBoxFrame(self, self.values)
        elif self.ctype == 7:
            # because we need a nice generic configurable checkbox list...
            StringBoxFrame(self, self.values, "Add Session Name", "Session")
        elif self.ctype == 8:
            StringBoxFrame(self, self.values, "Add Series Name", "Series")


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
            tip.SetTitleFont(wx.Font(15, wx.FONTFAMILY_SWISS, \
                                         wx.FONTSTYLE_NORMAL, \
                                         wx.FONTWEIGHT_BOLD, False))
            # Set the colour for the balloon title
            tip.SetTitleColour(wx.BLACK)
            # Leave the message font as default
            tip.SetMessageFont(wx.Font(9, wx.FONTFAMILY_SWISS, \
                                       wx.FONTSTYLE_NORMAL, \
                                       wx.FONTWEIGHT_BOLD, False))
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

    def set_design_matrix(self, design_matrix_obj):
        self.dmatrix_obj = design_matrix_obj


    #def get_listbox_selections(self):
    #    return self.listbox_selections



class ParametersCheckBox(wx.Frame):

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, \
                              title = "Select EV as measures from CPAC", \
                              size = (300,130))
        sizer = wx.BoxSizer(wx.VERTICAL)

        panel = wx.Panel(self)
        self.ctrl = wx.CheckListBox(panel, id = wx.ID_ANY,
                                    choices = ['MeanFD', 'MeanFD_Jenkinson', \
                                                   'MeanDVARS'])
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



class TextBoxCombo(wx.combo.ComboCtrl):
    def __init__(self, *args, **kw):
        wx.combo.ComboCtrl.__init__(self, *args, **kw)
        bmp = wx.Bitmap(p.resource_filename('CPAC', \
                                          'GUI/resources/images/plus12.jpg'),\
                                          wx.BITMAP_TYPE_ANY)
        self.SetButtonBitmaps(bmp, False)

    # Overridden from ComboCtrl, called when the combo button is clicked
    def OnButtonClick(self):
        ParametersCheckBox(self)



class CheckBoxGrid(wx.Panel):

    def __init__(self, parent, idx, selections, values, size):

        wx.Panel.__init__(self, parent, id=idx, size=size)

        all_sizer = wx.BoxSizer(wx.HORIZONTAL)
        window_sizer = wx.BoxSizer(wx.HORIZONTAL)
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)

        self.scrollWin = wx.ScrolledWindow(self, pos=(0,25), size=(600,205), \
                             style=wx.VSCROLL)
        self.scrollWin.SetBackgroundColour(wx.WHITE)

        self.selections = selections
        self.values = values

        # y-axis position of each new row (each new entry) in the checkbox
        # grid - this will change as entries are added
        self.y_pos = 0

        self.num_entries = 0

        # this sets the column names at the top of the checkbox grid
        #x_pos = 300
        x_pos = 5
        self.x_pos_increments = []

        for selection in self.selections:

            text = wx.StaticText(self, label=selection, pos=(x_pos, 0))

            spacing = 50 + text.GetSize().width - 20
            # spacing = 50 + (len(selection)*2) # adapt the spacing to
            #                                   # label length
            x_pos += spacing

            self.x_pos_increments.append(spacing)


        self.idx = 100

        self.grid_sizer = wx.GridBagSizer(wx.VERTICAL)

        self.scrollWin.SetSizer(self.grid_sizer)
        self.scrollWin.SetScrollRate(10,10)
        self.scrollWin.EnableScrolling(True,True)

        self.row_panel = wx.Panel(self.scrollWin,wx.HORIZONTAL)
        self.row_panel.SetBackgroundColour(wx.WHITE)

        bmp = wx.Bitmap(p.resource_filename('CPAC', \
                            'GUI/resources/images/plus12.jpg'), \
                            wx.BITMAP_TYPE_ANY)
        self.button = wx.BitmapButton(self, -1, bmp, \
                                      size=(bmp.GetWidth(), bmp.GetHeight()),\
                                      pos=(605,25))
        self.button.Bind(wx.EVT_BUTTON, self.onButtonClick)

        self.cbValuesDict = {}

        # this is for saving checkbox states during operation
        # for returning previous values when boxes come out from
        # being grayed out
        self.tempChoiceDict = {}

        # dictionary of 3 lists, each list is a list of chars
        # of either '1' or '0', a list for include, categorical,
        # demean
        self.choiceDict = {}

        # this dictionary holds all of the GUI controls associated with each
        # entry (row) in the checkbox grid. this will be used to appropriately
        # delete rows the user selects to delete
        self.entry_controls = {}


    def onButtonClick(self, event):

        FilepathBoxFrame(self)


    def onMinusClick(self, event, entry):

        for control in self.entry_controls[entry]:
            self.entry_controls[entry][control].Destroy()

        del self.entry_controls[entry]

        self.cbValuesDict = {}

        self.onCheck_UpdateValue()

        if entry in self.choiceDict:
            del self.choiceDict[entry]

        for row in self.entry_controls:
            for control in self.entry_controls[row]:
                self.entry_controls[row][control].Destroy()

        self.entry_controls = {}

        #self.idx = 100
        self.y_pos = 0

        self.onReload_set_selections(self.choiceDict)


    def add_checkbox_grid_value(self, entry):

        # set up the label of each header item
        #wx.StaticText(self.row_panel, label=entry, pos=(5,self.y_pos))

        #x_pos = 300
        x_pos = 5

        if entry not in self.entry_controls.keys():
            self.entry_controls[entry] = {}
        else:
            errmsg = "This file path has already been entered."
            errCon = wx.MessageDialog(self, errmsg, "ROI Path Exists",
                        wx.OK | wx.ICON_ERROR)
            errCon.ShowModal()
            errCon.Destroy()
            return -1

        for selection, spacing in zip(self.selections, self.x_pos_increments):

            # checkbox for possible selection
            self.cb = wx.CheckBox(self.row_panel, id=self.idx+1, \
                                      pos=(x_pos,self.y_pos))

            self.cb.SetValue(False)

            self.cbValuesDict[self.idx+1] = [entry, selection, False]

            self.cb.Bind(wx.EVT_CHECKBOX, \
                             lambda event: self.onCheck_UpdateValue(event))

            self.entry_controls[entry][selection] = self.cb

            # increment IDs
            self.idx += 1

            x_pos += spacing


        bmp = wx.Bitmap(p.resource_filename('CPAC', \
                        'GUI/resources/images/minus9.jpg'), \
                        wx.BITMAP_TYPE_ANY)
        self.button = wx.BitmapButton(self.row_panel, -1, bmp, \
                                    #size=(bmp.GetWidth(), bmp.GetHeight()),\
                                      pos=(x_pos,self.y_pos))
        self.button.Bind(wx.EVT_BUTTON, \
                             lambda event: self.onMinusClick(event, entry))

        self.path_text = wx.StaticText(self.row_panel, label=entry, \
                                           pos=(x_pos+30,self.y_pos))

        self.entry_controls[entry]["minus_button"] = self.button
        self.entry_controls[entry]["filepath_label"] = self.path_text


        # just a nice amount to space the rows out by
        self.y_pos += 30

        self.num_entries += 1

        # add the panel that contains all of the rows (labels and checkboxes)
        # to the grid sizer. the grid sizer is necessary for wxPython to know
        # when to provide a scrollbar in the scrollWin object
        try:
            self.grid_sizer.Add(self.row_panel, pos=(0,0))
        except:
            pass

        w,h = self.grid_sizer.GetMinSize()
        self.scrollWin.SetVirtualSize((w,h))


    '''
    def setup_checkbox_grid_values_for_gpa(self, value_list):

        # this is for when all of the checkbox grid row entries get populated
        # all at once from a pre-created list (i.e. group analysis model
        # builder)

        j = 0
        self.idx = 100
        self.includeCBList = []
        self.categoricalCBList = []
        self.demeanCBList = []

        self.includeCBValues = []
        self.categoricalCBValues = []
        self.demeanCBValues = []

        #self.maxIDNum = (len(value_list)*2)+101

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

            self.cbValuesDict[self.idx+1] = [name, 'categorical', False]

            self.cb.Bind(wx.EVT_CHECKBOX, \
                             lambda event: self.onCheck_UpdateValue(event))

            # Demean checkbox for header item
            self.cb = wx.CheckBox(row_panel, id=self.idx+2, pos=(400,j))
            self.cb.SetValue(False)
            self.demeanCBList.append(self.cb)


            self.cbValuesDict[self.idx+2] = [name, 'demean', False]

            self.cb.Bind(wx.EVT_CHECKBOX, \
                             lambda event: self.onCheck_UpdateValue(event))


            # just a nice amount to space the checkboxes out by
            j += 30

            # increment IDs
            self.idx += 2


        # automatically include some of the pre-calculated measures from
        # individual-level analysis as labels in the Model Setup checkbox
        # to remind users that they can include these into the design formula

        meanFDP_label = wx.StaticText(row_panel, \
                                      label='MeanFD_Power (demeaned)', \
                                      pos=(5,j))

        meanFDJ_label = wx.StaticText(row_panel, \
                                      label='MeanFD_Jenkinson (demeaned)', \
                                      pos=(5,j+30))

        measure_mean_label = wx.StaticText(row_panel, \
                                           label='Measure_Mean (demeaned)', \
                                           pos=(5,j+60))

        custom_roi_mean_label = wx.StaticText(row_panel, \
                                          label='Custom_ROI_Mean (demeaned)',\
                                          pos=(5,j+90))

        # add the panel that contains all of the rows (labels and checkboxes)
        # to the grid sizer. the grid sizer is necessary for wxPython to know
        # when to provide a scrollbar in the scrollWin object
        self.grid_sizer.Add(row_panel, pos=(0,0))

        w,h = self.grid_sizer.GetMinSize()
        self.scrollWin.SetVirtualSize((w,h))
    '''


    def onCheck_UpdateValue(self, event=None):

        # somehow take in the self.cbValuesDict[idx] (name, column, value)
        # and then GetValue from all idx, and update value for that idx

        for selection in self.selections:

            if selection not in self.choiceDict.keys():
                self.choiceDict[selection] = []

            for filepath in self.entry_controls.keys():

                # self.entry_controls[filepath] is another dictionary

                # and this is a checkbox
                cb_val = self.entry_controls[filepath][selection].GetValue()

                if cb_val == True:

                    if filepath not in self.choiceDict.keys():
                        self.choiceDict[filepath] = []

                    if selection not in self.choiceDict[filepath]:
                        self.choiceDict[filepath].append(selection)

                else:

                    if filepath in self.choiceDict.keys():

                        if selection in self.choiceDict[filepath]:
                            self.choiceDict[filepath].remove(selection)


        for roi_path in self.choiceDict.keys():

            if len(self.choiceDict[roi_path]) == 0:
                del self.choiceDict[roi_path]

        # after this is done, self.choiceDict will be a dictionary with a
        # format similar to:
        #     {"Avg": ["/path/to/ROI1.nii.gz","/path/to/ROI2.nii.gz"],
        #      "Voxel": ["/path/to/ROI1.nii.gz","/path/to/ROI3.nii.gz"]}


    def onReload_set_selections(self, choice_dict):

        for entry in choice_dict.keys():

            # re-populate the box with entries (but not selections)
            self.add_checkbox_grid_value(entry)

        # re-populate selections
        for entry in self.entry_controls:

            for ctrl in self.entry_controls[entry]:

                if ctrl in choice_dict[entry]:
                    cb = self.entry_controls[entry][ctrl]
                    cb.SetValue(True)

        # update the choiceDict
        self.onCheck_UpdateValue()


    def GetGridSelection(self):
        return self.choiceDict



class GPAModelCheckBoxGrid(wx.Panel):

    def __init__(self, parent, idx, values, size):

        wx.Panel.__init__(self, parent, id=idx, size=size)

        self.scrollWin = wx.ScrolledWindow(self, pos=(0,25), size=(450,205), \
                                               style=wx.VSCROLL)
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

        # dictionary of checkbox controls
        self.cb_dict = {}

        # dictionary of 3 lists, each list is a list of chars
        # of either '1' or '0', a list for include, categorical,
        # demean
        self.choiceDict = {}


        #set_checkbox_grid_values(self.values)


        #self.cbDict['include'] = self.includeCBList
        self.cbDict['categorical'] = self.categoricalCBList
        self.cbDict['demean'] = self.demeanCBList


    def set_checkbox_grid_values(self, value_list):

        j = 0
        self.idx = 100

        self.categoricalCBList = []
        self.demeanCBList = []

        self.categoricalCBValues = []
        self.demeanCBValues = []

        #self.maxIDNum = (len(value_list)*2)+101

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

            self.cb.Bind(wx.EVT_CHECKBOX, \
                             lambda event: self.onCheck_UpdateValue(event))

            if name not in self.cb_dict.keys():
                self.cb_dict[name] = {}

            self.cb_dict[name]["categorical"] = self.cb



            # Demean checkbox for header item
            self.cb = wx.CheckBox(row_panel, id=self.idx+2, pos=(400,j))
            self.cb.SetValue(False)
            self.demeanCBList.append(self.cb)


            self.cbValuesDict[self.idx+2] = [name, 'demean', False]

            self.cb.Bind(wx.EVT_CHECKBOX, \
                             lambda event: self.onCheck_UpdateValue(event))


            self.cb_dict[name]["demean"] = self.cb


            # just a nice amount to space the checkboxes out by
            j += 30

            # increment IDs
            self.idx += 2



        # automatically include some of the pre-calculated measures from
        # individual-level analysis as labels in the Model Setup checkbox
        # to remind users that they can include these into the design formula

        meanFDP_label = wx.StaticText(row_panel, \
                                          label='MeanFD_Power (demeaned)', \
                                          pos=(5,j))

        meanFDJ_label = wx.StaticText(row_panel, \
                                         label='MeanFD_Jenkinson (demeaned)',\
                                         pos=(5,j+30))

        measure_mean_label = wx.StaticText(row_panel, \
                                             label='Measure_Mean (demeaned)',\
                                             pos=(5,j+60))

        custom_roi_mean_label = wx.StaticText(row_panel, \
                                          label='Custom_ROI_Mean (demeaned)',\
                                          pos=(5,j+90))

        # add the panel that contains all of the rows (labels and checkboxes)
        # to the grid sizer. the grid sizer is necessary for wxPython to know
        # when to provide a scrollbar in the scrollWin object
        self.grid_sizer.Add(row_panel, pos=(0,0))

        w,h = self.grid_sizer.GetMinSize()
        self.scrollWin.SetVirtualSize((w,h))


    def onReload_set_selections(self, ev_selections):

        self.choiceDict["categorical"] = []
        self.choiceDict["demean"] = []

        for ev_name in self.cb_dict:

            if "categorical" in ev_selections.keys():

                if ev_name in ev_selections["categorical"]:

                    self.cb_dict[ev_name]["categorical"].SetValue(True)

                    self.choiceDict["categorical"].append(ev_name)

            if "demean" in ev_selections.keys():

                if ev_name in ev_selections["demean"]:

                    self.cb_dict[ev_name]["demean"].SetValue(True)

                    self.choiceDict["demean"].append(ev_name)


    def onCheck_UpdateValue(self, event):#, idNum):

        # somehow take in the self.cbValuesDict[idx] (name, column, value)
        # and then GetValue from all idx, and update value for that idx

        # then have another function which returns this entire dict

        for ev_name in self.cb_dict:

            for selection in self.cb_dict[ev_name]:

                cb = self.cb_dict[ev_name][selection]
                cb_val = cb.GetValue()

                if selection == "categorical":

                    if cb_val == True:

                        if "categorical" not in self.choiceDict.keys():
                            self.choiceDict["categorical"] = []

                        self.choiceDict["categorical"].append(ev_name)

                        if (ev_name,"demean") not in self.tempChoiceDict.keys():
                            self.tempChoiceDict[(ev_name, "demean")] = \
                                self.cb_dict[ev_name]["demean"].GetValue()

                        # set demean to N/A
                        self.cb_dict[ev_name]["demean"].Set3StateValue(2)

                    else:

                        if (ev_name,"demean") in self.tempChoiceDict.keys():
                            val = self.tempChoiceDict[(ev_name,"demean")]
                            self.cb_dict[ev_name]["demean"].SetValue(val)
                            del self.tempChoiceDict[(ev_name,"demean")]

                elif selection == "demean":

                    if cb_val == True:

                        if "demean" not in self.choiceDict.keys():
                            self.choiceDict["demean"] = []

                        self.choiceDict["demean"].append(ev_name)


    def GetGridSelection(self):
        return self.choiceDict
