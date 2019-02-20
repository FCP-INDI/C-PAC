from .custom_control import FileSelectorCombo, DirSelectorCombo, ListBoxCombo, TextBoxCombo, CheckBoxGrid, GPAModelCheckBoxGrid
from ..utils.constants import control as control_types
from ..utils.constants import dtype as data_types
from wx.lib import masked
from wx.lib.masked import NumCtrl
from wx.lib.intctrl import IntCtrl
import wx.lib.intctrl
import pkg_resources as p


class GenericClass(wx.ScrolledWindow):
    
    def __init__(self, parent, title="", static=True):
        self.parent = parent
        self.title = title
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.flexSizer = wx.FlexGridSizer(cols=2, hgap=15, vgap=15)
        self.flexSizer.AddGrowableCol(1)

        self.ctrl_list = []

        if static:
            self.__add_static()

        self.__set_scroll()
        self.switch = None
    
    def __set_scroll(self):
        maxWidth = 1000
        maxHeight = 1000
        
        self.parent.SetVirtualSize((maxWidth, maxHeight))
        self.parent.SetScrollRate(20,20)
        
    def add_static(self):
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        t = wx.StaticText(self.parent, -1, self.title)
        t.SetFont(wx.Font(18, wx.SWISS, wx.NORMAL, wx.BOLD))
        img_path = p.resource_filename('CPAC', 'GUI/resources/images/cpac_logo_2.jpg')
        img = wx.Image(img_path, wx.BITMAP_TYPE_JPEG).ConvertToBitmap()

        bmp = wx.StaticBitmap(self.parent, -1, img)
        
        hbox.Add(bmp)
        hbox.Add(t, 0, wx.TOP, 8)
        self.mainSizer.Add(hbox, proportion=0, flag=wx.ALL, border=5)
        self.mainSizer.Add(wx.StaticLine(self.parent), 0, wx.EXPAND|wx.TOP)
        
    __add_static = add_static

    def add_pheno_load_panel(self, sizer):
        
        buffer = wx.StaticText(self.parent, label="\t\t\t\t\t\t")

        self.flexSizer.Add(buffer)
        self.flexSizer.Add(sizer)

    def add(self, label, control, name="", type=0, 
            comment="", values="", style=0, size= wx.DefaultSize, 
            validator=wx.DefaultValidator, wkf_switch=False,
            validation_req=True, combo_type=None, selections=None):

        label = label.strip()
            
        label_text = wx.StaticText(self.parent, -1, label)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        help_img_path = p.resource_filename('CPAC', 'GUI/resources/images/help.png')
        help_img = wx.Image(help_img_path, wx.BITMAP_TYPE_ANY).ConvertToBitmap()

        button = wx.BitmapButton(self.parent, id=-1, bitmap=help_img,
                                 size=(help_img.GetWidth()+5, help_img.GetHeight()+5))

        button.Bind(wx.EVT_BUTTON, lambda event: \
                         self.OnShowDoc(event, comment))
        hbox.Add(label_text, flag=wx.CENTER | wx.ALL)
        hbox.AddSpacer(5)
        hbox.Add(button, flag=wx.CENTER | wx.ALL)
        self.flexSizer.Add(hbox)
        
        if control not in [v for k, v in control_types.__dict__.items() if not k.startswith('__')]:
            self.ctrl_list.append(control)
            self.type = control.__class__.__name__
            self.flexSizer.Add(control.get_ctrl(), flag=wx.CENTER | wx.EXPAND | wx.ALL)

        else:
            ctrl = Control(self.parent, name=name, 
                           type=control, datatype=type, 
                           values=values, size=size, 
                           validator=validator,
                           wkf_switch=wkf_switch,
                           help=comment, pretty_name=label,
                           validation_req=validation_req,
                           combo_type=combo_type,
                           selections=selections)
            
            self.ctrl_list.append(ctrl)
        
            if wkf_switch:
                self.switch = ctrl

            if control == control_types.CHOICE_BOX:
                self.parent.Bind(wx.EVT_CHOICE, lambda event: self.EvtChoice(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.TEXT_BOX:
                self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterBox(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.COMBO_BOX:
                self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterCombo(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.EXPAND | wx.ALL)
            elif control == control_types.INT_CTRL:
                self.parent.Bind(wx.lib.intctrl.EVT_INT, lambda event: self.TxtEnterBox(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.FLOAT_CTRL:
                self.parent.Bind(masked.EVT_NUM, lambda event: self.TxtEnterBox(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.DIR_COMBO_BOX:
                self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterCombo(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.EXPAND | wx.ALL)
            elif control == control_types.CHECKLIST_BOX:
                self.parent.Bind(wx.EVT_CHECKLISTBOX, lambda event: self.EvtCheckListBox(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.LISTBOX_COMBO:
                self.parent.Bind(wx.EVT_CHECKLISTBOX, lambda event: self.EvtListBoxCombo(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.TEXTBOX_COMBO:
                self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterCombo(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.EXPAND | wx.ALL)
            elif control == control_types.CHECKBOX_GRID:
                self.parent.Bind(wx.EVT_CHECKBOX, lambda event: self.EvtCheckBoxGrid(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            elif control == control_types.GPA_CHECKBOX_GRID:
                self.parent.Bind(wx.EVT_CHECKBOX, lambda event: self.EvtCheckBoxGrid(event, ctrl), id=ctrl.get_id())
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)
            else:
                self.flexSizer.Add(ctrl.get_ctrl(), flag=wx.CENTER | wx.ALL)

    def EvtChoice(self, event, ctrl):
        if type(event.GetString()) == unicode:
            value = event.GetString().encode('ascii', 'ignore')
        ctrl.set_selection(value)
        self.parent.Refresh()
        
    def TxtEnterBox(self, event, ctrl):
        ctrl.get_ctrl().SetBackgroundColour("white")
        ctrl.set_selection(ctrl.get_ctrl().GetValue())
        
    def TxtEnterCombo(self, event, ctrl):
        ctrl.text_ctrl.SetBackgroundColour("white")
        ctrl.set_selection(ctrl.text_ctrl.GetValue())

    def EvtCheckListBox(self, event, ctrl):
        index = event.GetSelection()
        label = ctrl.get_ctrl().GetString(index)
        if ctrl.get_ctrl().IsChecked(index):
            ctrl.set_selection(label, index)
        else:
            ctrl.set_selection(label,index, True)
    
    def EvtListBoxCombo(self, event, ctrl):
        
        index = event.GetSelection()
        label = ctrl.get_ctrl().GetListBoxCtrl().GetString(index)
        if ctrl.get_ctrl().GetListBoxCtrl().IsChecked(index):
            ctrl.set_selection(label, index)
        else:
            ctrl.set_selection(label,index, True)

    '''
    NEEDS DEV!
    '''
    def EvtCheckBoxGrid(self, event, ctrl):
        index = event.GetSelection()
        label = ctrl.get_ctrl().GetString(index)
        if ctrl.get_ctrl().IsChecked(index):
            ctrl.set_selection(label, index)
        else:
            ctrl.set_selection(label,index, True)
    
    def OnShowDoc(self, event, comment):
            wx.TipWindow(self.parent, comment, 500)

    def set_sizer(self):
        self.mainSizer.Add(self.flexSizer, 1, wx.EXPAND | wx.ALL, 15)
        self.parent.SetSizer(self.mainSizer)

    def get_ctrl_list(self):
        return self.ctrl_list

    def get_switch(self):
        return self.switch
        

class Control(wx.Control):

    def __init__(self, parent, name, type, datatype, values,
                 style=0, size=wx.DefaultSize, 
                 validator=wx.DefaultValidator, wkf_switch=False,
                 help="", pretty_name=None,
                 validation_req=True, combo_type=None, selections=None):
        
        self.name = name
        self.default_values = values
        self.type = type
        self.datatype = datatype
        self.wfk_switch = wkf_switch
        self.help = help
        self.id = None
        self.validation = validation_req
        self.pretty_name = pretty_name

        if type == control_types.CHOICE_BOX:
            self.ctrl = wx.Choice(parent, wx.ID_ANY,
                                  style=style, size=wx.DefaultSize, 
                                  validator=wx.DefaultValidator,
                                  choices=values)
            self.selection = values[0]
        
        elif type == control_types.TEXT_BOX:
            self.ctrl = wx.TextCtrl(parent, id=wx.ID_ANY, 
                                    value=values, style=style,
                                    validator=validator, 
                                    size=size)
            self.selection = self.ctrl.GetValue()
        
        elif type == control_types.COMBO_BOX:
            self.ctrl = FileSelectorCombo(parent, id=wx.ID_ANY, 
                                          size=size, style=style,
                                          validator=validator,
                                          value=values)
            
            self.text_ctrl = self.ctrl.GetTextCtrl()
            self.selection = self.text_ctrl.GetValue()
        
        elif type == control_types.INT_CTRL:
            self.ctrl = IntCtrl(parent, id=wx.ID_ANY,
                                size=size, validator=validator,
                                style=style, value=values,
                                max=10000)
            
            self.selection = self.ctrl.GetValue()
            
        elif type == control_types.FLOAT_CTRL:
            self.ctrl = NumCtrl(parent, id=wx.ID_ANY,
                                size=size, validator=validator,
                                style=style, value=values,
                                integerWidth=2, fractionWidth=3, 
                                allowNegative=False)
            
            self.selection = self.ctrl.GetValue()
            
        elif type == control_types.DIR_COMBO_BOX:
            self.ctrl = DirSelectorCombo(parent, 
                                         id=wx.ID_ANY, 
                                         size=size, 
                                         validator=validator,
                                         value=values)  
            self.text_ctrl = self.ctrl.GetTextCtrl()
            self.selection = self.text_ctrl.GetValue()
            
        elif type == control_types.CHECKLIST_BOX:
            self.ctrl = wx.CheckListBox(parent, id=wx.ID_ANY,
                                        size=size, validator=validator,
                                        style=style, choices=values)
            if datatype == data_types.LBOOL:
                self.selection = {}
                for val in values:
                    self.set_selection(val, values.index(val), True)
            else:
                self.selection = []
                
        elif type == control_types.LISTBOX_COMBO:
            self.ctrl = ListBoxCombo(parent, size=size, 
                                     validator=validator,style=style, 
                                     values=values, combo_type=combo_type)
            self.listbox_ctrl = self.ctrl.GetListBoxCtrl()
            self.id = self.listbox_ctrl.GetId()
            self.selection = []
            
            if (combo_type == 4) or (combo_type == 5):

                if values:

                    # values.keys() is a list of contrast options the user
                    # typed into the contrasts box. this only exists if the
                    # user entered contrasts, then saved them or went back to
                    # the first model builder window, then returned
                    for val in values:

                        # this re-checks the user's past contrast option
                        # SELECTIONS (which ones were checked) in the listbox
                        self.selection.append(val)
                            
           
            self.options = self.ctrl.get_listbox_options()

            
        elif type == control_types.TEXTBOX_COMBO:
            self.ctrl = TextBoxCombo(parent, id=wx.ID_ANY, 
                                     size=size, style=style,
                                     validator=validator,
                                     value=values)
            
            self.text_ctrl = self.ctrl.GetTextCtrl()
            self.selection = self.text_ctrl.GetValue()

         
        elif type == control_types.CHECKBOX_GRID:
            self.ctrl = CheckBoxGrid(parent, idx=wx.ID_ANY,
                                     selections=selections,
                                     values=values,
                                     size=wx.DefaultSize)
            
            self.default_values = selections
            self.selection = self.ctrl.GetGridSelection()

            add_string = "\n\nAvailable analyses: %s.\nDenote which " \
                         "analyses to run for each ROI path by listing the " \
                         "names above. For example, if you wish to run %s " \
                         "and %s, you would enter: '/path/to/ROI.nii.gz': " \
                         "%s, %s" % (selections, selections[0], \
                         selections[2], selections[0], selections[2])

            self.help = self.help + add_string


        elif type == control_types.GPA_CHECKBOX_GRID:
            self.ctrl = GPAModelCheckBoxGrid(parent, idx=wx.ID_ANY,
                                             values=values,
                                             size=wx.DefaultSize)
            
            self.selection = self.ctrl.GetGridSelection()
                
        self.set_id()

    def get_listbox_options(self):
        if self.get_type() == control_types.LISTBOX_COMBO:
            return self.options

    # this takes the list of available contrast names from modelDesign_window
    # and sends it to the 'Add Contrast' dialog box, so that it may do
    # validation immediately when the user enters contrast strings
    def set_available_contrasts(self, avail_cons):
        if self.get_type() == control_types.LISTBOX_COMBO:
            self.ctrl.set_available_contrasts(avail_cons)

    def set_design_matrix(self, design_matrix):
        if self.get_type() == control_types.LISTBOX_COMBO:
            self.ctrl.set_design_matrix(design_matrix)

    def set_id(self):
        if not self.id:
            self.id = self.ctrl.GetId()
    
    def get_id(self):
        return self.id
        
    def get_ctrl(self):
        return self.ctrl
    
    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def get_type(self):
        return self.type
    
    def get_datatype(self):
        return self.datatype
    
    def get_values(self):
        return self.default_values

    def get_validation(self):
        return self.validation
    
    def set_selection(self, value, index=0, remove=False,
                      convert_to_string=True):
        
        if isinstance(self.selection, list):
            if convert_to_string:
                value = str(value)

            # here, 'value' is the single element selected in the control by
            # the user, and 'self.selection' is the list containing these
            # selections
            if remove:
                self.selection.remove(value)
            else:
                self.selection.append(value)

        elif self.get_type() == control_types.CHECKBOX_GRID:
            self.ctrl.onReload_set_selections(value)

        elif self.get_type() == control_types.GPA_CHECKBOX_GRID:
            self.ctrl.onReload_set_selections(value)    

        elif isinstance(self.selection, dict):
            if remove:
                self.selection[value, index]= False
            else:
                self.selection[value, index]= True
            
        else:
            self.selection = value

        if self.get_type() == control_types.LISTBOX_COMBO:
            self.listbox_selections = self.selection

    def get_selection(self):
        return self.selection

    def set_new_selection(self, selection):
        self.selection = selection

    def get_switch(self):
        return self.wfk_switch
        
    def set_value(self, val):
        import ast

        if val in [None, "", "None", "none"]:
            val = self.get_values()
        else:
            if self.get_type() == control_types.LISTBOX_COMBO:
                listbox = self.ctrl.GetListBoxCtrl()
                for v in val:
                    if v:                       
                        listbox.Insert(v, 0)
                        listbox.Check(0)
                        self.set_selection(v)

            elif self.get_type() == control_types.CHECKLIST_BOX:
                # if the control is a checkbox, handle appropriately
                # this is for the derivative list in the group analysis GUI
                if "[" in val and "]" in val:
                    val = val.replace("[", "")
                    val = val.replace("]", "")
                    val = val.replace("'", "")
                    val = val.split(", ")

                try:
                    self.ctrl.SetCheckedStrings(val)
                except AssertionError:
                    pass

                strings = self.ctrl.GetCheckedStrings()
                sample_list = self.get_values()
                for s in strings:
                    self.set_selection(s, sample_list.index(s))

            elif self.get_type() == control_types.CHECKBOX_GRID:
                if "None" not in val:
                    self.ctrl.onReload_set_selections(val)

            elif self.get_type() == control_types.GPA_CHECKBOX_GRID:
                self.ctrl.set_checkbox_grid_values(val)

            else:
                if self.get_type() == control_types.CHOICE_BOX:
                    if isinstance(val, list):
                        val = val[0]
                    self.ctrl.SetStringSelection(val)
                elif self.get_type() in [control_types.INT_CTRL,
                                         control_types.FLOAT_CTRL]:
                    if str(val) != 'None':
                        try:
                            val = ast.literal_eval(val)
                        except:
                            pass
                    self.ctrl.SetValue(val)
                elif self.get_type() in [control_types.COMBO_BOX,
                                         control_types.DIR_COMBO_BOX,
                                         control_types.TEXTBOX_COMBO]:
                    self.text_ctrl.SetValue(val)
                else:
                    self.ctrl.SetValue(val)
                    
                self.set_selection(val) 
    
    def get_help(self):
        return self.help

    def get_pretty_name(self):
        try:
            page_title = self.get_ctrl().GetParent().page.title
        except:
            page_title = None
            
        control_name = self.name
        if type(self.pretty_name) == str:
            control_name = self.pretty_name

        return control_name \
            if not page_title \
            else "%s: %s" % (page_title, control_name)
