import wx
from .custom_control import FileSelectorCombo, DirSelectorCombo, ListBoxCombo, TextBoxCombo
from wx.lib import masked
from wx.lib.masked import NumCtrl
from wx.lib.intctrl import IntCtrl
import wx.lib.intctrl
import pkg_resources as p

class GenericClass(wx.ScrolledWindow):
    
    def __init__(self,parent,title=""):
        self.parent = parent
        self.title = title
        self.flexSizer = wx.FlexGridSizer(cols=2, hgap=15, vgap=15)   
        self.flexSizer.AddGrowableCol(1)
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.ctrl_list = []
        self.__add_static()
        self.__set_scroll()
        self.switch = None
    
    def __set_scroll(self):
        maxWidth = 1000
        maxHeight = 1000
        
        self.parent.SetVirtualSize((maxWidth, maxHeight))
        self.parent.SetScrollRate(20,20)
        
    def add_static(self):
        hbox= wx.BoxSizer(wx.HORIZONTAL)
        t = wx.StaticText(self.parent, -1, self.title)
        t.SetFont(wx.Font(18, wx.SWISS, wx.NORMAL, wx.BOLD))
        img_path = p.resource_filename('CPAC', 'GUI/resources/images/cpac_logo_2.jpg')
        img = wx.Image(img_path, wx.BITMAP_TYPE_JPEG).ConvertToBitmap()
        #img = wx.Image('images/cpac_logo2.jpg', wx.BITMAP_TYPE_JPEG).ConvertToBitmap()
        bmp = wx.StaticBitmap(self.parent, -1, img)
        
        hbox.Add(bmp)
        hbox.Add(t, 0, wx.TOP, 8)
        self.mainSizer.Add(hbox, proportion=0, flag=wx.ALL, border=5)
        self.mainSizer.Add(wx.StaticLine(self.parent), 0, wx.EXPAND|wx.TOP)
        
    __add_static = add_static
        
        
    def add(self, label, control, name, type = 0, 
            comment="", values="", style=0, size= wx.DefaultSize, 
            validator=wx.DefaultValidator, wkf_switch= False,
            validation_req = True, combo_type = None):
        
            
        label = wx.StaticText(self.parent, -1, label)
        hbox= wx.BoxSizer(wx.HORIZONTAL)
        img_path = p.resource_filename('CPAC', 'GUI/resources/images/help.png')
        image1 = wx.Image(img_path, wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        #image1 = wx.Image("images/help.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        button = wx.BitmapButton(self.parent, id=-1, bitmap=image1,
                                 pos=(10, 20), size = (image1.GetWidth()+5, image1.GetHeight()+5))
        #button = wx.Button(self.parent, id = wx.ID_HELP)
    
        button.Bind(wx.EVT_BUTTON, lambda event: \
                         self.OnShowDoc(event, comment))
        hbox.Add(label)
        hbox.Add(button)
        self.flexSizer.Add(hbox)
        
        ctrl = Control(self.parent, name=name, 
                      type=control, datatype=type, 
                      values=values, size=size, 
                      validator=validator,
                      wkf_switch = wkf_switch,
                      help = comment, 
                      validation_req = validation_req,
                      combo_type = combo_type)
        
        self.ctrl_list.append(ctrl)
        
        if wkf_switch:
            self.switch = ctrl
        
        if control == 0:
            self.parent.Bind(wx.EVT_CHOICE, lambda event: self.EvtChoice(event,ctrl), id =ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl(), proportion=0)    
        elif control == 1:
            self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterBox(event,ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl())
        elif control ==2 :
            self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterCombo(event,ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl(), flag = wx.EXPAND)
        elif control ==3:
            self.parent.Bind(wx.lib.intctrl.EVT_INT, lambda event: self.TxtEnterBox(event,ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl())
        elif control ==4:
            self.parent.Bind(masked.EVT_NUM, lambda event: self.TxtEnterBox(event,ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl())
        elif control ==5 :
            self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterCombo(event,ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl(), flag = wx.EXPAND)
        elif control ==6:
            self.parent.Bind(wx.EVT_CHECKLISTBOX,lambda event: self.EvtCheckListBox(event, ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl())
        elif control ==7:
            self.parent.Bind(wx.EVT_CHECKLISTBOX,lambda event: self.EvtListBoxCombo(event, ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl())
        elif control == 8 :
            self.parent.Bind(wx.EVT_TEXT, lambda event: self.TxtEnterCombo(event,ctrl), id = ctrl.get_id())
            self.flexSizer.Add(ctrl.get_ctrl(), flag = wx.EXPAND)


    def EvtChoice(self, event, ctrl):
        
        if type(event.GetString()) == unicode:
            value = event.GetString().encode('ascii', 'ignore')
        
        ctrl.set_selection(value)
            
        self.parent.Refresh()
        
    def TxtEnterBox(self, event, ctrl):
        ctrl.get_ctrl().SetBackgroundColour("white")
        #print "inside ctrl -->", ctrl
        #print "type event.GetString() -->", type(event.GetString())
        #print "ctrl.get_ctrl().GetValue() -->", ctrl.get_ctrl().GetValue()
        ctrl.set_selection(ctrl.get_ctrl().GetValue())
        
    def TxtEnterCombo(self, event, ctrl):
        ctrl.text_ctrl.SetBackgroundColour("white")
        #print "type ctrl.text_ctrl.GetValue() -->", type(ctrl.text_ctrl.GetValue())
        #print "ctrl.text_ctrl.GetValue() -->", ctrl.text_ctrl.GetValue()
        ctrl.set_selection(ctrl.text_ctrl.GetValue())
        
    
    def EvtCheckListBox(self, event, ctrl):
        index = event.GetSelection()
        label = ctrl.get_ctrl().GetString(index)
        if ctrl.get_ctrl().IsChecked(index):
            #print "label selected -->", label, index
            ctrl.set_selection(label, index)
        else:
            #print "label to be removed -->", label, index
            ctrl.set_selection(label,index, True)
    
    def EvtListBoxCombo(self, event, ctrl):
        
        index = event.GetSelection()
        label = ctrl.get_ctrl().GetListBoxCtrl().GetString(index)
        #print "EvtListBoxCombo label -->", label
        if ctrl.get_ctrl().GetListBoxCtrl().IsChecked(index):
            #print "label selected -->", label
            ctrl.set_selection(label, index)
        else:
            #print "label to be removed -->", label
            ctrl.set_selection(label,index, True)
        print ctrl.get_selection()
        
        
    
    def OnShowDoc(self, event, comment):
            wx.TipWindow(self.parent, comment, 500)


    def set_sizer(self):
        
        self.mainSizer.Add(self.flexSizer,1,wx.EXPAND|wx.ALL,15)
        self.parent.SetSizer(self.mainSizer)
        
    
    def get_ctrl_list(self):
        return self.ctrl_list

     
    def get_switch(self):    
        return self.switch

        

class Control(wx.Control):
    def __init__(self, parent, name, type, datatype, values,
                 style=0, size= wx.DefaultSize, 
                 validator=wx.DefaultValidator, 
                 wkf_switch=False, help="", validation_req = True, 
                 combo_type = None):
        
        self.name=name
        self.default_values = values
        self.type= type
        self.datatype = datatype
        self.wfk_switch = wkf_switch
        self.help = help
        self.id = None
        self.validation = validation_req
       
        if type ==0:
            self.ctrl= wx.Choice(parent, wx.ID_ANY,
                               style=style, size= wx.DefaultSize, 
                               validator=wx.DefaultValidator,
                               choices=values)
            self.selection = values[0]        
        
        elif type ==1:
            self.ctrl= wx.TextCtrl(parent, id = wx.ID_ANY, 
                               value = values, style= style,
                               validator = validator, 
                               size = size)
            self.selection = self.ctrl.GetValue()
        
        elif type ==2:
            self.ctrl= FileSelectorCombo(parent, id= wx.ID_ANY, 
                                         size=size, style=style,
                                         validator = validator,
                                         value = values)  
            
            self.text_ctrl = self.ctrl.GetTextCtrl()
            self.selection = self.text_ctrl.GetValue()
        
        elif type ==3:
            self.ctrl= IntCtrl(parent, id = wx.ID_ANY,
                               size = size, validator = validator,
                               style = style, value = values,
                               max = 10000)
            
            self.selection = self.ctrl.GetValue()
            
        elif type ==4:
            self.ctrl = NumCtrl(parent, id = wx.ID_ANY,
                                size = size, validator = validator,
                                style= style, value= values,
                                integerWidth=2, fractionWidth = 3, 
                                allowNegative=False)
            
            self.selection = self.ctrl.GetValue()
            
        elif type ==5:
            self.ctrl= DirSelectorCombo(parent, 
                                        id= wx.ID_ANY, 
                                        size=size, 
                                        validator = validator,
                                        value = values)  
            self.text_ctrl = self.ctrl.GetTextCtrl()
            self.selection = self.text_ctrl.GetValue()
            
        elif type ==6:
            self.ctrl = wx.CheckListBox(parent, id = wx.ID_ANY,
                                        size = size, validator = validator,
                                        style= style, choices = values)
            if datatype ==3:
                self.selection = {}
                for val in values:
                    self.set_selection(val, values.index(val), True)
                    
            else:
                self.selection = []
                
        elif type ==7:
            self.ctrl = ListBoxCombo(parent, size = size, 
                                     validator = validator,style= style, 
                                     values = values, combo_type= combo_type)
            self.listbox_ctrl = self.ctrl.GetListBoxCtrl()
            self.id = self.listbox_ctrl.GetId()
            self.selection = []
            
        elif type == 8:
            self.ctrl= TextBoxCombo(parent, id= wx.ID_ANY, 
                                         size=size, style=style,
                                         validator = validator,
                                         value = values)  
            
            self.text_ctrl = self.ctrl.GetTextCtrl()
            self.selection = self.text_ctrl.GetValue()
                
        self.set_id()
        
    def set_id(self):
        if self.id==None:
            self.id = self.ctrl.GetId()
    
    def get_id (self): 
        return self.id
        
    def get_ctrl(self):
        return self.ctrl
    
    def get_name(self):
        return self.name

    def get_type(self):
        return self.type
    
    def get_datatype(self):
        return self.datatype
    
    def get_values(self):
        return self.default_values

    def get_validation(self):
        return self.validation
    
    def set_selection(self, value, index = 0, remove=False, convert_to_string = True):
        
        if isinstance(self.selection, list):
            if convert_to_string:
                value = str(value)
                
            if remove:
                self.selection.remove(value)
            else:
                self.selection.append(value)
                
        elif isinstance(self.selection, dict):
            if remove:
                self.selection[value, index]= False
            else:
                self.selection[value, index]= True
            
        else:
            self.selection = value

    def get_selection(self):
        return self.selection

    def get_switch(self):
        return self.wfk_switch
        
    def set_value(self, val):
        import ast
        #print "self.get_name(), self.get_type() , val -->", self.get_name(), self.get_type(), val
        if val == None or val =="":
            val = self.get_values()
        else:
            if self.get_type()==7:
                listbox = self.ctrl.GetListBoxCtrl()
                for v in val:
                    if v:                       
                        listbox.Insert(v,0)
                        listbox.Check(0)
                        self.set_selection(v)         
            elif self.get_type()==6:
                #self.ctrl.SetChecked(val)
                self.ctrl.SetCheckedStrings(val)
                strings = self.ctrl.GetCheckedStrings() 
                sample_list = self.get_values()
                for s in strings:
                    self.set_selection(s, sample_list.index(s))
            else:
                if self.get_type() == 0:
                    if isinstance(val, list):
                        val = val[0]
                    #print "combo box value, default values---------------> ", val, self.get_values()
                    self.ctrl.SetStringSelection(val)
                elif self.get_type() ==3 or self.get_type()==4:
                    if str(val) != 'None':
                        val = ast.literal_eval(val)
                    self.ctrl.SetValue(val)
                elif self.get_type() ==2 or self.get_type()==5 or self.get_type()==8:
                    self.text_ctrl.SetValue(val)
                else:
                    self.ctrl.SetValue(val)
                    
                self.set_selection(val) 
    
    def get_help(self):
        return self.help
        
