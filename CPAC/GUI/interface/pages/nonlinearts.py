import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class NonLinearTimeSeriesAnalysis(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        #self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/nuisance.html'))
            
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/nuisance.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/nuisance.html')
#            else:
#                self.LoadFile('html/nuisance.html')
#        except:
#            self.LoadFile('html/nuisance.html')
            
            
    def get_counter(self):
        return self.counter
            
class InformationTheory(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        
        self.counter = counter

        
        self.page = GenericClass(self, "Information Theory Calculation Options")
        
        self.page.add(label="Run Information Theory Measures", 
                 control=control.CHOICE_BOX, 
                 name='runIT', 
                 type=dtype.LSTR, 
                 comment="Run Information Theory Measures", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
                 
        self.page.add(label="Voxelwise / ROI extraction", 
                 control=control.CHOICE_BOX, 
                 name='voxel_roi', 
                 type=dtype.LSTR, 
                 comment="Run Information Theory Measures voxelwise or after ROI timeseries extraction", 
                 values=["Voxelwise","ROI","Voxelwise/ROI"],
                 wkf_switch = True)         
        
        self.page.add(label="fMRI image", 
                     control=control.COMBO_BOX, 
                     name='input_image', 
                     type=dtype.STR, 
                     comment="fMRI image for calculation")
       
        self.page.add(label="Parcellation Mask", 
                     control=control.COMBO_BOX, 
                     name='input_mask', 
                     type=dtype.STR, 
                     comment="Parcellation Mask if you want to calculate")


        self.page.add(label = "Measures:",
                      #control = control.CHECKLISTBOX_COMBO,
                      control = control.LISTBOX_COMBO,
                      name = "Measures",
                      type = dtype.LDICT,
                      values = ['Entropy', 'Conditional Entropy','Mutual Information','Transfer Entropy','Entropy Correlation Coefficient'],
                      comment = "Select which IT measures to apply:\n"\
                                "ent = Entropy\n"\
                                 "condent = Conditional Entropy\n"\
                                 "mi = Mutual Information\n"\
                                 "te = Transfer Entropy\n"\
                                 "ecc = Entropy Correlation Coefficient\n",
                     size = (300,120),
                     combo_type =1)
                    
        self.page.add(label= "CompCor Components ",
                      control = control.TEXT_BOX,
                      name = "nComponents",
                      type = dtype.LNUM,
                      values = "5",
                      validator = CharValidator("no-alpha"),
                      comment = "Number of Principle Components to calculate when running CompCor. We recommend 5 or 6.")

        self.page.add(label="Output Options ",
                      control=control.CHECKLIST_BOX,
                      name="measure_options",
                      type=dtype.LBOOL,
                      values=['CSV', 'NUMPY','NIFTI'],
                      comment="By default, results are written as NIFTI files. Additional output formats are as a .csv spreadsheet or a Numpy array.")


        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        
class Causality(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Causality")
        
        self.page.add(label="Run Causality ", 
                 control=control.CHOICE_BOX, 
                 name='runCausality', 
                 type=dtype.LSTR, 
                 comment="Granger Causality", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
