import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p

class VMHC(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib.request import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        
        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/vmhc.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/vmhc.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/vmhc.html')
#            else:
#                self.LoadFile('html/vmhc.html')
#        except:
#            self.LoadFile('html/vmhc.html')
            
            
    def get_counter(self):
        return self.counter
            
class VMHCSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Voxel-mirrored Homotopic Connectivity (VMHC) Options")
        
        self.page.add(label="Calculate VMHC ", 
                 control=control.CHOICE_BOX, 
                 name='runVMHC', 
                 type=dtype.LSTR, 
                 comment="Calculate Voxel-mirrored Homotopic Connectivity (VMHC) for all voxels.", 
                 values=["Off","On"],
                 wkf_switch = True)
        
        self.page.add(label="Symmetric Template (Brain Only) ", 
         control=control.COMBO_BOX, 
         name='template_symmetric_brain_only', 
         type=dtype.STR, 
         values = "$FSLDIR/data/standard/MNI152_T1_${resolution_for_anat}_brain_symmetric.nii.gz",
         comment="Included as part of the 'Image Resource Files' package available on the Install page of the User Guide.\n\nIt is not necessary to change this path unless you intend to use a non-standard symmetric template.")
        
        self.page.add(label="Symmetric Template (With Skull) ", 
         control=control.COMBO_BOX, 
         name='template_symmetric_skull', 
         type=dtype.STR, 
         values = "$FSLDIR/data/standard/MNI152_T1_${resolution_for_anat}_symmetric.nii.gz",
         comment="Included as part of the 'Image Resource Files' package available on the Install page of the User Guide.\n\nIt is not necessary to change this path unless you intend to use a non-standard symmetric template.")

        self.page.add(label="Dilated Symmetric Brain Mask ", 
         control=control.COMBO_BOX, 
         name='dilated_symmetric_brain_mask', 
         type=dtype.STR, 
         values = "$FSLDIR/data/standard/MNI152_T1_${resolution_for_anat}_brain_mask_symmetric_dil.nii.gz",
         comment="Included as part of the 'Image Resource Files' package available on the Install page of the User Guide.\n\nIt is not necessary to change this path unless you intend to use a non-standard symmetric template.")
        
        self.page.add(label="FLIRT Configuration File ", 
         control=control.COMBO_BOX, 
         name='configFileTwomm', 
         type=dtype.STR, 
         values = "$FSLDIR/etc/flirtsch/T1_2_MNI152_2mm.cnf",
         comment="Included as part of the 'Image Resource Files' package available on the Install page of the User Guide.\n\nIt is not necessary to change this path unless you intend to use a non-standard symmetric template.")
        
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
