import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import os
import pkg_resources as p

class AnatomicalPreprocessing(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/anat.html'))
        
#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/anat.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/anat.html')
#            else:
#                self.LoadFile('html/anat.html')
#        except:
#            self.LoadFile('html/anat.html')
            
            
    def get_counter(self):
        return self.counter
            

class Segmentation(wx.ScrolledWindow):

    def __init__(self, parent, counter =0):
        wx.ScrolledWindow.__init__(self, parent)
        import os
        
        self.counter = counter
        
        fsl = os.environ.get('FSLDIR')
        if not fsl:
            fsl = "$FSLDIR"
                
        self.page = GenericClass(self, "Automatic Tissue Segmentation ")
        
        self.page.add(label="Run Tissue Segmentation ", 
                 control=control.CHOICE_BOX, 
                 name='runSegmentationPreprocessing', 
                 type=dtype.LSTR, 
                 comment="Automatically segment anatomical images into white matter, gray matter, and CSF based on prior probability maps.", 
                 values=["On","Off","On/Off"],
                 wkf_switch = True)

        self.page.add(label= "White Matter Probability Threshold ",
                 control=control.TEXT_BOX, 
                 name='whiteMatterThreshold', 
                 type=dtype.LNUM, 
                 values= "0.96",
                 validator = CharValidator("no-alpha"),
                 comment="Only voxels with a White Matter probability greater than this value will be classified as White Matter.\n\nCan be a single value or a list of values separated by commas.")
        
        self.page.add(label = "Gray Matter Probability Threshold ",
                 control =control.TEXT_BOX,
                 name = 'grayMatterThreshold',
                 type =dtype.LNUM,
                 values= "0.7",
                 validator = CharValidator("no-alpha"),
                 comment= "Only voxels with a Gray Matter probability greater than this value will be classified as Gray Matter.\n\nCan be a single value or a list of values separated by commas.")

        self.page.add(label= "CSF Probability Threshold ",
                 control=control.TEXT_BOX, 
                 name='cerebralSpinalFluidThreshold', 
                 type=dtype.LNUM, 
                 values = "0.96",
                 validator = CharValidator("no-alpha"),
                 comment="Only voxels with a CSF probability greater than this value will be classified as CSF.\n\nCan be a single value or a list of values separated by commas.")
        
        self.page.add(label= "Priors Directory ",
                 control=control.DIR_COMBO_BOX, 
                 name='priors_path', 
                 type=dtype.STR, 
                 values= os.path.join(fsl, 'data/standard/tissuepriors/2mm'),
                 comment="Full path to a directory containing binarized prior probability maps.\n\nThese maps are included as part of the 'Image Resource Files' package available on the Install page of the User Guide.\n\nIt is not necessary to change this path unless you intend to use non-standard priors.")

        self.page.add(label= "White Matter Prior Probability Map ",
                 control=control.COMBO_BOX, 
                 name='PRIORS_WHITE', 
                 type=dtype.STR, 
                 values = '$priors_path/avg152T1_white_bin.nii.gz',
                 comment="Full path to a binarized White Matter prior probability map.\n\nIt is not necessary to change this path unless you intend to use non-standard priors.")
        
        self.page.add(label= "Gray Matter Prior Probability Map ",
                 control=control.COMBO_BOX, 
                 name='PRIORS_GRAY', 
                 type=dtype.STR, 
                 values = '$priors_path/avg152T1_gray_bin.nii.gz',
                 comment="Full path to a binarized Gray Matter prior probability map.\n\nIt is not necessary to change this path unless you intend to use non-standard priors.")
        
        self.page.add(label= "CSF Prior Probability Map ",
                 control=control.COMBO_BOX, 
                 name='PRIORS_CSF', 
                 type=dtype.STR, 
                 values = '$priors_path/avg152T1_csf_bin.nii.gz',
                 comment="Full path to a binarized CSF prior probability map.\n\nIt is not necessary to change this path unless you intend to use non-standard priors.")        
                
        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter
    
class Registration(wx.ScrolledWindow):
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
        
        self.counter = counter
                
        self.page = GenericClass(self, "Anatomical Registration")
        
        fsl = os.environ.get('FSLDIR')
        if not fsl:
            fsl = "$FSLDIR"
        
        self.page.add(label="Run Anatomical Registration ", 
                     control=control.CHOICE_BOX, 
                     name='runRegistrationPreprocessing', 
                     type=dtype.LSTR, 
                     comment="Register anatomical images to a template.", 
                     values=["On","Off","On/Off"],
                     wkf_switch = True)
        
        self.page.add(label="Anatomical Template Resolution ", 
                      control=control.CHOICE_BOX, 
                      name='resolution_for_anat', 
                      type=dtype.STR, 
                      values = ["1mm", "2mm", "3mm"],
                      comment="The resolution to which anatomical images should be transformed during registration.\n\nThis is the resolution at which processed anatomical files will be output.")
        
        self.page.add(label="Anatomical Template (Brain Only) ", 
                     control=control.COMBO_BOX, 
                     name='template_brain_only_for_anat', 
                     type=dtype.STR, 
                     values = str(os.path.join(fsl, "data/standard/MNI152_T1_${resolution_for_anat}_brain.nii.gz")),
                     comment="Template to be used during registration.\n\nIt is not necessary to change this path unless you intend to use a non-standard template.")

        self.page.add(label="Anatomical Template (With Skull) ", 
                     control=control.COMBO_BOX, 
                     name='template_skull_for_anat', 
                     type=dtype.STR, 
                     values =  str(os.path.join(fsl, "data/standard/MNI152_T1_${resolution_for_anat}.nii.gz")),
                     comment="Template to be used during registration.\n\nIt is not necessary to change this path unless you intend to use a non-standard template.")

        self.page.add(label="Anatomical to Template Registration Method ", 
                     control=control.CHOICE_BOX, 
                     name='regOption', 
                     type=dtype.LSTR, 
                     comment="Use either ANTS or FSL (FLIRT and FNIRT) as your anatomical registration method.", 
                     values=["ANTS","FSL","ANTS & FSL"],
                     wkf_switch = True)

        self.page.add(label="FSL FNIRT Configuration File (FSL only) ", 
                     control=control.COMBO_BOX, 
                     name='fnirtConfig', 
                     type=dtype.STR, 
                     values =  str(os.path.join("T1_2_MNI152_2mm")),
                     comment="Configuration file to be used by FSL to set FNIRT parameters.\n\nIt is not necessary to change this path unless you intend to use custom FNIRT parameters or a non-standard template.")

        self.page.add(label="FSL FNIRT Reference Mask (FSL only) ", 
                     control=control.COMBO_BOX, 
                     name='ref_mask', 
                     type=dtype.STR, 
                     values =  str(os.path.join(fsl, "data/standard/MNI152_T1_${resolution_for_anat}_brain_mask_dil.nii.gz")),
                     comment="Configuration file to be used by FSL to set FNIRT parameters.\n\nIt is not necessary to change this path unless you intend to use custom FNIRT parameters or a non-standard template.")

        self.page.add(label="Use skull-on image to calculate transform? (ANTS only) ", 
                     control=control.CHOICE_BOX, 
                     name='regWithSkull', 
                     type=dtype.LSTR, 
                     comment="Register skull-on anatomical image to a template.", 
                     values=["Off","On"],
                     wkf_switch = True)

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
        return self.counter
                
