import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class Nuisance(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(self, parent, style= wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()
        
        self.counter = counter
        self.LoadPage(p.resource_filename('CPAC', 'GUI/resources/html/nuisance.html'))         
            
    def get_counter(self):
        return self.counter
            


class NuisanceRegression(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        import os
        
        self.counter = counter

        fsl = os.environ.get('FSLDIR')
        if not fsl:
            fsl = "$FSLDIR"
        
        self.page = GenericClass(self, "Nuisance Signal Regression Options")
        
        self.page.add(label="Run Nuisance Signal Regression ", 
                 control=control.CHOICE_BOX, 
                 name='runNuisance', 
                 type=dtype.LSTR, 
                 comment="Run Nuisance Signal Regression", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        
        self.page.add(label="Lateral Ventricles Mask (Standard Space) ", 
                     control=control.COMBO_BOX, 
                     name='lateral_ventricles_mask', 
                     type=dtype.STR, 
                     values = os.path.join(fsl, "data/atlases/HarvardOxford/HarvardOxford-lateral-ventricles-thr25-2mm.nii.gz"),
                     comment="Standard Lateral Ventricles Binary Mask")

        self.page.add(label="Anaticor ",
                      control=control.CHOICE_BOX,
                      name='include_anaticor',
                      type=dtype.LSTR,
                      comment="Include Anaticor as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label= "      Anaticor Radius (in mm) ",
                      control = control.TEXT_BOX,
                      name = "anaticor_radius",
                      type = dtype.LNUM,
                      values = "0",
                      validator = CharValidator("no-alpha"),
                      comment = "Radius (in mm) for Anaticor.")

        self.page.add(label="aCompCor ",
                      control=control.CHOICE_BOX,
                      name='include_aCompCor',
                      type=dtype.LSTR,
                      comment="Include aCompCor as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Components ",
                      control=control.TEXT_BOX,
                      name="aCompCor_num_pcs",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Number of components to retain in aCompCor.")

        self.page.add(label="      Tissues ",
                      control=control.CHOICE_BOX,
                      name = "aCompCor_tissues",
                      type=dtype.LSTR,
                      values = ['White Matter', 'CSF', 'White Matter + CSF'],
                      comment="")

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='aCompCor_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])
        
        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='aCompCor_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])
        
        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='aCompCor_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="tCompCor ",
                      control=control.CHOICE_BOX,
                      name='include_tCompCor',
                      type=dtype.LSTR,
                      comment="Include tCompCor as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Number of Components ",
                      control=control.TEXT_BOX,
                      name="tCompCor_num_pcs",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Number of components to retain in tCompCor.")

        self.page.add(label="      Threshold ",
                      control=control.TEXT_BOX,
                      name="tCompCor_threshold",
                      type=dtype.LNUM,
                      values="0.00",
                      validator=CharValidator("no-alpha"),
                      comment="Cutoff as raw variance value.\n\nA floating "
                              "point number followed by SD (ex. 1.5SD) = "
                              "mean + a multiple of the SD.\n\nA floating "
                              "point number followed by PCT (ex. 2PCT) = "
                              "percentile from the top (ex is top 2%).")

        self.page.add(label="      By Slice ",
                      control=control.CHOICE_BOX,
                      name='tCompCor_byslice',
                      type=dtype.LSTR,
                      comment="Whether or not the threshold criterion should "
                              "be applied by slice or across the entire "
                              "volume, makes most sense for SD or PCT.",
                      values=["False", "True"])

        self.page.add(label="      Tissues ",
                      control=control.CHOICE_BOX,
                      name="tCompCor_tissues",
                      type=dtype.LSTR,
                      values=['White Matter', 'CSF', 'White Matter + CSF'],
                      comment="")

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='tCompCor_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='tCompCor_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='tCompCor_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="White Matter ",
                      control=control.CHOICE_BOX,
                      name='include_white_matter',
                      type=dtype.LSTR,
                      comment="Include White Matter as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Summary Method ",
                      control=control.CHOICE_BOX,
                      name="white_matter_summary",
                      type=dtype.LSTR,
                      values=['PCA', 'Mean', 'NormMean', 'DetrendNormMean'],
                      comment="")
        
        self.page.add(label="      Number of Components ",
                      control=control.TEXT_BOX,
                      name="white_matter_num_pcs",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Number of components to retain in White Matter.")

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='white_matter_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='white_matter_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='white_matter_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="Ventricles ",
                      control=control.CHOICE_BOX,
                      name='include_ventricles',
                      type=dtype.LSTR,
                      comment="Include White Matter as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Summary Method ",
                      control=control.CHOICE_BOX,
                      name="ventricles_summary",
                      type=dtype.LSTR,
                      values=['PCA', 'Mean', 'NormMean', 'DetrendNormMean'],
                      comment="")

        self.page.add(label="      Number of Components ",
                      control=control.TEXT_BOX,
                      name="ventricles_num_pcs",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Number of components to retain in White Matter.")

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='ventricles_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='ventricles_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='ventricles_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="Grey Matter ",
                      control=control.CHOICE_BOX,
                      name='include_grey_matter',
                      type=dtype.LSTR,
                      comment="Include White Matter as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Summary Method ",
                      control=control.CHOICE_BOX,
                      name="grey_matter_summary",
                      type=dtype.LSTR,
                      values=['PCA', 'Mean', 'NormMean', 'DetrendNormMean'],
                      comment="")

        self.page.add(label="      Number of Components ",
                      control=control.TEXT_BOX,
                      name="grey_matter_num_pcs",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Number of components to retain in White Matter.")

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='grey_matter_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='grey_matter_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='grey_matter_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="Global Signal ",
                      control=control.CHOICE_BOX,
                      name='include_global_signal',
                      type=dtype.LSTR,
                      comment="Include White Matter as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Summary Method ",
                      control=control.CHOICE_BOX,
                      name="global_signal_summary",
                      type=dtype.LSTR,
                      values=['PCA', 'Mean', 'NormMean', 'DetrendNormMean'],
                      comment="")

        self.page.add(label="      Number of Components ",
                      control=control.TEXT_BOX,
                      name="global_signal_num_pcs",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Number of components to retain in White Matter.")

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='global_signal_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='global_signal_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='global_signal_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="Motion ",
                      control=control.CHOICE_BOX,
                      name='include_motion',
                      type=dtype.LSTR,
                      comment="Include Motion as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Include Delayed ",
                      control=control.CHOICE_BOX,
                      name='motion_delayed',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Squared ",
                      control=control.CHOICE_BOX,
                      name='motion_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Include Delayed Squared ",
                      control=control.CHOICE_BOX,
                      name='motion_delayed_squared',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="Censor ",
                      control=control.CHOICE_BOX,
                      name='include_censor',
                      type=dtype.LSTR,
                      comment="Include Censoring as a nuisance regressor.",
                      values=["Off", "On"])

        self.page.add(label="      Threshold Metric ",
                      control=control.CHOICE_BOX,
                      name="censor_thresh_metric",
                      type=dtype.LSTR,
                      values=['RMSD', 'DVARS', 'RMSD+DVARS'],
                      comment="")

        self.page.add(label="      Threshold ",
                      control=control.TEXT_BOX,
                      name="censor_threshold",
                      type=dtype.LNUM,
                      values="0.00",
                      validator=CharValidator("no-alpha"),
                      comment="Threshold to be applied to the metric.")

        self.page.add(label="      Number of Previous TRs to Remove ",
                      control=control.CHOICE_BOX,
                      name='censor_num_prev_trs',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Number of Subsequent TRs to Remove ",
                      control=control.CHOICE_BOX,
                      name='censor_num_subseq_trs',
                      type=dtype.LSTR,
                      comment="???????????",
                      values=["False", "True"])

        self.page.add(label="      Method ",
                      control=control.CHOICE_BOX,
                      name="censor_method",
                      type=dtype.LSTR,
                      values=['Kill', 'Zero', 'Interpolate',
                              'Spike Regression'],
                      comment="")

        self.page.add(label="PolyOrt ",
                      control=control.CHOICE_BOX,
                      name='include_PolyOrt',
                      type=dtype.LSTR,
                      comment="",
                      values=["Off", "On"])

        self.page.add(label="      Polynomial Degree ",
                      control=control.TEXT_BOX,
                      name="PolyOrtDegree",
                      type=dtype.LNUM,
                      values="0",
                      validator=CharValidator("no-alpha"),
                      comment="Polynomial degree up to which will be removed."
                              "\n\nex. 2 means constant + linear + quadratic,"
                              " probably the most that will be needed, "
                              "especially if band pass filtering is used.")

        self.page.add(label="Bandpass Filtering ",
                      control=control.CHOICE_BOX,
                      name='includeBandPass',
                      type=dtype.LSTR,
                      comment="",
                      values=["Off", "On"])

        self.page.add(label="      Bottom Frequency ",
                      control=control.TEXT_BOX,
                      name="bandPassBottomFrequency",
                      type=dtype.LNUM,
                      values="0.00",
                      validator=CharValidator("no-alpha"),
                      comment="Frequency in Hertz of the highpass part of "
                              "the passband, frequencies below this will be "
                              "removed.")

        self.page.add(label="      Top Frequency ",
                      control=control.TEXT_BOX,
                      name="bandPassTopFrequency",
                      type=dtype.LNUM,
                      values="0.00",
                      validator=CharValidator("no-alpha"),
                      comment="Frequency in Hertz of the lowpass part of "
                              "the passband, frequencies above this will be "
                              "removed.")

        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter
        


class MedianAngleCorrection(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Median Angle Correction Options")
        
        self.page.add(label="Run Median Angle Correction ", 
                 control=control.CHOICE_BOX, 
                 name='runMedianAngleCorrection', 
                 type=dtype.LSTR, 
                 comment="Correct for the global signal using Median Angle Correction.", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        
        self.page.add(label= "Target Angle (degrees) ",
                      control = control.TEXT_BOX,
                      name = "targetAngleDeg",
                      type = dtype.LNUM,
                      values = "90",
                      validator = CharValidator("no-alpha"),
                      comment = "Target angle used during Median Angle Correction.")
        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter



class FilteringSettings(wx.ScrolledWindow):
    
    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)
                
        self.counter = counter
        
        self.page = GenericClass(self, "Temporal Filtering Options")
        
        self.page.add(label="Run Temporal Filtering ", 
                 control=control.CHOICE_BOX, 
                 name='runFrequencyFiltering', 
                 type=dtype.LSTR, 
                 comment="Apply a temporal band-pass filter to functional data.", 
                 values=["Off","On","On/Off"],
                 wkf_switch = True)
        
        self.page.add(label = "Band-Pass Filters ",
                      control = control.LISTBOX_COMBO,
                      name = "nuisanceBandpassFreq",
                      type = dtype.LOFL,
                      values = [0.01, 0.1],
                      comment = "Define one or more band-pass filters by clicking the + button.",
                     size = (200,100),
                     combo_type = 2)

        
        self.page.set_sizer()
        parent.get_page_list().append(self)
        
    def get_counter(self):
            return self.counter



class Scrubbing(wx.ScrolledWindow):

    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter

        self.page = GenericClass(self, "Scrubbing Options")

        self.page.add(label="Run Scrubbing ",
                      control=control.CHOICE_BOX,
                      name='runScrubbing',
                      type=dtype.LSTR,
                      comment="Remove volumes exhibiting excessive motion.",
                      values=["Off", "On", "On/Off"],
                      wkf_switch=True)

        self.page.add(label="Framewise Displacement (FD) Threshold (mm) ",
                      control=control.TEXT_BOX,
                      name='scrubbingThreshold',
                      type=dtype.LNUM,
                      values="0.2",
                      validator=CharValidator("no-alpha"),
                      comment="Specify the maximum acceptable Framewise Displacement (FD) in millimeters.\n\nAny volume exhibiting FD greater than this value will be removed.",
                      size=(100, -1))

        self.page.add(label="Number of Preceeding Volumes to Remove ",
                      control=control.INT_CTRL,
                      name='numRemovePrecedingFrames',
                      type=dtype.NUM,
                      comment="Number of volumes to remove preceeding a volume with excessive FD.",
                      values=1)

        self.page.add(label="Number of Subsequent Volumes to Remove ",
                      control=control.INT_CTRL,
                      name='numRemoveSubsequentFrames',
                      type=dtype.NUM,
                      comment="Number of volumes to remove subsequent to a volume with excessive FD.",
                      values=2)

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
        return self.counter
