import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p


class AfterWarping(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        wx.html.HtmlWindow.__init__(
            self,
            parent,
            style=wx.html.HW_SCROLLBAR_AUTO
        )
        self.SetStandardFonts()

        self.counter = counter
        self.LoadFile(
            p.resource_filename('CPAC', 'GUI/resources/html/after_warp.html')
        )          

    def get_counter(self):
        return self.counter


class AfterWarpingOptions(wx.ScrolledWindow):

    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter

        self.page = GenericClass(self, "After Warping Options")

        self.page.add(label="Run Smoothing ",
                      control=control.CHOICE_BOX,
                      name='run_smoothing',
                      type=dtype.LSTR,
                      comment="Smooth the derivative outputs.\n\nOn - Run "
                              "smoothing and output only the smoothed "
                              "outputs.\nOn/Off - Run smoothing and output "
                              "both the smoothed and non-smoothed outputs.\n"
                              "Off - Don't run smoothing.",
                      values=["On", "On/Off", "Off"])

        self.page.add(label="Smoothing Kernel FWHM (in mm) ",
                      control=control.TEXT_BOX,
                      name='fwhm',
                      type=dtype.LNUM,
                      values= "4",
                      validator = CharValidator("no-alpha"),
                      comment="Full Width at Half Maximum of the Gaussian "
                              "kernel used during spatial smoothing.\n\nCan "
                              "be a single value or multiple values "
                              "separated by commas.\n\nNote that spatial "
                              "smoothing is run as the last step in the "
                              "individual-level analysis pipeline, such "
                              "that all derivatives are output both "
                              "smoothed and unsmoothed.")

        self.page.add(label="Smooth Before/After z-Scoring ",
                      control=control.CHOICE_BOX,
                      name='smoothing_order',
                      type=dtype.LSTR,
                      comment="Choose whether to smooth outputs before or "
                              "after z-scoring.",
                      values=["Before", "After"])

        self.page.add(label="z-score Standardize Derivatives ",
                      control=control.CHOICE_BOX,
                      name='runZScoring',
                      type=dtype.LSTR,
                      comment="z-score standardize the derivatives. This is "
                              "required for group-level analysis.\n\n"
                              "On - Run z-scoring and output only the "
                              "z-scored outputs.\nOn/Off - Run z-scoring and "
                              "output both the z-scored and raw score "
                              "versions of the outputs.\nOff - Don't run "
                              "z-scoring.",
                      values=["On", "On/Off", "Off"])

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
            return self.counter
