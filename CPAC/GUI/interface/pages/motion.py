import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import os
import pkg_resources as p


class Motion(wx.html.HtmlWindow):

    def __init__(self, parent, counter=0):
        from urllib2 import urlopen
        wx.html.HtmlWindow.__init__(
            self, parent, style=wx.html.HW_SCROLLBAR_AUTO)
        self.SetStandardFonts()

        self.counter = counter

        self.LoadFile(p.resource_filename('CPAC', 'GUI/resources/html/motion.html'))

#        try:
#            code = urlopen("http://fcp-indi.github.io/docs/user/compute_config.html").code
#            if (code / 100 < 4):
#                self.LoadPage('http://fcp-indi.github.io/docs/user/compute_config.html')
#            else:
#                self.LoadFile('html/settings.html')
#        except:
#            self.LoadFile('html/settings.html')

    def get_counter(self):
        return self.counter


class MotionOptions(wx.ScrolledWindow):

    def __init__(self, parent, counter=0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter

        self.page = GenericClass(self, "Motion Correction Options")

        self.page.add(label="Use Friston 24-Parameter Model ",
                      control=control.CHOICE_BOX,
                      name='runFristonModel',
                      type=dtype.LSTR,
                      comment="Use the Friston 24-Parameter Model during volume realignment.\n\nIf this option is turned off, only 6 parameters will be used.\n\nThese parameters will also be output as a spreadsheet.",
                      values=["On", "Off","On/Off"])

        self.page.add(label="Calculate Motion Statistics ",
                      control=control.CHOICE_BOX,
                      name='runGenerateMotionStatistics',
                      type=dtype.LSTR,
                      comment="Calculate motion statistics including Framewise Displacement (FD) and DVARS.\n\nRequired to run Scrubbing.\n\nThese parameters will also be output as a spreadsheet.",
                      values=["On", "Off","On/Off"])

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
