import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import os
import pkg_resources as p


class Motion(wx.html.HtmlWindow):

    def __init__(self, parent, counter=0):
        wx.html.HtmlWindow.__init__(
            self,
            parent,
            style=wx.html.HW_SCROLLBAR_AUTO
        )
        self.SetStandardFonts()

        self.counter = counter

        self.LoadFile(
            p.resource_filename('CPAC', 'GUI/resources/html/nuisance.html')
        )

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

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
            return self.counter
