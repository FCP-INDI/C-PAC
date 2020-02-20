import wx
import wx.html
from ..utils.generic_class import GenericClass
from ..utils.constants import control, dtype
from ..utils.validator import CharValidator
import pkg_resources as p

class ReHo(wx.html.HtmlWindow):

    def __init__(self, parent, counter  = 0):
        wx.html.HtmlWindow.__init__(
            self,
            parent,
            style=wx.html.HW_SCROLLBAR_AUTO
        )
        self.SetStandardFonts()

        self.counter = counter
        self.LoadFile(
            p.resource_filename('CPAC', 'GUI/resources/html/reho.html')
        )

    def get_counter(self):
        return self.counter

class ReHoSettings(wx.ScrolledWindow):

    def __init__(self, parent, counter = 0):
        wx.ScrolledWindow.__init__(self, parent)

        self.counter = counter

        self.page = GenericClass(self, "Regional Homogeneity (ReHo) Options")

        self.page.add(label="Calculate Regional Homogeneity (ReHo) ",
                 control=control.CHOICE_BOX,
                 name='runReHo',
                 type=dtype.LSTR,
                 comment="Calculate Regional Homogeneity (ReHo) for all voxels.",
                 values=["Off","On"],
                 wkf_switch = True)


        self.page.add(label="Voxel Cluster Size ",
                     control=control.CHOICE_BOX,
                     name='clusterSize',
                     type=dtype.NUM,
                     comment="Number of neighboring voxels used when calculating ReHo\n"
                             "7 (Faces)\n"
                             "19 (Faces + Edges)\n"
                             "27 (Faces + Edges + Corners)",
                     values=["7","19", "27"])

        self.page.set_sizer()
        parent.get_page_list().append(self)

    def get_counter(self):
            return self.counter
