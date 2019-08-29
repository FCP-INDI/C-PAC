import os
import wx
import sys
import pkg_resources as p
from xvfbwrapper import Xvfb

class Screenshot(object):

    window = None

    def __init__(self, output, size=(1200, 600), wait=1000):
        self.output = output
        self.size = size
        self.wait = wait

    def __enter__(self):
        self.vdisplay = Xvfb(width=self.size[0], height=self.size[1])
        self.vdisplay.start()
        self.app = wx.App(False)
        wx.UIActionSimulator().MouseMove(0, 0)
        return self

    def set_window(self, window):
        self.window = window

    def __exit__(self, type, value, traceback):
        def _screenshot():
            if not self.window: 
                return 
            rect = self.window.GetRect()
            if sys.platform == 'linux2':
                client_x, client_y = self.window.ClientToScreen((0, 0))
                border_width = client_x - rect.x
                title_bar_height = client_y - rect.y
                rect.width += (border_width * 2)
                rect.height += title_bar_height + border_width
            dcScreen = wx.ScreenDC()
            bmp = wx.EmptyBitmap(rect.width, rect.height)
            memDC = wx.MemoryDC()
            memDC.SelectObject(bmp)
            memDC.Blit(0, 0, rect.width, rect.height, dcScreen, rect.x, rect.y)
            memDC.SelectObject(wx.NullBitmap)
            img = bmp.ConvertToImage()

            try:
                os.makedirs(os.path.dirname(self.output))
            except:
                pass

            img.SaveFile(self.output, wx.BITMAP_TYPE_PNG)

        wx.FutureCall(self.wait, _screenshot)
        if self.window is not None:
            wx.FutureCall(self.wait + 1000, self.window.Destroy)

        self.app.MainLoop()
        self.vdisplay.stop()


def screenshot_page(page, output, size=(1200, 600), wait=1000, ind=True):

    from CPAC.GUI.interface.windows import config_window

    with Screenshot(output, size) as ss:

        config_mainframe = config_window.MainFrame(
            None, "load",
            path=p.resource_filename('CPAC', 'resources/configs/pipeline_config_template.yml'),
            size=(size), ind=ind
        )
        config_mainframe.Show(True)
        config_mainframe.OpenPage(page)

        ss.set_window(config_mainframe)
