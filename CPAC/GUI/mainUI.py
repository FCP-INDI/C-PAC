#!/usr/bin/env python

try:
    import wx
    from .interface.windows.main_window import ListBox
except ImportError as e:
    print('CPAC was imported and can be used, but wxPython was not '\
          'found. In order to use the GUI, wxPython must be installed.')
    print('Error: %s' % e)


def run():
    app = wx.App()
    ListBox(None, -1, 'Configure & Run CPAC')
    app.MainLoop()
    

    
if __name__ == "__main__":
    run()

