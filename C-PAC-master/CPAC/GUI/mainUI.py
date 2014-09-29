#!/usr/bin/env python

import wx
from interface.windows.main_window import ListBox


def run():
    app = wx.App()
    ListBox(None, -1, 'Configure & Run CPAC')
    app.MainLoop()
    

    
if __name__ == "__main__":
    run()

