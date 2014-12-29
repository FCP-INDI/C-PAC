# CPAC/tests/unit/GUI/interface/windows/unit_test.py
#
'''
This module performs testing on the functions in CPAC/GUI/interface/windows
'''

# Import packages
import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs': [np.get_include()]})

def test_config_window(sublist_path):
    '''
    This function tests the various functions and instance methods
    within the config_window.py module
    
    Parameters
    ----------
    sublist_path : string
        filepath to a C-PAC subject list
    
    Returns
    -------
    test_flg : boolean
        a flag indicating whether all of the tests performed on
        config_window passed or not
    '''
    
    # Import packages
    from CPAC.GUI.interface.windows import config_window
    import os
    import wx
    import yaml
    
    # Init variables
    test_flg = False
    wx_app = wx.App()
    config_main_frame = config_window.MainFrame(None)
    sublist = yaml.load(open(os.path.abspath(sublist_path),'r'))
    
    # Test instance methods
    # test_sublist
    sublist_flg = config_main_frame.test_sublist(sublist)
    
    if sublist_flg:
        test_flg = True
    
    # Return the pass status of the test
    return test_flg
