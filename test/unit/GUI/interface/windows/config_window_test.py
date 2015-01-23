# test/unit/GUI/interface/windows/test_config_window.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in CPAC/GUI/interface/windows
'''

# Import packages
import unittest
import wx
from test import RESOURCE_DIR


# Test case for the config_window.py module
class ConfigWindowTestCase(unittest.TestCase):
    '''
    This class is a test case for the config_window.py module
    
    Inherits
    --------
    unittest.TestCase class
    
    Attributes (class):
    -------------------
    see unittest.TestCase documentation
    
    Attributes (instance):
    ----------------------
    config_mainframe : CPAC.GUI.interface.windows.config_window.MainFrame 
                       instance
    '''
    
    # Method to setup the MainFrame object
    def setUp(self):
        '''
        Method to instantiate the MainFrame object as an instance attribute
        '''
        
        # Import packages
        from CPAC.GUI.interface.windows import config_window
        
        # Init variables
        wx_app = wx.App()
        self.config_mainframe = config_window.MainFrame(None)
    
    # Method to test the subject list
    def test_config_sublist(self):
        '''
        This method tests the various functions and instance methods
        within the config_window.py module
        '''
        
        # Import packages
        import os
        import yaml
        
        # Init variables
        sublist_path = os.path.join(RESOURCE_DIR, 'CPAC_subject_list.yml')
        sublist = yaml.load(open(sublist_path, 'r'))
        
        # Test subject list
        pass_flg = self.config_mainframe.test_sublist(sublist)
        
        # Assert that pass flag is true
        self.assertEqual(pass_flg, True, 'Subject failed')


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()