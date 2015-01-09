# CPAC/tests/unit/GUI/interface/windows/test_module.py
#
# Contributing authors (please append): Daniel Clark

'''
This module performs testing on the functions in CPAC/GUI/interface/windows
'''

# Import packages
import unittest
import wx
from CPAC.tests import RESOURCE_DIR


# Test case for the config_window.py module
class ConfigWindowTestCase(unittest.TestCase):
    '''
    This class is a test case for the config_window.py module
    
    Inherits
    --------
    unittest.TestCase class
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
        self.config_mainframe_instance = config_window.MainFrame(None)
    
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
        pass_flg = self.config_mainframe_instance.test_sublist(sublist)
        
        # Assert that pass flag is true
        self.assertEqual(pass_flg, True, 'Subject failed')


# Test case for the dataconfig_window.py module
class DataConfigWindowTestCase(unittest.TestCase):
    '''
    This class is a test case for the dataconfig_window.py module
    
    Inherits
    --------
    unittest.TestCase class
    '''
    
    # Method to setup the DataConfig object
    def setUp(self):
        '''
        Method to instantiate the DataConfig object as an instance attribute
        '''
        
        # Import packages
        from CPAC.GUI.interface.windows import dataconfig_window
        
        # Init variables
        wx_app = wx.App()
        self.data_config_instance = dataconfig_window.DataConfig(None)
    
    # Method to tear down the TestCase
    def tearDown(self):
        '''
        Method to close the DataConfig instance after testing
        '''
        
        # Use the instance's cancel function with the button event to close it
        self.data_config_instance.cancel(wx.EVT_BUTTON)
    
    # Method to test if the window launches
    def test_dataconfig_init(self):
        '''
        Method to test if the data config window launches
        '''
        
        # Init variables
        err_msg = 'Data config window never launched or is not of expected '\
                  'instance type'
        
        # Test to make sure it comes up
        self.assertIsInstance(self.data_config_instance, wx.Frame, msg=err_msg)


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()