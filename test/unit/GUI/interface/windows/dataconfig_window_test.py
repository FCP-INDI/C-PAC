# test/unit/GUI/interface/windows/test_dataconfig_window.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in CPAC/GUI/interface/windows
'''

# Import packages
import unittest
import wx
from test import RESOURCE_DIR, DATA_CONFIG

# Test case for the dataconfig_window.py module
class DataConfigWindowTestCase(unittest.TestCase):
    '''
    This class is a test case for the dataconfig_window.py module
    
    Inherits
    --------
    unittest.TestCase class
    
    Attributes (class):
    -------------------
    see unittest.TestCase documentation
    
    Attributes (instance):
    ----------------------
    data_config : CPAC.GUI.interface.windows.dataconfig_window.DataConfig
                  instance
    '''

    # Method to setup the DataConfig object
    def setUp(self):
        '''
        Method to instantiate the DataConfig object as an instance attribute
        '''

        # Import packages
        from CPAC.GUI.interface.windows import dataconfig_window

        # Init variables
        self.wx_app = wx.App()
        self.data_config = dataconfig_window.DataConfig(None)

    # Method to tear down the TestCase
    def tearDown(self):
        '''
        Method to close the DataConfig instance after testing
        '''

        # Use the instance's cancel function with the button event to close it
        self.data_config.cancel(wx.EVT_BUTTON)

    # Method to test if the window launches
    def test_dataconfig_init(self):
        '''
        Method to test if the data config window launches
        '''

        # Init variables
        err_msg = 'Data config window never launched or is not of expected '\
                  'instance type'

        # Test to make sure it comes up
        self.assertIsInstance(self.data_config, wx.Frame, msg=err_msg)

    # Method to test the run method
    def test_dataconfig_run(self):
        self.data_config.run(DATA_CONFIG)

# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()