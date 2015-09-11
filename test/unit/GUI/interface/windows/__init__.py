# test/unit/GUI/interface/windows/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/GUI/interface/windows package
'''

# Import module TestCases
from test.unit.GUI.interface.windows.config_window_test \
    import ConfigWindowTestCase
from test.unit.GUI.interface.windows.dataconfig_window_test \
    import DataConfigWindowTestCase


# Setup module environment
__all__ = ['ConfigWindowTestCase',
           'DataConfigWindowTestCase']
