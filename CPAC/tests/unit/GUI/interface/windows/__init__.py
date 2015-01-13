# CPAC/tests/unit/GUI/interface/windows/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/GUI/interface/windows package
'''

# Import module TestCases
from CPAC.tests.unit.GUI.interface.windows.test_config_window \
    import ConfigWindowTestCase
from CPAC.tests.unit.GUI.interface.windows.test_dataconfig_window \
    import DataConfigWindowTestCase


# Setup module environment
__all__ = ['ConfigWindowTestCase',
           'DataConfigWindowTestCase']
