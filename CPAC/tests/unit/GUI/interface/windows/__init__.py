# CPAC/tests/unit/GUI/interface/windows/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/GUI/interface/windows module
'''

# Import module TestCases
from test_module import ConfigWindowTestCase,\
                        DataConfigWindowTestCase

# Setup module environment
__all__ = ['ConfigWindowTestCase',
           'DataConfigWindowTestCase']
