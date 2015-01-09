# CPAC/tests/__init__.py
#
# Contributing authors (please append): Daniel Clark

'''
This module performs testing on the functions in CPAC
'''
# Import packages
import CPAC
import os
import unittest

# Init globals
cpac_base = CPAC.__file__
RESOURCE_DIR = '/'.join(cpac_base.split('/')[:-1]) + '/tests/resources'

# Import submodule TestCases
from CPAC.tests.unit.GUI.interface.windows import ConfigWindowTestCase,\
                                                  DataConfigWindowTestCase

# Import all TestCases to environment
__all__ = ['ConfigWindowTestCase',
           'DataConfigWindowTestCase']