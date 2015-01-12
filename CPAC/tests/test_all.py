# CPAC/tests/test_all.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This script runs all of the TestCases of all of the submodules within
the CPAC/tests folder and returns the results
'''

# Import packages
import unittest

# Import submodules
import CPAC.tests.unit.GUI.interface.windows as windows

# Init full test suite
full_suite = unittest.TestSuite()

# Add TestCases, by module, to the full TestSuite
full_suite.addTest(unittest.TestLoader().loadTestsFromModule(windows))

# Init the full suite runner and run the TestSuite
suite_runner = unittest.TextTestRunner()
suite_runner.run(full_suite)

