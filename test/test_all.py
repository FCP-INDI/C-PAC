# test/test_all.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This script runs all of the TestCases of all of the submodules within
the C-PAC/test package and returns the results
'''

# Import packages
import unittest

# Import submodules
import unit.AWS as aws
import unit.GUI.interface.windows as windows
import unit.network_centrality as centrality
import unit.pipeline as pipeline

# Init full test suite
full_suite = unittest.TestSuite()

# Add TestCases, by module, to the full TestSuite
full_suite.addTest(unittest.TestLoader().loadTestsFromModule(aws))
full_suite.addTest(unittest.TestLoader().loadTestsFromModule(windows))
full_suite.addTest(unittest.TestLoader().loadTestsFromModule(centrality))
full_suite.addTest(unittest.TestLoader().loadTestsFromModule(pipeline))

# Init the full suite runner and run the TestSuite
suite_runner = unittest.TextTestRunner()
suite_runner.run(full_suite)

