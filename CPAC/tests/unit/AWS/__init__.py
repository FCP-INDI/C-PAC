# CPAC/tests/unit/AWS/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/AWS package
'''

# Import module TestCases
from CPAC.tests.unit.AWS.test_fetch_creds import FetchCredsTestCase

# Setup module environment
__all__ = ['FetchCredsTestCase']