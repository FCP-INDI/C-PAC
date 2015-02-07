# test/unit/AWS/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/AWS package
'''

# Import module TestCases
from test.unit.AWS.fetch_creds_test import FetchCredsTestCase


# Setup module environment
__all__ = ['FetchCredsTestCase']