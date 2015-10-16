# test/unit/network_centrality/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/network_centrality package
'''

# Import module TestCases
from test.unit.network_centrality.resting_state_centrality_test import \
    CentralityWorkflowTestCase

# Setup module environment
__all__ = ['CentralityWorkflowTestCase']