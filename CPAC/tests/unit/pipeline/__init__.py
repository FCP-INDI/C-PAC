# CPAC/tests/unit/pipeline/__init__.py
#
# Contributing authors (please append):
# Daniel Clark
#
'''
This module performs testing via unittest.TestCases of the functionality of the
CPAC/pipeline package
'''

# Import module TestCases
from CPAC.tests.unit.pipeline.test_cpac_pipeline import CPACPipelineRunTestCase

# Setup module environment
__all__ = ['CPACPipelineRunTestCase']