# CPAC/tests/__init__.py
#
# Contributing authors (please append):
# Daniel Clark

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
# File resources
PIPELINE_CONFIG = os.path.join(RESOURCE_DIR, 'pipeline_config_test.yml')
SUBJECT_LIST = os.path.join(RESOURCE_DIR, 'CPAC_sublist_test.yml')
STRAT_FILE = os.path.join(RESOURCE_DIR, 'strategies.obj')
# AWS resources
AWS_CREDS = os.path.join(RESOURCE_DIR, 'aws_creds.csv')
DB_CREDS = os.path.join(RESOURCE_DIR, 'db_creds.csv')
BUCKET_NAME = 'fcp-indi'