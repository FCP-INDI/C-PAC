# test/__init__.py
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

# Init globals/constants for the test package
CPAC_DIR = os.path.abspath(os.path.dirname(CPAC.__file__)).replace('CPAC','')
TEST_DIR = os.path.join(CPAC_DIR, 'test')
RESOURCE_DIR = os.path.join(TEST_DIR, 'resources')
TEMPLATE_DIR = os.path.join(TEST_DIR, 'templates')
SETTINGS_DIR = os.path.join(RESOURCE_DIR, 'settings')
CONFIG_DIR = os.path.join(SETTINGS_DIR, 'configs')

# Create template data config and store in RESOURCE_DIR
DATA_CONFIG_TEMPLATE = os.path.join(TEMPLATE_DIR, 'data_config_template.yml')
DATA_CONFIG = os.path.join(CONFIG_DIR, 'data_config_test.yml')
tmp_f = open(DATA_CONFIG_TEMPLATE, 'r')
res_f = open(DATA_CONFIG, 'w')
for line in tmp_f:
    res_f.write(line.replace('RESOURCE_DIR', RESOURCE_DIR))
tmp_f.close()
res_f.close()

# File resources
PIPELINE_CONFIG = os.path.join(CONFIG_DIR, 'pipeline_config_test.yml')
SUBJECT_LIST = os.path.join(CONFIG_DIR, 'CPAC_subject_list_test.yml')
STRAT_FILE = os.path.join(CONFIG_DIR, 'strategies.obj')

# AWS resources
CREDS_DIR = os.path.join(SETTINGS_DIR, 'creds')
AWS_CREDS = os.path.join(CREDS_DIR, 'aws_creds.csv')
DB_CREDS = os.path.join(CREDS_DIR, 'db_creds.csv')
BUCKET_NAME = 'fcp-indi'