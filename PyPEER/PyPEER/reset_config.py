#!/usr/bin/env python
"""
Script used on the command line to reset configuration file

Authors:
    - Jake Son, 2017-2018  (jake.son@childmind.org)  http://jakeson.me

"""

from .peer_func import *
from pprint import pprint


def update_config():

    project_dir, data_dir, stimulus_path = scaffolding()
    os.chdir(project_dir)

    configs = load_config()

    print('\n\nThis is your current config file for reference:')
    print('====================================================\n')
    pprint(configs)
    print('\n')

    print('Update the config file:')
    print('====================================================\n')

    updated_configs = set_parameters(configs, new=True)

    print('\n\nThis is your new config file:')
    print('====================================================\n')
    pprint(updated_configs)
    print('\n')


if __name__ == "__main__":

    update_config()
