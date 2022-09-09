"""choose_participant_subset

Choose 1 sub per dataset or all subjects, only on main image

Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>."""
import random
import sys
from os.path import dirname
import json
import yaml

TEST_DATA = f'{dirname(dirname(__file__))}/test_data.yml'
with open(TEST_DATA, 'r', encoding='utf-8') as tds:
    TEST_DATA = yaml.safe_load(tds)


def all_participants():
    '''Return all subjects per dataset per species'''
    return {f'{species}_participants': [str(subject) for dataset in
            TEST_DATA['labels']['participant'][species] for subject in
            TEST_DATA['labels']['participant'][species][dataset]
            ] for species in TEST_DATA['labels']['participant']}


def random_participants():
    '''Return one subject per dataset per species'''
    return {f'{species}_participants': [str(random.choice(
        TEST_DATA['labels']['participant'][species][dataset]
    )) for dataset in TEST_DATA['labels']['participant'][species]
                                       ] for species in
        TEST_DATA['labels']['participant']}


if __name__ == '__main__':
    '''Return JSON to pass to `workflow_dispatch` inputs'''
    args = sys.argv[1:]
    if args:
        if args[0] == 'random':
            print("'" + json.dumps({'variant': [''], **random_participants()}
                                   ) + "'")
        elif args[0] == 'all':
            print("'" + json.dumps({'variant': [''], **all_participants()}) +
                  "'")
