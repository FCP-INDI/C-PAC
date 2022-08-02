'''Tests for C-PAC YAML

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
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.'''
import os
import tempfile
import yaml

from CPAC.utils.configuration import DEFAULT_PIPELINE_FILE
from CPAC.utils.yaml_template import create_yaml_from_template


def test_yaml_template():
    '''Test YAML pipeline template'''
    # Create temporary file
    config_file = tempfile.mkstemp(suffix='test_yaml_template')[1]

    # Create a new YAML configuration file based on the default pipeline
    # YAML file.
    with open(DEFAULT_PIPELINE_FILE, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    new_config = create_yaml_from_template(config, DEFAULT_PIPELINE_FILE, True)
    with open(config_file, 'w', encoding='utf-8') as f:
        f.write(new_config)

    # Verify that the output has preserved blank lines, comments
    with open(config_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Delete temporary file
    os.remove(config_file)

    # Assert first line declares YAML version
    assert lines[0].startswith('%YAML')
    assert lines[1] == '---\n'

    # Assert first lines starts with a comment
    assert lines[2].startswith('# ')

    # Assert that regressors configuration written in block style
    assert (line for line in lines if line.strip() == 'Regressors:')
    assert (line for line in lines if line.strip() == '- Bandpass:')

    # Assert that short lists of values are written in flow style
    # (e.g., "GM_label: [3, 42]")
    assert (line for line in lines if
            '[' in line and ',' in line and ']' in line and
            not line.lstrip().startswith('#'))
