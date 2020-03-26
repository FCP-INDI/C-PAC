from pathlib import Path

CPAC_DIR = str(Path(__file__).parent.parent.parent)

import sys

sys.path.append(CPAC_DIR)

import click
import os
import pytest

def test_build_data_config(cli_runner):
    from CPAC.__main__ import utils as CPAC_main_utils
    data_dir = os.path.join(CPAC_DIR, 'dev', 'circleci_data')
    
    os.chdir(data_dir)
    test_yaml = os.path.join(data_dir, "data_settings.yml")
    _delete_test_yaml(test_yaml)
    result = cli_runner.invoke(
        CPAC_main_utils.commands['data_config'].commands[
            'new_settings_template'
        ]
    )
    
    assert result.exit_code == 0
    assert result.output.startswith(
        "\nGenerated a default data_settings YAML file for editing"
    )
    assert os.path.exists(test_yaml)
    _delete_test_yaml(test_yaml)
    
def _delete_test_yaml(test_yaml):
    if os.path.exists(test_yaml):
        os.remove(test_yaml)