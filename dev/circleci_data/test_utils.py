from pathlib import Path

CPAC_DIR = str(Path(__file__).parent.parent.parent)

import sys

sys.path.append(CPAC_DIR)

import os

from CPAC.__main__ import utils as CPAC_main_utils

DATA_DIR = os.path.join(CPAC_DIR, 'dev', 'circleci_data')

def test_build_data_config(cli_runner):
    os.chdir(DATA_DIR)
    test_yaml = os.path.join(DATA_DIR, "data_settings.yml")
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
    
    
def test_new_settings_template(cli_runner):
    os.chdir(CPAC_DIR)
    
    result = cli_runner.invoke(
        CPAC_main_utils.commands['data_config'].commands['build'],
        [os.path.join(
            DATA_DIR,
            "data_settings_bids_examples_ds051_default_BIDS.yml"
        )]
    )
    
    participant_yaml = os.path.join(DATA_DIR, "data_config_ds051.yml")
    group_yaml = os.path.join(DATA_DIR, "group_analysis_participants_ds051.txt")
    
    assert result.exit_code == 0
    assert result.output.startswith("\nGenerating data configuration file..")
    assert os.path.exists(participant_yaml)
    assert os.path.exists(group_yaml)
    _delete_test_yaml(participant_yaml)
    _delete_test_yaml(group_yaml)
    
def _delete_test_yaml(test_yaml):
    if os.path.exists(test_yaml):
        os.remove(test_yaml)