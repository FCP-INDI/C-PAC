"""
Tests for surface configuration
"""
import os
import pkg_resources as p
import pytest
import yaml

from CPAC.pipeline.cpac_pipeline import run_workflow
from CPAC.utils.configuration import Configuration

@pytest.mark.skip(reason='timing out for unrelated reasons')
@pytest.mark.timeout(60)
def test_duplicate_freesurfer(tmp_path):
    """The pipeline should build fast if freesurfer is not self-duplicating"""
    config = Configuration(yaml.safe_load('FROM: abcd-options'))
    with open(p.resource_filename(
        "CPAC",
        os.path.join(
            "resources",
            "configs",
            "data_config_S3-BIDS-ABIDE.yml"
        )
    ), 'r') as data_config:
        sub_dict = yaml.safe_load(data_config)[0]
    for directory in ['output', 'working', 'log', 'crash_log']:
        directory_key = ['pipeline_setup', f'{directory}_directory', 'path']
        config[directory_key] = os.path.join(
            tmp_path, config[directory_key].lstrip('/'))
    run_workflow(sub_dict, config, False, test_config=True)
