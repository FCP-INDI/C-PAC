"""
Tests for surface configuration
"""
import os
import pkg_resources as p
import pytest
import yaml

from CPAC.pipeline.cpac_pipeline import run_workflow
from CPAC.utils.configuration import Configuration


@pytest.mark.timeout(60)
def test_duplicate_freesurfer():
    """The pipeline should build fast if freesurfer is not self-duplicating"""
    c = Configuration(yaml.safe_load('FROM: abcd-options'))
    sub_dict = yaml.safe_load(open(p.resource_filename(
        "CPAC",
        os.path.join(
            "resources",
            "configs",
            "data_config_S3-BIDS-ABIDE.yml"
        )
    ), 'r'))[0]

    run_workflow(sub_dict, c, False, test_config=True)
