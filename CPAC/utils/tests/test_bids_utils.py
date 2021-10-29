import os
import pytest
from CPAC.utils.bids_utils import create_cpac_data_config


def create_sample_bids_structure(root_dir):
    """Function to create temporary synthetic BIDS data for testing parsing"""
    def _prefix_entities(paths, path):
        return f'sub-{paths[path]["sub"]}_ses-{paths[path]["ses"]}'
    
    suffixes = {
        'anat': ['T1w.nii.gz', 'acq-VNavNorm_T1w.nii.gz'],
        'fmap': [],
        'func': [
            f'task-{task}_run-{run}_bold.nii.gz' for
           run in ['1', '2'] for task in ['rest']
        ]
    }

    paths = {
        os.path.join(root_dir, subpath): {
          'sub': sub,
          'ses': ses,
          'data_type': data_type
        } for data_type in ['anat', 'fmap', 'func'] for
        ses in ['1', '2'] for
        sub in ['0001', '0002'] for subpath in [
            f'sub-{sub}/ses-{ses}/{data_type}'
        ]
    }

    for path in paths:
        os.makedirs(path, exist_ok=True)
        for suffix in suffixes[paths[path]['data_type']]:
            open(os.path.join(path, '_'.join([
                _prefix_entities(paths, path),
                suffix
            ])), 'w')


@pytest.mark.parametrize('only_one_anat', [True, False])
def test_create_cpac_data_config_only_one_anat(tmp_path, only_one_anat):
    """Function to test 'only_one_anat' parameter of
    'create_cpac_data_config' function"""
    create_sample_bids_structure(tmp_path)
    assert isinstance(
        create_cpac_data_config(
            str(tmp_path), only_one_anat=only_one_anat
        )[0]['anat']['T1w'], str if only_one_anat else list)
