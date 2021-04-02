import os
import pytest
from CPAC.pipeline.cpac_pipeline import initialize_nipype_wf, \
    load_cpac_pipe_config
from CPAC.pipeline.engine import ResourcePool, ingress_raw_anat_data, \
    ingress_raw_func_data, ingress_pipeconfig_paths  # noqa F401

from CPAC.utils.bids_utils import create_cpac_data_config


@pytest.mark.skip(reason='needs refactoring (ingress_raw_data split into '
                  'anat and func)')
def test_ingress_raw_data(pipe_config, bids_dir, test_dir):

    sub_data_dct = create_cpac_data_config(bids_dir,
                                           skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup['output_directory']['path'] = \
        os.path.join(test_dir, 'out')
    cfg.pipeline_setup['working_directory']['path'] = \
        os.path.join(test_dir, 'work')

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    part_id = sub_data_dct['subject_id']
    ses_id = sub_data_dct['unique_id']

    unique_id = f'{part_id}_{ses_id}'

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    wf, rpool, diff, blip, fmap_rp_list = ingress_raw_data(wf, rpool, cfg,   # noqa F821
                                                           sub_data_dct,
                                                           unique_id,
                                                           part_id, ses_id)

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()


@pytest.mark.skip(reason='not a pytest test')
def test_ingress_pipeconfig_data(pipe_config, bids_dir, test_dir):

    sub_data_dct = create_cpac_data_config(bids_dir,
                                           skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup['output_directory']['path'] = \
        os.path.join(test_dir, 'out')
    cfg.pipeline_setup['working_directory']['path'] = \
        os.path.join(test_dir, 'work')

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    part_id = sub_data_dct['subject_id']
    ses_id = sub_data_dct['unique_id']

    unique_id = f'{part_id}_{ses_id}'

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    rpool = ingress_pipeconfig_paths(cfg, rpool, sub_data_dct, unique_id)

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()


cfg = "/code/default_pipeline.yml"
bids_dir = "/Users/steven.giavasis/data/HBN-SI_dataset/rawdata"
test_dir = "/test_dir"


if __name__ == '__main__':
    # test_ingress_raw_data(cfg, bids_dir, test_dir)
    test_ingress_pipeconfig_data(cfg, bids_dir, test_dir)
