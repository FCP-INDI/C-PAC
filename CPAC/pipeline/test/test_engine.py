
import os
from CPAC.pipeline.cpac_pipeline import initialize_nipype_wf, \
    load_cpac_pipe_config
from CPAC.pipeline.engine import ResourcePool, ingress_raw_anat_data, ingress_raw_func_data, \
    ingress_pipeconfig_paths

from CPAC.utils.bids_utils import create_cpac_data_config


def test_ingress_func_raw_data(pipe_config, bids_dir, test_dir):

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

    wf, rpool, diff, blip, fmap_rp_list = ingress_raw_func_data(wf, rpool, cfg,
                                                                sub_data_dct,
                                                                unique_id,
                                                                part_id, ses_id)

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()

def test_ingress_anat_raw_data(pipe_config, bids_dir, test_dir):

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

    rpool = ingress_raw_anat_data(wf, rpool, cfg,
                                     sub_data_dct,
                                     unique_id,
                                     part_id, ses_id)

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()

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


# cfg = "/code/default_pipeline.yml"
# bids_dir = "/Users/steven.giavasis/data/HBN-SI_dataset/rawdata"
# test_dir = "/test_dir"

cfg = "/Users/hecheng.jin/GitHub/DevBranch/dev/docker_data/default_pipeline.yml"
bids_dir = '/Users/hecheng.jin/Monkey/monkey_data_oxford/site-ucdavis'
test_dir = "/Users/hecheng.jin/GitHub/Test/T2ingress"

# test_ingress_func_raw_data(cfg, bids_dir, test_dir)
test_ingress_anat_raw_data(cfg, bids_dir, test_dir)
# test_ingress_pipeconfig_data(cfg, bids_dir, test_dir)