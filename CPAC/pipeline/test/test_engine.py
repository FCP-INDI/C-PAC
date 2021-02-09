
import os
import yaml
from CPAC.utils.configuration import Configuration
from CPAC.pipeline.cpac_pipeline import initialize_nipype_wf
from CPAC.pipeline.engine import ResourcePool, ingress_raw_data

from CPAC.utils.bids_utils import create_cpac_data_config


def test_ingress_raw_data(pipe_config, bids_dir, test_dir):

    sub_data_dct = create_cpac_data_config(bids_dir,
                                           skip_bids_validator=True)[0]

    # Load in pipeline config file
    config_file = os.path.realpath(pipe_config)
    try:
        if not os.path.exists(config_file):
            raise IOError
        else:
            cfg = Configuration(yaml.safe_load(open(config_file, 'r')))
    except IOError:
        print("config file %s doesn't exist" % config_file)
        raise
    except yaml.parser.ParserError as e:
        error_detail = "\"%s\" at line %d" % (
            e.problem,
            e.problem_mark.line
        )
        raise Exception(
            "Error parsing config file: {0}\n\n"
            "Error details:\n"
            "    {1}"
            "\n\n".format(config_file, error_detail)
        )
    except Exception as e:
        raise Exception(
            "Error parsing config file: {0}\n\n"
            "Error details:\n"
            "    {1}"
            "\n\n".format(config_file, e)
        )
    cfg.pipeline_setup['output_directory']['path'] = \
        os.path.join(test_dir, 'out')
    cfg.pipeline_setup['working_directory']['path'] = \
        os.path.join(test_dir, 'work')

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    part_id = sub_data_dct['subject_id']
    ses_id = sub_data_dct['unique_id']

    unique_id = f'{part_id}_{ses_id}'

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    wf, rpool, diff, blip, fmap_rp_list = ingress_raw_data(wf, rpool, cfg,
                                                           sub_data_dct,
                                                           unique_id,
                                                           part_id, ses_id)
    print(rpool.get_entire_rpool())
    print(diff)
    print('---')
    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()


cfg = "/code/default_pipeline.yml"
bids_dir = "/Users/steven.giavasis/data/HBN-SI_dataset/rawdata"
test_dir = "/test_dir"

test_ingress_raw_data(cfg, bids_dir, test_dir)