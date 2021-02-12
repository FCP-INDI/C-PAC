
import os

from CPAC.pipeline.cpac_runner import run_T1w_longitudinal
from CPAC.pipeline.cpac_pipeline import load_cpac_pipe_config
from CPAC.utils.bids_utils import create_cpac_data_config


def test_run_T1w_longitudinal(bids_dir, cfg, test_dir, part_id):

    sub_data_list = create_cpac_data_config(bids_dir,
                                            participant_label=part_id,
                                            skip_bids_validator=True)
    cfg = load_cpac_pipe_config(cfg)

    cfg.pipeline_setup['output_directory']['path'] = \
        os.path.join(test_dir, 'out')
    cfg.pipeline_setup['working_directory']['path'] = \
        os.path.join(test_dir, 'work')

    run_T1w_longitudinal(sub_data_list, cfg)


cfg = "/code/default_pipeline.yml"
bids_dir = "/Users/steven.giavasis/data/neurodata_hnu"
test_dir = "/test_dir"
part_id = "0025427"

test_run_T1w_longitudinal(bids_dir, cfg, test_dir, part_id)