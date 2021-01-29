import os
import pytest

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.qc.pipeline import create_qc_workflow
from CPAC.qc.utils import generate_qc_pages
from CPAC.utils import Configuration, Strategy, Outputs

from CPAC.pipeline.cpac_group_runner import gather_outputs


def file_node(path):
    input_node = pe.Node(
        util.IdentityInterface(fields=['file']), name='inputspec'
    )
    input_node.inputs.file = path
    return input_node, 'file'


@pytest.mark.skip(reason='needs refactoring')
def test_qc():

    outputs = Outputs()

    c = Configuration({
        "workingDirectory": "",
        "crashLogDirectory": "",
        "outputDirectory":""
    })

    workflow = pe.Workflow(name='workflow_name')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    strat_initial = Strategy()
    strat_list = [strat_initial]

    output_df_dct = gather_outputs(
        c.outputDirectory,
        [
            "functional_brain_mask",
            "functional_to_anat_linear_xfm",
            "anatomical_brain",
            "anatomical_reorient",
            "mean_functional_in_anat",
            "motion_params",
            "frame_wise_displacement_power",
            "frame_wise_displacement_jenkinson",
        ],
        None,
        get_motion=False,
        get_raw_score=False,
        get_func=True,
        derivatives=[
            "functional_brain_mask",
            "functional_to_anat_linear_xfm",
            "anatomical_brain",
            "anatomical_reorient",
            "mean_functional_in_anat",
            "motion_params",
            "frame_wise_displacement_power",
            "frame_wise_displacement_jenkinson",
        ],
        exts=['nii', 'nii.gz', '1D', 'mat', 'txt']
    )

    for (resource, _), df in output_df_dct.items():
        strat_initial.update_resource_pool({
            resource: file_node(df.Filepath[0])
        })

    qc_montage_id_a, qc_montage_id_s, qc_hist_id, qc_plot_id = \
        create_qc_workflow(workflow, c, strat_list, outputs.qc)
