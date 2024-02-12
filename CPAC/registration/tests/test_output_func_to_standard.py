import os
import nibabel as nb
import numpy as np
import pytest

from CPAC.registration import output_func_to_standard
from .mocks import configuration_strategy_mock
import nipype.interfaces.afni.utils as afni_utils
import nipype.interfaces.utility as util
import CPAC.utils.test_init as test_utils


@pytest.mark.skip(reason='needs refactoring')
def test_output_func_to_standard_ANTS():

    test_name = 'test_output_func_to_standard_ANTS'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock(method='ANTS')
    num_strat = 0

    # build the workflow
    workflow = pe.Workflow(name='test_output_func_to_standard_ANTS')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    output_func_to_standard(workflow,
            'mean_functional',
            'template_brain_for_func_preproc',
            'mean_functional_to_standard',
            strat, num_strat, c, input_image_type='func_derivative')

    out1_name = os.path.join(c.workingDirectory, test_name, 
            'apply_ants_warp_mean_functional_to_standard_0',
            'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_antswarp.nii.gz')

    in_node, in_file = strat['mean_functional']
    output_func_to_standard(workflow,
            (in_node, in_file),
            'template_brain_for_func_preproc',
            'mean_functional_to_standard_node',
            strat, num_strat, c, input_image_type='func_derivative')

    out2_name = os.path.join(c.workingDirectory, test_name, 
            'apply_ants_warp_mean_functional_to_standard_node_0',
            'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_antswarp.nii.gz')

    workflow.run()

    assert(test_utils.pearson_correlation(out1_name, out2_name) > 0.99)


@pytest.mark.skip(reason='needs refactoring')
def test_output_func_to_standard_FSL_linear():

    test_name = 'test_output_func_to_standard_FSL_linear'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock(method='FSL')
    strat.append_name('anat_mni_flirt_register_0')
    num_strat = 0

    # build the workflow
    workflow = pe.Workflow(name='test_output_func_to_standard_FSL_linear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    output_func_to_standard(workflow,
            'mean_functional',
            'template_brain_for_func_preproc',
            'mean_functional_to_standard',
            strat, num_strat, c, input_image_type='func_derivative')

    out1_name = os.path.join(c.workingDirectory, test_name, 
            'func_mni_fsl_warp_mean_functional_to_standard_0',
            'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_warp.nii.gz')

    func_node, func_out = strat['mean_functional']
    output_func_to_standard(workflow,
            (func_node, func_out),
            'template_brain_for_func_preproc',
            'mean_functional_to_standard_node',
            strat, num_strat, c, input_image_type='func_derivative')

    out2_name = os.path.join(c.workingDirectory, test_name, 
            'func_mni_fsl_warp_mean_functional_to_standard_node_0',
            'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_warp.nii.gz')

    workflow.run()

    assert(test_utils.pearson_correlation(out1_name, out2_name) > 0.99)


@pytest.mark.skip(reason='needs refactoring')
def test_output_func_to_standard_FSL_nonlinear():

    test_name = 'test_output_func_to_standard_FSL_nonlinear'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock(method='FSL')
    strat.append_name('anat_mni_fnirt_register_0')
    num_strat = 0

    # build the workflow
    workflow = pe.Workflow(name='test_output_func_to_standard_FSL_nonlinear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    output_func_to_standard(workflow,
            'mean_functional',
            'template_brain_for_func_preproc',
            'mean_functional_to_standard',
            strat, num_strat, c, input_image_type='func_derivative')

    out1_name = os.path.join(c.workingDirectory, test_name, 
            'func_mni_fsl_warp_mean_functional_to_standard_0',
            'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_warp.nii.gz')

    node, out_file = strat['mean_functional']
    output_func_to_standard(workflow,
            (node, out_file),
            'template_brain_for_func_preproc',
            'mean_functional_to_standard_node',
            strat, num_strat, c, input_image_type='func_derivative')

    out2_name = os.path.join(c.workingDirectory, test_name, 
            'func_mni_fsl_warp_mean_functional_to_standard_node_0',
            'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_warp.nii.gz')

    workflow.run()

    assert(test_utils.pearson_correlation(out1_name, out2_name) > .99)
