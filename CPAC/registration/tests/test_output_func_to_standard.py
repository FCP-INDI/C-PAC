import os
import nipype.pipeline.engine as pe
from CPAC.registration import output_func_to_standard
from mocks import configuration_strategy_mock

def test_output_func_to_standard_ANTS():

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

    workflow.run()

def test_output_func_to_standard_FSL_linear():

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

def test_output_func_to_standard_FSL_nonlinear():

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
