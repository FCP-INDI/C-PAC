import os
import pytest
from .mocks import configuration_strategy_mock
from CPAC.pipeline import nipype_pipeline_engine as pe
from ..output_func_to_standard import fsl_apply_transform_func_to_mni


@pytest.mark.skip(reason='needs refactoring')
def test_fsl_apply_transform_func_to_mni_nonlinear():

    c, strat = configuration_strategy_mock(method='FSL')

    strat.append_name('anat_mni_fnirt_register_0')

    # build the workflow
    workflow = pe.Workflow(name='test_fsl_apply_transform_func_to_mni_nonlinear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = fsl_apply_transform_func_to_mni(workflow, 'mean_functional_to_standard',
            'mean_functional', 'template_brain_for_func_preproc', 0, strat,
            c.funcRegFSLinterpolation)

    workflow.run()


@pytest.mark.skip(reason='needs refactoring')
def test_fsl_apply_transform_func_to_mni_nonlinear_mapnode():

    c, strat = configuration_strategy_mock(method='FSL')

    strat.append_name('anat_mni_fnirt_register_0')

    # build the workflow
    workflow = pe.Workflow(name='test_fsl_apply_transform_func_to_mni_nonlinear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = fsl_apply_transform_func_to_mni(workflow, 'dr_tempreg_to_standard',
            'dr_tempreg_maps_files', 'template_brain_for_func_preproc', 0, strat,
            c.funcRegFSLinterpolation, map_node=True)

    workflow.run()


@pytest.mark.skip(reason='needs refactoring')
def test_fsl_apply_transform_func_to_mni_linear():

    c, strat = configuration_strategy_mock(method='FSL')

    strat.append_name('anat_mni_flirt_register_0')

    # build the workflow
    workflow = pe.Workflow(name='test_fsl_apply_transform_func_to_mni_linear')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = fsl_apply_transform_func_to_mni(workflow, 'mean_functional_to_standard',
            'mean_functional', 'template_brain_for_func_preproc', 0, strat,
            c.funcRegFSLinterpolation)

    workflow.run()


@pytest.mark.skip(reason='needs refactoring')
def test_fsl_apply_transform_func_to_mni_linear_mapnode():

    c, strat = configuration_strategy_mock(method='FSL')

    strat.append_name('anat_mni_flirt_register_0')

    # build the workflow
    workflow = pe.Workflow(name='test_fsl_apply_transform_func_to_mni_linear_mapnode')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = fsl_apply_transform_func_to_mni(workflow, 'dr_tempreg_to_standard',
            'dr_tempreg_maps_files', 'template_brain_for_func_preproc', 0, strat,
            c.funcRegFSLinterpolation, map_node=True)

    workflow.run()
