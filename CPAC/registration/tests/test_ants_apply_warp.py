import os
import nipype.pipeline.engine as pe
from ..output_func_to_standard import ants_apply_warps_func_mni
from mocks import configuration_strategy_mock

def test_ants_apply_warp_func_mni():

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock()

    # build the workflow
    workflow = pe.Workflow(name='test_ants_apply_warps_func_mni')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = ants_apply_warps_func_mni(workflow, strat, 0, 8, 
            'mean_functional', 'mean_functional', 'mean_functional_to_standard', 
            c.funcRegANTSinterpolation, 'template_brain_for_func_preproc', 0,
            distcor=False)

    workflow = ants_apply_warps_func_mni(workflow, strat, 0, 8, 
            'mean_functional_to_standard', 'mean_functional', 'mean_functional_standard_to_original', 
            c.funcRegANTSinterpolation, 'template_brain_for_func_preproc', 0,
            distcor=False, inverse=True)

    workflow.run()

def test_ants_apply_warp_func_mni_mapnode():

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock()

    # build the workflow
    workflow = pe.Workflow(name='test_ants_apply_warps_func_mni')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = ants_apply_warps_func_mni(workflow, strat, 0, 8, 
            'dr_tempreg_maps_files', 'mean_functional', 'dr_tempreg_maps_to_standard', 
            c.funcRegANTSinterpolation, 'template_brain_for_func_preproc', 0,
            distcor=False, map_node=True)

    workflow = ants_apply_warps_func_mni(workflow, strat, 0, 8, 
            'dr_tempreg_maps_to_standard', 'mean_functional', 'dr_tempreg_maps_standard_to_original', 
            c.funcRegANTSinterpolation, 'template_brain_for_func_preproc', 0,
            distcor=False, map_node=True, inverse=True)


    workflow.run()

