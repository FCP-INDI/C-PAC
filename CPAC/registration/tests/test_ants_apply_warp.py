import os
from CPAC.pipeline import nipype_pipeline_engine as pe
from ..output_func_to_standard import ants_apply_warps_func_mni
from .mocks import configuration_strategy_mock
import CPAC.utils.test_init as test_utils
import nibabel as nb
import pytest

@pytest.mark.skip(reason="no way of currently testing this")
def test_ants_apply_warp_func_mni():

    test_name = 'test_ants_apply_warps_func_mni'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock()
    num_strat = 0

    node, out = strat['mean_functional']
    mean_functional = node.inputs.file

    # build the workflow
    workflow = pe.Workflow(name='test_ants_apply_warps_func_mni')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = ants_apply_warps_func_mni(workflow,
            'mean_functional_to_standard', 
            'mean_functional',
            'template_brain_for_func_preproc',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=False,
            inverse=False,
            input_image_type=0,
            num_ants_cores=1)

    workflow = ants_apply_warps_func_mni(workflow,
            'mean_functional_standard_to_original', 
            'mean_functional_to_standard',
            'mean_functional',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=False,
            inverse=True,
            input_image_type=0,
            num_ants_cores=1)

    workflow.run()

    mean_functional_after_transform=os.path.join(c.workingDirectory, test_name,
        'apply_ants_warp_mean_functional_standard_to_original_inverse_0',
        'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_antswarp_antswarp.nii.gz')

    assert(test_utils.pearson_correlation(mean_functional, mean_functional_after_transform) > .99)


@pytest.mark.skip(reason="no way of currently testing this")
def test_ants_apply_warps_func_mni_mapnode():

    test_name = 'test_ants_apply_warps_func_mni_mapnode'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock()
    num_strat = 0

    node, out = strat['dr_tempreg_maps_files']
    dr_spatmaps = node.inputs.file

    # build the workflow
    workflow = pe.Workflow(name='test_ants_apply_warps_func_mni_mapnode')
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = ants_apply_warps_func_mni(workflow,
            'dr_tempreg_maps_to_standard', 
            'dr_tempreg_maps_files',
            'template_brain_for_func_preproc',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=True,
            inverse=False,
            input_image_type=0,
            num_ants_cores=1)

    workflow = ants_apply_warps_func_mni(workflow,
            'dr_tempreg_maps_standard_to_original', 
            'dr_tempreg_maps_to_standard',
            'mean_functional',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=True,
            inverse=True,
            input_image_type=0,
            num_ants_cores=8)

    workflow.run()

    dr_spatmaps_after_transform=[os.path.join(c.workingDirectory, test_name,
        'apply_ants_warp_dr_tempreg_maps_standard_to_original_mapnode_inverse_0',
        'mapflow', '_apply_ants_warp_dr_tempreg_maps_standard_to_original_mapnode_inverse_0{0}'.format(n),
        'temp_reg_map_000{0}_antswarp_antswarp.nii.gz'.format(n))
        for n in range(0,10)]

    
    test_results = [ test_utils.pearson_correlation(orig_file, xformed_file) > 0.99 \
        for orig_file, xformed_file in zip(dr_spatmaps, dr_spatmaps_after_transform)]

    assert all(test_results)


@pytest.mark.skip(reason='needs refactoring')
def test_ants_apply_warp_func_mni_symm():

    test_name = 'test_ants_apply_warps_func_mni_symm'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock()
    num_strat = 0

    node, out = strat['mean_functional']
    mean_functional = node.inputs.file

    # build the workflow
    workflow = pe.Workflow(name=test_name)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = ants_apply_warps_func_mni(workflow,
            'mean_functional_to_standard_symm', 
            'mean_functional',
            'template_brain_for_func_preproc',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=False,
            inverse=False,
            symmetry='symmetric',
            input_image_type=0,
            num_ants_cores=8)

    workflow = ants_apply_warps_func_mni(workflow,
            'mean_functional_standard_to_original_symm', 
            'mean_functional_to_standard_symm',
            'mean_functional',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=False,
            inverse=True,
            symmetry='symmetric',
            input_image_type=0,
            num_ants_cores=1)

    retval = workflow.run()

    mean_functional_after_transform=os.path.join(c.workingDirectory, test_name,
        'apply_ants_warp_mean_functional_standard_to_original_symm_inverse_0',
        'sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_tstat_antswarp_antswarp.nii.gz')

    assert(test_utils.pearson_correlation(mean_functional, mean_functional_after_transform) > .93)


@pytest.mark.skip(reason='needs refactoring')
def test_ants_apply_warps_func_mni_mapnode_symm():

    test_name = 'test_ants_apply_warps_func_mni_mapnode_symm'

    # get the config and strat for the mock
    c, strat = configuration_strategy_mock()
    num_strat = 0

    node, out = strat['dr_tempreg_maps_files']
    dr_spatmaps = node.inputs.file

    # build the workflow
    workflow = pe.Workflow(name=test_name)
    workflow.base_dir = c.workingDirectory
    workflow.config['execution'] = {
        'hash_method': 'timestamp',
        'crashdump_dir': os.path.abspath(c.crashLogDirectory)
    }

    workflow = ants_apply_warps_func_mni(workflow,
            'dr_tempreg_maps_to_standard_symm', 
            'dr_tempreg_maps_files',
            'template_brain_for_func_preproc',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=True,
            inverse=False,
            symmetry='symmetric',
            input_image_type=0,
            num_ants_cores=8)

    workflow = ants_apply_warps_func_mni(workflow,
            'dr_tempreg_maps_standard_symm_to_original', 
            'dr_tempreg_maps_to_standard_symm',
            'mean_functional',
            num_strat, strat, 
            interpolation_method = c.funcRegANTSinterpolation,
            distcor=False,
            map_node=True,
            inverse=True,
            symmetry='symmetric',
            input_image_type=0,
            num_ants_cores=8)

    workflow.run()

    dr_spatmaps_after_transform=[os.path.join(c.workingDirectory, test_name,
        'apply_ants_warp_dr_tempreg_maps_standard_symm_to_original_mapnode_inverse_0',
        'mapflow', '_apply_ants_warp_dr_tempreg_maps_standard_symm_to_original_mapnode_inverse_0{0}'.format(n),
        'temp_reg_map_000{0}_antswarp_antswarp.nii.gz'.format(n))
        for n in range(0,10)]

    r = [test_utils.pearson_correlation(orig_file, xformed_file) \
            for orig_file, xformed_file in zip(dr_spatmaps, dr_spatmaps_after_transform)]

    print(r)
    test_results = [ r_value > 0.93 for r_value in r ] 

    assert all(test_results)

