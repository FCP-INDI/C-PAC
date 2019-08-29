import os
import yaml

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.utils.interfaces.function import Function

'''
tC-1.5PCT-PC5S-SDB
aC-WC-2mmE-PC5-SDB
WM-2mmE-PC5-SDB
CSF-2mmE-M-SDB
GM-2mmE-DNM-SDB
G-PC5-SDB
M-SDB
C-S-FD1.5SD-D1.5SD
PR-2
BP-T0.01-B0.1
'''

selector = """
Regressors:

 - PolyOrt:
     degree: 2
  
   tCompCor:
     summary:
       method: PC
       components: 5
     threshold: 1.5SD
     by_slice: true
 
   aCompCor:
     summary:
       method: PC
       components: 5
     tissues:
       - WhiteMatter
       - CerebrospinalFluid
     extraction_resolution: 2
   
   WhiteMatter:
     summary:
       method: PC
       components: 5
     extraction_resolution: 2

   CerebrospinalFluid:
     summary:
       method: PC
       components: 5
     extraction_resolution: 2
     erode_mask: true
     
   GlobalSignal:
     summary: Mean
     include_delayed: True
     include_squared: True
     include_delayed_squared: True

   Motion:
     include_delayed: True
     include_squared: True
     include_delayed_squared: True

   Censor:
     method: Interpolate
     thresholds:
       - type: FD_J
         value: 0.5
       - type: DVARS
         value: 17
         
"""

def identity_input(name, field, val):

    pool_input = pe.Node(util.IdentityInterface(
        fields=[field]),
        name=name
    )

    pool_input.inputs.set({ field: val })

    return pool_input, field



def test_nuisance_workflow_type1():

    return

    selector_test = yaml.load(selector)['Regressors']

    nuisance_regression_workflow, pipeline_resource_pool = \
        create_nuisance_workflow(
            nuisance_selectors=selector_test,
            use_ants=True,
            name='nuisance'
        )

    nuisance_regression_workflow.write_graph(graph2use='orig', simple_form=False)

    preprocessed = '/home/anibalsolon/cpac_tests/adhd/working/resting_preproc_sub-3899622_ses-1'

    nuisance_regression_workflow.inputs.inputspec.functional_file_path = preprocessed + '/func_preproc_automask_0/_scan_rest_run-1/func_normalize/sub-3899622_ses-1_task-rest_run-1_bold_calc_tshift_resample_volreg_calc_maths.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.functional_brain_mask_file_path = preprocessed + '/func_preproc_automask_0/_scan_rest_run-1/func_mask_normalize/sub-3899622_ses-1_task-rest_run-1_bold_calc_tshift_resample_volreg_calc_maths_maths.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.wm_mask_file_path = preprocessed + '/seg_preproc_0/WM/WM_mask/segment_seg_2_maths.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.csf_mask_file_path = preprocessed + '/seg_preproc_0/GM/GM_mask/segment_seg_1_maths.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.gm_mask_file_path = preprocessed + '/seg_preproc_0/GM/GM_mask/segment_seg_1_maths.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.lat_ventricles_mask_file_path = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_VentricleMask.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.mni_to_anat_linear_xfm_file_path = ''
    nuisance_regression_workflow.inputs.inputspec.func_to_anat_linear_xfm_file_path = preprocessed + '/anat_mni_ants_register_0/func_to_anat_bbreg_0/_scan_rest_run-1/bbreg_func_to_anat/sub-3899622_ses-1_task-rest_run-1_bold_calc_tshift_resample_volreg_calc_tstat_flirt.mat'
    nuisance_regression_workflow.inputs.inputspec.anat_to_mni_initial_xfm_file_path = preprocessed + '/anat_mni_ants_register_0/calc_ants_warp/transform0DerivedInitialMovingTranslation.mat'
    nuisance_regression_workflow.inputs.inputspec.anat_to_mni_rigid_xfm_file_path = preprocessed + '/anat_mni_ants_register_0/calc_ants_warp/transform1Rigid.mat'
    nuisance_regression_workflow.inputs.inputspec.anat_to_mni_affine_xfm_file_path = preprocessed + '/anat_mni_ants_register_0/calc_ants_warp/transform2Affine.mat'
    nuisance_regression_workflow.inputs.inputspec.motion_parameters_file_path = preprocessed + '/func_preproc_automask_0/_scan_rest_run-1/func_motion_correct_A/sub-3899622_ses-1_task-rest_run-1_bold_calc_tshift_resample.1D'
    nuisance_regression_workflow.inputs.inputspec.dvars_file_path = preprocessed + '/gen_motion_stats_0/_scan_rest_run-1/cal_DVARS/DVARS.1D'
    nuisance_regression_workflow.inputs.inputspec.fd_file_path =    preprocessed + '/gen_motion_stats_0/_scan_rest_run-1/calculate_FDJ/FD_J.1D'
    nuisance_regression_workflow.inputs.inputspec.brain_template_file_path = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
    nuisance_regression_workflow.inputs.inputspec.iterables = (
        ('selector', selector_test)
    )

    nuisance_regression_workflow.base_dir = '/tmp/working_dir'

    try:
        import shutil
        shutil.rmtree(nuisance_regression_workflow.base_dir)
    except:
        pass

    result_value = nuisance_regression_workflow.run()

    print("result {0}".format(result_value))

    assert 0 == 0


def summarize_timeseries(selector):

    return 'residual_output', 'regressors_output'



class NuisanceRegressor(object):

    def __init__(self, selectors):
        self.selectors = selectors

    def __str__(self):
        return "Stevenson"

    def _derivative_params(self, selector):
        nr_repr = ''
        if selector['include_squared']:
            nr_repr += 'S'
        if selector['include_delayed']:
            nr_repr += 'D'
        if selector['include_delayed_squared']:
            nr_repr += 'B'
        return nr_repr

    def _summary_params(self, selector):
        summ = selector['summary']
        nr_repr = ''
        if selector['include_squared']:
            nr_repr += 'S'
        if selector['include_delayed']:
            nr_repr += 'D'
        if selector['include_delayed_squared']:
            nr_repr += 'B'
        return nr_repr

    def __repr__(self):
        regs = [
            'GreyMatter',
            'WhiteMatter',
            'CerebrospinalFluid',
            'tCompCor',
            'aCompCor',
            'GlobalSignal',
            'Motion',
            'PolyOrt',
            'Censor',
            'Bandpass',
        ]

        nr_repr = ""
  
        for r in regs:
            if r not in self.selectors:
                continue
            
            nr_repr += "a"

        # tC-1.5PCT-PC5S-SDB
        # aC-WC-2mmE-PC5-SDB
        # WM-2mmE-PC5-SDB
        # CSF-2mmE-M-SDB
        # GM-2mmE-DNM-SDB
        # G-PC5-SDB
        # M-SDB
        # C-S-FD1.5SD-D1.5SD
        # PR-2
        # BP-T0.01-B0.1
        return "Stevenfather"



def test_iterable_selector():

    selector_test = yaml.load(selector)['Regressors']

    nuisance_wf = pe.Workflow(name='iterable_selector')
    nuisance_wf.base_dir = '/tmp/iter_working_dir'

    try:
        import shutil
        shutil.rmtree(nuisance_wf.base_dir)
    except:
        pass


    inputspec = pe.Node(util.IdentityInterface(fields=[
        'selector'
    ]), name='inputspec')

    summarize_timeseries_node = pe.Node(
        Function(
            input_names=[
                'selector'
            ],
            output_names=['residual_file_path',
                          'regressors_file_path'],
            function=summarize_timeseries,
            as_module=True,
        ),
        name='summarize_timeseries'
    )

    outputspec = pe.Node(util.IdentityInterface(fields=['residual_file_path',
                                                        'regressors_file_path']),
                         name='outputspec')


    nuisance_wf.connect(inputspec, 'selector', summarize_timeseries_node, 'selector')
    nuisance_wf.connect(summarize_timeseries_node, 'residual_file_path', outputspec, 'residual_file_path')
    nuisance_wf.connect(summarize_timeseries_node, 'regressors_file_path', outputspec, 'regressors_file_path')
    
    nuisance_wf.get_node('inputspec').iterables = (
        ('selector', [NuisanceRegressor(s) for s in selector_test])
    )

    nuisance_wf.run()