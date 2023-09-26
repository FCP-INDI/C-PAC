import glob
from os import environ
import pytest
import yaml

from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util

from CPAC.nuisance import create_nuisance_regression_workflow
from CPAC.nuisance.utils import NuisanceRegressor

selector_config = """
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
     include_backdiff: True
     include_backdiff_squared: True

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


@pytest.mark.skip(reason='needs refactoring')
def test_nuisance_workflow_type1():

    base_dir = '/tmp/nuisance_working_dir'

    try:
        import shutil
        shutil.rmtree(base_dir)
    except:
        pass

    selector_test = yaml.safe_load(selector_config)['Regressors']

    for selector in selector_test:

        nuisance_regression_workflow = \
            create_nuisance_regression_workflow(
                nuisance_selectors=NuisanceRegressor(selector),
                use_ants=True,
                name='nuisance'
            )

        nuisance_regression_workflow.write_graph(graph2use='orig', simple_form=False)

        preprocessed = '/cc_dev/cpac_working/old_compcor'

        nuisance_regression_workflow.inputs.inputspec.anatomical_file_path = glob.glob(preprocessed + '/anat_preproc_already_0/anat_reorient/*_resample.nii.gz')[0]
        nuisance_regression_workflow.inputs.inputspec.functional_file_path = glob.glob(preprocessed + '/func_preproc_automask_0/_scan_*/func_normalize/*_calc_tshift_resample_volreg_calc_maths.nii.gz')[0]
        nuisance_regression_workflow.inputs.inputspec.functional_brain_mask_file_path = glob.glob(preprocessed + '/func_preproc_automask_0/_scan_*/func_mask_normalize/*_calc_tshift_resample_volreg_calc_maths_maths.nii.gz')[0]
        nuisance_regression_workflow.inputs.inputspec.wm_mask_file_path = glob.glob(preprocessed + '/seg_preproc_0/WM/WM_mask/segment_seg_2_maths.nii.gz')[0]
        nuisance_regression_workflow.inputs.inputspec.csf_mask_file_path = glob.glob(preprocessed + '/seg_preproc_0/GM/GM_mask/segment_seg_1_maths.nii.gz')[0]
        nuisance_regression_workflow.inputs.inputspec.gm_mask_file_path = glob.glob(preprocessed + '/seg_preproc_0/GM/GM_mask/segment_seg_1_maths.nii.gz')[0]
        nuisance_regression_workflow.inputs.inputspec.lat_ventricles_mask_file_path = glob.glob(
            f'{environ.get("FSLDIR")}/data/standard/'
            'MNI152_T1_2mm_VentricleMask.nii.gz')[0]
        # nuisance_regression_workflow.inputs.inputspec.mni_to_anat_linear_xfm_file_path = glob.glob('')[0]
        nuisance_regression_workflow.inputs.inputspec.func_to_anat_linear_xfm_file_path = glob.glob(preprocessed + '/func_to_anat_bbreg_0/_scan_*/bbreg_func_to_anat/*_calc_tshift_resample_volreg_calc_tstat_flirt.mat')[0]
        nuisance_regression_workflow.inputs.inputspec.anat_to_mni_initial_xfm_file_path = glob.glob(preprocessed + '/anat_mni_ants_register_0/calc_ants_warp/transform0DerivedInitialMovingTranslation.mat')[0]
        nuisance_regression_workflow.inputs.inputspec.anat_to_mni_rigid_xfm_file_path = glob.glob(preprocessed + '/anat_mni_ants_register_0/calc_ants_warp/transform1Rigid.mat')[0]
        nuisance_regression_workflow.inputs.inputspec.anat_to_mni_affine_xfm_file_path = glob.glob(preprocessed + '/anat_mni_ants_register_0/calc_ants_warp/transform2Affine.mat')[0]
        nuisance_regression_workflow.inputs.inputspec.motion_parameters_file_path = glob.glob(preprocessed + '/func_preproc_automask_0/_scan_*/func_motion_correct_A/*_calc_tshift_resample.1D')[0]
        nuisance_regression_workflow.inputs.inputspec.dvars_file_path = glob.glob(preprocessed + '/gen_motion_stats_0/_scan_*/cal_DVARS/DVARS.*')[0]
        nuisance_regression_workflow.inputs.inputspec.fd_j_file_path = glob.glob(preprocessed + '/gen_motion_stats_0/_scan_*/calculate_FDJ/FD_J.1D')[0]
        nuisance_regression_workflow.get_node('inputspec').iterables = ([
            ('selector', [NuisanceRegressor(selector)]),
        ])

        nuisance_regression_workflow.base_dir = base_dir

        nuisance_regression_workflow.run()
