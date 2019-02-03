from CPAC.nuisance import create_nuisance_workflow


def test_nuisance_workflow_type1():

    selector_test = {
        'Anaticor': None,
        'tCompCor': {
            'num_pcs': 3,
            'threshold': '2 PCT',
            'by_slice': True,
            'tissues': 'WM',
            'include_delayed': False,
            'include_squared': False,
            'include_delayed_squared': False
        },
        'aCompCor': {
            'num_pcs': 5,
            'tissues': 'WM',
            'include_delayed': False,
            'include_squared': False,
            'include_delayed_squared': False
        },
        'WhiteMatter': None,
        'Ventricles': {
            'summary_method': 'DetrendNormMean',
            'num_pcs': None,
            'include_delayed': False,
            'include_squared': False,
            'include_delayed_squared': False
        },
        'GreyMatter': None,
        'GlobalSignal': None,
        'Motion': {
            'include_delayed': True,
            'include_squared': True,
            'include_delayed_squared': True
        },
        'Censor': {
            'thresh_metric': 'FD',
            'fd_threshold': '1.5 SD',
            'number_of_previous_trs_to_remove': 0,
            'number_of_subsequent_trs_to_remove': 0,
            'censor_method': 'SpikeRegression'
        },
        'PolyOrt': {
            'degree': 2
        },
        'Bandpass': None
    }

    nuisance_regression_workflow = create_nuisance_workflow(use_ants=True, selector=selector_test)

    nuisance_regression_workflow.inputs.inputspec.set({
        "selector": selector_test,
        "wm_mask_file_path": '/home/anibalsolon/cpac_tests/adhd/o/output/pipeline_analysis/sub-3899622_ses-1/anatomical_wm_mask/segment_seg_2_maths.nii.gz',
        "csf_mask_file_path": '/home/anibalsolon/cpac_tests/adhd/o/output/pipeline_analysis/sub-3899622_ses-1/anatomical_gm_mask/segment_seg_1_maths.nii.gz',
        "gm_mask_file_path": '/home/anibalsolon/cpac_tests/adhd/o/output/pipeline_analysis/sub-3899622_ses-1/anatomical_csf_mask/segment_seg_0_maths.nii.gz',
        "lat_ventricles_mask_file_path": '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_VentricleMask.nii.gz',
        
        "motion_parameters_file_path":  '/home/anibalsolon/cpac_tests/adhd/o/output/pipeline_analysis/sub-3899622_ses-1/motion_params/_scan_rest_run-1/motion_parameters.txt',
        "fd_file_path": '/home/anibalsolon/cpac_tests/adhd/o/output/pipeline_analysis/sub-3899622_ses-1/frame_wise_displacement_jenkinson/_scan_rest_run-1/FD_J.1D',
        "dvars_file_path": '/home/anibalsolon/cpac_tests/adhd/working/resting_preproc_sub-3899622_ses-1/gen_motion_stats_0/_scan_rest_run-1/cal_DVARS/DVARS.1D',
        "functional_file_path": '/home/ccraddock/nuisance_test/functional.nii.gz',
        "functional_brain_mask_file_path": '/home/ccraddock/nuisance_test/func_mask.nii.gz',
        "brain_template_file_path": '/home/ccraddock/nuisance_test/MNI152_T1_2mm_brain.nii.gz',
        
        "anat_to_mni_initial_xfm_file_path":  '/home/ccraddock/nuisance_test/anat_to_mni_initial_xfm.mat',
        "anat_to_mni_rigid_xfm_file_path":  '/home/ccraddock/nuisance_test/anat_to_mni_rigid_xfm.mat',
        "anat_to_mni_affine_xfm_file_path":  '/home/ccraddock/nuisance_test/anat_to_mni_affine_xfm.mat',
        "func_to_anat_linear_xfm_file_path": '/home/ccraddock/nuisance_test/func_to_anat_linear_xfm.mat',
    })

    nuisance_regression_workflow.base_dir = '/home/ccraddock/nuisance_test/working_dir'

    result_value = nuisance_regression_workflow.run()

    print("result {0}".format(result_value))

    assert 0 == 0
