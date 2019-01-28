from CPAC.nuisance import create_nuisance_workflow


def test_nuisance_workflow_type1():

    """
    test_selector = {'Anaticor' : None | {radius = <radius in mm>},
        'aCompCor' : None | {num_pcs = <number of components to retain>,
                            tissues = 'WM' | 'CSF' | 'WM+CSF',
                            include_delayed = True | False,
                            include_squared = True | False,
                            include_delayed_squared = True | False},
        'WhiteMatter' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                       num_pcs = <number of components to retain>,
                       include_delayed = True | False,
                       include_squared = True | False,
                       include_delayed_squared = True | False},
        'Ventricles' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                       num_pcs = <number of components to retain>,
                       include_delayed = True | False,
                       include_squared = True | False,
                       include_delayed_squared = True | False},
        'GreyMatter' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                       num_pcs = <number of components to retain>,
                       include_delayed = True | False,
                       include_squared = True | False,
                       include_delayed_squared = True | False},
        'GlobalSignal' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                           num_pcs = <number of components to retain>,
                           include_delayed = True | False,
                           include_squared = True | False,
                           include_delayed_squared = True | False},
        'Motion' : None | {include_delayed = True | False,
                           include_squared = True | False,
                           include_delayed_squared = True | False},
        'Censor' : None | { thresh_metric = 'RMSD','DVARS', or 'RMSD+DVARS',
                            threshold = <threshold to be applied to metric, if using
                              RMSD+DVARS, this should be a tuple (RMSD thresh, DVARS thresh)>,
                            number_of_previous_trs_to_remove = True | False,
                            number_of_subsequent_trs_to_remove = True | False,
                            method = 'Kill', 'Zero', 'Interpolate', 'SpikeRegression'},
        'PolyOrt' : None | { degree = <polynomial degree up to which will be removed, e.g. 2 means
                                       constant + linear + quadratic, practically that is probably,
                                       the most that will be need esp. if band pass filtering>},
        'Bandpass' : None | { bottom_frequency = <frequency in hertz of the highpass part of the pass
                                                  band, frequencies below this will be removed>,
                              top_frequency = <frequency in hertz of the lowpass part of the pass
                                               band, frequencies above this will be removed>},
        }

    """

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
        "wm_mask_file_path": '/home/ccraddock/nuisance_test/wm_mask.nii.gz',
        "csf_mask_file_path": '/home/ccraddock/nuisance_test/csf_mask.nii.gz',
        "gm_mask_file_path": '/home/ccraddock/nuisance_test/gm_mask.nii.gz',
        "lat_ventricles_mask_file_path": '/home/ccraddock/nuisance_test/MNI152_T1_2mm_VentricleMask.nii.gz',
        
        "motion_parameters_file_path":  '/home/ccraddock/nuisance_test/motion_parameters.1D',
        "fd_file_path": '/home/ccraddock/nuisance_test/fd.1D',
        "dvars_file_path": '/home/ccraddock/nuisance_test/dvars.1d',
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
