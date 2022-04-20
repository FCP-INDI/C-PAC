'''Validation schema for C-PAC pipeline configurations'''
# pylint: disable=too-many-lines
from itertools import chain, permutations

import numpy as np
from voluptuous import All, ALLOW_EXTRA, Any, In, Length, Match, Optional, \
                       Range, Required, Schema
from voluptuous.error import Error as VoluptuousError
from voluptuous.validators import Maybe

from CPAC.utils import coerce_to_list
from CPAC.utils.utils import lookup_nested_value, set_nested_value
from .backwards_compatibility import backwards_compatible
from .constants import ALWAYS_LISTS, ANTs_PARAMETERS, forkable, MUTEX, \
                       Number,  permutation_message, RESOLUTION_REGEX, \
                       TOGGLED_OPTIONS, VALID_OPTIONS
from .exceptions import handle_custom_error

latest_schema = Schema({
    'FROM': Maybe(str),
    'pipeline_setup': {
        'pipeline_name': All(str, Length(min=1)),
        'output_directory': {
            'path': str,
            'source_outputs_dir': Maybe(str),
            'pull_source_once': bool,
            'write_func_outputs': bool,
            'write_debugging_outputs': bool,
            'output_tree': str,
            'quality_control': {
                'generate_quality_control_images': bool,
                'generate_xcpqc_files': bool}},
        'working_directory': {'path': str, 'remove_working_dir': bool},
        'log_directory': {'run_logging': bool, 'path': str},
        'crash_log_directory': {'path': Maybe(str)},
        'system_config': {
            'FSLDIR': Maybe(str),
            'on_grid': {
                'run': bool,
                'resource_manager': Maybe(str),
                'SGE': {
                    'parallel_environment': Maybe(str),
                    'queue': Maybe(str)}},
            'maximum_memory_per_participant': Number,
            'raise_insufficient': bool,
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_OMP_threads': int,
            'num_participants_at_once': int,
            'random_seed': Maybe(Any(
                'random',
                All(int, Range(min=1, max=np.iinfo(np.int32).max)))),
            'observed_usage': {'callback_log': Maybe(str), 'buffer': Number}},
        'Amazon-AWS': {'aws_output_bucket_credentials': Maybe(str),
                       's3_encryption': bool},
        'Debugging': {'verbose': bool}},
    'anatomical_preproc': {
        'run': bool,
        'run_t2': bool,
        'non_local_means_filtering': {'run': forkable,
                                      'noise_model': Maybe(str)},
        'n4_bias_field_correction': {'run': forkable, 'shrink_factor': int},
        't1t2_bias_field_correction': Required(
            # require 'T1w_brain_ACPC_template' if 'acpc_target' is 'brain'
            Any({'run': False, 'BiasFieldSmoothingSigma': Maybe(int)},
                {'run': True, 'BiasFieldSmoothingSigma': Maybe(int)})),
        'acpc_alignment': Required(
            # require 'T1w_brain_ACPC_template' and
            # 'T2w_brain_ACPC_template' if 'acpc_target' is 'brain'
            Any({'run': False,
                 'run_before_preproc': Maybe(bool),
                 'brain_size': Maybe(int),
                 'FOV_crop': Maybe(In({'robustfov', 'flirt'})),
                 'acpc_target': Maybe(In(VALID_OPTIONS['acpc']['target'])),
                 'align_brain_mask': Maybe(bool),
                 'T1w_ACPC_template': Maybe(str),
                 'T1w_brain_ACPC_template': Maybe(str),
                 'T2w_ACPC_template': Maybe(str),
                 'T2w_brain_ACPC_template': Maybe(str)},
                {'run': True,
                 'run_before_preproc': bool,
                 'brain_size': int,
                 'FOV_crop': In({'robustfov', 'flirt'}),
                 'acpc_target': VALID_OPTIONS['acpc']['target'][1],
                 'align_brain_mask': Maybe(bool),
                 'T1w_ACPC_template': str,
                 'T1w_brain_ACPC_template': Maybe(str),
                 'T2w_ACPC_template': Maybe(str),
                 'T2w_brain_ACPC_template': Maybe(str)},
                {'run': True,
                 'run_before_preproc': bool,
                 'brain_size': int,
                 'FOV_crop': In({'robustfov', 'flirt'}),
                 'acpc_target': VALID_OPTIONS['acpc']['target'][0],
                 'align_brain_mask': Maybe(bool),
                 'T1w_ACPC_template': str,
                 'T1w_brain_ACPC_template': str,
                 'T2w_ACPC_template': Maybe(str),
                 'T2w_brain_ACPC_template': Maybe(str)}),
            msg='\'brain\' requires \'T1w_brain_ACPC_template\' and '
                '\'T2w_brain_ACPC_template\' to '
                'be populated if \'run\' is not set to Off'),
        'brain_extraction': {
            'run': bool,
            'using': [In(VALID_OPTIONS['brain_extraction']['using'])],
            'AFNI-3dSkullStrip': {
                'mask_vol': bool,
                'shrink_factor': Number,
                'var_shrink_fac': bool,
                'shrink_factor_bot_lim': Number,
                'avoid_vent': bool,
                'n_iterations': int,
                'pushout': bool,
                'touchup': bool,
                'fill_hole': int,
                'NN_smooth': int,
                'smooth_final': int,
                'avoid_eyes': bool,
                'use_edge': bool,
                'exp_frac': Number,
                'push_to_edge': bool,
                'use_skull': bool,
                'perc_int': Number,
                'max_inter_iter': int,
                'fac': Number,
                'blur_fwhm': Number,
                'monkey': bool},
            'FSL-FNIRT': {
                'interpolation': In({'trilinear', 'sinc', 'spline'})},
            'FSL-BET': {
                'frac': Number,
                'mask_boolean': bool,
                'mesh_boolean': bool,
                'outline': bool,
                'padding': bool,
                'radius': int,
                'reduce_bias': bool,
                'remove_eyes': bool,
                'robust': bool,
                'skull': bool,
                'surfaces': bool,
                'threshold': bool,
                'vertical_gradient': Range(min=-1, max=1)},
            'UNet': {'unet_model': str},
            'niworkflows-ants': {
                'template_path': str,
                'mask_path': str,
                'regmask_path': str},
            'FreeSurfer-BET': {'T1w_brain_template_mask_ccs': str}}},
    'segmentation': {
        'run': bool,
        'tissue_segmentation': {
            'using': [In(
                {'FSL-FAST', 'FreeSurfer', 'ANTs_Prior_Based',
                 'Template_Based'})],
            'FSL-FAST': {
                'thresholding': {
                    'use': In({'Auto', 'Custom'}),
                    'Custom': {
                        'CSF_threshold_value': float,
                        'WM_threshold_value': float,
                        'GM_threshold_value': float}},
                'use_priors': {
                    'run': bool,
                    'priors_path': str,
                    'WM_path': str,
                    'GM_path': str,
                    'CSF_path': str}},
            'FreeSurfer': {
                'erode': int,
                'CSF_label': [int],
                'GM_label': [int],
                'WM_label': [int]},
            'ANTs_Prior_Based': {
                'run': forkable,
                'template_brain_list': [str],
                'template_segmentation_list': [str],
                'CSF_label': [int],
                'GM_label': [int],
                'WM_label': [int]},
            'Template_Based': {
                'run': forkable,
                'template_for_segmentation': [In(
                    VALID_OPTIONS['segmentation']['template'])],
                'WHITE': str,
                'GRAY': str,
                'CSF': str}}},
    'registration_workflows': {
        'anatomical_registration': {
            'run': bool,
            'resolution_for_anat': All(str, Match(RESOLUTION_REGEX)),
            'T1w_brain_template': str,
            'T1w_template': str,
            'T1w_brain_template_mask': Maybe(str),
            'reg_with_skull': bool,
            'registration': {
                'using': [In({'ANTS', 'FSL', 'FSL-linear'})],
                'ANTs': {
                    'use_lesion_mask': bool,
                    'T1_registration': Maybe(ANTs_PARAMETERS),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    })},
                'FSL-FNIRT': {
                    'fnirt_config': Maybe(str),
                    'ref_resolution': All(str, Match(RESOLUTION_REGEX)),
                    'FNIRT_T1w_brain_template': Maybe(str),
                    'FNIRT_T1w_template': Maybe(str),
                    'interpolation': In({'trilinear', 'sinc', 'spline'}),
                    'identity_matrix': str,
                    'ref_mask': str,
                    'ref_mask_res-2': str,
                    'T1w_template_res-2': str}},
            'overwrite_transform': {'run': bool, 'using': In({'FSL'})}},
        'functional_registration': {
            'coregistration': {
                'run': bool,
                'reference': In({'brain', 'restore-brain'}),
                'interpolation': In({'trilinear', 'sinc', 'spline'}),
                'using': str,
                'input': str,
                'cost': str,
                'dof': int,
                'arguments': Maybe(str),
                'func_input_prep': {
                    'reg_with_skull': bool,
                    'input': [In({
                        'Mean_Functional', 'Selected_Functional_Volume',
                        'fmriprep_reference'})],
                    'Mean Functional': {'n4_correct_func': bool},
                    'Selected Functional Volume': {
                        'func_reg_input_volume': int}},
                'boundary_based_registration': {
                    'run': forkable,
                    'bbr_schedule': str,
                    'bbr_wm_map': In({'probability_map',
                                      'partial_volume_map'}),
                    'bbr_wm_mask_args': str,
                    'reference': In({'whole-head', 'brain'})}},
            'EPI_registration': {
                'run': bool,
                'using': [In({'ANTS', 'FSL', 'FSL-linear'})],
                'EPI_template': str,
                'EPI_template_mask': Maybe(str),
                'ANTs': {
                    'parameters': Maybe(ANTs_PARAMETERS),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    })},
                'FSL-FNIRT': {
                    'fnirt_config': Maybe(str),
                    'interpolation': In({'trilinear', 'sinc', 'spline'}),
                    'identity_matrix': str}},
            'func_registration_to_template': {
                'run': bool,
                'run_EPI': bool,
                'output_resolution': {
                    'func_preproc_outputs': All(
                        str, Match(RESOLUTION_REGEX)),
                    'func_derivative_outputs': All(
                        str, Match(RESOLUTION_REGEX))},
                'target_template': {
                    'using': [In({'T1_template', 'EPI_template'})],
                    'T1_template': {
                        'T1w_brain_template_funcreg': str,
                        'T1w_template_funcreg': Maybe(str),
                        'T1w_brain_template_mask_funcreg': Maybe(str),
                        'T1w_template_for_resample': Maybe(str)},
                    'EPI_template': {
                        'EPI_template_funcreg': str,
                        'EPI_template_mask_funcreg': Maybe(str),
                        'EPI_template_for_resample': Maybe(str)}},
                'ANTs_pipelines': {
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'})},
                'FNIRT_pipelines': {
                    'interpolation': In({'trilinear', 'sinc', 'spline'}),
                    'identity_matrix': str},
                'apply_transform': {
                    'using': In({'default', 'abcd', 'single_step_resampling',
                                 'dcan_nhp'})}}}},
    'surface_analysis': {
        'freesurfer': {
            'run': bool,
            'reconall_args': Maybe(str),
            'freesurfer_dir': Maybe(str)},
        'post_freesurfer': {
            'run': bool,
            'surf_atlas_dir': Maybe(str),
            'gray_ordinates_dir': Maybe(str),
            'gray_ordinates_res': Maybe(int),
            'high_res_mesh': Maybe(int),
            'low_res_mesh': Maybe(int),
            'subcortical_gray_labels': Maybe(str),
            'freesurfer_labels': Maybe(str),
            'fmri_res': Maybe(int),
            'smooth_fwhm': Maybe(int)}},
    'longitudinal_template_generation': {
        'run': bool,
        'average_method': In({'median', 'mean', 'std'}),
        'dof': In({12, 9, 7, 6}),
        'interp': In({'trilinear', 'nearestneighbour', 'sinc', 'spline'}),
        'cost': In({
            'corratio', 'mutualinfo', 'normmi', 'normcorr', 'leastsq',
            'labeldiff', 'bbr'}),
        'thread_pool': int,
        'convergence_threshold': Number},
    'functional_preproc': {
        'run': bool,
        'truncation': {'start_tr': int, 'stop_tr': Maybe(Any(int, 'End'))},
        'scaling': {'run': bool, 'scaling_factor': Number},
        'despiking': {'run': forkable},
        'slice_timing_correction': {
            'run': forkable,
            'tpattern': Maybe(str),
            'tzero': Maybe(int)},
        'motion_estimates_and_correction': {
            'run': bool,
            'motion_estimates': {
                'calculate_motion_first': bool,
                'calculate_motion_after': bool},
            'motion_correction': {
                'using': [In({'3dvolreg', 'mcflirt'})],
                'AFNI-3dvolreg': {
                    'functional_volreg_twopass': bool},
                'motion_correction_reference': [In({
                    'mean', 'median', 'selected_volume', 'fmriprep_reference'
                })],
                'motion_correction_reference_volume': int},
            'motion_estimate_filter': {
                'run': forkable,
                'filters': Any(
                    None,  # no filters
                    VALID_OPTIONS['motion_estimate_filter'],  # one filter
                    [VALID_OPTIONS['motion_estimate_filter']]  # filter series
                )}},
        'distortion_correction': {
            'run': forkable,
            'using': [In(['PhaseDiff', 'Blip', 'Blip-FSL-TOPUP'])],
            'PhaseDiff': {
                'fmap_skullstrip_option': In(['BET', 'AFNI']),
                'fmap_skullstrip_BET_frac': float,
                'fmap_skullstrip_AFNI_threshold': float},
            'Blip-FSL-TOPUP': {
                'warpres': int,
                'subsamp': int,
                'fwhm': int,
                'miter': int,
                'lambda': int,
                'ssqlambda': int,
                'regmod': In({'bending_energy', 'membrane_energy'}),
                'estmov': int,
                'minmet': int,
                'splineorder': int,
                'numprec': str,
                'interp': In({'spline', 'linear'}),
                'scale': int,
                'regrid': int}},
        'func_masking': {
            'using': [In(
                ['AFNI', 'FSL', 'FSL_AFNI', 'Anatomical_Refined',
                 'Anatomical_Based', 'Anatomical_Resampled',
                 'CCS_Anatomical_Refined'])],
            # handle validating mutually-exclusive booleans for FSL-BET
            # functional_mean_boolean must be True if one of the mutually-
            # exclusive options are
            # see MUTEX definition for more definition
            'FSL-BET': Any(*(
                # exactly one mutually exclusive option on
                [{k: d[k] for d in r for k in d} for r in [[
                    {**MUTEX['FSL-BET']['rem'],
                     'functional_mean_boolean': True,
                     k1: True,
                     k2: False} for k2 in MUTEX['FSL-BET']['mutex'] if k2 != k1
                ] for k1 in MUTEX['FSL-BET']['mutex']]] + [{
                    # no mutually-exclusive options on
                    **MUTEX['FSL-BET']['rem'],
                    'functional_mean_boolean': bool,
                    **{k: False for k in MUTEX['FSL-BET']['mutex']}}])),
            'FSL_AFNI': {
                'bold_ref': str,
                'brain_mask': str,
                'brain_probseg': str},
            'Anatomical_Refined': {
                'anatomical_mask_dilation': bool},
            'apply_func_mask_in_native_space': bool},
        'generate_func_mean': {
            'run': bool},
        'normalize_func': {
            'run': bool}},
    'nuisance_corrections': {
        '1-ICA-AROMA': {
            'run': forkable,
            'denoising_type': In({'aggr', 'nonaggr'})},
        '2-nuisance_regression': {
            'run': forkable,
            'create_regressors': bool,
            'Regressors': Maybe([Schema({
                'Name': Required(str),
                'Censor': {
                    'method': str,
                    'thresholds': [{
                        'type': str,
                        'value': float}],
                    'number_of_previous_trs_to_censor': Maybe(int),
                    'number_of_subsequent_trs_to_censor': Maybe(int)},
                'Motion': {
                    'include_delayed': bool,
                    'include_squared': bool,
                    'include_delayed_squared': bool},
                'aCompCor': VALID_OPTIONS['Regressors']['CompCor'],
                'tCompCor': VALID_OPTIONS['Regressors']['CompCor'],
                'CerebrospinalFluid': VALID_OPTIONS[
                    'Regressors'
                ]['segmentation'],
                'WhiteMatter': VALID_OPTIONS[
                    'Regressors'
                ]['segmentation'],
                'GreyMatter': VALID_OPTIONS[
                    'Regressors'
                ]['segmentation'],
                'GlobalSignal': {'summary': str},
                'PolyOrt': {'degree': int},
                'Bandpass': {
                    'bottom_frequency': float,
                    'top_frequency': float,
                    'method': str,
                }  # how to check if [0] is > than [1]?
            }, extra=ALLOW_EXTRA)]),
            'lateral_ventricles_mask': Maybe(str),
            'bandpass_filtering_order': Maybe(
                In({'After', 'Before'})),
            'regressor_masks': {
                'erode_anatomical_brain_mask': {
                    'run': bool,
                    'brain_mask_erosion_prop': Number,
                    'brain_mask_erosion_mm': Number,
                    'brain_erosion_mm': Number
                },
                'erode_csf': {
                    'run': bool,
                    'csf_erosion_prop': Number,
                    'csf_mask_erosion_mm': Number,
                    'csf_erosion_mm': Number,
                },
                'erode_wm': {
                    'run': bool,
                    'wm_erosion_prop': Number,
                    'wm_mask_erosion_mm': Number,
                    'wm_erosion_mm': Number,
                },
                'erode_gm': {
                    'run': bool,
                    'gm_erosion_prop': Number,
                    'gm_mask_erosion_mm': Number,
                    'gm_erosion_mm': Number}}}},
    'amplitude_low_frequency_fluctuation': {
        'run': bool,
        'highpass_cutoff': [float],
        'lowpass_cutoff': [float]},
    'voxel_mirrored_homotopic_connectivity': {
        'run': bool,
        'symmetric_registration': {
            'T1w_brain_template_symmetric': str,
            'T1w_brain_template_symmetric_for_resample': str,
            'T1w_template_symmetric': str,
            'T1w_template_symmetric_for_resample': str,
            'dilated_symmetric_brain_mask': str,
            'dilated_symmetric_brain_mask_for_resample': str}},
    'regional_homogeneity': {
        'run': bool,
        'cluster_size': In({7, 19, 27})},
    'post_processing': {
        'spatial_smoothing': {
            'output': [In({'smoothed', 'nonsmoothed'})],
            'smoothing_method': [In({'FSL', 'AFNI'})],
            'fwhm': [int]},
        'z-scoring': {
            'output': [In({'z-scored', 'raw'})]}},
    'timeseries_extraction': {
        'run': bool,
        Optional('roi_paths_fully_specified'): bool,
        'tse_roi_paths': Optional(
            Maybe({
                str: In({', '.join(options) for options in
                         list(chain.from_iterable([list(
                             permutations(
                                 VALID_OPTIONS['timeseries']['roi_paths'],
                                 number_of)
                         ) for number_of in range(1, 6)]))})}),
            msg=permutation_message('tse_roi_paths',
                                    VALID_OPTIONS['timeseries']['roi_paths'])),
        'realignment': In({'ROI_to_func', 'func_to_ROI'}),
        'connectivity_matrix': {
            option: Maybe([In(VALID_OPTIONS['connectivity_matrix'][option])])
            for option in ['using', 'measure']}},
    'seed_based_correlation_analysis': {
        'run': bool,
        Optional('roi_paths_fully_specified'): bool,
        'sca_roi_paths': Optional(
            Maybe({
                str: In({', '.join(options) for options in
                         list(chain.from_iterable([list(
                              permutations(VALID_OPTIONS['sca']['roi_paths'],
                                           number_of)
                              ) for number_of in range(1, 4)]))})}),
            msg=permutation_message('sca_roi_paths',
                                    VALID_OPTIONS['sca']['roi_paths'])),
        'norm_timeseries_for_DR': bool},
    'network_centrality': {
        'run': bool,
        'memory_allocation': Number,
        'template_specification_file': str,
        'degree_centrality': {
            'weight_options': [In(
                VALID_OPTIONS['centrality']['weight_options'])],
            'correlation_threshold_option': In(
                VALID_OPTIONS['centrality']['threshold_options']),
            'correlation_threshold': Range(min=-1, max=1)},
        'eigenvector_centrality': {
            'weight_options': [In(
                VALID_OPTIONS['centrality']['weight_options'])],
            'correlation_threshold_option': In(
                VALID_OPTIONS['centrality']['threshold_options']),
            'correlation_threshold': Range(min=-1, max=1)},
        'local_functional_connectivity_density': {
            'weight_options': [In(
                VALID_OPTIONS['centrality']['weight_options'])],
            'correlation_threshold_option': In([
                o for o in VALID_OPTIONS['centrality']['threshold_options'] if
                o != 'Sparsity threshold']),
            'correlation_threshold': Range(min=-1, max=1)}},
    'PyPEER': {
        'run': bool,
        'eye_scan_names': [str],
        'data_scan_names': [str],
        'eye_mask_path': str,
        'stimulus_path': Maybe(str),
        'minimal_nuisance_correction': {
            'peer_gsr': bool,
            'peer_scrub': bool,
            'scrub_thresh': float}}
})


def schema(config_dict):
    '''Function to test the schema validity of a given config dict.

    Parameters
    ----------
    config_dict : dict

    Returns
    -------
    voluptuous.schema_builder.Schema
        validated configuration schema
    '''
    for keys in ALWAYS_LISTS:
        config_dict = coerce_to_list(config_dict, keys)
    config_dict = backwards_compatible(config_dict)
    for option in TOGGLED_OPTIONS:
        # Allow options to be mutually incompatible if those options are off
        switch = lookup_nested_value(config_dict, option['switch'])
        if True not in switch:
            latest_schema.schema = set_nested_value(latest_schema.schema,
                                                    option['key'],
                                                    option['Off'])
    try:
        return latest_schema(config_dict)
    except VoluptuousError as voluptuous_error:
        return handle_custom_error(voluptuous_error, config_dict)


schema.schema = latest_schema.schema
