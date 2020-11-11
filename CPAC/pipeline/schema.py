from itertools import chain, permutations
from voluptuous import Schema, Required, All, Any, Length, Range, Match, In, \
                       ALLOW_EXTRA
from voluptuous.validators import Maybe

centrality_options = {
    'method_options': ['degree_centrality', 'eigenvector_centrality',
                       'local_functional_connectivity_density'],
    'threshold_options': ['Significance threshold', 'Sparsity threshold',
                          'Correlation threshold'],
    'weight_options': ['Binarized', 'Weighted']
}

schema = Schema({
    Required('pipeline_setup'): {
        Required('pipeline_name'): All(str, Length(min=1)),
        Required('working_directory'): {
            Required('path'): str,
        },
        # Required('crash_log_directory'): str,
        Required('log_directory'): {
            Required('path'): str,
        },
        Required('output_directory'): {
            Required('path'): str,
        },
        Required('system_config'): {
            'maximum_memory_per_participant': Any(float, int),
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_participants_at_once': int
        },
    },
    # 
    # 'FSLDIR': str,
    # 'runOnGrid': bool,
    # 'resourceManager': str,
    # 'parallelEnvironment': str,
    # 'queue': str,
    # 
    # 'awsOutputBucketCredentials': str,
    # 's3Encryption': bool,  # check/normalize
    # 
    # 'maximumMemoryPerParticipant': float,
    # 'maxCoresPerParticipant': All(int, Range(min=1)),
    # 'numParticipantsAtOnce': All(int, Range(min=1)),
    # 'num_ants_threads': All(int, Range(min=1)),
    # 
    # 'write_func_outputs': bool,
    # 'write_debugging_outputs': bool,  # check/normalize
    # 'generateQualityControlImages': bool, # check/normalize
    # 'removeWorkingDir': bool,
    # 'run_logging': bool,
    # 'reGenerateOutputs': bool, # check/normalize
    # 'runSymbolicLinks': bool, # check/normalize
    # 
    # 'resolution_for_anat': All(str, Match(r'^[0-9]+mm$')),
    # 'template_brain_only_for_anat': str,
    # 'template_skull_for_anat': str,
    # 'template_symmetric_brain_only': str,
    # 'template_symmetric_skull': str,
    # 'dilated_symmetric_brain_mask': str,
    # 
    # 'already_skullstripped': bool,
    # 'skullstrip_option': In(['AFNI', 'BET']),
    # 'skullstrip_shrink_factor': float,
    # 'skullstrip_var_shrink_fac': bool,
    # 'skullstrip_shrink_factor_bot_lim':float,
    # 'skullstrip_avoid_vent': bool,
    # 'skullstrip_n_iterations': int,
    # 'skullstrip_pushout': bool,
    # 'skullstrip_touchup': bool,
    # 'skullstrip_fill_hole': int,
    # 'skullstrip_NN_smooth': int,
    # 'skullstrip_smooth_final': int,
    # 'skullstrip_avoid_eyes': bool,
    # 'skullstrip_use_edge': bool,
    # 'skullstrip_exp_frac': float,
    # 'skullstrip_push_to_edge': bool,  # check/normalize
    # 'skullstrip_use_skull': bool,  # check/normalize
    # 'skullstrip_perc_int': int,
    # 'skullstrip_max_inter_iter': int,
    # 'skullstrip_fac': int,
    # 'skullstrip_blur_fwhm': int,
    # 'bet_frac': float,
    # 'bet_mask_boolean': bool,  # check/normalize
    # 'bet_mesh_boolean': bool,  # check/normalize
    # 'bet_outline': bool,  # check/normalize
    # 'bet_padding': bool,  # check/normalize
    # 'bet_radius': int,
    # 'bet_reduce_bias': bool,  # check/normalize
    # 'bet_remove_eyes': bool,  # check/normalize
    # 'bet_robust': bool,  # check/normalize
    # 'bet_skull': bool,  # check/normalize
    # 'bet_surfaces': bool,  # check/normalize
    # 'bet_threshold': bool,  # check/normalize
    # 'bet_vertical_gradient': float,
    Required('anatomical_preproc'): {
        Required('registration_workflow'): {
            Required('registration'): {
                Required('using'): [In({'ANTS', 'FSL'})],
                'ANTs': {
                    'EPI_registration': Any(
                        None, 'None', dict, [dict]
                    ),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    }),
                    'use_lesion_mask': bool
                },
                'FSL-FNIRT': {
                    'interpolation': In({
                        'trilinear', 'sinc', 'spline'
                    }),
                },
            },
            'reg_with_skull': bool,
        },
        Required('segmentation_workflow'): {
            '1-segmentation': {
                'ANTs_Prior_Based': {
                    Required('run'): [bool],
                },
                'Template_Based': {
                    'template_for_segmentation': [
                        In({'EPI Template', 'T1 Template'})
                    ],
                },
                'WHITE': str,
                'GRAY': str,
                'CSF': str,
            },
            '2-use_priors': {
                'run': bool,
                'priors_path': str,
                'WM_path': str,
                'GM_path': str,
                'CSF_path': str
            },
            '3-custom_thresholding': {
                'run': bool,
                'CSF_threshold_value': float,
                'WM_threshold_value': float,
                'GM_threshold_value': float
            },
            Required('4-erosion'): {
                Required('erode_anatomical_brain_mask'): {
                    Required('run'): bool,
                    'brain_mask_erosion_prop': Any(float, int),
                    'brain_mask_erosion_mm': Any(float, int),
                    'brain_erosion_mm': Any(float, int)
                },
                'erode_csf': {
                    'run': bool,
                    'csf_erosion_prop': Any(float, int),
                    'csf_mask_erosion_mm': Any(float, int),
                    'csf_erosion_mm': Any(float, int),
                },
                'erode_wm': {
                    'run': bool,
                    'wm_erosion_prop': Any(float, int),
                    'wm_mask_erosion_mm': Any(float, int),
                    'wm_erosion_mm': Any(float, int),
                },
                'erode_gm': {
                    'run': bool,
                    'gm_erosion_prop': Any(float, int),
                    'gm_mask_erosion_mm': Any(float, int),
                    'gm_erosion_mm': Any(float, int),
                }
            }
        }
    },
    # 
    # 'fnirtConfig': str,
    # 'ref_mask': str,
    # 'regWithSkull': [bool], # check/normalize
    # 
    # 'runSegmentationPreprocessing': [bool], # check/normalize
    # 'priors_path': str,
    # 
    # 'slice_timing_correction': [bool], # check/normalize
    # 'TR': Any(None, float),
    # 'slice_timing_pattern': Any(str, int), # check for nifti header option
    # 'startIdx': All(int, Range(min=0)),
    # 'stopIdx': Any(None, All(int, Range(min=1))),
    # 
    # 'runEPI_DistCorr': [bool], # check/normalize
    # 'fmap_distcorr_skullstrip': [In(['BET', 'AFNI'])],
    # 'fmap_distcorr_frac': float, # check if it needs to be a list
    # 'fmap_distcorr_threshold': float,
    # 'fmap_distcorr_deltaTE': float, # check if it needs to be a list
    # 'fmap_distcorr_dwell_time': float, # check if it needs to be a list
    # 'fmap_distcorr_dwell_asym_ratio': float, # check if it needs to be a list
    # 'fmap_distcorr_pedir': In(["x", "y", "z", "-x", "-y", "-z"]),
    # 
    # 'runBBReg': [bool], # check/normalize
    # 'boundaryBasedRegistrationSchedule': str,
    # 
    # 'func_reg_input': In(['Mean Functional', 'Selected Functional Volume']),
    # 'func_reg_input_volume': All(int, Range(min=0)),
    # 'functionalMasking': [In(['3dAutoMask', 'BET'])],
    # 
    # 'runRegisterFuncToMNI': [bool], # check/normalize
    # 'resolution_for_func_preproc': All(str, Match(r'^[0-9]+mm$')),
    Required('functional_registration'): {
        Required('1-coregistration'): {
            Required('run'): [bool],
            'func_input_prep': {
                Required('input'): [In({
                    'Mean Functional', 'Selected Functional Volume'
                })]
            },
            'boundary_based_registration': {
                Required('run'): [bool],
                'bbr_schedule': str
            }
        },
        Required('2-func_registration_to_template'): {
            Required('target_template'): {
                Required('using'): [In({'T1_template', 'EPI_template'})]
            },
            'output_resolution': {
                'func_derivative_outputs': All(str, Match(r'^[0-9]+mm$')),
            }
        }
    },
    # 'template_brain_only_for_func': str,
    # 'template_skull_for_func': str,
    # 'identityMatrix': str,
    # 'configFileTwomm': str,
    Required('nuisance_corrections'): {
        Required('1-ICA-AROMA'): {
            Required('run'): [bool],
            'denoising_type': In({'aggr', 'nonaggr'}),
        },
        Required('2-nuisance_regression'): {
            Required('run'): [bool],
            'Regressors': [{
                'Motion': {
                    'include_delayed': bool,
                    'include_squared': bool,
                    'include_delayed_squared': bool
                },
                'aCompCor': {
                    'summary': {
                        'method': str,
                        'components': int
                    },
                    'tissues': [str],
                    'extraction_resolution': int
                },
                'CerebrospinalFluid': {
                    'summary': str,
                    'extraction_resolution': int,
                    'erode_mask': bool
                },
                'GlobalSignal': {'summary': str},
                'PolyOrt': {'degree': int},
                'Bandpass': {
                    'bottom_frequency': float,
                    'top_frequency': float
                }  # how to check if [0] is > than [1]?
            }],
            'lateral_ventricles_mask': str,
            Required('bandpass_filtering_order'): In({'After', 'Before'})
        },
        '3-median-angle-correction': {
            Required('run'): [bool],
            Required('target_angle_deg'): Any(int, float),
        }
    },
    Required('amplitude_low_frequency_fluctuation'): {
        Required('run'): bool,
        'highpass_cutoff': [float],
        'lowpass_cutoff': [float],
    },
    Required('voxel_mirrored_homotopic_connectivity'): {
        Required('run'): bool,
    },
    Required('regional_homogeneity'): {
        Required('run'): bool,
        'clusterSize': In({7, 19, 27}),
    },
    # 
    # 'nComponents': int, # check if list
    # 'runFristonModel': [bool], # check/normalize
    # 'runMotionSpike': [Any(None, In(['despiking', 'scrubbing']))], # check/normalize, check None
    # 
    # 'fdCalc': In(['power', 'jenkinson']), # check if it needs to be a list
    # 'spikeThreshold': float, # check if it needs to be a list
    # 'numRemovePrecedingFrames': int,
    # 'numRemoveSubsequentFrames': int,
    'timeseries_extraction': {
        Required('run'): bool,
        Required('tse_roi_paths'): Any(None, {
            str: In({', '.join([
                option for option in options
            ]) for options in list(chain.from_iterable([list(
                permutations({'Avg', 'Voxel', 'SpatialReg', 'PearsonCorr',
                              'PartialCorr'}, number_of)
            ) for number_of in range(1, 6)]))}),
        }),
        'realignment': In({'ROI_to_func', 'func_to_ROI'}),
        'roi_tse_outputs': Any(None, [In({None, 'csv', 'numpy'})]),
    },

    Required('seed_based_correlation_analysis'): {
        Required('run'): bool,
        'sca_roi_paths': Any(None, {
            str: In({', '.join([
                option for option in options
            ]) for options in list(chain.from_iterable([list(
                permutations({'Avg', 'DualReg', 'MultReg'}, number_of)
            ) for number_of in range(1, 4)]))})
        }),
    },
    # 'sca_roi_paths': Any(None, {
    #     str: [In(['average', 'dual_regression', 'multiple_regression'])],
    #     # normalize before running thrugh schema
    # }),
    # 'mrsNorm': bool,
    # 
    # 
    # 'run_smoothing': [bool], # check/normalize
    # 'fwhm': float,
    # 'smoothing_order': In(['after', 'before']),
    # 
    # 'runZScoring': [bool], # check/normalize
    # 
    # 'run_fsl_feat': bool, # check/normalize
    # 'numGPAModelsAtOnce': int,
    # 'modelConfigs': [str],
    # 
    # 'run_basc': bool,
    # 'basc_resolution': All(str, Match(r'^[0-9]+mm$')),
    # 'basc_proc': int,
    # 'basc_memory': float,
    # 'basc_roi_mask_file': str,
    # 'basc_cross_cluster_mask_file': str,
    # 'basc_similarity_metric_list': [In(['correlation', 'euclidean', 'cityblock', 'cosine'])],
    # 'basc_timeseries_bootstrap_list': int,
    # 'basc_dataset_bootstrap_list': int,
    # 'basc_n_clusters_list': int,
    # 'basc_affinity_thresh': float,
    # 'basc_output_sizes': int,
    # 'basc_cross_cluster': bool,
    # 'basc_blocklength_list': int,
    # 'basc_group_dim_reduce': bool,
    # 'basc_inclusion': str,
    # 'basc_pipeline': str,
    # 'basc_scan_inclusion': str,
    # 
    # 'runMDMR': bool,
    # 'mdmr_inclusion': str,
    # 'mdmr_roi_file': str,
    # 'mdmr_regressor_file': str,
    # 'mdmr_regressor_participant_column': str,
    # 'mdmr_regressor_columns': str,
    # 'mdmr_permutations': int,
    # 'mdmr_parallel_nodes': int,
    # 
    # 'runISC': bool,
    # 'runISFC': bool,
    # 'isc_voxelwise': bool,
    # 'isc_roiwise': bool,
    # 'isc_permutations': int,
    Required('network_centrality'): {
        Required('run'): [bool],
        'memory_allocation': Any(float, int),
        'template_specification_file': str,
        'degree_centrality': {
            'weight_options': [Maybe(In(
                centrality_options['weight_options']
            ))],
            'correlation_threshold_option': In(
                centrality_options['threshold_options']),
            'correlation_threshold': Range(min=-1, max=1)
        },
        'eigenvector_centrality': {
            'weight_options': [Maybe(In(
                centrality_options['weight_options']
            ))],
            'correlation_threshold_option': In(
                centrality_options['threshold_options']
            ),
            'correlation_threshold': Range(min=-1, max=1)
        },
        'local_functional_connectivity_density': {
            'weight_options': [Maybe(In(
                centrality_options['weight_options']
            ))],
            'correlation_threshold_option': In([
                o for o in centrality_options['threshold_options'] if
                o != 'Sparsity threshold'
            ]),
            'correlation_threshold': Range(min=-1, max=1)
        },
    },
    Required('PyPEER'): {
        Required('run'): [bool],
    },
}, extra=ALLOW_EXTRA)
