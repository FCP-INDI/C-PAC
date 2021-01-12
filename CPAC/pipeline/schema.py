from itertools import chain, permutations
from voluptuous import All, Any, In, Length, Match, Range, Schema
from voluptuous.validators import Maybe

valid_options = {
    'centrality': {
       'method_options': ['degree_centrality', 'eigenvector_centrality',
                          'local_functional_connectivity_density'],
       'threshold_options': ['Significance threshold', 'Sparsity threshold',
                             'Correlation threshold'],
       'weight_options': ['Binarized', 'Weighted']
    },
    'segmentation': {
        'using': ['FSL-FAST', 'ANTs_Prior_Based', 'Template_Based'],
        'template': ['EPI Template', 'T1 Template']
    }
}
mutex = {  # mutually exclusive booleans
    'FSL-BET': {
        # exactly zero or one of each of the following can be True for FSL-BET
        'mutex': ['reduce_bias', 'robust', 'padding', 'remove_eyes',
                  'surfaces'],
        # the remaining keys: validators for FSL-BET
        'rem': {
            'frac': float,
            'mesh_boolean': bool,
            'outline': bool,
            'radius': int,
            'skull': bool,
            'threshold': bool,
            'vertical_gradient': Range(min=-1, max=1, min_included=False,
                                       max_included=False),
        }
    }
}
resolution_regex = r'^(x*[0-9](\.[0-9]+)*mm)*$'

schema = Schema({
    'FROM': Maybe(str),
    'pipeline_setup': {
        'pipeline_name': All(str, Length(min=1)),
        'output_directory': {
            'path': str,
            'write_func_outputs': bool,
            'write_debugging_outputs': bool,
            'output_tree': str,
            'generate_quality_control_images': bool,
        },
        'working_directory': {
            'path': str,
            'remove_working_dir': bool,
        },
        'log_directory': {
            'run_logging': bool,
            'path': str,
        },
        'crash_log_directory': {
            'path': str,
        },
        'system_config': {
            'on_grid': {
                'run': bool,
                'resource_manager': Maybe(str),
                'SGE': {
                    'parallel_environment': Maybe(str),
                    'queue': Maybe(str),
                },
            },
            'maximum_memory_per_participant': Any(float, int),
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_participants_at_once': int
        },
        'Amazon-AWS': {
            'aws_output_bucket_credentials': Maybe(str),
            's3_encryption': bool,
        },
    },
    'FSLDIR': Maybe(str),
    'anatomical_preproc': {
        'run': bool,
        'non_local_means_filtering': bool,
        'n4_bias_field_correction': bool,
        'acpc_alignment': {
            'run': bool,
            'brain_size': int,
            'template_skull': str,
            'template_brain': Maybe(str),
        },
        'brain_extraction': {
            'already_skullstripped': bool,
            'extraction': {
                'using': [In({
                    '3dSkullStrip', 'BET', 'UNet', 'niworkflows-ants',
                    'FreeSurfer-ABCD'
                })],
                'AFNI-3dSkullStrip': {
                    'mask_vol': bool,
                    'shrink_factor': Any(float, int),
                    'var_shrink_fac': bool,
                    'shrink_factor_bot_lim': Any(float, int),
                    'avoid_vent': bool,
                    'n_iterations': int,
                    'pushout': bool,
                    'touchup': bool,
                    'fill_hole': int,
                    'NN_smooth': int,
                    'smooth_final': int,
                    'avoid_eyes': bool,
                    'use_edge': bool,
                    'exp_frac': Any(float, int),
                    'push_to_edge': bool,
                    'use_skull': bool,
                    'perc_int': Any(float, int),
                    'max_inter_iter': int,
                    'fac': Any(float, int),
                    'blur_fwhm': Any(float, int),
                    'monkey': bool,
                },
                'FSL-FNIRT': {
                    'interpolation': In({
                        'trilinear', 'sinc', 'spline'
                    }),
                },
                'FSL-BET': {
                    'frac': Any(float, int),
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
                    'vertical_gradient': Range(min=-1, max=1)
                },
                'UNet': {
                    'unet_model': str,
                },
                'niworkflows-ants': {
                    'template_path': str,
                    'mask_path': str,
                    'regmask_path': str,
                },
            },
            'reg_with_skull': bool,
        },
        'segmentation_workflow': {
            'run': [bool],
            '1-segmentation': {
                'using': [
                    In({'FSL-FAST', 'ANTs_Prior_Based', 'Template_Based'})
                ],
                'ANTs_Prior_Based': {
                    'run': Maybe([bool]),
                    'template_brain_list': [str],
                    'template_segmentation_list': [str],
                    'CSF_label': int,
                    'left_GM_label': int,
                    'right_GM_label': int,
                    'left_WM_label': int,
                    'right_WM_label': int,
                },
                'Template_Based': {
                    'run': [bool],
                    'template_for_segmentation': [In(
                        valid_options['segmentation']['template'])],
                    'WHITE': str,
                    'GRAY': str,
                    'CSF': str,
                },
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
            '4-erosion': {
                'erode_anatomical_brain_mask': {
                    'run': bool,
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
        },
        'registration_workflow': {
            'resolution_for_anat': All(str, Match(resolution_regex)),
            'template_brain_only_for_anat': str,
            'template_skull_for_anat': str,
            'reg_with_skull': bool,
            'registration': {
                'using': [In({'ANTS', 'FSL'})],
                'ANTs': {
                    'use_lesion_mask': bool,
                    'T1_registration': Maybe(Any(
                        'None', dict, [dict]
                    )),
                    'EPI_registration': Maybe(Any(
                        'None', dict, [dict]
                    )),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    }),
                },
                'FSL-FNIRT': {
                    'fsl_linear_reg_only': list,
                    'fnirt_config': str,
                    'ref_mask': str,
                    'interpolation': In({
                        'trilinear', 'sinc', 'spline'
                    }),
                },
            },
        },
    },
    'longitudinal_template_generation': {
        'run': bool,
        'average_method': In({'median', 'mean', 'std'}),
        'dof': In({12, 9, 7, 6}),
        'interp': In({'trilinear', 'nearestneighbour', 'sinc', 'spline'}),
        'cost': In({
            'corratio', 'mutualinfo', 'normmi', 'normcorr', 'leastsq',
            'labeldiff', 'bbr'}),
        'thread_pool': int,
        'convergence_threshold': Any(int, float),
    },
    'functional_preproc': {
        'run': bool,
        'truncation': {
            'start_tr': int,
            'stop_tr': Maybe(Any(int, 'End'))
        },
        'scaling': {
            'run': bool,
            'scaling_factor': Any(int, float)
        },
        'despiking': {
            'run': [bool]
        },
        'slice_timing_correction': {
            'run': [bool]
        },
        'motion_estimates_and_correction': {
            'calculate_motion_first': [bool],
            'motion_correction': {
                'using': [In({'3dvolreg', 'mcflirt'})],
                'AFNI-3dvolreg': {
                    'functional_volreg_twopass': bool,
                },
                'motion_correction_reference': [In({
                    'mean', 'median', 'selected volume'})],
                'motion_correction_reference_volume': int,
            },
            'motion_estimate_filter': {
                'run': [bool],
                'filter_type': In({'notch', 'lowpass'}),
                'filter_order': int,
                'breathing_rate_min': Maybe(Any(float, int)),
                'breathing_rate_max': Maybe(Any(float, int)),
                'center_frequency': Maybe(Any(float, int)),
                'filter_bandwidth': Maybe(Any(float, int)),
                'lowpass_cutoff': Maybe(Any(float, int)),
            },
        },
        'distortion_correction': {
            'run': [bool],
            'using': [In(['PhaseDiff', 'Blip'])],
            'PhaseDiff': {
                'fmap_skullstrip_option': In(['BET', 'AFNI']),
                'fmap_skullstrip_frac': float,
                'fmap_skullstrip_threshold': float,
            },
        },
        'func_masking': {
            'using': [In(['AFNI', 'FSL', 'FSL_AFNI', 'Anatomical_Refined'])],
            # handle validating mutually-exclusive booleans for FSL-BET
            # functional_mean_boolean must be True if one of the mutually-
            # exclusive options are
            # see mutex definition for more definition
            'FSL-BET': Any(*(
                # exactly one mutually exclusive option on
                [{k: d[k] for d in r for k in d} for r in [[
                    {
                        **mutex['FSL-BET']['rem'],
                        'functional_mean_boolean': True,
                        k1: True,
                        k2: False
                    } for k2 in mutex['FSL-BET']['mutex'] if k2 != k1
                ] for k1 in mutex['FSL-BET']['mutex']]] +
                # no mutually-exclusive options on
                [{
                    **mutex['FSL-BET']['rem'],
                    'functional_mean_boolean': bool,
                    **{k: False for k in mutex['FSL-BET']['mutex']}
                }])
            ),
            'Anatomical_Refined': {
                'anatomical_mask_dilation': bool,
            },
        },
    },
    'functional_registration': {
        '1-coregistration': {
            'run': [bool],
            'func_input_prep': {
                'input': [In({
                    'Mean Functional', 'Selected Functional Volume'
                })],
                'Mean Functional': {
                    'n4_correct_func': bool
                },
                'Selected Functional Volume': {
                    'func_reg_input_volume': int
                },
            },
            'boundary_based_registration': {
                'run': [bool],
                'bbr_schedule': str
            }
        },
        '2-func_registration_to_template': {
            'run': bool,
            'output_resolution': {
                'func_preproc_outputs': All(str, Match(resolution_regex)),
                'func_derivative_outputs': All(str, Match(resolution_regex)),
                'template_for_resample': str,
            },
            'target_template': {
                'using': [In({'T1_template', 'EPI_template'})],
                'T1_template': {
                    'template_brain': str,
                    'template_skull': str,
                },
                'EPI_template': {
                    'template_epi': str,
                },
            },
            'ANTs_pipelines': {
                'interpolation': In({
                    'Linear', 'BSpline', 'LanczosWindowedSinc'})
            },
            'FNIRT_pipelines': {
                'interpolation': In({'trilinear', 'sinc', 'spline'}),
                'identity_matrix': str
            },
        },
    },
    'nuisance_corrections': {
        '1-ICA-AROMA': {
            'run': [bool],
            'denoising_type': In({'aggr', 'nonaggr'}),
        },
        '2-nuisance_regression': {
            'run': [bool],
            'Regressors': Maybe([{
                'Censor': {
                    'method': str,
                    'thresholds': [{
                        'type': str,
                        'value': float,
                    }],
                    'number_of_previous_trs_to_censor': Maybe(int),
                    'number_of_subsequent_trs_to_censor': Maybe(int),
                },
                'Motion': {
                    'include_delayed': bool,
                    'include_squared': bool,
                    'include_delayed_squared': bool
                },
                'aCompCor': {
                    'degree': int,
                    'erode_mask_mm': bool,
                    'summary': {
                        'method': str,
                        'components': int,
                        'filter': str,
                    },
                    'threshold': str,
                    'tissues': [str],
                    'extraction_resolution': int
                },
                'tCompCor': {
                    'degree': int,
                    'erode_mask_mm': bool,
                    'summary': {
                        'method': str,
                        'components': int,
                        'filter': str,
                    },
                    'threshold': str,
                    'tissues': [str],
                    'extraction_resolution': int
                },
                'CerebrospinalFluid': {
                    'summary': Any(
                        str, {'components': int, 'method': str}
                    ),
                    'extraction_resolution': int,
                    'erode_mask': bool
                },
                'GlobalSignal': {'summary': str},
                'PolyOrt': {'degree': int},
                'WhiteMatter': {
                    'erode_mask': bool,
                    'extraction_resolution': Any(int, float, str),
                    'summary': Any(
                        str, {'components': int, 'method': str}
                    ),
                },
                'Bandpass': {
                    'bottom_frequency': float,
                    'top_frequency': float,
                    'method': str,
                }  # how to check if [0] is > than [1]?
            }]),
            'lateral_ventricles_mask': Maybe(str),
            'bandpass_filtering_order': Maybe(
                In({'After', 'Before'}))
        },
    },
    'amplitude_low_frequency_fluctuation': {
        'run': bool,
        'highpass_cutoff': [float],
        'lowpass_cutoff': [float],
    },
    'voxel_mirrored_homotopic_connectivity': {
        'run': bool,
        'symmetric_registration': {
            'template_symmetric_brain_only': str,
            'template_symmetric_brain_for_resample': str,
            'template_symmetric_skull': str,
            'template_symmetric_skull_for_resample': str,
            'dilated_symmetric_brain_mask': str,
            'dilated_symmetric_brain_mask_for_resample': str,
            'FNIRT_pipelines': {
                'config_file': str
            },
        },
    },
    'regional_homogeneity': {
        'run': bool,
        'cluster_size': In({7, 19, 27}),
    },
    'post_processing': {
        'spatial_smoothing': {
            'run': Any(bool, [bool]),
            'smoothing_method': In({'FSL', 'AFNI'}),
            'fwhm': [int],
            'smoothing_order': In({'Before', 'After'})
        },
        'z-scoring': {
            'run': Any(bool, [bool]),
        },
    },
    'timeseries_extraction': {
        'run': bool,
        'tse_roi_paths': Maybe({
            str: In({', '.join([
                option for option in options
            ]) for options in list(chain.from_iterable([list(
                permutations({'Avg', 'Voxel', 'SpatialReg', 'PearsonCorr',
                              'PartialCorr'}, number_of)
            ) for number_of in range(1, 6)]))}),
        }),
        'realignment': In({'ROI_to_func', 'func_to_ROI'}),
        'roi_tse_outputs': Maybe([In({None, 'csv', 'numpy'})]),
    },

    'seed_based_correlation_analysis': {
        'run': bool,
        'sca_roi_paths': Maybe({
            str: In({', '.join([
                option for option in options
            ]) for options in list(chain.from_iterable([list(
                permutations({'Avg', 'DualReg', 'MultReg'}, number_of)
            ) for number_of in range(1, 4)]))})
        }),
        'norm_timeseries_for_DR': bool,
    },
    'network_centrality': {
        'run': [bool],
        'memory_allocation': Any(float, int),
        'template_specification_file': str,
        'degree_centrality': {
            'weight_options': [In(
                valid_options['centrality']['weight_options']
            )],
            'correlation_threshold_option': In(
                valid_options['centrality']['threshold_options']),
            'correlation_threshold': Range(min=-1, max=1)
        },
        'eigenvector_centrality': {
            'weight_options': [In(
                valid_options['centrality']['weight_options']
            )],
            'correlation_threshold_option': In(
                valid_options['centrality']['threshold_options']
            ),
            'correlation_threshold': Range(min=-1, max=1)
        },
        'local_functional_connectivity_density': {
            'weight_options': [In(
                valid_options['centrality']['weight_options']
            )],
            'correlation_threshold_option': In([
                o for o in valid_options['centrality']['threshold_options'] if
                o != 'Sparsity threshold'
            ]),
            'correlation_threshold': Range(min=-1, max=1)
        },
    },
    'PyPEER': {
        'run': [bool],
        'eye_scan_names': [str],
        'data_scan_names': [str],
        'eye_mask_path': str,
        'stimulus_path': Maybe(str),
        'minimal_nuisance_correction': {
            'peer_gsr': bool,
            'peer_scrub': bool,
            'scrub_thresh': float,
        },
    },
})
