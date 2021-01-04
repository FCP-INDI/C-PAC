from itertools import chain, permutations
from voluptuous import Schema, Required, All, Any, Length, Range, Match, In, \
                       ALLOW_EXTRA
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

schema = Schema({
    'pipeline_setup': {
        Required('pipeline_name'): All(str, Length(min=1)),
        Required('output_directory'): {
            Required('path'): str,
            'write_func_outputs': bool,
            'write_debugging_outputs': bool,
            'output_tree': str,
            'generate_quality_control_images': bool,
        },
        Required('working_directory'): {
            Required('path'): str,
            'remove_working_dir': bool,
        },
        Required('log_directory'): {
            'run_logging': bool,
            Required('path'): str,
        },
        Required('crash_log_directory'): {
            Required('path'): str,
        },
        Required('system_config'): {
            'on_grid':{
                Required('run'): bool,
                'resource_manager': Any(None, str),
                'SGE': {
                    'parallel_environment': Any(None, str),
                    'queue': Any(None, str),
                },
            },
            'maximum_memory_per_participant': Any(float, int),
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_participants_at_once': int
        },
        'Amazon-AWS':{
            'aws_output_bucket_credentials': Any(None, str),
            's3_encryption': bool,
        },
        'pipeline_IMPORT': Any(None, str),
    },
    'FSLDIR': Any(None, str),

    'anatomical_preproc': {
        Required('run'): bool,
        Required('non_local_means_filtering'): bool,
        Required('n4_bias_field_correction'): bool,
        Required('acpc_alignment'): {
            'run': bool,
            'brain_size': int,
            'template_skull': str,
            'template_brain': Any(None, str),
        },
        Required('brain_extraction'):{
            'already_skullstripped': bool,
            'extraction': {
                Required('using'): [In({'3dSkullStrip', 'BET', 'UNet', 'niworkflows-ants'})],
                'AFNI-3dSkullStrip': {
                    'mask_vol': bool,
                    'shrink_factor': Any(float, int),
                    'var_shrink_fac': bool,
                    'shrink_factor_bot_lim':Any(float, int),
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
        },
        Required('segmentation_workflow'): {
            'run': bool,
            '1-segmentation': {
                'using': [In({'FSL-FAST', 'ANTs_Prior_Based', 'Template_Based'})],
                'ANTs_Prior_Based': {
                    Required('run'): [bool],
                    'template_brain_list': list,
                    'template_segmentation_list': list,
                    'CSF_label': int,
                    'left_GM_label': int,
                    'right_GM_label': int,
                    'left_WM_label': int,
                    'right_WM_label': int,
                },
                'Template_Based': {
                    Required('run'): [bool],
                    'template_for_segmentation': [
                        In(valid_options['segmentation']['template'])
                    ],
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
        },
        Required('registration_workflow'): {
            Required('resolution_for_anat'): All(str, Match(r'^[0-9]+mm$')),
            'template_brain_only_for_anat': str,
            'template_skull_for_anat': str,
            'reg_with_skull': bool,
            Required('registration'): {
                Required('using'): [In({'ANTS', 'FSL'})],
                'ANTs': {
                    'use_lesion_mask': bool,
                    'T1_registration': Any(
                        None, 'None', dict, [dict]
                    ),
                    'EPI_registration': Any(
                        None, 'None', dict, [dict]
                    ),
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
    'functional_registration': {
        Required('1-coregistration'): {
            Required('run'): [bool],
            'func_input_prep': {
                Required('input'): [In({
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
                Required('run'): [bool],
                'bbr_schedule': str
            }
        },
        Required('2-func_registration_to_template'): {
            Required('run'): bool,
            Required('output_resolution'):{
                'func_preproc_outputs': All(str, Match(r'^[0-9]+mm$')),
                'func_derivative_outputs': All(str, Match(r'^[0-9]+mm$')),
                'template_for_resample': str,
            },
            Required('target_template'): {
                Required('using'): [In({'T1_template', 'EPI_template'})],
                'T1_template': {
                    'template_brain': str,
                    'template_skull': str,
                },
                'EPI_template': {
                    'template_epi': str, 
                },
            },
            Required('ANTs_pipelines'): {
                'interpolation': In({'Linear', 'BSpline', 'LanczosWindowedSinc'})
            },
            Required('FNIRT_pipelines'): {
                'interpolation': In({'trilinear', 'sinc', 'spline'}),
                'identity_matrix': str
            },
        },
    },

    'nuisance_corrections': {
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
    },
    'amplitude_low_frequency_fluctuation': {
        Required('run'): bool,
        'highpass_cutoff': [float],
        'lowpass_cutoff': [float],
    },

    'voxel_mirrored_homotopic_connectivity': {
        Required('run'): bool,
    },

    'regional_homogeneity': {
        'run': bool,
        'cluster_size': In({7, 19, 27}),
    },
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

    'seed_based_correlation_analysis': {
        Required('run'): bool,
        'sca_roi_paths': Any(None, {
            str: In({', '.join([
                option for option in options
            ]) for options in list(chain.from_iterable([list(
                permutations({'Avg', 'DualReg', 'MultReg'}, number_of)
            ) for number_of in range(1, 4)]))})
        }),
        'norm_timeseries_for_DR': bool,
    },
    'network_centrality': {
        Required('run'): [bool],
        'memory_allocation': Any(float, int),
        'template_specification_file': str,
        'degree_centrality': {
            'weight_options': [Maybe(In(
                valid_options['centrality']['weight_options']
            ))],
            'correlation_threshold_option': In(
                valid_options['centrality']['threshold_options']),
            'correlation_threshold': Range(min=-1, max=1)
        },
        'eigenvector_centrality': {
            'weight_options': [Maybe(In(
                valid_options['centrality']['weight_options']
            ))],
            'correlation_threshold_option': In(
                valid_options['centrality']['threshold_options']
            ),
            'correlation_threshold': Range(min=-1, max=1)
        },
        'local_functional_connectivity_density': {
            'weight_options': [Maybe(In(
                valid_options['centrality']['weight_options']
            ))],
            'correlation_threshold_option': In([
                o for o in valid_options['centrality']['threshold_options'] if
                o != 'Sparsity threshold'
            ]),
            'correlation_threshold': Range(min=-1, max=1)
        },
    },
    'PyPEER': {
        Required('run'): [bool],
    },
}, extra=ALLOW_EXTRA)
