from itertools import chain, permutations
from voluptuous import All, ALLOW_EXTRA, Any, In, Length, Match, Optional, \
                       Range, Required, Schema
from voluptuous.validators import ExactSequence, Maybe
from CPAC import __version__

# 1 or more digits, optional decimal, 'e', optional '-', 1 or more digits
scientific_notation_str_regex = r'^([0-9]+(\.[0-9]*)*(e)-{0,1}[0-9]+)*$'

# (1 or more digits, optional decimal, 0 or more lowercase characters (units))
# ('x',
#  1 or more digits, optional decimal, 0 or more lowercase characters (units)
# ) 0 or more times
resolution_regex = r'^[0-9]+(\.[0-9]*){0,1}[a-z]*' \
                   r'(x[0-9]+(\.[0-9]*){0,1}[a-z]*)*$'

Number = Any(float, int, All(str, Match(scientific_notation_str_regex)))
forkable = Any(bool, [bool])
valid_options = {
    'acpc': {
        'target': ['brain', 'whole-head']
    },
    'boundary_based_registration': {
        'using': ['FSL', 'FreeSurfer']
    },
    'brain_extraction': {
        'using': ['3dSkullStrip', 'BET', 'UNet', 'niworkflows-ants',
                  'FreeSurfer-ABCD']
    },
    'centrality': {
       'method_options': ['degree_centrality', 'eigenvector_centrality',
                          'local_functional_connectivity_density'],
       'threshold_options': ['Significance threshold', 'Sparsity threshold',
                             'Correlation threshold'],
       'weight_options': ['Binarized', 'Weighted']
    },
    'sca': {
        'roi_paths': {'Avg', 'DualReg', 'MultReg'},
    },
    'segmentation': {
        'using': ['FSL-FAST', 'ANTs_Prior_Based', 'Template_Based'],
        'template': ['EPI_Template', 'T1_Template'],
    },
    'timeseries': {
        'roi_paths': {'Avg', 'Voxel', 'SpatialReg', 'PearsonCorr',
                      'PartialCorr'},
    },
    'Regressors': {
        'CompCor': {
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
        'segmentation': {
            'erode_mask': bool,
            'extraction_resolution': Any(
                int, float, All(str, Match(resolution_regex))
            ),
            'include_delayed': bool,
            'include_delayed_squared': bool,
            'include_squared': bool,
            'summary': Any(
                str, {'components': int, 'method': str}
            ),
        },
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
ANTs_parameter_transforms = {
    'gradientStep': Number,
    'metric': {
        'type': str,
        'metricWeight': int,
        'numberOfBins': int,
        'samplingStrategy': str,
        'samplingPercentage': Number,
        'radius': Number,
    },
    'convergence': {
        'iteration': All(str, Match(resolution_regex)),
        'convergenceThreshold': Number,
        'convergenceWindowSize': int,
    },
    'smoothing-sigmas': All(str, Match(resolution_regex)),
    'shrink-factors': All(str, Match(resolution_regex)),
    'use-histogram-matching': bool,
    'updateFieldVarianceInVoxelSpace': Number,
    'totalFieldVarianceInVoxelSpace': Number,
    'winsorize-image-intensities': {
        'lowerQuantile': float,
        'upperQuantile': float,
    },
}
ANTs_parameters = [Any(
    {
        'collapse-output-transforms': int
    }, {
        'dimensionality': int
    }, {
        'initial-moving-transform': {
            'initializationFeature': int,
        },
    }, {
        'transforms': [Any({
            'Rigid': ANTs_parameter_transforms,
        }, {
            'Affine': ANTs_parameter_transforms,
        }, {
            'SyN': ANTs_parameter_transforms,
        })],
    }, {
        'verbose': Any(bool, In({0, 1})),
    }, {
        'float': Any(bool, In({0, 1})),
    }, {
        'masks': {
            'fixed_image_mask': bool,
            'moving_image_mask': bool,
        },
    }, dict  # TODO: specify other valid ANTs parameters
)]
_url_version = 'nightly' if __version__.endswith(
    '-dev') else f'v{__version__.lstrip("v")}'


def permutation_message(key, options):
    '''Function to give a clean, human-readable error message for keys that accept permutation values

    Parameters
    ----------
    key: str

    options: list or set

    Returns
    -------
    msg: str'''  # noqa E501
    return f'''

\'{key}\' takes a dictionary with paths to region-of-interest (ROI)
 NIFTI files (.nii or .nii.gz) as keys and a comma separated string
 of analyses to run. For example, if you wish to run Avg and
 MultReg, you would enter:

    '/path/to/ROI.nii.gz': Avg, MultReg

Available analyses for \'{key}\' are {options}

'''


schema = Schema({
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
            'path': Maybe(str),
        },
        'system_config': {
            'FSLDIR': Maybe(str),
            'on_grid': {
                'run': bool,
                'resource_manager': Maybe(str),
                'SGE': {
                    'parallel_environment': Maybe(str),
                    'queue': Maybe(str),
                },
            },
            'maximum_memory_per_participant': Number,
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_OMP_threads': int,
            'num_participants_at_once': int
        },
        'Amazon-AWS': {
            'aws_output_bucket_credentials': Maybe(str),
            's3_encryption': bool,
        },
        'Debugging': {
            'verbose': bool
        },
    },
    'anatomical_preproc': {
        'run': bool,
        'non_local_means_filtering': forkable,
        'n4_bias_field_correction': forkable,
        'acpc_alignment': Required(
            # require 'T1w_brain_ACPC_template' if 'acpc_target' is 'brain'
            Any({
                'run': False,
                'run_before_preproc': Maybe(bool),
                'brain_size': Maybe(int),
                'acpc_target': Maybe(In(valid_options['acpc']['target'])),
                'T1w_ACPC_template': Maybe(str),
                'T1w_brain_ACPC_template': Maybe(str),
            }, {
                'run': True,
                'run_before_preproc': bool,
                'brain_size': int,
                'acpc_target': valid_options['acpc']['target'][1],
                'T1w_ACPC_template': str,
                'T1w_brain_ACPC_template': Maybe(str),
            }, {
                'run': True,
                'run_before_preproc': bool,
                'brain_size': int,
                'acpc_target': valid_options['acpc']['target'][0],
                'T1w_ACPC_template': str,
                'T1w_brain_ACPC_template': str,
            },),
            msg='\'brain\' requires \'T1w_brain_ACPC_template\' to '
                'be populated if \'run\' is not set to Off',
        ),
        'brain_extraction': {
            'using': [In(valid_options['brain_extraction']['using'])],
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
                'monkey': bool,
            },
            'FSL-FNIRT': {
                'interpolation': In({
                    'trilinear', 'sinc', 'spline'
                }),
            },
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
    'segmentation': {
        'run': bool,
        'tissue_segmentation': {
            'using': [In(
                {'FSL-FAST', 'FreeSurfer', 'ANTs_Prior_Based',
                 'Template_Based'}
            )],
            'FSL-FAST': {
                'thresholding': {
                    'use': In({'Auto', 'Custom'}),
                    'Custom': {
                        'CSF_threshold_value': float,
                        'WM_threshold_value': float,
                        'GM_threshold_value': float,
                    },
                },
                'use_priors': {
                    'run': bool,
                    'priors_path': str,
                    'WM_path': str,
                    'GM_path': str,
                    'CSF_path': str
                },
            },
            'Freesurfer': Maybe(dict),
            'ANTs_Prior_Based': {
                'run': forkable,
                'template_brain_list': [str],
                'template_segmentation_list': [str],
                'CSF_label': int,
                'left_GM_label': int,
                'right_GM_label': int,
                'left_WM_label': int,
                'right_WM_label': int,
            },
            'Template_Based': {
                'run': forkable,
                'template_for_segmentation': [In(
                    valid_options['segmentation']['template']
                )],
                'WHITE': str,
                'GRAY': str,
                'CSF': str,
            },
        },
    },
    'registration_workflows': {
        'anatomical_registration': {
            'run': bool,
            'resolution_for_anat': All(str, Match(resolution_regex)),
            'T1w_brain_template': str,
            'T1w_template': str,
            'T1w_brain_template_mask': Maybe(str),
            'reg_with_skull': bool,
            'registration': {
                'using': [In({'ANTS', 'FSL', 'FSL-linear'})],
                'ANTs': {
                    'use_lesion_mask': bool,
                    'T1_registration': Maybe(ANTs_parameters),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    }),
                },
                'FSL-FNIRT': {
                    'fnirt_config': Maybe(str),
                    'ref_mask': Maybe(str),
                    'interpolation': In({
                        'trilinear', 'sinc', 'spline'
                    }),
                    'identity_matrix': str,
                },
            },
        },
        'functional_registration': {
            'coregistration': {
                'run': bool,
                'func_input_prep': {
                    'input': [In({
                        'Mean_Functional', 'Selected_Functional_Volume'
                    })],
                    'Mean Functional': {
                        'n4_correct_func': bool
                    },
                    'Selected Functional Volume': {
                        'func_reg_input_volume': int
                    },
                },
                'boundary_based_registration': {
                    'using': Maybe(In(
                        valid_options['boundary_based_registration']['using']
                    )),
                    'run': forkable,
                    'bbr_schedule': str
                },
            },
            'EPI_registration': {
                'run': bool,
                'using': [In({'ANTS', 'FSL', 'FSL-linear'})],
                'EPI_template': str,
                'EPI_template_mask': Maybe(str),
                'ANTs': {
                    'parameters': Maybe(ANTs_parameters),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    }),
                },
                'FSL-FNIRT': {
                    'fnirt_config': Maybe(str),
                    'interpolation': In({'trilinear', 'sinc', 'spline'}),
                    'identity_matrix': str,
                },
            },
            'func_registration_to_template': {
                'run': bool,
                'output_resolution': {
                    'func_preproc_outputs': All(
                        str, Match(resolution_regex)),
                    'func_derivative_outputs': All(
                        str, Match(resolution_regex)
                    ),
                },
                'target_template': {
                    'using': [In({'T1_template', 'EPI_template'})],
                    'T1_template': {
                        'T1w_brain_template_funcreg': str,
                        'T1w_template_funcreg': Maybe(str),
                        'T1w_brain_template_mask_funcreg': Maybe(str),
                        'T1w_template_for_resample': Maybe(str),
                    },
                    'EPI_template': {
                        'EPI_template_funcreg': str,
                        'EPI_template_mask_funcreg': Maybe(str),
                        'EPI_template_for_resample': Maybe(str)
                    },
                },
                'ANTs_pipelines': {
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'})
                },
                'FNIRT_pipelines': {
                    'interpolation': In({'trilinear', 'sinc', 'spline'}),
                    'identity_matrix': str,
                },
            },
        },
    },
    'surface_analysis': {
        'run_freesurfer': bool,
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
        'convergence_threshold': Number,
    },
    'functional_preproc': {
        'run': bool,
        'truncation': {
            'start_tr': int,
            'stop_tr': Maybe(Any(int, 'End'))
        },
        'scaling': {
            'run': bool,
            'scaling_factor': Number
        },
        'despiking': {
            'run': forkable
        },
        'slice_timing_correction': {
            'run': forkable
        },
        'motion_estimates_and_correction': {
            'calculate_motion_first': bool,
            'motion_correction': {
                'using': [In({'3dvolreg', 'mcflirt'})],
                'AFNI-3dvolreg': {
                    'functional_volreg_twopass': bool,
                },
                'motion_correction_reference': [In({
                    'mean', 'median', 'selected volume'})],
                'motion_correction_reference_volume': int,
            },
            'motion_estimate_filter': Required(
                Any({  # no motion estimate filter
                    'run': Maybe(Any(
                        ExactSequence([False]), ExactSequence([]), False)),
                    'filter_type': Maybe(In({'notch', 'lowpass'})),
                    'filter_order': Maybe(int),
                    'breathing_rate_min': Maybe(Number),
                    'breathing_rate_max': Maybe(Number),
                    'center_frequency': Maybe(Number),
                    'filter_bandwidth': Maybe(Number),
                    'lowpass_cutoff': Maybe(Number),
                }, {  # notch filter with breathing_rate_* set
                    Required('run'): forkable,
                    Required('filter_type'): 'notch',
                    Required('filter_order'): int,
                    Required('breathing_rate_min'): Number,
                    'breathing_rate_max': Number,
                    'center_frequency': Maybe(Number),
                    'filter_bandwidth': Maybe(Number),
                    'lowpass_cutoff': Maybe(Number),
                }, {  # notch filter with manual parameters set
                    Required('run'): forkable,
                    Required('filter_type'): 'notch',
                    Required('filter_order'): int,
                    'breathing_rate_min': None,
                    'breathing_rate_max': None,
                    Required('center_frequency'): Number,
                    Required('filter_bandwidth'): Number,
                    'lowpass_cutoff': Maybe(Number),
                }, {  # lowpass filter with breathing_rate_min
                    Required('run'): forkable,
                    Required('filter_type'): 'lowpass',
                    Required('filter_order'): int,
                    Required('breathing_rate_min'): Number,
                    'breathing_rate_max': Maybe(Number),
                    'center_frequency': Maybe(Number),
                    'filter_bandwidth': Maybe(Number),
                    'lowpass_cutoff': Maybe(Number),
                }, {  # lowpass filter with lowpass_cutoff
                    Required('run'): forkable,
                    Required('filter_type'): 'lowpass',
                    Required('filter_order'): int,
                    Required('breathing_rate_min', default=None): None,
                    'breathing_rate_max': Maybe(Number),
                    'center_frequency': Maybe(Number),
                    'filter_bandwidth': Maybe(Number),
                    Required('lowpass_cutoff'): Number,
                },),
                msg='`motion_estimate_filter` configuration is invalid. See '
                    f'https://fcp-indi.github.io/docs/{_url_version}/user/'
                    'func#motion_estimate_filter_valid_options for details.\n',
            ),
        },
        'distortion_correction': {
            'run': forkable,
            'using': [In(['PhaseDiff', 'Blip'])],
            'PhaseDiff': {
                'fmap_skullstrip_option': In(['BET', 'AFNI']),
                'fmap_skullstrip_BET_frac': float,
                'fmap_skullstrip_AFNI_threshold': float,
            },
        },
        'func_masking': {
            'using': [In(
                ['AFNI', 'FSL', 'FSL_AFNI', 'Anatomical_Refined',
                 'Anatomical_Based']
            )],
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
    'nuisance_corrections': {
        '1-ICA-AROMA': {
            'run': forkable,
            'denoising_type': In({'aggr', 'nonaggr'}),
        },
        '2-nuisance_regression': {
            'run': forkable,
            'Regressors': Maybe([Schema({
                'Name': Required(str),
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
                'aCompCor': valid_options['Regressors']['CompCor'],
                'tCompCor': valid_options['Regressors']['CompCor'],
                'CerebrospinalFluid': valid_options[
                    'Regressors'
                ]['segmentation'],
                'WhiteMatter': valid_options[
                    'Regressors'
                ]['segmentation'],
                'GreyMatter': valid_options[
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
                    'gm_erosion_mm': Number,
                }
            },
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
            'T1w_brain_template_symmetric': str,
            'T1w_brain_template_symmetric_for_resample': str,
            'T1w_template_symmetric': str,
            'T1w_template_symmetric_for_resample': str,
            'dilated_symmetric_brain_mask': str,
            'dilated_symmetric_brain_mask_for_resample': str,
        },
    },
    'regional_homogeneity': {
        'run': bool,
        'cluster_size': In({7, 19, 27}),
    },
    'post_processing': {
        'spatial_smoothing': {
            'output': [In({'smoothed', 'nonsmoothed'})],
            'smoothing_method': [In({'FSL', 'AFNI'})],
            'fwhm': [int]
        },
        'z-scoring': {
            'output': [In({'z-scored', 'raw'})],
        },
    },
    'timeseries_extraction': {
        'run': bool,
        'tse_roi_paths': Optional(
            Maybe({
                str: In({', '.join([
                    option for option in options
                ]) for options in list(chain.from_iterable([list(
                    permutations(valid_options['timeseries']['roi_paths'],
                                 number_of)
                ) for number_of in range(1, 6)]))}),
            }),
            msg=permutation_message(
                'tse_roi_paths', valid_options['timeseries']['roi_paths'])
        ),
        'realignment': In({'ROI_to_func', 'func_to_ROI'}),
    },

    'seed_based_correlation_analysis': {
        'run': bool,
        'sca_roi_paths': Optional(
            Maybe({
                str: In({', '.join([
                    option for option in options
                ]) for options in list(chain.from_iterable([list(
                    permutations(valid_options['sca']['roi_paths'], number_of)
                ) for number_of in range(1, 4)]))})
            }),
            msg=permutation_message(
                'sca_roi_paths', valid_options['sca']['roi_paths'])
        ),
        'norm_timeseries_for_DR': bool,
    },
    'network_centrality': {
        'run': bool,
        'memory_allocation': Number,
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
        'run': bool,
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
