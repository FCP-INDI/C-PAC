# Copyright (C) 2022  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Validation schema for C-PAC pipeline configurations"""
# pylint: disable=too-many-lines
import re
from itertools import chain, permutations
from pathvalidate import sanitize_filename
from voluptuous import All, ALLOW_EXTRA, Any, Capitalize, Coerce, Equal, \
                       ExactSequence, ExclusiveInvalid, In, Length, Lower, \
                       Match, Maybe, Optional, Range, Required, Schema
from CPAC import docs_prefix
from CPAC.pipeline.random_state.seed import MAX_SEED
from CPAC.utils.datatypes import ListFromItem

# 1 or more digits, optional decimal, 'e', optional '-', 1 or more digits
scientific_notation_str_regex = r'^([0-9]+(\.[0-9]*)*(e)-{0,1}[0-9]+)*$'

# (1 or more digits, optional decimal, 0 or more lowercase characters (units))
# ('x',
#  1 or more digits, optional decimal, 0 or more lowercase characters (units)
# ) 0 or more times
resolution_regex = r'^[0-9]+(\.[0-9]*){0,1}[a-z]*' \
                   r'(x[0-9]+(\.[0-9]*){0,1}[a-z]*)*$'

Number = Any(float, int, All(str, Match(scientific_notation_str_regex)))
forkable = All(Coerce(ListFromItem), [bool], Length(max=2))
valid_options = {
    'acpc': {
        'target': ['brain', 'whole-head']
    },
    'brain_extraction': {
        'using': ['3dSkullStrip', 'BET', 'UNet', 'niworkflows-ants',
                  'FreeSurfer-BET-Tight', 'FreeSurfer-BET-Loose',
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
        'roi_paths': {'Avg', 'Voxel', 'SpatialReg'},
    },
    'connectivity_matrix': {
        'using': {'AFNI', 'Nilearn', 'ndmg'},
        'measure': {'Pearson', 'Partial', 'Spearman', 'MGC',
                    # 'TangentEmbed'  # "Skip tangent embedding for now"
                    },
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
                int, float, 'Functional', All(str, Match(resolution_regex))
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
            'functional_mean_thr': {
                'run': bool,
                'threshold_value': Maybe(int),
            },
            'functional_mean_bias_correction': bool,
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


def permutation_message(key, options):
    '''Function to give a clean, human-readable error message for keys that accept permutation values

    Parameters
    ----------
    key: str

    options: list or set

    Returns
    -------
    msg: str'''  # noqa: E501
    return f'''

\'{key}\' takes a dictionary with paths to region-of-interest (ROI)
 NIFTI files (.nii or .nii.gz) as keys and a comma separated string
 of analyses to run. For example, if you wish to run Avg and
 MultReg, you would enter:

    '/path/to/ROI.nii.gz': Avg, MultReg

Available analyses for \'{key}\' are {options}

'''


def sanitize(filename):
    '''Sanitize a filename and replace whitespaces with underscores'''
    return re.sub(r'\s+', '_', sanitize_filename(filename))


latest_schema = Schema({
    'FROM': Maybe(str),
    'pipeline_setup': {
        'pipeline_name': All(str, Length(min=1), sanitize),
        'output_directory': {
            'path': str,
            'source_outputs_dir': Maybe(str),
            'pull_source_once': bool,
            'write_func_outputs': bool,
            'write_debugging_outputs': bool,
            'output_tree': str,
            'quality_control': {
                'generate_quality_control_images': bool,
                'generate_xcpqc_files': bool,
            },
        },
        'working_directory': {
            'path': str,
            'remove_working_dir': bool,
        },
        'log_directory': {
            'run_logging': bool,
            'path': str,
            'graphviz': {
                'entire_workflow': {
                    'generate': bool,
                    'graph2use': Maybe(All(Coerce(ListFromItem),
                                           [All(Lower,
                                            In(('orig', 'hierarchical', 'flat',
                                                'exec', 'colored')))])),
                    'format': Maybe(All(Coerce(ListFromItem),
                                        [All(Lower, In(('png', 'svg')))])),
                    'simple_form': Maybe(bool)}}
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
            'raise_insufficient': bool,
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_OMP_threads': int,
            'num_participants_at_once': int,
            'random_seed': Maybe(Any(
                'random',
                All(int, Range(min=1, max=MAX_SEED)))),
            'observed_usage': {
                'callback_log': Maybe(str),
                'buffer': Number,
            },
        },
        'Amazon-AWS': {
            'aws_output_bucket_credentials': Maybe(str),
            's3_encryption': bool,
        },
        'Debugging': {
            'verbose': bool,
        },
    },
    'anatomical_preproc': {
        'run': bool,
        'run_t2': bool,
        'non_local_means_filtering': {
            'run': forkable,
            'noise_model': Maybe(str),
        },
        'n4_bias_field_correction': {
            'run': forkable,
            'shrink_factor': int,
        },
        't1t2_bias_field_correction': Required(
            # require 'T1w_brain_ACPC_template' if 'acpc_target' is 'brain'
            Any({
                'run': False,
                'BiasFieldSmoothingSigma': Maybe(int),
            }, {
                'run': True,
                'BiasFieldSmoothingSigma': Maybe(int),
            },),
        ),
        'acpc_alignment': Required(
            # require 'T1w_brain_ACPC_template' and
            # 'T2w_brain_ACPC_template' if 'acpc_target' is 'brain'
            Any({
                'run': False,
                'run_before_preproc': Maybe(bool),
                'brain_size': Maybe(int),
                'FOV_crop': Maybe(In({'robustfov', 'flirt'})),
                'acpc_target': Maybe(In(valid_options['acpc']['target'])),
                'align_brain_mask': Maybe(bool),
                'T1w_ACPC_template': Maybe(str),
                'T1w_brain_ACPC_template': Maybe(str),
                'T2w_ACPC_template': Maybe(str),
                'T2w_brain_ACPC_template': Maybe(str),
            }, {
                'run': True,
                'run_before_preproc': bool,
                'brain_size': int,
                'FOV_crop': In({'robustfov', 'flirt'}),
                'acpc_target': valid_options['acpc']['target'][1],
                'align_brain_mask': Maybe(bool),
                'T1w_ACPC_template': str,
                'T1w_brain_ACPC_template': Maybe(str),
                'T2w_ACPC_template': Maybe(str),
                'T2w_brain_ACPC_template': Maybe(str),
            }, {
                'run': True,
                'run_before_preproc': bool,
                'brain_size': int,
                'FOV_crop': In({'robustfov', 'flirt'}),
                'acpc_target': valid_options['acpc']['target'][0],
                'align_brain_mask': Maybe(bool),
                'T1w_ACPC_template': str,
                'T1w_brain_ACPC_template': str,
                'T2w_ACPC_template': Maybe(str),
                'T2w_brain_ACPC_template': Maybe(str),
            },),
            msg='\'brain\' requires \'T1w_brain_ACPC_template\' and '
                '\'T2w_brain_ACPC_template\' to '
                'be populated if \'run\' is not set to Off',
        ),
        'brain_extraction': {
            'run': bool,
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
                'unet_model': Maybe(str),
            },
            'niworkflows-ants': {
                'template_path': Maybe(str),
                'mask_path': Maybe(str),
                'regmask_path': Maybe(str),
            },
            'FreeSurfer-BET': {
                'T1w_brain_template_mask_ccs': Maybe(str)
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
                    'priors_path': Maybe(str),
                    'WM_path': Maybe(str),
                    'GM_path': Maybe(str),
                    'CSF_path': Maybe(str)
                },
            },
            'FreeSurfer': {
                'erode': Maybe(int),
                'CSF_label': Maybe([int]),
                'GM_label': Maybe([int]),
                'WM_label': Maybe([int]),
            },
            'ANTs_Prior_Based': {
                'run': forkable,
                'template_brain_list': Maybe(Any([str], [])),
                'template_segmentation_list': Maybe(Any([str], [])),
                'CSF_label': [int],
                'GM_label': [int],
                'WM_label': [int],
            },
            'Template_Based': {
                'run': forkable,
                'template_for_segmentation': [In(
                    valid_options['segmentation']['template']
                )],
                'WHITE': Maybe(str),
                'GRAY': Maybe(str),
                'CSF': Maybe(str),
            },
        },
    },
    'registration_workflows': {
        'quality_thresholds': {
            metric: Maybe(float) for
            metric in ('Dice', 'Jaccard', 'CrossCorr', 'Coverage')},
        'anatomical_registration': {
            'run': bool,
            'resolution_for_anat': All(str, Match(resolution_regex)),
            'T1w_brain_template': Maybe(str),
            'T1w_template': Maybe(str),
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
                    'ref_resolution': All(str, Match(resolution_regex)),
                    'FNIRT_T1w_brain_template': Maybe(str),
                    'FNIRT_T1w_template': Maybe(str),
                    'interpolation': In({
                        'trilinear', 'sinc', 'spline'
                    }),
                    'identity_matrix': Maybe(str),
                    'ref_mask': Maybe(str),
                    'ref_mask_res-2': Maybe(str),
                    'T1w_template_res-2': Maybe(str),
                },
            },
            'overwrite_transform': {
                'run': bool,
                'using': In({'FSL'}),
            },
        },
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
                        'fmriprep_reference'
                    })],
                    'Mean Functional': {
                        'n4_correct_func': bool
                    },
                    'Selected Functional Volume': {
                        'func_reg_input_volume': int
                    },
                },
                'boundary_based_registration': {
                    'run': All(Coerce(ListFromItem),
                               [Any(bool, All(Lower, Equal('fallback')))],
                               Length(max=3)),
                    'bbr_schedule': str,
                    'bbr_wm_map': In(('probability_map',
                                      'partial_volume_map')),
                    'bbr_wm_mask_args': str,
                    'reference': In(('whole-head', 'brain'))
                },
            },
            'EPI_registration': {
                'run': bool,
                'using': [In({'ANTS', 'FSL', 'FSL-linear'})],
                'EPI_template': Maybe(str),
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
                    'identity_matrix': Maybe(str),
                },
            },
            'func_registration_to_template': {
                'run': bool,
                'run_EPI': bool,
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
                        'T1w_brain_template_funcreg': Maybe(str),
                        'T1w_template_funcreg': Maybe(str),
                        'T1w_brain_template_mask_funcreg': Maybe(str),
                        'T1w_template_for_resample': Maybe(str),
                    },
                    'EPI_template': {
                        'EPI_template_funcreg': Maybe(str),
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
                    'identity_matrix': Maybe(str),
                },
                'apply_transform': {
                    'using': In({'default', 'abcd', 'dcan_nhp',
                                 'single_step_resampling_from_stc'}),
                },
            },
        },
    },
    'surface_analysis': {
        'freesurfer': {
            'run': bool,
            'reconall_args': Maybe(str),
            'freesurfer_dir': Maybe(str)
        },
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
            'smooth_fwhm': Maybe(int),
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
        'convergence_threshold': Number,
    },
    'functional_preproc': {
        'run': bool,
        'truncation': {
            'start_tr': int,
            'stop_tr': Maybe(Any(int, All(Capitalize, 'End')))
        },
        'scaling': {
            'run': bool,
            'scaling_factor': Number
        },
        'despiking': {
            'run': forkable
        },
        'slice_timing_correction': {
            'run': forkable,
            'tpattern': Maybe(str),
            'tzero': Maybe(int),
        },
        'motion_estimates_and_correction': {
            'run': bool,
            'motion_estimates': {
                'calculate_motion_first': bool,
                'calculate_motion_after': bool,
            },
            'motion_correction': {
                'using': [In({'3dvolreg', 'mcflirt'})],
                'AFNI-3dvolreg': {
                    'functional_volreg_twopass': bool,
                },
                'motion_correction_reference': [In({
                    'mean', 'median', 'selected_volume', 'fmriprep_reference'})],
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
                    f'{docs_prefix}/user/'
                    'func#motion_estimate_filter_valid_options for details.\n',
            ),
        },
        'distortion_correction': {
            'run': forkable,
            'using': [In(['PhaseDiff', 'Blip', 'Blip-FSL-TOPUP'])],
            'PhaseDiff': {
                'fmap_skullstrip_option': In(['BET', 'AFNI']),
                'fmap_skullstrip_BET_frac': float,
                'fmap_skullstrip_AFNI_threshold': float,
            },
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
                'regrid': int                
            }
        },
        'func_masking': {
            'using': [In(
                ['AFNI', 'FSL', 'FSL_AFNI', 'Anatomical_Refined',
                 'Anatomical_Based', 'Anatomical_Resampled',
                 'CCS_Anatomical_Refined']
            )],
            # handle validating mutually-exclusive booleans for FSL-BET
            # functional_mean_boolean must be True if one of the mutually-
            # exclusive options are
            # see mutex definition for more definition
            'FSL-BET': Maybe(Any(*(
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
                }]))
            ),
            'FSL_AFNI': {
                'bold_ref': Maybe(str),
                'brain_mask': Maybe(str),
                'brain_probseg': Maybe(str),
            },
            'Anatomical_Refined': {
                'anatomical_mask_dilation': Maybe(bool),
            },
            'apply_func_mask_in_native_space': bool,
        },
        'generate_func_mean': {
            'run': bool,
        },
        'normalize_func': {
            'run': bool,
        },
    },
    'nuisance_corrections': {
        '1-ICA-AROMA': {
            'run': forkable,
            'denoising_type': In({'aggr', 'nonaggr'}),
        },
        '2-nuisance_regression': {
            'run': forkable,
            'space': All(Coerce(ListFromItem),
                         [All(Lower, In({'native', 'template'}))]),
            'create_regressors': bool,
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
                    'brain_mask_erosion_prop': Maybe(Number),
                    'brain_mask_erosion_mm': Maybe(Number),
                    'brain_erosion_mm': Maybe(Number)
                },
                'erode_csf': {
                    'run': bool,
                    'csf_erosion_prop': Maybe(Number),
                    'csf_mask_erosion_mm': Maybe(Number),
                    'csf_erosion_mm': Maybe(Number),
                },
                'erode_wm': {
                    'run': bool,
                    'wm_erosion_prop': Maybe(Number),
                    'wm_mask_erosion_mm': Maybe(Number),
                    'wm_erosion_mm': Maybe(Number),
                },
                'erode_gm': {
                    'run': bool,
                    'gm_erosion_prop': Maybe(Number),
                    'gm_mask_erosion_mm': Maybe(Number),
                    'gm_erosion_mm': Maybe(Number),
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
            'T1w_brain_template_symmetric': Maybe(str),
            'T1w_brain_template_symmetric_for_resample': Maybe(str),
            'T1w_template_symmetric': Maybe(str),
            'T1w_template_symmetric_for_resample': Maybe(str),
            'dilated_symmetric_brain_mask': Maybe(str),
            'dilated_symmetric_brain_mask_for_resample': Maybe(str),
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
        Optional('roi_paths_fully_specified'): bool,
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
        'connectivity_matrix': {
            option: Maybe([In(valid_options['connectivity_matrix'][option])])
            for option in ['using', 'measure']
        },
    },
    'seed_based_correlation_analysis': {
        'run': bool,
        Optional('roi_paths_fully_specified'): bool,
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
        'template_specification_file': Maybe(str),
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
        'eye_scan_names': Maybe(Any([str], [])),
        'data_scan_names': Maybe(Any([str], [])),
        'eye_mask_path': Maybe(str),
        'stimulus_path': Maybe(str),
        'minimal_nuisance_correction': {
            'peer_gsr': bool,
            'peer_scrub': bool,
            'scrub_thresh': float,
        },
    },
})


def schema(config_dict):
    '''Validate a pipeline configuration against the latest validation schema
    by first applying backwards-compatibility patches, then applying
    Voluptuous validation, then handling complex configuration interaction
    checks before returning validated config_dict.

    Parameters
    ----------
    config_dict : dict

    Returns
    -------
    dict
    '''
    from CPAC.utils.utils import _changes_1_8_0_to_1_8_1
    partially_validated = latest_schema(_changes_1_8_0_to_1_8_1(config_dict))
    try:
        if (partially_validated['registration_workflows'][
            'functional_registration'
        ]['func_registration_to_template']['apply_transform'][
            'using'
        ] == 'single_step_resampling_from_stc' and partially_validated[
            'nuisance_corrections'
        ]['2-nuisance_regression']['space'] != ['template']):
            raise ExclusiveInvalid(
                '``single_step_resampling_from_stc`` requires template-space '
                'nuisance regression. Either set ``nuisance_corrections: '
                '2-nuisance_regression: space`` to ``template`` or choose a '
                'different option for ``registration_workflows: '
                'functional_registration: func_registration_to_template: '
                'apply_transform: using``')
    except KeyError:
        pass
    return partially_validated


schema.schema = latest_schema.schema
