# Copyright (C) 2022-2024  C-PAC Developers

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
import numpy as np
from pathvalidate import sanitize_filename
from subprocess import CalledProcessError
from voluptuous import All, ALLOW_EXTRA, Any, BooleanInvalid, Capitalize, \
                       Coerce, CoerceInvalid, ExclusiveInvalid, In, Length, \
                       LengthInvalid, Lower, Match, Maybe, MultipleInvalid, \
                       Optional, Range, Required, Schema, Title
from CPAC.utils.datatypes import ItemFromList, ListFromItem
from CPAC.utils.docs import DOCS_URL_PREFIX
from CPAC.utils.utils import YAML_BOOLS

# 1 or more digits, optional decimal, 'e', optional '-', 1 or more digits
SCIENTIFIC_NOTATION_STR_REGEX = r'^([0-9]+(\.[0-9]*)*(e)-{0,1}[0-9]+)*$'

# (1 or more digits, optional decimal, 0 or more lowercase characters (units))
# ('x',
#  1 or more digits, optional decimal, 0 or more lowercase characters (units)
# ) 0 or more times
RESOLUTION_REGEX = r'^[0-9]+(\.[0-9]*){0,1}[a-z]*' \
                   r'(x[0-9]+(\.[0-9]*){0,1}[a-z]*)*$'

Number = Any(float, int, All(str, Match(SCIENTIFIC_NOTATION_STR_REGEX)))


def str_to_bool1_1(x):  # pylint: disable=invalid-name
    '''Convert strings to Booleans for YAML1.1 syntax

    Ref https://yaml.org/type/bool.html

    Parameters
    ----------
    x : any

    Returns
    -------
    bool
    '''
    if isinstance(x, str):
        try:
            x = float(x)
            if x == 0:
                return False
        except ValueError:
            pass
        x = (True if str(x).lower() in YAML_BOOLS[True] else
             False if str(x).lower() in YAML_BOOLS[False] else x)
    if not isinstance(x, (bool, int)):
        raise BooleanInvalid('Type boolean value was expected, type '
                             f'{getattr(type(x), "__name__", str(type(x)))} '
                             f'value\n\n{x}\n\nwas provided')
    return bool(x)


bool1_1 = All(str_to_bool1_1, bool)
forkable = All(Coerce(ListFromItem), [bool1_1], Length(max=2))
valid_options = {
    'acpc': {
        'target': ['brain', 'whole-head']
    },
    'brain_extraction': {
        'using': ['3dSkullStrip', 'BET', 'UNet', 'niworkflows-ants',
                  'FreeSurfer-BET-Tight', 'FreeSurfer-BET-Loose',
                  'FreeSurfer-ABCD', 'FreeSurfer-Brainmask']
    },
    'centrality': {
       'method_options': ['degree_centrality', 'eigenvector_centrality',
                          'local_functional_connectivity_density'],
       'threshold_options': ['Significance threshold', 'Sparsity threshold',
                             'Correlation threshold'],
       'weight_options': ['Binarized', 'Weighted']
    },
    'motion_correction': ['3dvolreg', 'mcflirt'],
    'sca': {
        'roi_paths': ['Avg', 'DualReg', 'MultReg'],
    },
    'segmentation': {
        'using': ['FSL-FAST', 'ANTs_Prior_Based', 'Template_Based'],
        'template': ['EPI_Template', 'T1_Template'],
    },
    'timeseries': {
        'roi_paths': ['Avg', 'Voxel', 'SpatialReg'],
    },
    'connectivity_matrix': {
        'using': ['AFNI', 'Nilearn', 'ndmg'],
        'measure': ['Pearson', 'Partial', 'Spearman', 'MGC',
                    # 'TangentEmbed'  # "Skip tangent embedding for now"
        ],
    },
    'Regressors': {
        'CompCor': {
            'degree': int,
            'erode_mask_mm': bool1_1,
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
            'erode_mask': bool1_1,
            'extraction_resolution': Any(
                int, float, 'Functional', All(str, Match(RESOLUTION_REGEX))
            ),
            'include_delayed': bool1_1,
            'include_delayed_squared': bool1_1,
            'include_squared': bool1_1,
            'summary': Any(
                str, {'components': int, 'method': str}
            ),
        },
    },
    'target_space': ['Native', 'Template']
}
valid_options['space'] = list({option.lower() for option in
                               valid_options['target_space']})
mutex = {  # mutually exclusive booleans
    'FSL-BET': {
        # exactly zero or one of each of the following can be True for FSL-BET
        'mutex': ['reduce_bias', 'robust', 'padding', 'remove_eyes',
                  'surfaces'],
        # the remaining keys: validators for FSL-BET
        'rem': {
            'frac': float,
            'mesh_boolean': bool1_1,
            'outline': bool1_1,
            'radius': int,
            'skull': bool1_1,
            'threshold': bool1_1,
            'vertical_gradient': Range(min=-1, max=1, min_included=False,
                                       max_included=False),
            'functional_mean_thr': {
                'run': bool1_1,
                'threshold_value': Maybe(int),
            },
            'functional_mean_bias_correction': bool1_1,
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
        'iteration': All(str, Match(RESOLUTION_REGEX)),
        'convergenceThreshold': Number,
        'convergenceWindowSize': int,
    },
    'smoothing-sigmas': All(str, Match(RESOLUTION_REGEX)),
    'shrink-factors': All(str, Match(RESOLUTION_REGEX)),
    'use-histogram-matching': bool1_1,
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
        'verbose': Any(Coerce(int), In({0, 1})),
    }, {
        'float': Any(Coerce(int), In({0, 1})),
    }, {
        'masks': {
            'fixed_image_mask': bool1_1,
            'moving_image_mask': bool1_1,
        },
    }, dict  # TODO: specify other valid ANTs parameters
)]
motion_estimate_filter = Any({  # notch filter with breathing_rate_* set
        Required('filter_type'): 'notch',
        Required('filter_order'): int,
        Required('breathing_rate_min'): Number,
        'breathing_rate_max': Number,
        'center_frequency': Maybe(Number),
        'filter_bandwidth': Maybe(Number),
        'lowpass_cutoff': Maybe(Number),
        'Name': Maybe(str)
    }, {  # notch filter with manual parameters set
        Required('filter_type'): 'notch',
        Required('filter_order'): int,
        'breathing_rate_min': None,
        'breathing_rate_max': None,
        Required('center_frequency'): Number,
        Required('filter_bandwidth'): Number,
        'lowpass_cutoff': Maybe(Number),
        'Name': Maybe(str)
    }, {  # lowpass filter with breathing_rate_min
        Required('filter_type'): 'lowpass',
        Required('filter_order'): int,
        Required('breathing_rate_min'): Number,
        'breathing_rate_max': Maybe(Number),
        'center_frequency': Maybe(Number),
        'filter_bandwidth': Maybe(Number),
        'lowpass_cutoff': Maybe(Number),
        'Name': Maybe(str)
    }, {  # lowpass filter with lowpass_cutoff
        Required('filter_type'): 'lowpass',
        Required('filter_order'): int,
        Required('breathing_rate_min', default=None): None,
        'breathing_rate_max': Maybe(Number),
        'center_frequency': Maybe(Number),
        'filter_bandwidth': Maybe(Number),
        Required('lowpass_cutoff'): Number,
        'Name': Maybe(str)},
    msg='`motion_estimate_filter` configuration is invalid.\nSee '
        f'{DOCS_URL_PREFIX}/user/'
        'func#motion-estimate-filter-valid-options for details.\n')
target_space = All(Coerce(ListFromItem),
                   [All(Title, In(valid_options['target_space']))])


def name_motion_filter(mfilter, mfilters=None):
    '''Given a motion filter, create a short string for the filename

    Parameters
    ----------
    mfilter : dict

    mfliters : list or None

    Returns
    -------
    str

    Examples
    --------
    >>> name_motion_filter({'filter_type': 'notch', 'filter_order': 2,
    ...     'center_frequency': 0.31, 'filter_bandwidth': 0.12})
    'notch2fc0p31bw0p12'
    >>> name_motion_filter({'filter_type': 'notch', 'filter_order': 4,
    ...     'breathing_rate_min': 0.19, 'breathing_rate_max': 0.43})
    'notch4fl0p19fu0p43'
    >>> name_motion_filter({'filter_type': 'lowpass', 'filter_order': 4,
    ...     'lowpass_cutoff': .0032})
    'lowpass4fc0p0032'
    >>> name_motion_filter({'filter_type': 'lowpass', 'filter_order': 2,
    ...     'breathing_rate_min': 0.19})
    'lowpass2fl0p19'
    >>> name_motion_filter({'filter_type': 'lowpass', 'filter_order': 2,
    ...     'breathing_rate_min': 0.19}, [{'Name': 'lowpass2fl0p19'}])
    'lowpass2fl0p19dup1'
    >>> name_motion_filter({'filter_type': 'lowpass', 'filter_order': 2,
    ...     'breathing_rate_min': 0.19}, [{'Name': 'lowpass2fl0p19'},
    ...     {'Name': 'lowpass2fl0p19dup1'}])
    'lowpass2fl0p19dup2'
    '''
    if mfilters is None:
        mfilters = []
    if 'Name' in mfilter:
        name = mfilter['Name']
    else:
        if mfilter['filter_type'] == 'notch':
            if mfilter.get('breathing_rate_min'):
                range_str = (f'fl{mfilter["breathing_rate_min"]}'
                             f'fu{mfilter["breathing_rate_max"]}')
            else:
                range_str = (f'fc{mfilter["center_frequency"]}'
                             f'bw{mfilter["filter_bandwidth"]}')
        else:
            if mfilter.get('breathing_rate_min'):
                range_str = f'fl{mfilter["breathing_rate_min"]}'
            else:
                range_str = f'fc{mfilter["lowpass_cutoff"]}'
        range_str = range_str.replace('.', 'p')
        name = f'{mfilter["filter_type"]}{mfilter["filter_order"]}{range_str}'
    dupes = 'Name' not in mfilter and len([_ for _ in (_.get('Name', '') for
                                           _ in mfilters) if
                                           _.startswith(name)])
    if dupes:
        dup = re.search('(?=[A-Za-z0-9]*)(dup[0-9]*)', name)
        if dup:  # Don't chain 'dup' suffixes
            name = name.replace(dup.group(), f'dup{dupes}')
        else:
            name = f'{name}dup{dupes}'
    return name


def permutation_message(key, options):
    '''Function to give a clean, human-readable error message for keys
    that accept permutation values

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
    'skip env check': Maybe(bool),  # flag for skipping an environment check
    'pipeline_setup': {
        'pipeline_name': All(str, Length(min=1), sanitize),
        'output_directory': {
            'path': str,
            'source_outputs_dir': Maybe(str),
            'pull_source_once': bool1_1,
            'write_func_outputs': bool1_1,
            'write_debugging_outputs': bool1_1,
            'output_tree': str,
            'quality_control': {
                'generate_quality_control_images': bool1_1,
                'generate_xcpqc_files': bool1_1,
            },
            'user_defined': Maybe(str),
        },
        'working_directory': {
            'path': str,
            'remove_working_dir': bool1_1,
        },
        'log_directory': {
            'run_logging': bool1_1,
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
                    'simple_form': Maybe(bool)}},
        },
        'crash_log_directory': {
            'path': Maybe(str),
        },
        'system_config': {
            'fail_fast': bool1_1,
            'FSLDIR': Maybe(str),
            'on_grid': {
                'run': bool1_1,
                'resource_manager': Maybe(str),
                'SGE': {
                    'parallel_environment': Maybe(str),
                    'queue': Maybe(str),
                },
            },
            'maximum_memory_per_participant': Number,
            'raise_insufficient': bool1_1,
            'max_cores_per_participant': int,
            'num_ants_threads': int,
            'num_OMP_threads': int,
            'num_participants_at_once': int,
            'random_seed': Maybe(Any(
                'random',
                All(int, Range(min=1, max=np.iinfo(np.int32).max)))),
            'observed_usage': {
                'callback_log': Maybe(str),
                'buffer': Number,
            },
        },
        'Amazon-AWS': {
            'aws_output_bucket_credentials': Maybe(str),
            's3_encryption': bool1_1,
        },
        'Debugging': {
            'verbose': bool1_1,
        },
        'outdir_ingress': {
            'run': bool1_1,
            'Template': Maybe(str),
        },
    },
    'anatomical_preproc': {
        'run': bool1_1,
        'run_t2': bool1_1,
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
                'run_before_preproc': Maybe(bool1_1),
                'brain_size': Maybe(int),
                'FOV_crop': Maybe(In({'robustfov', 'flirt'})),
                'acpc_target': Maybe(In(valid_options['acpc']['target'])),
                'align_brain_mask': Maybe(bool1_1),
                'T1w_ACPC_template': Maybe(str),
                'T1w_brain_ACPC_template': Maybe(str),
                'T2w_ACPC_template': Maybe(str),
                'T2w_brain_ACPC_template': Maybe(str),
            }, {
                'run': True,
                'run_before_preproc': bool1_1,
                'brain_size': int,
                'FOV_crop': In({'robustfov', 'flirt'}),
                'acpc_target': valid_options['acpc']['target'][1],
                'align_brain_mask': Maybe(bool1_1),
                'T1w_ACPC_template': str,
                'T1w_brain_ACPC_template': Maybe(str),
                'T2w_ACPC_template': Maybe(str),
                'T2w_brain_ACPC_template': Maybe(str),
            }, {
                'run': True,
                'run_before_preproc': bool1_1,
                'brain_size': int,
                'FOV_crop': In({'robustfov', 'flirt'}),
                'acpc_target': valid_options['acpc']['target'][0],
                'align_brain_mask': Maybe(bool1_1),
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
            'run': bool1_1,
            'using': [In(valid_options['brain_extraction']['using'])],
            'AFNI-3dSkullStrip': {
                'mask_vol': bool1_1,
                'shrink_factor': Number,
                'var_shrink_fac': bool1_1,
                'shrink_factor_bot_lim': Number,
                'avoid_vent': bool1_1,
                'n_iterations': int,
                'pushout': bool1_1,
                'touchup': bool1_1,
                'fill_hole': int,
                'NN_smooth': int,
                'smooth_final': int,
                'avoid_eyes': bool1_1,
                'use_edge': bool1_1,
                'exp_frac': Number,
                'push_to_edge': bool1_1,
                'use_skull': bool1_1,
                'perc_int': Number,
                'max_inter_iter': int,
                'fac': Number,
                'blur_fwhm': Number,
                'monkey': bool1_1,
            },
            'FSL-FNIRT': {
                'interpolation': In({
                    'trilinear', 'sinc', 'spline'
                }),
            },
            'FSL-BET': {
                'frac': Number,
                'Robustfov': bool1_1,
                'mesh_boolean': bool1_1,
                'outline': bool1_1,
                'padding': bool1_1,
                'radius': int,
                'reduce_bias': bool1_1,
                'remove_eyes': bool1_1,
                'robust': bool1_1,
                'skull': bool1_1,
                'surfaces': bool1_1,
                'threshold': bool1_1,
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
        'run': bool1_1,
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
                    'run': bool1_1,
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
        'anatomical_registration': {
            'run': bool1_1,
            'resolution_for_anat': All(str, Match(RESOLUTION_REGEX)),
            'T1w_brain_template': Maybe(str),
            'T1w_template': Maybe(str),
            'T1w_brain_template_mask': Maybe(str),
            'reg_with_skull': bool1_1,
            'registration': {
                'using': [In({'ANTS', 'FSL', 'FSL-linear'})],
                'ANTs': {
                    'use_lesion_mask': bool1_1,
                    'T1_registration': Maybe(ANTs_parameters),
                    'interpolation': In({
                        'Linear', 'BSpline', 'LanczosWindowedSinc'
                    }),
                },
                'FSL-FNIRT': {
                    'fnirt_config': Maybe(str),
                    'ref_resolution': All(str, Match(RESOLUTION_REGEX)),
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
                'run': bool1_1,
                'using': In({'FSL'}),
            },
        },
        'functional_registration': {
            'coregistration': {
                'run': bool1_1,
                'reference': In({'brain', 'restore-brain'}),
                'interpolation': In({'trilinear', 'sinc', 'spline'}),
                'using': str,
                'input': str,
                'cost': str,
                'dof': int,
                'arguments': Maybe(str),
                'func_input_prep': {
                    'reg_with_skull': bool1_1,
                    'input': [In({
                        'Mean_Functional', 'Selected_Functional_Volume',
                        'fmriprep_reference'
                    })],
                    'Mean Functional': {
                        'n4_correct_func': bool1_1
                    },
                    'Selected Functional Volume': {
                        'func_reg_input_volume': int
                    },
                },
                'boundary_based_registration': {
                    'run': forkable,
                    'bbr_schedule': str,
                    'bbr_wm_map': In({'probability_map', 'partial_volume_map'}),
                    'bbr_wm_mask_args': str,
                    'reference': In({'whole-head', 'brain'})
                },
            },
            'EPI_registration': {
                'run': bool1_1,
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
                'run': bool1_1,
                'run_EPI': bool1_1,
                'output_resolution': {
                    'func_preproc_outputs': All(
                        str, Match(RESOLUTION_REGEX)),
                    'func_derivative_outputs': All(
                        str, Match(RESOLUTION_REGEX)
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
        'abcd_prefreesurfer_prep':{
            'run': bool1_1,
        },
        'freesurfer': {
            'run_reconall': bool1_1,
            'reconall_args': Maybe(str),
            # 'generate_masks': bool1_1,
            'ingress_reconall': bool1_1,
        },
        'post_freesurfer': {
            'run': bool1_1,
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
        'amplitude_low_frequency_fluctuation': {
            'run': bool1_1,
        },
        'regional_homogeneity': {
            'run': bool1_1,
        },
        'surface_connectivity': {
            'run': bool1_1,
            'surface_parcellation_template': Maybe(str),
        },
    },
    'longitudinal_template_generation': {
        'run': bool1_1,
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
        'run': bool1_1,
        'truncation': {
            'start_tr': int,
            'stop_tr': Maybe(Any(int, All(Capitalize, 'End')))
        },
        'update_header': {
            'run': bool1_1,
        },
        'scaling': {
            'run': bool1_1,
            'scaling_factor': Number
        },
        'despiking': {
            'run': forkable,
            'space': In({'native', 'template'})
        },
        'slice_timing_correction': {
            'run': forkable,
            'tpattern': Maybe(str),
            'tzero': Maybe(int),
        },
        'motion_estimates_and_correction': {
            'run': bool1_1,
            'motion_estimates': {
                'calculate_motion_first': bool1_1,
                'calculate_motion_after': bool1_1,
            },
            'motion_correction': {
                'using': Optional(All(Coerce(ListFromItem),
                                      Length(min=0, max=1,
                    msg='Forking is currently broken for this option. '
                        'Please use separate configs if you want to '
                        'use each of 3dvolreg and mcflirt. Follow '
                        'https://github.com/FCP-INDI/C-PAC/issues/1935 '
                        'to see when this issue is resolved.'),
                             [In(valid_options['motion_correction'])])),
                'AFNI-3dvolreg': {
                    'functional_volreg_twopass': bool1_1,
                },
                'motion_correction_reference': [In({
                    'mean', 'median', 'selected_volume',
                    'fmriprep_reference'})],
                'motion_correction_reference_volume': int,
            },
            'motion_estimate_filter': Required(
                Any({'run': forkable,
                     'filters': [motion_estimate_filter]},
                    {'run': All(forkable, [In([False], [])]),
                     'filters': Maybe(list)})
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
            'run': bool1_1,
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
                    'functional_mean_boolean': bool1_1,
                    **{k: False for k in mutex['FSL-BET']['mutex']}
                }]))
            ),
            'FSL_AFNI': {
                'bold_ref': Maybe(str),
                'brain_mask': Maybe(str),
                'brain_probseg': Maybe(str),
            },
            'Anatomical_Refined': {
                'anatomical_mask_dilation': Maybe(bool1_1),
            },
            'apply_func_mask_in_native_space': bool1_1,
        },
        'generate_func_mean': {
            'run': bool1_1,
        },
        'normalize_func': {
            'run': bool1_1,
        },
        'coreg_prep': {
            'run': bool1_1,
        },
    },
    'nuisance_corrections': {
        '1-ICA-AROMA': {
            'run': forkable,
            'denoising_type': In({'aggr', 'nonaggr'}),
        },
        '2-nuisance_regression': {
            'run': forkable,
            'space': All(Coerce(ItemFromList),
                         Lower, In({'native', 'template'})),
            'create_regressors': bool1_1,
            'ingress_regressors': {
                'run': bool1_1,
                'Regressors': {
                    'Name': Maybe(str),
                    'Columns': [str]},
            },
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
                    'include_delayed': bool1_1,
                    'include_squared': bool1_1,
                    'include_delayed_squared': bool1_1
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
                    'run': bool1_1,
                    'brain_mask_erosion_prop': Maybe(Number),
                    'brain_mask_erosion_mm': Maybe(Number),
                    'brain_erosion_mm': Maybe(Number)
                },
                'erode_csf': {
                    'run': bool1_1,
                    'csf_erosion_prop': Maybe(Number),
                    'csf_mask_erosion_mm': Maybe(Number),
                    'csf_erosion_mm': Maybe(Number),
                },
                'erode_wm': {
                    'run': bool1_1,
                    'wm_erosion_prop': Maybe(Number),
                    'wm_mask_erosion_mm': Maybe(Number),
                    'wm_erosion_mm': Maybe(Number),
                },
                'erode_gm': {
                    'run': bool1_1,
                    'gm_erosion_prop': Maybe(Number),
                    'gm_mask_erosion_mm': Maybe(Number),
                    'gm_erosion_mm': Maybe(Number),
                }
            },
        },
    },
    'amplitude_low_frequency_fluctuation': {
        'run': bool1_1,
        'target_space': target_space,
        'highpass_cutoff': [float],
        'lowpass_cutoff': [float],
    },
    'voxel_mirrored_homotopic_connectivity': {
        'run': bool1_1,
        'symmetric_registration': {
            'T1w_brain_template_symmetric': Maybe(str),
            'T1w_brain_template_symmetric_funcreg': Maybe(str),
            'T1w_brain_template_symmetric_for_resample': Maybe(str),
            'T1w_template_symmetric': Maybe(str),
            'T1w_template_symmetric_funcreg': Maybe(str),
            'T1w_template_symmetric_for_resample': Maybe(str),
            'dilated_symmetric_brain_mask': Maybe(str),
            'dilated_symmetric_brain_mask_for_resample': Maybe(str),
        },
    },
    'regional_homogeneity': {
        'run': bool1_1,
        'target_space': target_space,
        'cluster_size': In({7, 19, 27}),
    },
    'post_processing': {
        'spatial_smoothing': {
            'run': bool1_1,
            'output': [In({'smoothed', 'nonsmoothed'})],
            'smoothing_method': [In({'FSL', 'AFNI'})],
            'fwhm': [int]
        },
        'z-scoring': {
            'run': bool1_1,
            'output': [In({'z-scored', 'raw'})],
        },
    },
    'timeseries_extraction': {
        'run': bool1_1,
        Optional('roi_paths_fully_specified'): bool1_1,
        'tse_roi_paths': Optional(
            Maybe({
                str: In({', '.join(
                    list(options)
                ) for options in list(chain.from_iterable([list(
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
        'run': bool1_1,
        Optional('roi_paths_fully_specified'): bool1_1,
        'sca_roi_paths': Optional(
            Maybe({
                str: In({', '.join(list(
                    options
                )) for options in list(chain.from_iterable([list(
                    permutations(valid_options['sca']['roi_paths'], number_of)
                ) for number_of in range(1, 4)]))})
            }),
            msg=permutation_message(
                'sca_roi_paths', valid_options['sca']['roi_paths'])
        ),
        'norm_timeseries_for_DR': bool1_1,
    },
    'network_centrality': {
        'run': bool1_1,
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
        'run': bool1_1,
        'eye_scan_names': Maybe(Any([str], [])),
        'data_scan_names': Maybe(Any([str], [])),
        'eye_mask_path': Maybe(str),
        'stimulus_path': Maybe(str),
        'minimal_nuisance_correction': {
            'peer_gsr': bool1_1,
            'peer_scrub': bool1_1,
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
    try:
        partially_validated = latest_schema(
            _changes_1_8_0_to_1_8_1(config_dict))
    except MultipleInvalid as multiple_invalid:
        if (multiple_invalid.path == ['nuisance_corrections',
                                      '2-nuisance_regression', 'space'] and
                isinstance(multiple_invalid.errors[0], CoerceInvalid)):
            raise CoerceInvalid(
                'Nusiance regression space is not forkable. Please choose '
                f'only one of {valid_options["space"]}',
                path=multiple_invalid.path) from multiple_invalid
        raise multiple_invalid
    try:
        if (partially_validated['registration_workflows'][
            'functional_registration'
        ]['func_registration_to_template']['apply_transform'][
            'using'
        ] == 'single_step_resampling_from_stc'):
            or_else = ('or choose a different option for '
                       '``registration_workflows: functional_registration: '
                       'func_registration_to_template: apply_transform: '
                       'using``')
            if True in partially_validated['nuisance_corrections'][
                '2-nuisance_regression']['run'] and partially_validated[
                'nuisance_corrections'
            ]['2-nuisance_regression']['space'] != 'template':
                raise ExclusiveInvalid(
                    '``single_step_resampling_from_stc`` requires '
                    'template-space nuisance regression. Either set '
                    '``nuisance_corrections: 2-nuisance_regression: space`` '
                    f'to ``template`` {or_else}')
            if any(registration != 'ANTS' for registration in
                   partially_validated['registration_workflows'][
                       'anatomical_registration']['registration']['using']):
                raise ExclusiveInvalid(
                    '``single_step_resampling_from_stc`` requires '
                    'ANTS registration. Either set '
                    '``registration_workflows: anatomical_registration: '
                    f'registration: using`` to ``ANTS`` {or_else}')
    except KeyError:
        pass
    try:
        motion_filters = partially_validated['functional_preproc'][
            'motion_estimates_and_correction']['motion_estimate_filter']
        if True in motion_filters['run']:
            for motion_filter in motion_filters['filters']:
                motion_filter['Name'] = name_motion_filter(
                    motion_filter, motion_filters['filters'])
        else:
            motion_filters['filters'] = []
    except KeyError:
        pass
    try:
        # 'motion_correction.using' is only optional if 'run' is Off
        mec = partially_validated['functional_preproc'][
            'motion_estimates_and_correction']
        if mec['run']:
            try:
                # max should be len(valid_options['motion_correction'])
                # once #1935 is resolved
                Length(min=1, max=1)(mec['motion_correction']['using'])
            except LengthInvalid:
                mec_path = ['functional_preproc',
                            'motion_estimates_and_correction']
                raise LengthInvalid(  # pylint: disable=raise-missing-from
                    f'If data[{"][".join(map(repr, mec_path))}][\'run\'] is '
                    # length must be between 1 and
                    # len(valid_options['motion_correction']) once #1935 is
                    # resolved
                    'True, length of list must be exactly 1',
                    path=[*mec_path, 'motion_correction', 'using'])
    except KeyError:
        pass
    try:
        # Check for mutually exclusive options
        if (partially_validated['nuisance_corrections'][
            '2-nuisance_regression']['ingress_regressors']['run'] and
            partially_validated['nuisance_corrections'][
                '2-nuisance_regression']['create_regressors']):
            raise ExclusiveInvalid(
                "[!] Ingress_regressors and create_regressors can't both run! "
                " Try turning one option off.\n ")
    except KeyError:
        pass
    try:
        if not partially_validated.get("skip env check"
                                       ) and 'unet' in [using.lower() for using in
                partially_validated['anatomical_preproc'][
                    'brain_extraction']['using']]:
            try:
                from importlib import import_module
                import_module('CPAC.unet')
            except (CalledProcessError, ImportError, ModuleNotFoundError, OSError) as error:
                import site
                raise OSError(
                    'U-Net brain extraction requires torch to be installed, '
                    'but the installation path in this container is '
                    'read-only. Please bind a local writable path to '
                    f'"{site.USER_BASE}" in the container to use U-Net.'
                ) from error
    except KeyError:
        pass
    return partially_validated


schema.schema = latest_schema.schema
