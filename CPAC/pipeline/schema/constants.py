'''Constants for configuration validation'''
from voluptuous import All, Any, In, Match, Range
from voluptuous.validators import Maybe
from CPAC import __version__

# 1 or more digits, optional decimal, 'e', optional '-', 1 or more digits
SCIENTIFIC_NOTATION_STR_REGEX = r'^([0-9]+(\.[0-9]*)*(e)-{0,1}[0-9]+)*$'

# (1 or more digits, optional decimal, 0 or more lowercase characters (units))
# ('x',
#  1 or more digits, optional decimal, 0 or more lowercase characters (units)
# ) 0 or more times
RESOLUTION_REGEX = r'^[0-9]+(\.[0-9]*){0,1}[a-z]*' \
                   r'(x[0-9]+(\.[0-9]*){0,1}[a-z]*)*$'

Number = Any(float, int, All(str, Match(SCIENTIFIC_NOTATION_STR_REGEX)))
forkable = Any(bool, [bool])
VALID_OPTIONS = {
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
        'using': {'AFNI', 'Nilearn'},
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
                int, float, 'Functional', All(str, Match(RESOLUTION_REGEX))
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
MUTEX = {  # mutually exclusive booleans
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
ANTs_PARAMETER_TRANSFORMS = {
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
    'use-histogram-matching': bool,
    'updateFieldVarianceInVoxelSpace': Number,
    'totalFieldVarianceInVoxelSpace': Number,
    'winsorize-image-intensities': {
        'lowerQuantile': float,
        'upperQuantile': float,
    },
}
ANTs_PARAMETERS = [Any(  # pylint: disable=invalid-name
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
            'Rigid': ANTs_PARAMETER_TRANSFORMS,
        }, {
            'Affine': ANTs_PARAMETER_TRANSFORMS,
        }, {
            'SyN': ANTs_PARAMETER_TRANSFORMS,
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
URL_VERSION = 'nightly' if __version__.endswith(
    '-dev') else f'v{__version__.lstrip("v")}'


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
