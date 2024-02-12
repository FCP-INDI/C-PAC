# Copyright (C) 2012-2023  C-PAC Developers

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
import os
from pathlib import Path
from typing import Optional, Union
import nibabel as nib
from CPAC.pipeline.nipype_pipeline_engine import Node
from CPAC.pipeline.schema import valid_options
from CPAC.utils.docs import docstring_parameter
from CPAC.utils.interfaces.function import Function
from CPAC.utils.typing import ITERABLE, LIST


def convert_pvalue_to_r(datafile, p_value, two_tailed=False):
    '''
    Method to calculate correlation threshold from p_value

    Parameters
    ----------
    datafile : string
        filepath to dataset to extract number of time pts from
    p_value : float
        significance threshold p-value
    two_tailed : boolean (optional); default=False
        flag to indicate whether to calculate the two-tailed t-test
        threshold for the returned correlation value

    Returns
    -------
    r_value : float
        correlation threshold value
    '''

    import nibabel as nb
    import numpy as np
    import scipy.stats

    # Get two-tailed distribution
    if two_tailed:
        p_value = p_value / 2

    # Load in data and number of time pts
    img = nb.load(datafile).get_fdata()
    t_pts = img.shape[-1]

    # N-2 degrees of freedom with Pearson correlation (two sample means)
    deg_freedom = t_pts-2

    # Inverse Survival Function (Inverse of SF)
    # Note: survival function (SF) is also known as the complementary
    # cumulative distribution function (CCDF): F_(x) = p = P(X > x) = 1 - F(x)
    # The inverse will yield: x = F_^-1(p) = F_^-1(P(X > x))
    # where x is a value under the distribution of the random variable X
    # such that the probability of getting greater than x, is p
    t_value = scipy.stats.t.isf(p_value, deg_freedom)
    r_value = np.sqrt(t_value ** 2 / (deg_freedom + t_value ** 2))

    # Return correlation coefficient
    return r_value


def merge_lists(deg_list: Optional[LIST[str]] = None,
                eig_list: Optional[LIST[str]] = None,
                lfcd_list: Optional[LIST[str]] = None):
    '''Function to actually do the list merging.

    Parameters
    ----------
    deg_list : list of str, optional
        list of paths to degree centrality outputs

    eig_list : list of str, optional
        list of paths to eigenvector centrality outputs

    lfcd_list : list of str, optional
        list of paths to local functional connectivity density outputs

    Returns
    -------
    degree_weighted : str
        path to weighted degree centrality output

    degree_binarized : str
        path to binarized degree centrality output

    eigen_weighted : str
        path to weighted eigenvector centrality output

    eigen_binarized : str
        path to binarized eigenvector centrality output

    lfcd_weighted : str
        path to weighted local functional connectivity density output

    lfcd_binarized : str
        path to binarized local functional connectivity density output
    '''
    if deg_list is None:
        deg_list = []
    if eig_list is None:
        eig_list = []
    if lfcd_list is None:
        lfcd_list = []
    merged_list = []
    merged_list.extend(deg_list)
    merged_list.extend(eig_list)
    merged_list.extend(lfcd_list)

    # and we're hijacking this function to parse out the files instead!
    degree_weighted = None
    degree_binarized = None
    eigen_weighted = None
    eigen_binarized = None
    lfcd_weighted = None
    lfcd_binarized = None
    for path in merged_list:
        if 'degree' in path and 'Weighted' in path:
            degree_weighted = path
        elif 'degree' in path and 'Binarize' in path:
            degree_binarized = path
        elif 'eigen' in path and 'Weighted' in path:
            eigen_weighted = path
        elif 'eigen' in path and 'Binarize' in path:
            eigen_binarized = path
        elif 'local_functional' in path and 'Weighted' in path:
            lfcd_weighted = path
        elif 'local_functional' in path and 'Binarize' in path:
            lfcd_binarized = path

    return (degree_weighted, degree_binarized, eigen_weighted,
            eigen_binarized, lfcd_weighted, lfcd_binarized)


def create_merge_node(pipe_num: int) -> Node:
    '''Create a Function Node to merge lists for the centrality workflow

    Parameters
    ----------
    pipe_num : int

    Returns
    -------
    Node
        a Function Node to merge lists for the centrality workflow

    Notes
    -----
    Node Inputs::

        deg_list : list of strings
            list of paths to degree centrality outputs

        eig_list : list of strings
            list of paths to eigenvector centrality outputs

        lfcd_list : list of strings
            list of paths to local functional connectivity density outputs

    Node Outputs::

        degree_weighted : string
            path to weighted degree centrality output

        degree_binarized : string
            path to binarized degree centrality output

        eigen_weighted : string
            path to weighted eigenvector centrality output

        eigen_binarized : string
            path to binarized eigenvector centrality output

        lfcd_weighted : string
            path to weighted local functional connectivity density output

        lfcd_binarized : string
            path to binarized local functional connectivity density output
    '''
    return Node(Function(input_names=['deg_list', 'eig_list', 'lfcd_list'],
                         output_names=['degree_weighted',
                                       'degree_binarized',
                                       'eigen_weighted',
                                       'eigen_binarized',
                                       'lfcd_weighted',
                                       'lfcd_binarized'],
                         function=merge_lists, as_module=True),
                name=f'centrality_merge_node_{pipe_num}')


@Function.sig_imports(['from typing import Union', 'import os',
                       'from pathlib import Path', 'import nibabel as nib',
                       'from CPAC.pipeline.schema import valid_options',
                       'from CPAC.utils.docs import docstring_parameter',
                       'from CPAC.utils.typing import ITERABLE, LIST'])
@docstring_parameter(
    weight_options=tuple(valid_options['centrality']['weight_options']))
def sep_nifti_subbriks(nifti_file: Union[Path, str], out_names: ITERABLE[str]
                       ) -> LIST[str]:
    '''Separate sub-briks of niftis and save specified out

    Parameters
    ----------
    nifti_file : ~pathlib.Path or str
        path to NIfTI output of an AFNI centrality tool

    out_names : iterable of str
        an iterable of strings, each ending with one of {weight_options}

    Returns
    -------
    list of str
        paths to each of the specified outputs as its own file
    '''
    output_niftis = []
    weight_options = valid_options['centrality']['weight_options']
    selected_options = {_[::-1].split('_', 1)[0][::-1]: _ for _ in out_names}

    nii_img = nib.load(nifti_file)
    nii_arr = nii_img.get_fdata()
    nii_affine = nii_img.affine
    nii_dims = nii_arr.shape

    if nii_dims[-1] != len(weight_options):
        if len(nii_dims) == 3 and len(out_names) == 1:
            pass
        else:
            err_msg = 'out_names must have same number of elements as '\
                      'nifti sub-briks'
            raise Exception(err_msg)

    for brik, option in enumerate(weight_options):
        if option in selected_options:
            if len(nii_dims) == 3:
                brik_arr = nii_arr
            elif len(nii_dims) > 3:
                brik_arr = nii_arr[:, :, :, 0, brik]
            out_file = os.path.join(os.getcwd(),
                                    selected_options[option] + '.nii.gz')
            out_img = nib.Nifti1Image(brik_arr, nii_affine)
            out_img.to_filename(out_file)
            output_niftis.append(out_file)

    return output_niftis


@docstring_parameter(m_options=valid_options['centrality']['method_options'],
                     t_options=valid_options['centrality'][
                         'threshold_options'])
def check_centrality_params(method_option, threshold_option, threshold):
    '''
    Function to check the centrality parameters.

    Parameters
    ----------
    method_option : str or int
        one of {m_options} or index of option

    threshold_option : str
        one of {t_options} or index of option

    threshold: float

    Returns
    -------
    method_option : str
        one of {m_options}

    threshold_option : str
        one of {t_options}
    '''

    # Check method option
    if isinstance(method_option, int):
        if method_option < len(valid_options['centrality']['method_options']):
            method_option = valid_options[
                'centrality']['method_options'][method_option]
        else:
            raise MethodOptionError(method_option)
    elif not isinstance(method_option, str):
        raise TypeError('Method option must be a string, but type \'%s\' '
                        'provided' % type(method_option).__name__)

    # Check threshold option
    if type(threshold_option) is list:
        threshold_option = threshold_option[0]
    if type(threshold_option) is int:
        if threshold_option < len(
            valid_options['centrality']['threshold_options']
        ):
            threshold_option = valid_options[
                'centrality']['threshold_options'][threshold_option]
        else:
            raise ThresholdOptionError(threshold_option, method_option)
    elif type(threshold_option) is not str:
        raise TypeError('Threshold option must be a string, but type \'%s\' '
                        'provided' % type(threshold_option).__name__)

    # Format input strings
    method_option = method_option.lower().rstrip(' ')
    method_options_v1 = ['degree', 'eigenvector', 'lfcd']
    if method_option in method_options_v1:
        method_option = valid_options['centrality']['method_options'][
            method_options_v1.index(method_option)
        ]
    if ' ' not in threshold_option:
        threshold_option = ' '.join([threshold_option, 'threshold'])
    threshold_option = threshold_option.capitalize().rstrip(' ')

    # Check for strings properly formatted
    if method_option not in valid_options['centrality']['method_options']:
        raise MethodOptionError(method_option)

    # Check for strings properly formatted
    if threshold_option not in valid_options['centrality'][
        'threshold_options'
    ]:
        raise ThresholdOptionError(threshold_option, method_option)

    # Check for invalid combinations of method_option + threshold_option
    if (
        method_option == 'local_functional_connectivity_density' and
        threshold_option == 'Sparsity threshold'
    ):
        raise ThresholdOptionError(threshold_option, method_option)

    # If it's significance/sparsity thresholding, check for (0,1]
    if (
        threshold_option == 'Significance threshold' or
        threshold_option == 'Sparsity threshold'
    ):
        if threshold <= 0 or threshold > 1:
            raise ThresholdError(threshold_option, threshold)

    # If it's correlation, check for [-1,1]
    elif threshold_option == 'Correlation threshold':
        if threshold < -1 or threshold > 1:
            raise ThresholdError(threshold_option, threshold)
    else:
        raise ThresholdOptionError(threshold_option, method_option)

    # Return valid method and threshold options
    return method_option, threshold_option


class MethodOptionError(ValueError):
    """Raised when a selected centrality method option is not supported.
    """
    def __init__(self, method_option):
        self.method_option = method_option
        self.message = 'Method option \'%s\' not supported' % method_option
        super().__init__(self.message)


class ThresholdError(ValueError):
    """Raised when a selected threshold value is not supported for a
    selected threshold option.
    """
    def __init__(self, threshold_option, threshold):
        self.threshold_option = threshold_option
        self.threshold = threshold
        print(type(threshold))
        self.message = f'For \'{threshold_option}\', threshold value must be '
        if (
            threshold_option == 'Significance threshold' or
            threshold_option == 'Sparsity threshold'
        ):
            self.message += 'a positive number greater than 0 '
        elif threshold_option == 'Correlation threshold':
            self.message += 'greater than or equal to -1 '
        else:
            raise ThresholdOptionError(threshold_option)
        self.message += 'and less than or equal to 1.\n Currently it is set ' \
                        f'at {threshold}'
        super().__init__(self.message)


class ThresholdOptionError(ValueError):
    """Raised when a selected threshold option is not supported for a
    selected centrality measure.
    """
    def __init__(self, threshold_option, method_option=None):
        self.method_option = method_option
        self.threshold_option = threshold_option
        self.message = f'Threshold option \'{threshold_option}\' not supported'
        if self.method_option:
            self.message += ' for network centrality measure ' \
                            f'\'{method_option}\''
        self.message += '; fix this in the pipeline config'
        if (
            method_option == 'local_functional_connectivity_density' and
            threshold_option == 'Sparsity threshold'
        ):
            _valid_options = ' or '.join([
                f"'{t}'" for t in valid_options[
                    'centrality'
                ]['threshold_options'] if t != threshold_option
            ])
            self.message += f'. \'{method_option}\' must use {_valid_options}.'
        super().__init__(self.message)
