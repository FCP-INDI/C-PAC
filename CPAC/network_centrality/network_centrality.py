# Copyright (C) 2015-2023  C-PAC Developers

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
from pathlib import Path
from typing import Optional, Union
from nipype.interfaces.afni.preprocess import DegreeCentrality, LFCD
from nipype.pipeline.engine import Workflow
from CPAC.network_centrality.utils import ThresholdOptionError
from CPAC.pipeline.schema import valid_options
from CPAC.utils.docs import docstring_parameter
from CPAC.utils.interfaces.afni import AFNI_GTE_21_1_1, ECM
from CPAC.utils.typing import LIST


@docstring_parameter(m_options=valid_options['centrality']['method_options'],
                     t_options=valid_options['centrality'][
                        'threshold_options'],
                     w_options=valid_options['centrality']['weight_options'])
def create_centrality_wf(wf_name: str, method_option: str,
                         weight_options: LIST[str], threshold_option: str,
                         threshold: float, num_threads: Optional[int] = 1,
                         memory_gb: Optional[float] = 1.0,
                         base_dir: Optional[Union[Path, str]] = None
                         ) -> Workflow:
    """
    Function to create the afni-based centrality workflow.

    .. seealso::

        * :py:func:`~CPAC.network_centrality.pipeline.connect_centrality_workflow`
        * :py:func:`~CPAC.network_centrality.utils.create_merge_node`
        * :py:func:`~CPAC.network_centrality.utils.sep_nifti_subbriks`

    Parameters
    ----------
    wf_name : string
        the name of the workflow
    method_option : string
        one of {m_options}
    weight_options : list
        one or more of {w_options}
    threshold_option : string
        one of {t_options}
    threshold : float
        the threshold value for thresholding the similarity matrix
    num_threads : integer, optional
        the number of threads to utilize for centrality computation; default=1
    memory_gb : float,optional
        the amount of memory the centrality calculation will take (GB);
        default=1.0
    base_dir : path or str, optional
        the base directory for the workflow; default=None

    Returns
    -------
    centrality_wf : nipype Workflow
        the initialized nipype workflow for the afni centrality command

    Notes
    -----
    Workflow Inputs::

        inputspec.in_file : string
            path to input functional data NIfTI file

        inputspec.template : string
            path to input mask template NIfTI file

        inputspec.threshold : float
            threshold value for thresholding the similarity matrix

    Workflow Outputs::

        outputspec.outfile_list : list of strings
            list of paths to output files (binarized and weighted)
    """  # pylint: disable=line-too-long
    from CPAC.pipeline import nipype_pipeline_engine as pe
    from nipype.interfaces import utility as util
    from CPAC.network_centrality import utils
    from CPAC.utils.interfaces.function import Function

    test_thresh = threshold

    if threshold_option == 'Sparsity threshold':
        test_thresh = threshold / 100.0

    method_option, threshold_option = \
        utils.check_centrality_params(method_option, threshold_option,
                                      test_thresh)
    # Eigenvector centrality and AFNI ≥ 21.1.1?
    ecm_gte_21_1_01 = ((method_option == 'eigenvector_centrality') and
                       AFNI_GTE_21_1_1)
    out_names = tuple(f'{method_option}_{x}' for x in weight_options)
    if base_dir is None:
        centrality_wf = pe.Workflow(name=wf_name)
    else:
        centrality_wf = pe.Workflow(name=wf_name, base_dir=base_dir)

    input_node = pe.Node(util.IdentityInterface(fields=['in_file',
                                                        'template',
                                                        'threshold']),
                         name='inputspec')
    input_node.inputs.threshold = threshold
    output_node = pe.Node(util.IdentityInterface(fields=['outfile_list']),
                          name='outputspec')

    # Degree centrality
    if method_option == 'degree_centrality':
        afni_centrality_node = pe.Node(DegreeCentrality(environ={
                'OMP_NUM_THREADS': str(num_threads)
            }), name='afni_centrality', mem_gb=memory_gb)
        afni_centrality_node.inputs.out_file = \
            'degree_centrality_merged.nii.gz'

    # Eigenvector centrality
    elif method_option == 'eigenvector_centrality':
        if ecm_gte_21_1_01:
            afni_centrality_node = pe.MapNode(ECM(environ={
                'OMP_NUM_THREADS': str(num_threads)
            }), name='afni_centrality', mem_gb=memory_gb,
                iterfield=['do_binary', 'out_file'])
            afni_centrality_node.inputs.out_file = [
                f'eigenvector_centrality_{w_option}.nii.gz' for
                w_option in weight_options]
            afni_centrality_node.inputs.do_binary = [
                w_option == 'Binarized' for w_option in weight_options]
            centrality_wf.connect(afni_centrality_node, 'out_file',
                                  output_node, 'outfile_list')
        else:
            afni_centrality_node = pe.Node(ECM(environ={
                'OMP_NUM_THREADS': str(num_threads)
            }), name='afni_centrality', mem_gb=memory_gb)
            afni_centrality_node.inputs.out_file = \
                'eigenvector_centrality_merged.nii.gz'
        afni_centrality_node.inputs.memory = memory_gb  # 3dECM input only

    # lFCD
    elif method_option == 'local_functional_connectivity_density':
        afni_centrality_node = pe.Node(LFCD(environ={
                'OMP_NUM_THREADS': str(num_threads)
            }), name='afni_centrality', mem_gb=memory_gb)
        afni_centrality_node.inputs.out_file = 'lfcd_merged.nii.gz'

    if not ecm_gte_21_1_01:
        # Need to separate sub-briks except for 3dECM if AFNI > 21.1.01
        sep_subbriks_node = \
            pe.Node(Function(input_names=['nifti_file', 'out_names'],
                             output_names=['output_niftis'],
                             function=utils.sep_nifti_subbriks),
                    name='sep_nifti_subbriks')
        sep_subbriks_node.inputs.out_names = out_names
        centrality_wf.connect([(afni_centrality_node, sep_subbriks_node,
                                [('out_file', 'nifti_file')]),
                               (sep_subbriks_node, output_node,
                                [('output_niftis', 'outfile_list')])])

    afni_centrality_node.interface.num_threads = num_threads

    # Connect input image and mask template
    centrality_wf.connect([(input_node, afni_centrality_node,
                            [('in_file', 'in_file'),
                             ('template', 'mask')])])

    # If we're doing significance thresholding, convert to correlation
    if threshold_option == 'Significance threshold':
        # Check and (possibly) conver threshold
        convert_thr_node = pe.Node(
            Function(input_names=['datafile',
                                  'p_value',
                                  'two_tailed'],
                     output_names=['rvalue_threshold'],
                     function=utils.convert_pvalue_to_r),
            name='convert_threshold')
        # Wire workflow to connect in conversion node
        centrality_wf.connect([(input_node, convert_thr_node,
                                [('in_file', 'datafile'),
                                 ('threshold', 'p_value')]),
                               (convert_thr_node, afni_centrality_node,
                                [('rvalue_threshold', 'thresh')])])

    # Sparsity thresholding
    elif threshold_option == 'Sparsity threshold':
        # Check to make sure it's not lFCD
        if method_option == 'local_functional_connectivity_density':
            raise ThresholdOptionError(threshold_option, method_option)

        # Otherwise, connect threshold to sparsity input
        centrality_wf.connect(input_node, 'threshold',
                              afni_centrality_node, 'sparsity')

    # Correlation thresholding
    elif threshold_option == 'Correlation threshold':
        centrality_wf.connect(input_node, 'threshold',
                              afni_centrality_node, 'thresh')

    return centrality_wf
