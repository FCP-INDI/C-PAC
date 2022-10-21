from nipype.interfaces.afni.preprocess import DegreeCentrality, ECM, LFCD
from CPAC.pipeline.schema import valid_options
from CPAC.utils.docs import docstring_parameter
from CPAC.utils.interfaces.function import Function
import subprocess


@docstring_parameter(m_options=valid_options['centrality']['method_options'],
                     t_options=valid_options['centrality'][
                        'threshold_options'],
                     w_options=valid_options['centrality']['weight_options'])

def eig_centrality_binarized(memory,environ,in_file,mask,sparsity=None,thresh=None):
    out_file = 'eigenvector_centrality_merged.nii.gz'
    
    cmd = ['3dECM', '-mask', mask, '-memory', memory, '-prefix', out_file, \
        '-sparsity', sparsity, '-thresh', thresh, in_file, '-do_binary']

    subprocess.check_output(cmd)
    
    return out_file

def create_centrality_wf(wf_name, method_option, weight_options,
                         threshold_option, threshold, num_threads=1,
                         memory_gb=1.0):
    """
    Function to create the afni-based centrality workflow

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
    num_threads : integer (optional); default=1
        the number of threads to utilize for centrality computation
    memory_gb : float (optional); default=1.0
        the amount of memory the centrality calculation will take (GB)

    Returns
    -------
    centrality_wf : nipype Workflow
        the initialized nipype workflow for the afni centrality command
    """

    from CPAC.pipeline import nipype_pipeline_engine as pe
    import nipype.interfaces.utility as util
    import CPAC.network_centrality.utils as utils

    test_thresh = threshold

    if threshold_option == 'sparsity':
        test_thresh = threshold / 100.0

    method_option, threshold_option = \
        utils.check_centrality_params(method_option, threshold_option,
                                      test_thresh)

    centrality_wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['in_file',
                                                        'template',
                                                        'threshold']),
                         name='inputspec')

    input_node.inputs.threshold = threshold

    # Degree centrality
    if method_option == 'degree_centrality':
        afni_centrality_node = pe.Node(DegreeCentrality(environ={
                'OMP_NUM_THREADS': str(num_threads)
            }), name='afni_centrality', mem_gb=memory_gb)
        afni_centrality_node.inputs.out_file = \
            'degree_centrality_merged.nii.gz'

    # Eigenvector centrality
    elif method_option == 'eigenvector_centrality':
        if 'Binarized' in weight_options:
            eig_centrality_binarized_import = ['import subprocess']
            afni_centrality_node = pe.Node(Function(input_names=['memory',
                                              'environ',
                                              'in_file',
                                              'mask',
                                              'sparsity',
                                              'thresh'],
                                 output_names=['out_file'],
                                 function=eig_centrality_binarized),
                                imports=eig_centrality_binarized_import,
                   name=f'afni_centrality')
            afni_centrality_node.inputs.environ = {'OMP_NUM_THREADS': str(num_threads)}
            afni_centrality_node.inputs.memory = memory_gb
            afni_centrality_node.inputs.out_file = \
            'eigenvector_centrality_merged.nii.gz'
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

    out_names = tuple(f'{method_option}_{x}' for x in weight_options)

    afni_centrality_node.interface.num_threads = num_threads

    # Connect input image and mask tempalte
    centrality_wf.connect(input_node, 'in_file',
                          afni_centrality_node, 'in_file')
    centrality_wf.connect(input_node, 'template',
                          afni_centrality_node, 'mask')

    # If we're doing significan thresholding, convert to correlation
    if threshold_option == 'Significance threshold':
        # Check and (possibly) conver threshold
        convert_thr_node = pe.Node(
            util.Function(input_names=['datafile',
                                       'p_value',
                                       'two_tailed'],
                          output_names=['rvalue_threshold'],
                          function=utils.convert_pvalue_to_r),
            name='convert_threshold')
        # Wire workflow to connect in conversion node
        centrality_wf.connect(input_node, 'in_file',
                              convert_thr_node, 'datafile')
        centrality_wf.connect(input_node, 'threshold',
                              convert_thr_node, 'p_value')
        centrality_wf.connect(convert_thr_node, 'rvalue_threshold',
                              afni_centrality_node, 'thresh')

    # Sparsity thresholding
    elif threshold_option == 'Sparsity threshold':
        # Check to make sure it's not lFCD
        if method_option == 'local_functional_connectivity_density':
            raise Exception('Sparsity thresholding is not supported for lFCD')

        # Otherwise, connect threshold to sparsity input
        centrality_wf.connect(input_node, 'threshold',
                              afni_centrality_node, 'sparsity')

    # Correlation thresholding
    elif threshold_option == 'Correlation threshold':
        centrality_wf.connect(input_node, 'threshold',
                              afni_centrality_node, 'thresh')

    # Need to separate sub-briks
    sep_subbriks_node = \
        pe.Node(util.Function(input_names=['nifti_file', 'out_names'],
                              output_names=['output_niftis'],
                              function=utils.sep_nifti_subbriks),
                name='sep_nifti_subbriks')

    sep_subbriks_node.inputs.out_names = out_names

    centrality_wf.connect(afni_centrality_node, 'out_file',
                          sep_subbriks_node, 'nifti_file')

    output_node = pe.Node(util.IdentityInterface(fields=['outfile_list',
                                                         'oned_output']),
                          name='outputspec')

    centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                          output_node, 'outfile_list')

    return centrality_wf
