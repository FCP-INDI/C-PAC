# CPAC/network_centrality/network_centraliy.py
#
# Authors: Daniel Clark

'''
This module contains functions which build and return the network
centrality nipype workflow
'''


# Calculate eigenvector centrality from one_d file
def calc_eigen_from_1d(oned_file, num_threads, mask_file):
    '''
    '''

    # Import packages
    import os

    import nibabel as nib
    import numpy as np
    import scipy.sparse.linalg as linalg
    from nipype import logging

    from CPAC.network_centrality import utils

    # Init variables
    num_eigs = 1
    which_eigs = 'LM'
    max_iter = 1000
    logger = logging.getLogger('workflow')

    # See if mask was ROI atlas or mask file
    mask_img = nib.load(mask_file)
    mask_arr = mask_img.get_data()
    mask_affine = mask_img.get_affine()
    if len(np.unique(mask_arr)) > 2:
        template_type = 1
    else:
        template_type = 0

    # Temporarily set MKL threads to num_threads for this process only
    os.system('export MKL_NUM_THREADS=%d' % num_threads)

    # Get the similarity matrix from the 1D file
    bin_sim_matrix, wght_sim_matrix = \
        utils.parse_and_return_mats(oned_file, mask_arr)

    # Use scipy's sparse linalg library to get eigen-values/vectors
    bin_eig_val, bin_eig_vect = linalg.eigsh(bin_sim_matrix, k=num_eigs,
                                             which=which_eigs, maxiter=max_iter)
    wght_eig_val, wght_eig_vect = linalg.eigsh(wght_sim_matrix, k=num_eigs,
                                               which=which_eigs, maxiter=max_iter)

    # Create eigenvector tuple
    centrality_tuple = ('eigenvector_centrality_binarized', np.abs(bin_eig_val*bin_eig_vect))
    bin_outfile = utils.map_centrality_matrix(centrality_tuple, mask_affine,
                                              mask_arr, template_type)

    centrality_tuple = ('eigenvector_centrality_weighted', np.abs(wght_eig_val*wght_eig_vect))
    wght_outfile = utils.map_centrality_matrix(centrality_tuple, mask_affine,
                                               mask_arr, template_type)

    # Grab outfile paths
    eigen_outfiles = [bin_outfile, wght_outfile]

    # Record eigenvalues in logger just in case
    logger.info('Eigenvalues for %s are - bin: %.5f, wght: %.5f' \
                % (oned_file, bin_eig_val, wght_eig_val))

    # Return the eigenvector output file
    return eigen_outfiles


# Return the afni centrality/lfcd workflow
def create_afni_centrality_wf(wf_name, method_option, threshold, threshold_option,
                              num_threads=1, memory=1):
    '''
    '''

    # Import packages
    import nipype.pipeline.engine as pe
    import nipype.interfaces.afni as afni
    import nipype.interfaces.utility as util
    import CPAC.network_centrality.utils as utils

    # Check the centrality parameters
    method_option, threshold_option = \
        utils.check_centrality_params(method_option, threshold_option, threshold)

    # Init variables
    centrality_wf = pe.Workflow(name=wf_name)

    # Create inputspec node
    input_node = pe.Node(util.IdentityInterface(fields=['datafile',
                                                        'template',
                                                        'threshold']),
                         name='inputspec')

    # Input threshold
    input_node.inputs.threshold = threshold

    # Define main input/function node
    # If it's degree or eigenvector, initiate the degree centrality node
    if method_option == 'degree' or method_option == 'eigenvector':
        afni_centrality_node = \
            pe.Node(afni.DegreeCentrality(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                    name='afni_centrality')
        afni_centrality_node.inputs.out_file = 'degree_centrality.nii.gz'
    # Otherwise, if it's lFCD, initiate the lFCD node
    elif method_option == 'lfcd':
        afni_centrality_node = \
            pe.Node(afni.LFCD(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                    name='afni_centrality')
        afni_centrality_node.inputs.out_file = 'lfcd.nii.gz'

    # Limit its num_threads and memory via ResourceMultiProc plugin
    afni_centrality_node.interface.num_threads = num_threads
    afni_centrality_node.interface.estimated_memory = memory

    # Connect input image and mask tempalte
    centrality_wf.connect(input_node, 'datafile',
                          afni_centrality_node, 'in_file')
    centrality_wf.connect(input_node, 'template',
                          afni_centrality_node, 'mask')

    # If we're doing significan thresholding, convert to correlation
    if threshold_option == 'significance':
        # Check and (possibly) conver threshold
        convert_thr_node = pe.Node(util.Function(input_names=['datafile',
                                                              'p_value',
                                                              'two_tailed'],
                                                 output_names=['rvalue_threshold'],
                                                 function=utils.convert_pvalue_to_r),
                                   name='convert_threshold')
        # Wire workflow to connect in conversion node
        centrality_wf.connect(input_node, 'datafile',
                              convert_thr_node, 'datafile')
        centrality_wf.connect(input_node, 'threshold',
                              convert_thr_node, 'p_value')
        centrality_wf.connect(convert_thr_node, 'rvalue_threshold',
                              afni_centrality_node, 'thresh')

    # Sparsity thresholding
    elif threshold_option == 'sparsity':
        # Check to make sure it's not lFCD
        if method_option == 'lfcd':
            err_msg = 'Sparsity thresholding is not supported for lFCD'
            raise Exception(err_msg)
        # Otherwise, connect threshold to sparsity input
        centrality_wf.connect(input_node, 'threshold',
                              afni_centrality_node, 'sparsity')

    # Correlation thresholding
    elif threshold_option == 'correlation':
        centrality_wf.connect(input_node, 'threshold',
                              afni_centrality_node, 'thresh')

    # Define degree seperate bin/wght node
    sep_subbriks_node = \
        pe.Node(util.Function(input_names=['nifti_file', 'out_names'],
                              output_names=['output_niftis'],
                              function=utils.sep_nifti_subbriks),
                name='sep_nifti_subbriks')

    # Define outputs node
    output_node = pe.Node(util.IdentityInterface(fields=['outfile_list',
                                                         'oned_output']),
                          name='output_node')

    # Degree centrality
    if method_option == 'degree':
        # Connect the degree centrality output image to seperate subbriks node
        centrality_wf.connect(afni_centrality_node, 'out_file',
                              sep_subbriks_node, 'nifti_file')
        # Name output files
        sep_subbriks_node.inputs.out_names = ('degree_centrality_binarize',
                                              'degree_centrality_weighted')
        # Connect the degree centrality outputs to output_node
        centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                              output_node, 'outfile_list')

    # If running eigenvector centrality, insert additional node
    elif method_option == 'eigenvector':
        # Tell 3dDegreeCentrality to create 1D file
        afni_centrality_node.inputs.oned_file = 'similarity_matrix.1D'

        # Init run eigenvector centrality node
        run_eigen_node = \
            pe.Node(util.Function(input_names=['oned_file',
                                               'num_threads',
                                               'mask_file'],
                                  output_names=['eigen_outfiles'],
                                  function=calc_eigen_from_1d),
                    name='afni_eigen_centrality')
 
        # Limit its num_threads and memory via ResourceMultiProce plugin
        run_eigen_node.interface.num_threads = num_threads
        run_eigen_node.interface.estimated_memory = memory
        # MKL threads for scipy eigenvector built with Intel MKL
        run_eigen_node.inputs.num_threads = num_threads
 
        # Connect in the run eigenvector node to the workflow
        centrality_wf.connect(afni_centrality_node, 'oned_file',
                              run_eigen_node, 'oned_file')
        centrality_wf.connect(input_node, 'template',
                              run_eigen_node, 'mask_file')
        # Connect outputs
        centrality_wf.connect(run_eigen_node, 'eigen_outfiles',
                              output_node, 'outfile_list')
        centrality_wf.connect(afni_centrality_node, 'oned_file',
                              output_node, 'oned_output')

    # lFCD
    elif method_option == 'lfcd':
        # Connect the output image to seperate subbriks node
        centrality_wf.connect(afni_centrality_node, 'out_file',
                              sep_subbriks_node, 'nifti_file')
        # Name output files
        sep_subbriks_node.inputs.out_names = ('lfcd_binarize', 'lfcd_weighted')
        # Connect the outputs to output_node
        centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                              output_node, 'outfile_list')

    # Return the centrality workflow
    return centrality_wf
