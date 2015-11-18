# CPAC/network_centrality/network_centraliy.py
#
# Authors: Daniel Clark

'''
This module contains functions which build and return the network
centrality nipype workflow
'''

# Import packages
from nipype.interfaces.base import (CommandLine, CommandLineInputSpec,
                                    TraitedSpec)
from nipype.interfaces.base import traits, File

# Input spec class
class afniDegreeCentralityInputSpec(CommandLineInputSpec):
    '''
    '''

    # Define class variables
    prefix = traits.Str(exists=True, argstr='-prefix %s', position=0,
                        desc='Output file name prefix', mandatory=True)
    mask = File(argstr='-mask %s', exists=True, position=1,
                desc='Mask file to use on input data')
    thresh = traits.Float(argstr='-thresh %f', position=2,
                          desc='Threshold to exclude where corr <= thresh')
    sparsity = traits.Float(argstr='-sparsity %f', position=3,
                            desc='The percentage of correlations to keep')
    out_1d = traits.Str(argstr='-out1D %s', position=4,
                        desc='Filepath to output 1D file with similarity matrix')
    polort = traits.Int(argstr='-polort %d', position=5,
                        desc='')
    autoclip = traits.Bool(argstr='-autoclip', desc='Clip off low-intensity '\
                                                    'regions in the dataset')
    automask = traits.Bool(argstr='-automask', desc='Mask the dataset to target '\
                                                    'brain-only voxels')
    dataset = File(argstr='%s', exists=True, position=-1,
                   desc='Functional input dataset to use')


# Output spec class
class afniDegreeCentralityOutputSpec(TraitedSpec):
    '''
    '''

    # Define command outputs
    img_outfile = File(desc='The binarized and weighted degree centrality '\
                            'images stored in two sub-briks of a nifti',
                       exists=True)
    one_d_outfile = File(desc='The text output of the similarity matrix computed'\
                             'after thresholding with one-dimensional and '\
                             'ijk voxel indices, correlations, image extents, '\
                             'and affine matrix')


# Command line execution class
class afniDegreeCentrality(CommandLine):
    '''
    '''

    # Define command, input, and output spec
    _cmd = '3dDegreeCentrality'
    input_spec = afniDegreeCentralityInputSpec
    output_spec = afniDegreeCentralityOutputSpec

    # Gather generated outputs
    def _list_outputs(self):

        # Import packages
        import os

        # Get generated outputs dictionary and assign generated outputs
        # to out output spec
        outputs = self.output_spec().get()
        outputs['img_outfile'] = os.path.abspath(self.inputs.prefix)
        if self.inputs.out_1d:
            outputs['one_d_outfile'] = os.path.abspath(self.inputs.out_1d)

        # Return outputs
        return outputs


# Input spec class
class afniLFCDInputSpec(CommandLineInputSpec):
    '''
    '''

    # Define class variables
    prefix = traits.Str(exists=True, argstr='-prefix %s', position=0,
                        desc='Output file name prefix', mandatory=True)
    mask = File(argstr='-mask %s', exists=True, position=1,
                desc='Mask file to use on input data')
    thresh = traits.Float(argstr='-thresh %f', position=2,
                          desc='Threshold to exclude where corr <= thresh')
    out_1d = traits.Str(argstr='-out1D %s', position=4,
                        desc='Filepath to output 1D file with similarity matrix')
    polort = traits.Int(argstr='-polort %d', position=5,
                        desc='')
    autoclip = traits.Bool(argstr='-autoclip', desc='Clip off low-intensity '\
                                                    'regions in the dataset')
    automask = traits.Bool(argstr='-automask', desc='Mask the dataset to target '\
                                                    'brain-only voxels')
    dataset = File(argstr='%s', exists=True, position=-1,
                   desc='Functional input dataset to use')


# Output spec class
class afniLFCDOutputSpec(TraitedSpec):
    '''
    '''

    # Define command outputs
    img_outfile = File(desc='The binarized and weighted lFCD images stored in two '\
                            'sub-briks of a nifti',
                       exists=True)


# Command line execution class
class afniLFCD(CommandLine):
    '''
    '''

    # Define command, input, and output spec
    _cmd = '3dLFCD'
    input_spec = afniLFCDInputSpec
    output_spec = afniLFCDOutputSpec

    # Gather generated outputs
    def _list_outputs(self):

        # Import packages
        import os

        # Get generated outputs dictionary and assign generated outputs
        # to out output spec
        outputs = self.output_spec().get()
        outputs['img_outfile'] = os.path.abspath(self.inputs.prefix)

        # Return outputs
        return outputs


# Calculate eigenvector centrality from one_d file
def calc_eigen_from_1d(one_d_file, num_threads, mask_file):
    '''
    '''

    # Import packages
    import os
    import logging

    import nibabel as nib
    import numpy as np
    import scipy.sparse.linalg as linalg

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
        utils.parse_and_return_mats(one_d_file, mask_arr)

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
                % (one_d_file, bin_eig_val, wght_eig_val))

    # Return the eigenvector output file
    return eigen_outfiles


# Return the afni centrality/lfcd workflow
def create_afni_centrality_wf(wf_name='network_centrality',
                              method_option, threshold_option,
                              num_threads=1, memory=1):
    '''
    '''

    # Import packages
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import CPAC.network_centrality.utils as utils

    # Init variables
    centrality_wf = pe.Workflow(name=wf_name)

    # Create inputspec node
    input_node = pe.Node(util.IdentityInterface(fields=['datafile',
                                                        'template',
                                                        'threshold']),
                         name='inputspec')

    # Define main input/function node
    # If it's degree or eigenvector, initiate the degree centrality node
    if method_option == 'deg' or method_option == 'eig':
        afni_centrality_node = \
            pe.Node(afniDegreeCentrality(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                    name='afni_centrality')
        afni_centrality_node.inputs.prefix = 'degree_centrality.nii.gz'
    # Otherwise, if it's lFCD, initiate the lFCD node
    elif method_option == 'lfcd':
        afni_centrality_node = \
            pe.Node(afniLFCD(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                    name='afni_centrality')
        afni_centrality_node.inputs.prefix = 'lfcd.nii.gz'

    # Limit its num_threads and memory via ResourceMultiProc plugin
    afni_centrality_node.interface.num_threads = num_threads
    afni_centrality_node.interface.memory = memory

    # Connect input image and mask tempalte
    centrality_wf.connect(input_node, 'dataset_file',
                          afni_centrality_node, 'dataset')
    centrality_wf.connect(input_node, 'template',
                          afni_centrality_node, 'mask')

    # If we're doing significan thresholding, convert to correlation
    if threshold_option == 'pval':
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
    elif threshold_option == 'sparse':
        # Check to make sure it's not lFCD
        if method_option == 'lfcd':
            err_msg = 'Sparsity thresholding is not supported for lFCD'
            raise Exception(err_msg)
        # Otherwise, connect threshold to sparsity input
        centrality_wf.connect(input_node, 'theshold',
                              afni_centrality_node, 'sparsity')

    # Correlation thresholding
    elif threshold_option == 'rval':
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
                                                         'one_d_output']),
                          name='output_node')
    # Degree centrality
    if method_option == 'deg':
        # Connect the degree centrality output image to seperate subbriks node
        centrality_wf.connect(afni_centrality_node, 'img_outfile',
                              sep_subbriks_node, 'nifti_file')
        # Name output files
        sep_subbriks_node.inputs.out_names = ('degree_centrality_binarize',
                                              'degree_centrality_weighted')
        # Connect the degree centrality outputs to output_node
        centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                              output_node, 'outfile_list')
    # If running eigenvector centrality, insert additional node
    elif method_option == 'eig':
        # Tell 3dDegreeCentrality to create 1D file
        afni_centrality_node.inputs.out_1d = 'similarity_matrix.1D'

        # Init run eigenvector centrality node
        run_eigen_node = \
            pe.Node(util.Function(input_names=['one_d_file',
                                               'num_threads',
                                               'mask_file'],
                                  output_names=['eigen_outfiles'],
                                  function=calc_eigen_from_1d),
                    name='afni_eigen_centrality')

        # Limit its num_threads and memory via ResourceMultiProce plugin
        run_eigen_node.interface.num_threads = num_threads
        run_eigen_node.interface.memory = memory
        # MKL threads for scipy eigenvector built with Intel MKL
        run_eigen_node.inputs.num_threads = num_threads

        # Connect in the run eigenvector node to the workflow
        centrality_wf.connect(afni_centrality_node, 'one_d_outfile',
                              run_eigen_node, 'one_d_file')
        centrality_wf.connect(input_node, 'template',
                              run_eigen_node, 'mask_file')
        # Connect outputs
        centrality_wf.connect(run_eigen_node, 'eigen_outfiles',
                              output_node, 'outfile_list')
        centrality_wf.connect(afni_centrality_node, 'one_d_outfile',
                              output_node, 'one_d_output')
    # lFCD
    elif method_option == 'lfcd':
        # Connect the output image to seperate subbriks node
        centrality_wf.connect(afni_centrality_node, 'img_outfile',
                              sep_subbriks_node, 'nifti_file')
        # Name output files
        sep_subbriks_node.inputs.out_names = ('lfcd_binarize', 'lfcd_weighted')
        # Connect the outputs to output_node
        centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                              output_node, 'outfile_list')

    # Return the centrality workflow
    return centrality_wf
