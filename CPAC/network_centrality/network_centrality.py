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
    degree_outfile = File(desc='The binarized and weighted degree centrality '\
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
        outputs['degree_outfile'] = os.path.abspath(self.inputs.prefix)
        if self.inputs.out_1d:
            outputs['one_d_outfile'] = os.path.abspath(self.inputs.out_1d)

        # Return outputs
        return outputs


# Calculate eigenvector centrality from one_d file
def calc_eigen_from_1d(one_d_file, num_threads, mask_file):
    '''
    '''

    # Import packages
    import os

    import nibabel as nib
    import numpy as np
    import scipy.sparse.linalg as linalg

    from CPAC.network_centrality import utils

    # Init variables
    num_eigs = 1
    which_eigs = 'LM'
    max_iter = 1000

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

    # Return the eigenvector output file
    return eigen_outfiles


# Return the network centrality workflow
def create_network_centrality_wf(wf_name='network_centrality', num_threads=1,
                                 memory=1, run_eigen=False):
    '''
    '''

    # Import packages
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import CPAC.network_centrality.utils as utils

    # Init variables
    centrality_wf = pe.Workflow(name=wf_name)

    # Define main input/function node
    afni_centrality_node = \
        pe.Node(interface=\
                afniDegreeCentrality(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                name='afni_degree_centrality')

    # Limit its num_threads and memory via ResourceMultiProc plugin
    afni_centrality_node.interface.num_threads = num_threads
    afni_centrality_node.interface.memory = memory

    # Define degree seperate bin/wght node
    sep_subbriks_node = \
        pe.Node(interface=util.Function(input_names=['nifti_file', 'out_names'],
                                        output_names=['output_niftis'],
                                        function=utils.sep_nifti_subbriks),
                name='sep_nifti_subbriks')

    # Connect the degree centrality output image to seperate subbriks node
    centrality_wf.connect(afni_centrality_node, 'degree_outfile',
                          sep_subbriks_node, 'nifti_file')
    sep_subbriks_node.inputs.out_names = ('degree_centrality_binarize',
                                          'degree_centrality_weighted')

    # Define outputs node
    output_node = \
        pe.Node(interface=util.IdentityInterface(fields=['degree_outfile_list',
                                                         'eigen_outfile_list',
                                                         'one_d_output']),
                name='output_node')

    # If run_eigen is set, connect 1d output to run_eigen node
    if run_eigen:
        # Tell 3dDegreeCentrality to create 1D file
        afni_centrality_node.inputs.out_1d = 'similarity_matrix.1D'

        # Init run eigenvector centrality node
        run_eigen_node = \
            pe.Node(interface=util.Function(input_names=['one_d_file',
                                                         'num_threads',
                                                         'mask_file'],
                                            output_names=['eigen_outfiles'],
                                            function=calc_eigen_from_1d),
                    name='afni_eigen_centrality')
        # And pass in the number of threads for it to use
        run_eigen_node.inputs.num_threads = num_threads

        # Limit its num_threads and memory via ResourceMultiProce plugin
        run_eigen_node.interface.num_threads = num_threads
        run_eigen_node.interface.memory = memory

        # Connect in the run eigenvector node to the workflow 
        centrality_wf.connect(afni_centrality_node, 'one_d_outfile',
                              run_eigen_node, 'one_d_file')
        centrality_wf.connect(run_eigen_node, 'eigen_outfiles',
                              output_node, 'eigen_outfile_list')

    # Connect the degree centrality outputs to output_node
    centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                          output_node, 'degree_outfile_list')
    centrality_wf.connect(afni_centrality_node, 'one_d_outfile',
                          output_node, 'one_d_output')

    # Return the centrality workflow
    return centrality_wf
