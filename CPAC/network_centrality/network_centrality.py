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

# Input spec class
class afniDegreeCentralityInputSpec(CommandLineInputSpec):
    '''
    '''

    # Import packages
    from nipype.interfaces.base import traits, File

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


# Output spec class
class afniDegreeCentralityOutputSpec(TraitedSpec):
    '''
    '''

    # Import packages
    from nipype.interfaces.base import File

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
def parse_and_return_mats(one_d_file):
    '''
    '''

    # Import packages
    import re
    import numpy as np
    import scipy.sparse as sparse

    # Init variables
    # Capture all positive/negative floats/ints
    reg_pattern = r'[-+]?\d*\.\d+|[-+]?\d+'

    # Parse one_d file
    print 'Reading 1D file...'
    with open(one_d_file, 'r') as fopen:
        lines = fopen.readlines()

    # Store line and build a list of numbers
    affine_line = lines[2]
    affine_elements = re.findall(reg_pattern, affine_line)
    # Store into numpy array and reshape, grab only 3x3 and
    # cast as floating point
    affine_matrix = np.array(affine_elements)
    affine_matrix = affine_matrix.reshape((4,4))[:3,:3].astype('float32')

    # Store line and build a list of numbers
    extents_line = lines[4]
    extents_elements = re.findall(reg_pattern, extents_line)
    # Store as tuple of integers
    img_dims = tuple([int(el) for el in extents_elements])
    # Get the one-dim size of matrix
    one_d = np.prod(img_dims)

    # Parse out numbers
    print 'Parsing contents...'
    graph = [re.findall(reg_pattern, line) for line in lines[6:]]

    # Cast as numpy arrays and extract i, j, w
    print 'Creating arrays...'
    graph_arr = np.array(graph)
    i_array = graph_arr[:,0].astype('int32')
    j_array = graph_arr[:,1].astype('int32')
    w_array = graph_arr[:,-1].astype('float32')

    # Construct the sparse matrix
    print 'Constructing sparse matrix...'
    mat_upper_tri = sparse.coo_matrix((w_array, (i_array, j_array)),
                                      shape=(one_d, one_d))

    # Make symmetric
    similarity_matrix = mat_upper_tri + mat_upper_tri.T

    # Return the symmetric matrix
    return similarity_matrix, affine_matrix


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
    mask_img = nib.load(mask_file).get_data()
    if len(np.unique(mask_img) > 2):
        template_type = 1
    else:
        template_type = 0

    # Temporarily set MKL threads to num_threads for this process only
    os.system('export MKL_NUM_THREADS=%d' % num_threads)

    # Get the similarity matrix from the 1D file
    sim_matrix, affine_matrix = parse_and_return_mats(one_d_file)

    # Use scipy's sparse linalg library to get eigen-values/vectors
    eig_val, eig_vect = linalg.eigsh(sim_matrix, k=num_eigs, which=which_eigs,
                                     maxiter=max_iter)

    # Create eigenvector tuple
    centrality_tuple = ('eigenvector_centrality', np.abs(eig_vect))

    eigen_outfile = utils.map_centrality_matrix(centrality_tuple, affine_matrix,
                                                mask_file, template_type)

    # Return the eigenvector output file
    return eigen_outfile


# Return the network centrality workflow
def create_network_centrality_wf(wf_name='network_centrality', num_threads=1,
                                 memory=1, run_eigen=False):
    '''
    '''

    # Import packages
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    # Init variables
    centrality_wf = pe.Workflow(name='afni_centrality')

    # Define main input/function node
    afni_centrality_node = \
        pe.Node(interface=\
                afniDegreeCentrality(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                name='afni_degree_centrality')

    # Limit its num_threads and memory via ResourceMultiProc plugin
    afni_centrality_node.interface.num_threads = num_threads
    afni_centrality_node.interface.memory = memory

    # Define outputs node
    output_node = \
        pe.Node(interface=util.IdentityInterface(fields=['degree_output',
                                                         'eigen_output',
                                                         'one_d_output']),
                name='output_node')

    # If run_eigen is set, connect 1d output to run_eigen node
    if run_eigen:
        # Tell 3dDegreeCentrality to create 1D file
        afni_centrality_node.inputs.out_1d = 'similarity_matrix.1D'

        # Init run eigenvector centrality node
        run_eigen_node = \
            pe.Node(interface=util.Function(input_names=['one_d_file',
                                                         'num_threads'],
                                            output_names=['eigen_outfile'],
                                            function=calc_eigen_from_1d),
                    name='run_eigen_node')
        # And pass in the number of threads for it to use
        run_eigen_node.inputs.num_threads = num_threads

        # Limit its num_threads and memory via ResourceMultiProce plugin
        run_eigen_node.interface.num_threads = num_threads
        run_eigen_node.interface.memory = memory

        # Connect in the run eigenvector node to the workflow 
        centrality_wf.connect(afni_centrality_node, 'one_d_outfile',
                              run_eigen_node, 'one_d_file')
        centrality_wf.connect(run_eigen_node, 'eigen_outfile',
                              output_node, 'eigen_output')

    # Connect the degree centrality outputs to output_node
    centrality_wf.connect(afni_centrality_node, 'degree_outfile',
                          output_node, 'degree_output')
    centrality_wf.connect(afni_centrality_node, 'one_d_outfile',
                          output_node, 'one_d_output')
