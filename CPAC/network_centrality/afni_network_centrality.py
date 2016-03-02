# CPAC/network_centrality/afni_network_centrality.py
#
# Authors: Daniel Clark

'''
This module contains functions which build and return the network
centrality nipype workflow
'''

# Import packages
from nipype.interfaces.afni.base import AFNICommand, AFNICommandInputSpec,\
                                        AFNICommandOutputSpec
from nipype.interfaces.base import traits, File


class CentralityInputSpec(AFNICommandInputSpec):
    """
    inherits the out_file parameter from AFNICommandOutputSpec base class
    """

    mask = File(desc='mask file to mask input data',
                   argstr="-mask %s",
                   exists=True)

    thresh = traits.Float(desc='threshold to exclude connections where corr <= thresh',
                          argstr='-thresh %f')

    polort = traits.Int(desc='', argstr='-polort %d')

    autoclip = traits.Bool(desc='Clip off low-intensity regions in the dataset',
                           argstr='-autoclip')

    automask = traits.Bool(desc='Mask the dataset to target brain-only voxels',
                           argstr='-automask')


class DegreeCentralityInputSpec(CentralityInputSpec):
    """
    inherits the out_file parameter from AFNICommandOutputSpec base class
    """

    in_file = File(desc='input file to 3dDegreeCentrality',
                   argstr='%s',
                   position=-1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)

    pearson = traits.Bool(desc='Correlation is normal pearson (default)',
                          argstr='-pearson')

    sparsity = traits.Float(desc='only take the top percent of connections',
                            argstr='-sparsity %f')

    oned_file = traits.Str(desc='output filepath to text dump of correlation matrix',
                           argstr='-out1D %s', mandatory=False)


class DegreeCentralityOutputSpec(AFNICommandOutputSpec):
    """
    inherits the out_file parameter from AFNICommandOutputSpec base class
    """

    oned_file = File(desc='The text output of the similarity matrix computed'\
                          'after thresholding with one-dimensional and '\
                          'ijk voxel indices, correlations, image extents, '\
                          'and affine matrix')


class DegreeCentrality(AFNICommand):
    """Performs degree centrality on a dataset using a given maskfile
    via 3dDegreeCentrality
    For complete details, see the `3dDegreeCentrality Documentation.
    <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDegreeCentrality.html>
    Examples
    ========
    >>> from nipype.interfaces import afni as afni
    >>> degree = afni.DegreeCentrality()
    >>> degree.inputs.in_file = 'func_preproc.nii'
    >>> degree.inputs.mask = 'mask.nii'
    >>> degree.inputs.sparsity = 1 # keep the top one percent of connections
    >>> degree.cmdline
    '3dDegreeCentrality -sparsity 1 -mask mask.nii func_preproc.nii'
    >>> res = degree.run() # doctest: +SKIP
    """

    _cmd = '3dDegreeCentrality'
    input_spec = DegreeCentralityInputSpec
    output_spec = DegreeCentralityOutputSpec

    # Re-define generated inputs
    def _list_outputs(self):
        # Import packages
        import os

        # Update outputs dictionary if oned file is defined
        outputs = super(DegreeCentrality, self)._list_outputs()
        if self.inputs.oned_file:
            outputs['oned_file'] = os.path.abspath(self.inputs.oned_file)

        return outputs


class ECMInputSpec(CentralityInputSpec):
    """
    inherits the out_file parameter from AFNICommandOutputSpec base class
    """

    in_file = File(desc='input file to 3dECM',
                   argstr='%s',
                   position=-1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)

    full = traits.Bool(desc='use full power method for thresholded and/or '\
                            'binary eigenvector centrality',
                       argstr='-full')

    fecm = traits.Bool(desc='use shortcut method for computing fast '\
                            'eigenvector centrality mapping; used during non-'\
                            'thresholded and/or non-binary centrality',
                       argstr='-fecm')

    binary = traits.Bool(desc='calculate binarized centrlaity rather than '\
                              'weighted; requires threshold; forces scale=1.0 '\
                              'and shift=0.0',
                         argstr='-binary')

    sparsity = traits.Float(desc='only take the top percent of connections',
                            argstr='-sparsity %f')

    shift = traits.Float(desc='shifting value to enforce non-negativity among '\
                              'correlations; should be >= 0',
                         argstr='-shift %f')

    scale = traits.Float(desc='scale correlations after shift; should be >= 0',
                         argstr='-scale %f')

    eps = traits.Float(desc='the stopping criterion for power iteration; '\
                            'default=0.1',
                       argstr='-eps %f')

    max_iter = traits.Int(desc='maximum number of iterations for power '\
                               'iteration; default=1000',
                          argstr='-max_iter %d')

    memory = traits.Float(desc='Memory limit (GB)- if 3dECM uses more than '\
                               'this the program will error out instead of '\
                               'risk crashing the system; default=2',
                          argstr='-memory %f')


class ECM(AFNICommand):
    """Performs eigenvector centrality on a dataset using a given maskfile
    via 3dECM
    For complete details, see the `3dECM Documentation.
    <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dECM.html>
    Examples
    ========
    >>> from nipype.interfaces import afni as afni
    >>> eig = afni.ECM()
    >>> eig.inputs.in_file = 'func_preproc.nii'
    >>> eig.inputs.mask = 'mask.nii'
    >>> eig.inputs.sparsity = 1 # keep the top one percent of connections
    >>> eig.cmdline
    '3dECM -sparsity 1 -mask mask.nii func_preproc.nii'
    >>> res = eig.run() # doctest: +SKIP
    """

    _cmd = '3dECM'
    input_spec = ECMInputSpec
    output_spec = AFNICommandOutputSpec


class LFCDInputSpec(CentralityInputSpec):
    """
    inherits the out_file parameter from AFNICommandOutputSpec base class
    """

    in_file = File(desc='input file to 3dLFCD',
                   argstr='%s',
                   position=-1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)


class LFCD(AFNICommand):
    """Performs degree centrality on a dataset using a given maskfile
    via 3dLFCD
    For complete details, see the `3dLFCD Documentation.
    <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dLFCD.html>
    Examples
    ========
    >>> from nipype.interfaces import afni as afni
    >>> lfcd = afni.LFCD()
    >>> lfcd.inputs.in_file = 'func_preproc.nii'
    >>> lfcd.inputs.mask = 'mask.nii'
    >>> lfcd.inputs.threshold = .8 # keep all connections with corr >= 0.8
    >>> lfcd.cmdline
    '3dLFCD -threshold 0.8 -mask mask.nii func_preproc.nii'
    >>> res = lfcd.run() # doctest: +SKIP
    """

    _cmd = '3dLFCD'
    input_spec = LFCDInputSpec
    output_spec = AFNICommandOutputSpec


# Return the afni centrality/lfcd workflow
def create_afni_centrality_wf(wf_name, method_option, threshold_option,
                              threshold, num_threads=1, memory_gb=1):
    '''
    '''

    # Import packages
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import CPAC.network_centrality.utils as utils

    # Check the centrality parameters
    method_option, threshold_option = \
        utils.check_centrality_params(method_option, threshold_option, threshold)

    # Init variables
    centrality_wf = pe.Workflow(name=wf_name)

    # Create inputspec node
    input_node = pe.Node(util.IdentityInterface(fields=['in_file',
                                                        'template',
                                                        'threshold']),
                         name='inputspec')

    # Input threshold
    input_node.inputs.threshold = threshold

    # Define main input/function node
    # Degree centrality
    if method_option == 'degree':
        afni_centrality_node = \
            pe.Node(DegreeCentrality(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                    name=method_option)
        afni_centrality_node.inputs.out_file = 'degree_centrality_merged.nii.gz'
        out_names = ('degree_centrality_binarize', 'degree_centrality_weighted')
    # Eigenvector centrality
    elif method_option == 'eigenvector':
        # For now, two nodes for bin, wght, respectibvely
        afni_centrality_node = \
        pe.Node(ECM(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                name=method_option)
        afni_centrality_node.inputs.out_file = 'eigenvector_centrality_weighted.nii.gz'
        #afni_centrality_node.inputs.binary = True
    # lFCD
    elif method_option == 'lfcd':
        afni_centrality_node = \
            pe.Node(LFCD(environ={'OMP_NUM_THREADS' : str(num_threads)}),
                    name=method_option)
        afni_centrality_node.inputs.out_file = 'lfcd_merged.nii.gz'
        out_names = ('lfcd_binarize', 'lfcd_weighted')

    # Limit its num_threads and memory via ResourceMultiProc plugin
    afni_centrality_node.interface.num_threads = num_threads
    afni_centrality_node.interface.estimated_memory_gb = memory_gb

    # Connect input image and mask tempalte
    centrality_wf.connect(input_node, 'in_file',
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

    # Define outputs node
    output_node = pe.Node(util.IdentityInterface(fields=['outfile_list',
                                                         'oned_output']),
                          name='outputspec')

    # If we're doing degree or lfcd, need to seprate sub-briks
    if method_option == 'degree' or method_option == 'lfcd':
        # Define seperate bin/wght node
        sep_subbriks_node = \
            pe.Node(util.Function(input_names=['nifti_file', 'out_names'],
                                  output_names=['output_niftis'],
                                  function=utils.sep_nifti_subbriks),
                    name='sep_nifti_subbriks')
        sep_subbriks_node.inputs.out_names = out_names

        # Connect the degree centrality output image to seperate subbriks node
        centrality_wf.connect(afni_centrality_node, 'out_file',
                              sep_subbriks_node, 'nifti_file')
        centrality_wf.connect(sep_subbriks_node, 'output_niftis',
                              output_node, 'outfile_list')
    else:
        centrality_wf.connect(afni_centrality_node, 'out_file',
                              output_node, 'outfile_list')


    # Return the centrality workflow
    return centrality_wf
