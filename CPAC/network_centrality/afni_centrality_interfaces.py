# CPAC/network_centrality/afni_centrality_interfaces.py
#

'''
This module contains the nipype afni interfaces for the centrality
functions: 3dDegreeCentrality, 3dECM, and 3dLFCD in case the nipype
version installed does not contain those interfaces
'''

# Import packages
from nipype.interfaces.afni.base import (AFNICommand, AFNICommandInputSpec,
                                         AFNICommandOutputSpec)
from nipype.interfaces.base import (traits, File)


class CentralityInputSpec(AFNICommandInputSpec):
    """Common input spec class for all centrality-related commmands
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
    """DegreeCentrality inputspec
    """

    in_file = File(desc='input file to 3dDegreeCentrality',
                   argstr='%s',
                   position=-1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)

    sparsity = traits.Float(desc='only take the top percent of connections',
                            argstr='-sparsity %f')

    oned_file = traits.Str(desc='output filepath to text dump of correlation matrix',
                           argstr='-out1D %s', mandatory=False)


class DegreeCentralityOutputSpec(AFNICommandOutputSpec):
    """DegreeCentrality outputspec
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
    >>> degree.inputs.out_file = 'out.nii'
    >>> degree.cmdline
    '3dDegreeCentrality -sparsity 1 -mask mask.nii -prefix out.nii func_preproc.nii'
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
    """ECM inputspec
    """

    in_file = File(desc='input file to 3dECM',
                   argstr='%s',
                   position=-1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)

    sparsity = traits.Float(desc='only take the top percent of connections',
                            argstr='-sparsity %f')

    full = traits.Bool(desc='Full power method; enables thresholding; '\
                            'automatically selected if -thresh or -sparsity '\
                            'are set',
                       argstr='-full')

    fecm = traits.Bool(desc='Fast centrality method; substantial speed '\
                            'increase but cannot accomodate thresholding; '\
                            'automatically selected if -thresh or -sparsity '\
                            'are not set',
                       argstr='-fecm')

    shift = traits.Float(desc='shift correlation coefficients in similarity '\
                              'matrix to enforce non-negativity, s >= 0.0; '\
                              'default = 0.0 for -full, 1.0 for -fecm',
                            argstr='-shift %f')

    scale = traits.Float(desc='scale correlation coefficients in similarity '\
                              'matrix to after shifting, x >= 0.0; '\
                              'default = 1.0 for -full, 0.5 for -fecm',
                            argstr='-scale %f')

    eps = traits.Float(desc='sets the stopping criterion for the power '\
                            'iteration; l2|v_old - v_new| < eps*|v_old|; '\
                            'default = 0.001',
                            argstr='-eps %f')

    max_iter = traits.Int(desc='sets the maximum number of iterations to use '\
                               'in the power iteration; default = 1000',
                          argstr='-max_iter %d')

    memory = traits.Float(desc='Limit memory consumption on system by setting '\
                               'the amount of GB to limit the algorithm to; '\
                               'default = 2GB',
                          argstr='-memory %f')


class ECM(AFNICommand):
    """Performs degree centrality on a dataset using a given maskfile
    via the 3dLFCD command

    For complete details, see the `3dECM Documentation.
    <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dECM.html>

    Examples
    ========

    >>> from nipype.interfaces import afni as afni
    >>> ecm = afni.ECM()
    >>> ecm.inputs.in_file = 'func_preproc.nii'
    >>> ecm.inputs.mask = 'mask.nii'
    >>> ecm.inputs.sparsity = 0.1 # keep top 0.1% of connections
    >>> ecm.inputs.out_file = 'out.nii'
    >>> ecm.cmdline
    '3dECM -sparsity 0.1 -mask mask.nii -prefix out.nii func_preproc.nii'
    >>> res = ecm.run() # doctest: +SKIP
    """

    _cmd = '3dECM'
    input_spec = ECMInputSpec
    output_spec = AFNICommandOutputSpec


class LFCDInputSpec(CentralityInputSpec):
    """LFCD inputspec
    """

    in_file = File(desc='input file to 3dLFCD',
                   argstr='%s',
                   position=-1,
                   mandatory=True,
                   exists=True,
                   copyfile=False)


class LFCD(AFNICommand):
    """Performs degree centrality on a dataset using a given maskfile
    via the 3dLFCD command

    For complete details, see the `3dLFCD Documentation.
    <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dLFCD.html>

    Examples
    ========

    >>> from nipype.interfaces import afni as afni
    >>> lfcd = afni.LFCD()
    >>> lfcd.inputs.in_file = 'func_preproc.nii'
    >>> lfcd.inputs.mask = 'mask.nii'
    >>> lfcd.inputs.thresh = 0.8 # keep all connections with corr >= 0.8
    >>> lfcd.inputs.out_file = 'out.nii'
    >>> lfcd.cmdline
    '3dLFCD -thresh 0.8 -mask mask.nii -prefix out.nii func_preproc.nii'
    >>> res = lfcd.run() # doctest: +SKIP
    """

    _cmd = '3dLFCD'
    input_spec = LFCDInputSpec
    output_spec = AFNICommandOutputSpec