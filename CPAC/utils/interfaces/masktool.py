from builtins import str, bytes
import inspect

from nipype import logging
from nipype.interfaces.base import (CommandLineInputSpec, CommandLine, Directory, TraitedSpec,
                    traits, isdefined, File, InputMultiObject, InputMultiPath,
                    Undefined, Str)

from nipype.interfaces.afni.base import (AFNICommandBase, AFNICommand, AFNICommandInputSpec,
                   AFNICommandOutputSpec, AFNIPythonCommandInputSpec,
                   AFNIPythonCommand)
from nipype.interfaces.io import IOBase, add_traits
from nipype.utils.filemanip import ensure_list
from nipype.utils.functions import getsource, create_function_from_source

iflogger = logging.getLogger('nipype.interface')


class MaskToolInputSpec(AFNICommandInputSpec):
    in_files = InputMultiPath(
        File(desc='input file to 3dmerge', exists=True),
        argstr='-input %s',
        position=-1,
        mandatory=True,
        copyfile=False)
    out_file = File(
        name_template='%s_mask',
        desc='output image file name',
        argstr='-prefix %s',
        name_source='in_files')
    count = traits.Bool(
        desc='Instead of created a binary 0/1 mask dataset, create one with '
        'counts of voxel overlap, i.e., each voxel will contain the '
        'number of masks that it is set in.',
        argstr='-count',
        position=2)
    datum = traits.Enum(
        'byte',
        'short',
        'float',
        argstr='-datum %s',
        desc='specify data type for output. Valid types are \'byte\', '
        '\'short\' and \'float\'.')
    dilate_inputs = Str(
        desc='Use this option to dilate and/or erode datasets as they are '
        'read. ex. \'5 -5\' to dilate and erode 5 times',
        argstr='-dilate_inputs %s')
    dilate_results = Str(
        desc='dilate and/or erode combined mask at the given levels.',
        argstr='-dilate_results %s')
    frac = traits.Float(
        desc='When combining masks (across datasets and sub-bricks), use '
        'this option to restrict the result to a certain fraction of the '
        'set of volumes',
        argstr='-frac %s')
    inter = traits.Bool(
        desc='intersection, this means -frac 1.0', argstr='-inter')
    union = traits.Bool(desc='union, this means -frac 0', argstr='-union')
    fill_holes = traits.Bool(
        desc='This option can be used to fill holes in the resulting mask, '
        'i.e. after all other processing has been done.',
        argstr='-fill_holes')
    fill_dirs = Str(
        desc='fill holes only in the given directions. This option is for use '
        'with -fill holes. should be a single string that specifies '
        '1-3 of the axes using {x,y,z} labels (i.e. dataset axis order), '
        'or using the labels in {R,L,A,P,I,S}.',
        argstr='-fill_dirs %s',
        requires=['fill_holes'])
    verbose = traits.Int(
        desc='specify verbosity level, for 0 to 3', argstr='-verb %s')


class MaskToolOutputSpec(TraitedSpec):
    out_file = File(desc='mask file', exists=True)


class MaskTool(AFNICommand):
    """3dmask_tool - for combining/dilating/eroding/filling masks
    For complete details, see the `3dmask_tool Documentation.
    <https://afni.nimh.nih.gov/pub../pub/dist/doc/program_help/3dmask_tool.html>`_
    Examples
    ========
    >>> from nipype.interfaces import afni
    >>> masktool = afni.MaskTool()
    >>> masktool.inputs.in_file = 'functional.nii'  # doctest: +SKIP
    >>> masktool.inputs.outputtype = 'NIFTI'
    >>> masktool.cmdline  # doctest: +SKIP
    '3dmask_tool -prefix functional_mask.nii -input functional.nii'
    >>> res = automask.run()  # doctest: +SKIP
    """

    _cmd = '3dmask_tool'
    input_spec = MaskToolInputSpec
    output_spec = MaskToolOutputSpec
