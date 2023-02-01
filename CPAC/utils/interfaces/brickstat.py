import os

from builtins import str, bytes
import inspect

from nipype import logging
from nipype.interfaces.base import (CommandLineInputSpec, CommandLine, Directory, TraitedSpec,
                    traits, isdefined, File, InputMultiObject, InputMultiPath,
                    Undefined, Str)

from nipype.utils.filemanip import (load_json, save_json, split_filename)

from nipype.interfaces.afni.base import (AFNICommandBase, AFNICommand, AFNICommandInputSpec,
                   AFNICommandOutputSpec, AFNIPythonCommandInputSpec,
                   AFNIPythonCommand)
from nipype.interfaces.io import IOBase, add_traits
from nipype.utils.filemanip import ensure_list
from nipype.utils.functions import getsource, create_function_from_source

iflogger = logging.getLogger('nipype.interface')


class BrickStatInputSpec(CommandLineInputSpec):
    in_file = File(
        desc='input file to 3dmaskave',
        argstr='%s',
        position=-1,
        mandatory=True,
        exists=True)
    mask = File(
        desc='-mask dset = use dset as mask to include/exclude voxels',
        argstr='-mask %s',
        position=2,
        exists=True)
    min = traits.Bool(
        desc='print the minimum value in dataset', argstr='-min', position=1)
    slow = traits.Bool(
        desc='read the whole dataset to find the min and max values',
        argstr='-slow')
    max = traits.Bool(
        desc='print the maximum value in the dataset', argstr='-max')
    mean = traits.Bool(
        desc='print the mean value in the dataset', argstr='-mean')
    sum = traits.Bool(
        desc='print the sum of values in the dataset', argstr='-sum')
    var = traits.Bool(desc='print the variance in the dataset', argstr='-var')
    percentile = traits.Tuple(
        traits.Float,
        traits.Float,
        traits.Float,
        desc='p0 ps p1 write the percentile values starting '
        'at p0% and ending at p1% at a step of ps%. '
        'only one sub-brick is accepted.',
        argstr='-percentile %.3f %.3f %.3f')


class BrickStatOutputSpec(TraitedSpec):
    min_val = traits.Float(desc='output')


class BrickStat(AFNICommandBase):
    """Computes maximum and/or minimum voxel values of an input dataset.
    TODO Add optional arguments.

    For complete details, see the `3dBrickStat Documentation.
    <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dBrickStat.html>`_

    Examples
    ========

    >>> from nipype.interfaces import afni
    >>> brickstat = afni.BrickStat()
    >>> brickstat.inputs.in_file = 'functional.nii'  # doctest: +SKIP
    >>> brickstat.inputs.mask = 'skeleton_mask.nii.gz'  # doctest: +SKIP
    >>> brickstat.inputs.min = True
    >>> brickstat.cmdline  # doctest: +SKIP
    '3dBrickStat -min -mask skeleton_mask.nii.gz functional.nii'
    >>> res = brickstat.run()  # doctest: +SKIP

    """
    _cmd = '3dBrickStat'
    input_spec = BrickStatInputSpec
    output_spec = BrickStatOutputSpec

    def aggregate_outputs(self, runtime=None, needed_outputs=None):

        outputs = self._outputs()

        outfile = os.path.join(os.getcwd(), 'stat_result.json')

        if runtime is None:
            try:
                min_val = load_json(outfile)['stat']
            except IOError:
                return self.run().outputs
        else:
            min_val = []
            for line in runtime.stdout.split('\n'):
                if line:
                    values = line.split()
                    if len(values) > 1:
                        min_val.append([float(val) for val in values])
                    else:
                        min_val.extend([float(val) for val in values])

            if len(min_val) == 1:
                min_val = min_val[0]
            save_json(outfile, dict(stat=min_val))
        if type(min_val) == list:
            min_val = min_val[-1]
        outputs.min_val = min_val

        return outputs