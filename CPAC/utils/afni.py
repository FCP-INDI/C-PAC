
from __future__ import (print_function, division, unicode_literals,
                        absolute_import)
from builtins import str, bytes

import os
import os.path as op
import re
import numpy as np

from nipype.interfaces.base import (
    CommandLineInputSpec, CommandLine, Directory, TraitedSpec,
    traits, isdefined, File, InputMultiPath, Undefined, Str
)

from nipype.interfaces.afni.base import (
    AFNICommand,
    AFNICommandInputSpec,
    AFNICommandOutputSpec
)


class Edge3InputSpec(AFNICommandInputSpec):
    in_file = File(
        desc='input file to 3dedge3',
        argstr='-input %s',
        position=0,
        mandatory=True,
        exists=True,
        copyfile=False)
    out_file = File(
        desc='output image file name', position=-1, argstr='-prefix %s',
        name_source=['in_file'], name_template='%s_edge3')
    datum = traits.Enum(
        'byte',
        'short',
        'float',
        argstr='-datum %s',
        desc='specify data type for output. Valid types are \'byte\', '
        '\'short\' and \'float\'.')
    fscale = traits.Bool(
        desc='Force scaling of the output to the maximum integer range.',
        argstr='-fscale',
        xor=['gscale', 'nscale', 'scale_floats'])
    gscale = traits.Bool(
        desc='Same as \'-fscale\', but also forces each output sub-brick to '
        'to get the same scaling factor.',
        argstr='-gscale',
        xor=['fscale', 'nscale', 'scale_floats'])
    nscale = traits.Bool(
        desc='Don\'t do any scaling on output to byte or short datasets.',
        argstr='-nscale',
        xor=['fscale', 'gscale', 'scale_floats'])
    scale_floats = traits.Float(
        desc='Multiply input by VAL, but only if the input datum is '
        'float. This is needed when the input dataset '
        'has a small range, like 0 to 2.0 for instance. '
        'With such a range, very few edges are detected due to '
        'what I suspect to be truncation problems. '
        'Multiplying such a dataset by 10000 fixes the problem '
        'and the scaling is undone at the output.',
        argstr='-scale_floats %f',
        xor=['fscale', 'gscale', 'nscale'])
    verbose = traits.Bool(
        desc='Print out some information along the way.', argstr='-verbose')


class Edge3(AFNICommand):
    """Does 3D Edge detection using the library 3DEdge
    by Gregoire Malandain (gregoire.malandain@sophia.inria.fr).
    For complete details, see the `3dedge3 Documentation.
    <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dedge3.html>`_
    references_ = [{'entry': BibTeX('@article{Deriche1987,'
                                    'author={R. Deriche},'
                                    'title={Optimal edge detection using recursive filtering},'
                                    'journal={International Journal of Computer Vision},'
                                    'volume={2},',
                                    'pages={167-187},'
                                    'year={1987},'
                                    '}'),
                    'tags': ['method'],
                    },
                   {'entry': BibTeX('@article{MongaDericheMalandainCocquerez1991,'
                                    'author={O. Monga, R. Deriche, G. Malandain, J.P. Cocquerez},'
                                    'title={Recursive filtering and edge tracking: two primary tools for 3D edge detection},'
                                    'journal={Image and vision computing},'
                                    'volume={9},',
                                    'pages={203-214},'
                                    'year={1991},'
                                    '}'),
                    'tags': ['method'],
                    },
                   ]
    Examples
    ========
    >>> from nipype.interfaces import afni
    >>> edge3 = afni.Edge3()
    >>> edge3.inputs.in_file = 'functional.nii'
    >>> edge3.inputs.out_file = 'edges.nii'
    >>> edge3.inputs.datum = 'byte'
    >>> edge3.cmdline
    '3dedge3 -input functional.nii -datum byte -prefix edges.nii'
    >>> res = edge3.run()  # doctest: +SKIP
    """

    _cmd = '3dedge3'
    input_spec = Edge3InputSpec
    output_spec = AFNICommandOutputSpec
