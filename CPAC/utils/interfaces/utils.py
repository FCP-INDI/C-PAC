# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

# This functionality is adapted from poldracklab/niworkflows:
# https://github.com/poldracklab/niworkflows/blob/master/niworkflows/interfaces/utils.py  
# https://fmriprep.readthedocs.io/
# https://poldracklab.stanford.edu/
# We are temporarily maintaining our own copy for more granular control.

"""
Utilities

"""


import os
import re
import json
import shutil
import numpy as np
import nibabel as nb
import nilearn.image as nli
# from textwrap import indent
from collections import OrderedDict

import scipy.ndimage as nd
from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.io import add_traits
from nipype.interfaces.base import (
    traits, isdefined, File, InputMultiPath,
    TraitedSpec, BaseInterfaceInputSpec, SimpleInterface,
    DynamicTraitedSpec
)
from CPAC.info import __version__
import copy

LOG = logging.getLogger('nipype.interface')


class CopyXFormInputSpec(DynamicTraitedSpec, BaseInterfaceInputSpec):
    hdr_file = File(exists=True, mandatory=True, desc='the file we get the header from')


class CopyXForm(SimpleInterface):
    """
    Copy the x-form matrices from `hdr_file` to `out_file`.
    """
    input_spec = CopyXFormInputSpec
    output_spec = DynamicTraitedSpec

    def __init__(self, fields=None, **inputs):
        self._fields = fields or ['in_file']
        if isinstance(self._fields, str):
            self._fields = [self._fields]

        super(CopyXForm, self).__init__(**inputs)

        add_traits(self.inputs, self._fields)
        for f in set(self._fields).intersection(list(inputs.keys())):
            setattr(self.inputs, f, inputs[f])

    def _outputs(self):
        base = super(CopyXForm, self)._outputs()
        if self._fields:
            fields = copy.copy(self._fields)
            if 'in_file' in fields:
                idx = fields.index('in_file')
                fields.pop(idx)
                fields.insert(idx, 'out_file')

            base = add_traits(base, fields)
        return base

    def _run_interface(self, runtime):
        for f in self._fields:
            in_files = getattr(self.inputs, f)
            self._results[f] = []
            if not isinstance(in_files, list):  # if isinstance(in_files, str):
                in_files = [in_files]
            for in_file in in_files:
                out_name = fname_presuffix(
                    in_file, suffix='_xform', newpath=runtime.cwd)
                # Copy and replace header
                shutil.copy(in_file, out_name)
                _copyxform(self.inputs.hdr_file, out_name,
                           message='CopyXForm (niworkflows v%s)' % __version__)
                self._results[f].append(out_name)

            # Flatten out one-element lists
            if len(self._results[f]) == 1:
                self._results[f] = self._results[f][0]

        default = self._results.pop('in_file', None)
        if default:
            self._results['out_file'] = default
        return runtime


def _copyxform(ref_image, out_image, message=None):
    # Read in reference and output
    # Use mmap=False because we will be overwriting the output image
    resampled = nb.load(out_image, mmap=False)
    orig = nb.load(ref_image)

    if not np.allclose(orig.affine, resampled.affine):
        LOG.debug(
            'Affines of input and reference images do not match, '
            'FMRIPREP will set the reference image headers. '
            'Please, check that the x-form matrices of the input dataset'
            'are correct and manually verify the alignment of results.')

    # Copy xform infos
    qform, qform_code = orig.header.get_qform(coded=True)
    sform, sform_code = orig.header.get_sform(coded=True)
    header = resampled.header.copy()
    header.set_qform(qform, int(qform_code))
    header.set_sform(sform, int(sform_code))
    header['descrip'] = 'xform matrices modified by %s.' % (message or '(unknown)')

    newimg = resampled.__class__(resampled.get_fdata(), orig.affine, header)
    newimg.to_filename(out_image)
