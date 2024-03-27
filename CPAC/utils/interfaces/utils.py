# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Imports local names from C-PAC instead of from niworkflows
#     * Replaces `niworkflows.interfaces.utils.CopyXForm._outputs._fields.copy` with `copy.copy`
#     * Changes `if isinstance(in_files, str)` to `if not isinstance(in_files, list)`
#     * Removes classes and functions not used by C-PAC
#     * Removes unused imports
#     * Changes deprecated `get_data` to `get_fdata`
#     * Docstrings updated accordingly

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#    Copyright 2020 The NiPreps Developers
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

# Modifications copyright (C) 2019 - 2024  C-PAC Developers
# This file is part of C-PAC.
"""Utilities for Nipype translation of ANTs workflows.

This functionality is adapted from poldracklab/niworkflows:
https://github.com/poldracklab/niworkflows/blob/994dd2dc/niworkflows/interfaces/utils.py
https://fmriprep.readthedocs.io/
https://poldracklab.stanford.edu/
We are temporarily maintaining our own copy for more granular control.
"""

import copy
import shutil

import numpy as np
import nibabel as nib
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    DynamicTraitedSpec,
    File,
    SimpleInterface,
)
from nipype.interfaces.io import add_traits
from nipype.utils.filemanip import fname_presuffix

from CPAC.info import __version__

LOG = logging.getLogger("nipype.interface")


class CopyXFormInputSpec(DynamicTraitedSpec, BaseInterfaceInputSpec):
    """InputSpec for CopyXForm."""

    hdr_file = File(exists=True, mandatory=True, desc="the file we get the header from")


class CopyXForm(SimpleInterface):
    """Copy the x-form matrices from `hdr_file` to `out_file`."""

    input_spec = CopyXFormInputSpec
    output_spec = DynamicTraitedSpec

    def __init__(self, fields=None, **inputs):
        self._fields = fields or ["in_file"]
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
            if "in_file" in fields:
                idx = fields.index("in_file")
                fields.pop(idx)
                fields.insert(idx, "out_file")

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
                    in_file, suffix="_xform", newpath=runtime.cwd
                )
                # Copy and replace header
                shutil.copy(in_file, out_name)
                _copyxform(
                    self.inputs.hdr_file,
                    out_name,
                    message="CopyXForm (niworkflows v%s)" % __version__,
                )
                self._results[f].append(out_name)

            # Flatten out one-element lists
            if len(self._results[f]) == 1:
                self._results[f] = self._results[f][0]

        default = self._results.pop("in_file", None)
        if default:
            self._results["out_file"] = default
        return runtime


def _copyxform(ref_image, out_image, message=None):
    # Read in reference and output
    # Use mmap=False because we will be overwriting the output image
    resampled = nib.load(out_image, mmap=False)
    orig = nib.load(ref_image)

    if not np.allclose(orig.affine, resampled.affine):
        LOG.debug(
            "Affines of input and reference images do not match, "
            "FMRIPREP will set the reference image headers. "
            "Please, check that the x-form matrices of the input dataset"
            "are correct and manually verify the alignment of results."
        )

    # Copy xform infos
    qform, qform_code = orig.header.get_qform(coded=True)
    sform, sform_code = orig.header.get_sform(coded=True)
    header = resampled.header.copy()
    header.set_qform(qform, int(qform_code))
    header.set_sform(sform, int(sform_code))
    header["descrip"] = "xform matrices modified by %s." % (message or "(unknown)")

    newimg = resampled.__class__(resampled.get_fdata(), orig.affine, header)
    newimg.to_filename(out_image)
