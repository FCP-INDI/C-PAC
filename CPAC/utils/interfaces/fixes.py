# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Imports local names from C-PAC instead of from niworkflows
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
"""
Fixes for Nipype translation of ANTs workflows.

This functionality is adapted from poldracklab/niworkflows:
https://github.com/poldracklab/niworkflows/blob/994dd2dc/niworkflows/interfaces/fixes.py
https://fmriprep.readthedocs.io/
https://poldracklab.stanford.edu/
We are temporarily maintaining our own copy for more granular control.
"""

import os

from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants.resampling import ApplyTransforms

from CPAC.info import __version__
from CPAC.utils.interfaces.utils import _copyxform


class FixHeaderApplyTransforms(ApplyTransforms):
    """Replace nipype.interfaces.ants.resampling.ApplyTransforms.

    Fixes the resampled image header to match the xform of the reference image.
    """

    def _run_interface(self, runtime, correct_return_codes=(0,)):
        # Run normally
        runtime = super(FixHeaderApplyTransforms, self)._run_interface(
            runtime, correct_return_codes
        )

        _copyxform(
            self.inputs.reference_image,
            os.path.abspath(self._gen_filename("output_image")),
            message="%s (niworkflows v%s)" % (self.__class__.__name__, __version__),
        )
        return runtime


class FixHeaderRegistration(Registration):
    """Replace nipype.interfaces.ants.registration.Registration.

    Fixes the resampled image header to match the xform of the reference image.
    """

    def _run_interface(self, runtime, correct_return_codes=(0,)):
        # Run normally
        runtime = super(FixHeaderRegistration, self)._run_interface(
            runtime, correct_return_codes
        )

        # Forward transform
        out_file = self._get_outputfilenames(inverse=False)
        if out_file is not None and out_file:
            _copyxform(
                self.inputs.fixed_image[0],
                os.path.abspath(out_file),
                message="%s (niworkflows v%s)" % (self.__class__.__name__, __version__),
            )

        # Inverse transform
        out_file = self._get_outputfilenames(inverse=True)
        if out_file is not None and out_file:
            _copyxform(
                self.inputs.moving_image[0],
                os.path.abspath(out_file),
                message="%s (niworkflows v%s)" % (self.__class__.__name__, __version__),
            )

        return runtime
