#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Removes classes unused by C-PAC
#     * Removes niworkflows doctests
#     * Adds `PrintHeader`, `SetDirectionByMatrix` and their respective (InputSpec and OutputSpec)s
#     * Made private classes public
#     * Removed original example code
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
"""Nipype interfaces for ANTs commands.

Some of this functionality is adapted from nipreps/niworkflows:
- https://github.com/nipreps/niworkflows/blob/994dd2dc/niworkflows/interfaces/ants.py
- https://fmriprep.readthedocs.io/
- https://poldracklab.stanford.edu/
We are temporarily maintaining our own copy for more granular control.
"""

import os

from nipype.interfaces import base
from nipype.interfaces.ants.base import ANTSCommand, ANTSCommandInputSpec
from nipype.interfaces.base import isdefined, traits


class ImageMathInputSpec(ANTSCommandInputSpec):
    """InputSpec for ImageMath."""

    dimension = traits.Int(
        3, usedefault=True, position=1, argstr="%d", desc="dimension of output image"
    )
    output_image = base.File(
        position=2,
        argstr="%s",
        name_source=["op1"],
        name_template="%s_maths",
        desc="output image file",
        keep_extension=True,
    )
    operation = base.Str(
        mandatory=True, position=3, argstr="%s", desc="operations and intputs"
    )
    op1 = base.File(
        exists=True, mandatory=True, position=-2, argstr="%s", desc="first operator"
    )
    op2 = traits.Either(
        base.File(exists=True),
        base.Str,
        position=-1,
        argstr="%s",
        desc="second operator",
    )


class ImageMathOuputSpec(base.TraitedSpec):
    """OutputSpec for ImageMath."""

    output_image = base.File(exists=True, desc="output image file")


class ImageMath(ANTSCommand):
    """
    Operations over images.

    Example:
    --------

    """

    _cmd = "ImageMath"
    input_spec = ImageMathInputSpec
    output_spec = ImageMathOuputSpec


class PrintHeaderInputSpec(ANTSCommandInputSpec):
    """InputSpec for ``PrintHeader``.

    See `PrintHeader: DESCRIPTION <https://manpages.debian.org/testing/ants/PrintHeader.1.en.html#DESCRIPTION>`_ for ``what_information`` values.
    """  # pylint: disable=line-too-long

    image = base.File(
        position=2,
        argstr="%s",
        name_source=["image"],
        desc="image to read header from",
        exists=True,
        mandatory=True,
    )

    what_information = traits.Int(
        position=3, argstr="%i", name="what_information", desc="read what from header"
    )


class PrintHeaderOutputSpec(base.TraitedSpec):
    """OutputSpec for ``PrintHeader``."""

    header = traits.String(name="header")


class PrintHeader(ANTSCommand):
    """Print image header information."""

    _cmd = "PrintHeader"
    # pylint: disable=protected-access
    _gen_filename = base.StdOutCommandLine._gen_filename
    input_spec = PrintHeaderInputSpec
    output_spec = PrintHeaderOutputSpec
    _terminal_output = "stream"

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        """Populate outputs."""
        outputs = super().aggregate_outputs(runtime, needed_outputs)
        outputs.trait_set(header=runtime.stdout)
        self.output_spec().trait_set(header=runtime.stdout)
        return outputs

    def _list_outputs(self):
        return self._outputs().get()


class ResampleImageBySpacingInputSpec(ANTSCommandInputSpec):
    """InputSpec for ResampleImageBySpacing."""

    dimension = traits.Int(
        3, usedefault=True, position=1, argstr="%d", desc="dimension of output image"
    )
    input_image = base.File(
        exists=True, mandatory=True, position=2, argstr="%s", desc="input image file"
    )
    output_image = base.File(
        position=3,
        argstr="%s",
        name_source=["input_image"],
        name_template="%s_resampled",
        desc="output image file",
        keep_extension=True,
    )
    out_spacing = traits.Either(
        traits.List(traits.Float, minlen=2, maxlen=3),
        traits.Tuple(traits.Float, traits.Float, traits.Float),
        traits.Tuple(traits.Float, traits.Float),
        position=4,
        argstr="%s",
        mandatory=True,
        desc="output spacing",
    )
    apply_smoothing = traits.Bool(
        False, argstr="%d", position=5, desc="smooth before resampling"
    )
    addvox = traits.Int(
        argstr="%d",
        position=6,
        requires=["apply_smoothing"],
        desc="addvox pads each dimension by addvox",
    )
    nn_interp = traits.Bool(
        argstr="%d", desc="nn interpolation", position=-1, requires=["addvox"]
    )


class ResampleImageBySpacingOutputSpec(base.TraitedSpec):
    """OutputSpec for ResampleImageBySpacing."""

    output_image = traits.File(exists=True, desc="resampled file")


class ResampleImageBySpacing(ANTSCommand):
    """Resamples an image with a given spacing."""

    _cmd = "ResampleImageBySpacing"
    input_spec = ResampleImageBySpacingInputSpec
    output_spec = ResampleImageBySpacingOutputSpec

    def _format_arg(self, name, trait_spec, value):
        if name == "out_spacing":
            if len(value) != self.inputs.dimension:
                msg = "out_spacing dimensions should match dimension"
                raise ValueError(msg)

            value = " ".join(["%d" % d for d in value])

        return super(ResampleImageBySpacing, self)._format_arg(name, trait_spec, value)


class SetDirectionByMatrixInputSpec(ANTSCommandInputSpec):
    """InputSpec for ``SetDirectionByMatrix``."""

    infile = base.File(
        position=2,
        argstr="%s",
        name_source=["infile"],
        desc="image to copy header to",
        exists=True,
        mandatory=True,
    )
    outfile = base.File(
        position=3,
        argstr="%s",
        name_sources=["infile", "outfile"],
        desc="output image file",
        usedefault=True,
    )
    direction = traits.String(argstr="%s", position=4, desc="dimensions, x-delimited")


class SetDirectionByMatrixOutputSpec(base.TraitedSpec):
    """OutputSpec for ``SetDirectionByMatrix``."""

    outfile = base.File(exists=True, desc="output image file")


class SetDirectionByMatrix(ANTSCommand):
    """Set image header information from a matrix of dimensions."""

    _cmd = "SetDirectionByMatrix"
    # pylint: disable=protected-access
    _gen_filename = base.StdOutCommandLine._gen_filename
    input_spec = SetDirectionByMatrixInputSpec
    output_spec = SetDirectionByMatrixOutputSpec

    def _format_arg(self, name, trait_spec, value):
        if name == "direction":
            return value.replace("x", " ")
        return super()._format_arg(name, trait_spec, value)

    def _list_outputs(self):
        return {"outfile": self.inputs.outfile}


class ThresholdImageInputSpec(ANTSCommandInputSpec):
    """InputSpec for ThresholdImage."""

    dimension = traits.Int(
        3, usedefault=True, position=1, argstr="%d", desc="dimension of output image"
    )
    input_image = base.File(
        exists=True, mandatory=True, position=2, argstr="%s", desc="input image file"
    )
    output_image = base.File(
        position=3,
        argstr="%s",
        name_source=["input_image"],
        name_template="%s_resampled",
        desc="output image file",
        keep_extension=True,
    )

    mode = traits.Enum(
        "Otsu",
        "Kmeans",
        argstr="%s",
        position=4,
        requires=["num_thresholds"],
        xor=["th_low", "th_high"],
        desc="whether to run Otsu / Kmeans thresholding",
    )
    num_thresholds = traits.Int(position=5, argstr="%d", desc="number of thresholds")
    input_mask = base.File(
        exists=True,
        requires=["num_thresholds"],
        argstr="%s",
        desc="input mask for Otsu, Kmeans",
    )

    th_low = traits.Float(position=4, argstr="%f", xor=["mode"], desc="lower threshold")
    th_high = traits.Float(
        position=5, argstr="%f", xor=["mode"], desc="upper threshold"
    )
    inside_value = traits.Float(
        1, position=6, argstr="%f", requires=["th_low"], desc="inside value"
    )
    outside_value = traits.Float(
        0, position=7, argstr="%f", requires=["th_low"], desc="outside value"
    )


class ThresholdImageOutputSpec(base.TraitedSpec):
    """OutputSpec for ThresholdImage."""

    output_image = traits.File(exists=True, desc="resampled file")


class ThresholdImage(ANTSCommand):
    """Apply thresholds on images."""

    _cmd = "ThresholdImage"
    input_spec = ThresholdImageInputSpec
    output_spec = ThresholdImageOutputSpec


class AIInputSpec(ANTSCommandInputSpec):
    """InputSpec for AffineInitializer."""

    dimension = traits.Int(
        3, usedefault=True, argstr="-d %d", desc="dimension of output image"
    )
    verbose = traits.Bool(
        False, usedefault=True, argstr="-v %d", desc="enable verbosity"
    )

    fixed_image = traits.File(
        exists=True,
        mandatory=True,
        desc="Image to which the moving_image should be transformed",
    )
    moving_image = traits.File(
        exists=True,
        mandatory=True,
        desc="Image that will be transformed to fixed_image",
    )

    fixed_image_mask = traits.File(exists=True, argstr="-x %s", desc="fixed mage mask")
    moving_image_mask = traits.File(
        exists=True, requires=["fixed_image_mask"], desc="moving mage mask"
    )

    metric_trait = (
        traits.Enum("Mattes", "GC", "MI"),
        traits.Int(32),
        traits.Enum("Regular", "Random", "None"),
        traits.Range(value=0.2, low=0.0, high=1.0),
    )
    metric = traits.Tuple(
        *metric_trait, argstr="-m %s", mandatory=True, desc="the metric(s) to use."
    )

    transform = traits.Tuple(
        traits.Enum("Affine", "Rigid", "Similarity"),
        traits.Range(value=0.1, low=0.0, exclude_low=True),
        argstr="-t %s[%f]",
        usedefault=True,
        desc="Several transform options are available",
    )

    principal_axes = traits.Bool(
        False,
        usedefault=True,
        argstr="-p %d",
        xor=["blobs"],
        desc="align using principal axes",
    )
    search_factor = traits.Tuple(
        traits.Float(20),
        traits.Range(value=0.12, low=0.0, high=1.0),
        usedefault=True,
        argstr="-s [%f,%f]",
        desc="search factor",
    )

    search_grid = traits.Either(
        traits.Tuple(
            traits.Float, traits.Tuple(traits.Float, traits.Float, traits.Float)
        ),
        traits.Tuple(traits.Float, traits.Tuple(traits.Float, traits.Float)),
        argstr="-g %s",
        desc="Translation search grid in mm",
    )

    convergence = traits.Tuple(
        traits.Range(low=1, high=10000, value=10),
        traits.Float(1e-6),
        traits.Range(low=1, high=100, value=10),
        usedefault=True,
        argstr="-c [%d,%f,%d]",
        desc="convergence",
    )

    output_transform = traits.File(
        "initialization.mat", usedefault=True, argstr="-o %s", desc="output file name"
    )


class AIOuputSpec(base.TraitedSpec):
    """OutputSpec for AffineInitializer."""

    output_transform = traits.File(exists=True, desc="output file name")


class AI(ANTSCommand):
    """Replaces ``AffineInitializer``."""

    _cmd = "antsAI"
    input_spec = AIInputSpec
    output_spec = AIOuputSpec

    def _run_interface(self, runtime, correct_return_codes=(0,)):
        runtime = super(AI, self)._run_interface(runtime, correct_return_codes)

        setattr(
            self,
            "_output",
            {
                "output_transform": os.path.join(
                    runtime.cwd, os.path.basename(self.inputs.output_transform)
                )
            },
        )
        return runtime

    def _format_arg(self, opt, spec, val):
        if opt == "metric":
            val = "%s[{fixed_image},{moving_image},%d,%s,%f]" % val
            val = val.format(
                fixed_image=self.inputs.fixed_image,
                moving_image=self.inputs.moving_image,
            )
            return spec.argstr % val

        if opt == "search_grid":
            val1 = "x".join(["%f" % v for v in val[1]])
            fmtval = "[%s]" % ",".join([str(val[0]), val1])
            return spec.argstr % fmtval

        if opt == "fixed_image_mask":
            if isdefined(self.inputs.moving_image_mask):
                return spec.argstr % ("[%s,%s]" % (val, self.inputs.moving_image_mask))

        return super(AI, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        return getattr(self, "_output")
