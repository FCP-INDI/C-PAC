# -*- coding: utf-8 -*-
# Copyright (C) 2024-2025  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Create longitudinal template using ``mri_robust_template``."""

import os
from typing import cast, Literal

from nipype.interfaces.base import (
    File,
    InputMultiPath,
    isdefined,
    OutputMultiPath,
    traits,
)
from nipype.interfaces.freesurfer import longitudinal
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nipype.interfaces.freesurfer.utils import LTAConvert

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.configuration import Configuration


class RobustTemplateInputSpec(longitudinal.RobustTemplateInputSpec):  # noqa: D101
    affine = traits.Bool(default_value=False, desc="compute 12 DOF registration")
    mapmov = traits.Either(
        InputMultiPath(File(exists=False)),
        traits.Bool,
        argstr="--mapmov %s",
        desc="output images: map and resample each input to template",
    )
    maxit = traits.Int(
        argstr="--maxit %d",
        mandatory=False,
        desc="iterate max # times (if #tp>2 default 6, else 5 for 2tp reg.)",
    )


class RobustTemplateOutputSpec(longitudinal.RobustTemplateOutputSpec):  # noqa: D101
    mapmov = OutputMultiPath(
        File(),
        desc="each input mapped and resampled to longitudinal template",
    )


class RobustTemplate(longitudinal.RobustTemplate):  # noqa: D101
    # STATEMENT OF CHANGES:
    #     This class is derived from sources licensed under the Apache-2.0 terms,
    #     and this class has been changed.

    # CHANGES:
    #     * Added handling for `affine`, `mapmov` and `maxit`.
    #     * Renamed transform outputs.

    # ORIGINAL WORK'S ATTRIBUTION NOTICE:
    #    Copyright (c) 2009-2016, Nipype developers

    #    Licensed under the Apache License, Version 2.0 (the "License");
    #    you may not use this file except in compliance with the License.
    #    You may obtain a copy of the License at

    #        http://www.apache.org/licenses/LICENSE-2.0

    #    Unless required by applicable law or agreed to in writing, software
    #    distributed under the License is distributed on an "AS IS" BASIS,
    #    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    #    See the License for the specific language governing permissions and
    #    limitations under the License.

    #    Prior to release 0.12, Nipype was licensed under a BSD license.

    # Modifications copyright (C) 2024  C-PAC Developers
    input_spec = RobustTemplateInputSpec
    output_spec = RobustTemplateOutputSpec

    def _format_arg(self, name, spec, value):
        if name == "average_metric":
            # return enumeration value
            return spec.argstr % {"mean": 0, "median": 1}[value]
        if name in ("mapmov", "transform_outputs", "scaled_intensity_outputs"):
            value = self._list_outputs()[name]
        return super()._format_arg(name, spec, value)

    def _list_outputs(self):
        """:py:meth:`~nipype.interfaces.freesurfer.RobustTemplate._list_outputs` + `mapmov`."""
        outputs = self.output_spec().get()
        outputs["out_file"] = os.path.abspath(self.inputs.out_file)
        n_files = len(self.inputs.in_files)
        fmt = "{}{:02d}.{}" if n_files > 9 else "{}{:d}.{}"  # noqa: PLR2004
        for key, prefix, ext in [
            ("transform_outputs", "space-longitudinal", "lta"),
            ("scaled_intensity_outputs", "is", "txt"),
            ("mapmov", "space-longitudinal", "nii.gz"),
        ]:
            if isdefined(getattr(self.inputs, key)):
                fnames = getattr(self.inputs, key)
                if fnames is True:
                    fnames = [fmt.format(prefix, i + 1, ext) for i in range(n_files)]
                outputs[key] = [os.path.abspath(x) for x in fnames]
        return outputs


def mri_robust_template(
    name: str, cfg: Configuration, num_sessions: int
) -> pe.Workflow:
    """Return a subworkflow to run `mri_robust_template` with common options.

    Converts transform files to FSL format.
    """
    wf = pe.Workflow(name=name)
    node = pe.Node(
        RobustTemplate(
            affine=cfg["longitudinal_template_generation", "dof"] == 12,  # noqa: PLR2004
            average_metric=cfg["longitudinal_template_generation", "average_method"],
            auto_detect_sensitivity=True,
            mapmov=True,
            out_file=f"{name}.mgz",
            transform_outputs=True,
        ),
        name="mri_robust_template",
    )
    max_iter = cast(
        int | Literal["default"], cfg["longitudinal_template_generation", "max_iter"]
    )
    if isinstance(max_iter, int):
        node.set_input("maxit", max_iter)

    nifti_template = pe.Node(MRIConvert(out_type="niigz"), name="NIfTI-template")
    wf.connect(node, "out_file", nifti_template, "in_file")

    nifti_outputs = pe.MapNode(MRIConvert(), name="NIfTI-mapmov", iterfield=["in_file"])
    wf.connect(node, "mapmov", nifti_outputs, "in_file")
    reorient_outputs = cfg.orientation_node(
        "reorient_longitudinal_template", pe.MapNode
    )
    wf.connect(nifti_outputs, "out_file", reorient_outputs, "in_file")
    reorient_outputs.set_input(
        "out_file", [f"space-longitudinal{i + 1}.nii.gz" for i in range(num_sessions)]
    )

    convert = pe.MapNode(
        LTAConvert(), name="convert-to-FSL", iterfield=["in_lta", "out_fsl"]
    )
    wf.connect(node, "transform_outputs", convert, "in_lta")
    convert.set_input(
        "out_fsl", [f"space-longitudinal{i + 1}.mat" for i in range(num_sessions)]
    )

    return wf
