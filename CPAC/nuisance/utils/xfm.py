# Copyright (C) 2025  C-PAC Developers

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
"""Transformation utilities for nuisance regression."""

from typing import cast, Literal

from nipype.pipeline.engine import Workflow

from CPAC.pipeline.engine import ResourcePool
from CPAC.registration.registration import apply_transform
from CPAC.utils.configuration import Configuration


def transform_bold_mask_to_native(
    wf: Workflow,
    strat_pool: ResourcePool,
    cfg: Configuration,
    pipe_num: int,
    reg_tool: Literal["ants", "fsl"],
) -> tuple[Workflow, str]:
    """Transform a template-space BOLD mask to native space."""
    num_cpus = cast(
        int, cfg["pipeline_setup", "system_config", "max_cores_per_participant"]
    )
    num_ants_cores = cast(
        int, cfg["pipeline_setup", "system_config", "num_ants_threads"]
    )
    apply_xfm = apply_transform(
        f"xfm_from-template_to-bold_mask_{pipe_num}",
        reg_tool,
        time_series=True,
        num_cpus=num_cpus,
        num_ants_cores=num_ants_cores,
    )
    apply_xfm.inputs.inputspec.interpolation = cfg[
        "registration_workflows",
        "functional_registration",
        "func_registration_to_template",
        f"{'ANTs' if reg_tool == 'ants' else 'FNIRT'}_pipelines",
        "interpolation",
    ]
    sbref = strat_pool.node_data("sbref")
    bold_mask = strat_pool.node_data(
        ["space-template_desc-bold_mask", "space-template_desc-brain_mask"]
    )
    xfm = strat_pool.node_data("from-template_to-bold_mode-image_xfm")
    wf.connect(
        [
            (bold_mask.node, apply_xfm, [(bold_mask.out, "inputspec.input_image")]),
            (sbref.node, apply_xfm, [(sbref.out, "inputspec.reference")]),
            (xfm.node, apply_xfm, [(xfm.out, "inputspec.transform")]),
        ]
    )

    return apply_xfm, "outputspec.output_image"
