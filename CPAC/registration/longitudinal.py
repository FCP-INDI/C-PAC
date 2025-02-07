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
# pylint: disable=too-many-lines,ungrouped-imports,wrong-import-order
"""Longitudial registration workflows and utilities."""

from nipype.interfaces.utility import IdentityInterface
from nipype.pipeline.engine import Node, Workflow

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.engine import ResourcePool
from CPAC.registration.utils import apply_transform, compose_ants_warp


def t1w_to_longitudinal_to_template_ants(
    wf: Workflow, strat_pool: ResourcePool, interpolation: str
) -> dict[str, tuple[Node, str]]:
    """Combine and apply T1w transforms from native to longitudinal to template."""
    t1_to_long_to_template = pe.Node(
        IdentityInterface(
            fields=[
                "input_brain",
                "reference_brain",
                "input_head",
                "reference_head",
                "input_mask",
                "reference_mask",
                "transform",
                "interpolation",
            ]
        ),
        name="t1_to_long_to_template_inputspec",
    )
    list_of_xfms = pe.Node(
        IdentityInterface(fields=["T1w_to_longitudinal", "longitudinal_to_template"]),
        "list_of_xfms",
    )
    wf.connect(
        *strat_pool.node_data("from-T1w_to-longitudinal_mode-image_xfm"),
        list_of_xfms,
        "T1w_to_longitudinal",
    )
    wf.connect(
        *strat_pool.node_data("from-longitudinal_to-template_mode-image_xfm"),
        list_of_xfms,
        "longitudinal_to_template",
    )
    composite_xfm = compose_ants_warp(
        wf=wf,
        name="T1wLinearTemplate_xfm",
        input_node=t1_to_long_to_template,
        warp_from="input_brain",
        warp_to="reference_brain",
        inputs=[(list_of_xfms, ["longitudinal_to_template", "T1w_to_longitudinal"])],
    )
    list_of_invxfms = pe.Node(
        IdentityInterface(fields=["longitudinal_to_T1w", "template_to_longitudinal"]),
        "list_of_inv_xfms",
    )
    wf.connect(
        *strat_pool.node_data("from-longitudinal_to-T1w_mode-image_xfm"),
        list_of_invxfms,
        "longitudinal_to_T1w",
    )
    wf.connect(
        *strat_pool.node_data("from-template_to-template_mode-image_xfm"),
        list_of_invxfms,
        "template_to_longitudinal",
    )
    inverse_composite_xfm = compose_ants_warp(
        wf=wf,
        name="T1wLinearTemplateInv_xfm",
        input_node=t1_to_long_to_template,
        warp_from="reference_brain",
        warp_to="input_brain",
        inputs=[(list_of_invxfms, ["longitudinal_to_T1w", "template_to_longitudinal"])],
    )
    preproc = apply_transform("warp_t1_to_longitudinal_to_template", "ants")
    preproc.inputs.inputspec.interpolation = interpolation
    wf.connect(
        [
            (
                t1_to_long_to_template,
                preproc,
                [
                    ("input_brain", "inputspec.input_image"),
                    ("reference_brain", "reference"),
                ],
            ),
            (composite_xfm, preproc, [("output_image", "transform")]),
        ]
    )
    outputs: dict[str, tuple[Node | Workflow, str]] = {
        "space-template_desc-preproc_T1w": (preproc, "output_image"),
        # "from-T1w_to-template_mode-image_desc-linear_xfm": None,
        # "from-template_to-T1w_mode-image_desc-linear_xfm": None,
        # "from-T1w_to-template_mode-image_desc-nonlinear_xfm": None,
        # "from-template_to-T1w_mode-image_desc-nonlinear_xfm": None,
        "from-T1w_to-template_mode-image_xfm": (composite_xfm, "output_image"),
        "from-template_to-T1w_mode-image_xfm": (inverse_composite_xfm, "output_image"),
    }

    return outputs


t1w_to_longitudinal_to_template = {"ants": t1w_to_longitudinal_to_template_ants}
