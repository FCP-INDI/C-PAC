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
from CPAC.registration.utils import (
    apply_transform,
    CommonRegistrationInputs,
    compose_ants_warp,
    RegistrationTemplates,
)


def get_common_reg_inputs(
    strat_pool: ResourcePool,
    registration_templates: RegistrationTemplates = RegistrationTemplates(),
) -> CommonRegistrationInputs:
    """Get common longitudinal registration inputs."""
    orig = "longitudinal"
    has_longitudinal = strat_pool.check_rpool(
        "longitudinal-template_space-longitudinal_desc-brain_T1w"
    )
    is_longitudinal = False
    if has_longitudinal and not strat_pool.check_rpool("desc-preproc_T1w"):
        orig = "longitudinal"
        has_longitudinal = False
        is_longitudinal = True
        input_brain = strat_pool.node_data(
            "longitudinal-template_space-longitudinal_desc-brain_T1w"
        )
        input_head = strat_pool.node_data(
            "longitudinal-template_space-longitudinal_desc-head_T1w"
        )
        reference_mask = strat_pool.node_data(
            "longitudinal-template_space-longitudinal_desc-brain_mask"
        )
        lesion_mask = None
    else:
        orig = "T1w"
        input_brain = strat_pool.node_data("desc-preproc_T1w")
        input_head = strat_pool.node_data(
            [  # TODO: check the order of T1w
                "desc-restore_T1w",
                "desc-head_T1w",
                "desc-preproc_T1w",
            ]
        )
        reference_mask = (
            strat_pool.node_data(registration_templates.reference_mask)
            if strat_pool.check_rpool(registration_templates.reference_mask)
            else None
        )
        lesion_mask = (
            strat_pool.node_data("label-lesion_mask")
            if strat_pool.check_rpool("label-lesion_mask")
            else None
        )
    if has_longitudinal:
        t1w_brain_template = strat_pool.node_data(
            "longitudinal-template_space-longitudinal_desc-brain_T1w"
        )
        t1w_template = strat_pool.node_data(
            "longitudinal-template_space-longitudinal_desc-head_T1w"
        )
    else:
        t1w_brain_template = strat_pool.node_data(
            registration_templates.reference_brain
        )
        t1w_template = strat_pool.node_data(registration_templates.reference_head)
    if is_longitudinal:
        brain_mask = strat_pool.node_data(
            "longitudinal-template_space-longitudinal_desc-brain_mask"
        )
    else:
        brain_mask = strat_pool.node_data(
            [
                "space-T1w_desc-brain_mask",
                "space-T1w_desc-acpcbrain_mask",
            ]
        )
    return CommonRegistrationInputs(
        orig,
        has_longitudinal,
        input_brain,
        input_head,
        reference_mask,
        lesion_mask,
        t1w_brain_template,
        t1w_template,
        brain_mask,
    )


def t1w_to_longitudinal_to_template_ants(
    wf: Workflow, strat_pool: ResourcePool, template: str = "template"
) -> dict[str, tuple[Node | Workflow, str]]:
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
            ],
            interpolation=strat_pool.ants_interp,
        ),
        name=f"t1_to_long_to_{template}_inputspec",
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
        *strat_pool.node_data(f"from-longitudinal_to-{template}_mode-image_xfm"),
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
        *strat_pool.node_data(f"from-{template}_to-longitudinal_mode-image_xfm"),
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
    preproc = apply_transform(f"warp_t1_to_longitudinal_to_{template}", "ants")
    preproc.inputs.inputspec.interpolation = strat_pool.ants_interp
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
        f"space-{template}_desc-preproc_T1w": (preproc, "output_image"),
        # "from-T1w_to-template_mode-image_desc-linear_xfm": None,
        # "from-template_to-T1w_mode-image_desc-linear_xfm": None,
        # "from-T1w_to-template_mode-image_desc-nonlinear_xfm": None,
        # "from-template_to-T1w_mode-image_desc-nonlinear_xfm": None,
        f"from-T1w_to-{template}_mode-image_xfm": (composite_xfm, "output_image"),
        f"from-{template}_to-T1w_mode-image_xfm": (
            inverse_composite_xfm,
            "output_image",
        ),
    }

    return outputs


t1w_to_longitudinal_to_template = {"ants": t1w_to_longitudinal_to_template_ants}
