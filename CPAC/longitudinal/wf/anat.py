# -*- coding: utf-8 -*-
# Copyright (C) 2020-2024  C-PAC Developers

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
"""Longitudinal workflows for anatomical data."""

from typing import cast, Optional

from nipype import config as nipype_config
from nipype.interfaces import fsl
from nipype.interfaces.utility import Merge

from CPAC.longitudinal.preproc import subject_specific_template
from CPAC.longitudinal.robust_template import mri_robust_template
from CPAC.longitudinal.wf.utils import (
    check_creds_path,
    cross_graph_connections,
    cross_pool_resources,
    select_session_node,
)
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.cpac_pipeline import (
    build_anat_preproc_stack,
    build_segmentation_stack,
    build_T1w_registration_stack,
    connect_pipeline,
    initialize_nipype_wf,
)
from CPAC.pipeline.engine import ingress_output_dir, initiate_rpool, ResourcePool
from CPAC.pipeline.nodeblock import nodeblock, NODEBLOCK_RETURN
from CPAC.registration.registration import apply_transform
from CPAC.utils.configuration import Configuration
from CPAC.utils.utils import check_prov_for_regtool


@nodeblock(
    name="mask_T1w_longitudinal_template",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=["desc-brain_T1w"],
    outputs=["space-T1w_desc-brain_mask"],
)
def mask_T1w_longitudinal_template(
    wf: pe.Workflow, cfg, strat_pool, pipe_num, opt=None
) -> NODEBLOCK_RETURN:
    """Create a native-space brain mask for longitudinal template generation."""
    brain_mask = pe.Node(
        interface=fsl.maths.MathsCommand(),
        name=f"longitudinal_anatomical_brain_mask_{pipe_num}",
    )
    brain_mask.inputs.args = "-bin"

    node, out = strat_pool.get_data("desc-brain_T1w")
    wf.connect(node, out, brain_mask, "in_file")

    outputs = {"space-T1w_desc-brain_mask": (brain_mask, "out_file")}

    return wf, outputs


def pick_map(
    file_list: list[list[str]] | list[str], index: str, file_type: str
) -> Optional[str]:
    """Choose a file from a list of files."""
    if isinstance(file_list, list):
        if len(file_list) == 1 and isinstance(file_list[0], list):
            file_list = file_list[0]
        for file_name in file_list:
            assert isinstance(file_name, str)
            if file_name.endswith(f"{file_type}_{index}.nii.gz"):
                return file_name
    return None


@nodeblock(
    name="mask_longitudinal_T1w_brain",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=["space-longitudinal_desc-brain_T1w"],
    outputs=["space-longitudinal_desc-brain_mask"],
)
def mask_longitudinal_T1w_brain(
    wf, cfg, strat_pool, pipe_num, opt=None
) -> NODEBLOCK_RETURN:
    """Create brain mask for longitudinal T1w image."""
    brain_mask = pe.Node(
        interface=fsl.maths.MathsCommand(),
        name=f"longitudinal_T1w_brain_mask_{pipe_num}",
    )
    brain_mask.inputs.args = "-bin"

    node, out = strat_pool.get_data("space-longitudinal_desc-brain_T1w")
    wf.connect(node, out, brain_mask, "in_file")

    outputs = {"space-longitudinal_desc-brain_mask": (brain_mask, "out_file")}

    return (wf, outputs)


@nodeblock(
    name="warp_longitudinal_T1w_to_template",
    config=["longitudinal_template_generation"],
    switch=["run"],
    option_key="using",
    option_val="C-PAC legacy",
    inputs=[
        (
            "space-longitudinal_desc-brain_T1w",
            "from-longitudinal_to-template_mode-image_xfm",
        ),
        "T1w-brain-template",
    ],
    outputs=["space-template_desc-brain_T1w"],
)
def warp_longitudinal_T1w_to_template(
    wf, cfg, strat_pool, pipe_num, opt=None
) -> NODEBLOCK_RETURN:
    """Transform longitudinal T1w images to template space."""
    xfm_prov = strat_pool.get_cpac_provenance(
        "from-longitudinal_to-template_mode-image_xfm"
    )
    reg_tool = check_prov_for_regtool(xfm_prov)
    num_cpus = cfg.pipeline_setup["system_config"]["max_cores_per_participant"]

    num_ants_cores = cfg.pipeline_setup["system_config"]["num_ants_threads"]

    apply_xfm = apply_transform(
        f"warp_longitudinal_to_T1template_{pipe_num}",
        reg_tool,
        time_series=False,
        num_cpus=num_cpus,
        num_ants_cores=num_ants_cores,
    )

    if reg_tool == "ants":
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            "anatomical_registration"
        ]["registration"]["ANTs"]["interpolation"]
    elif reg_tool == "fsl":
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            "anatomical_registration"
        ]["registration"]["FSL-FNIRT"]["interpolation"]

    node, out = strat_pool.get_data("space-longitudinal_desc-brain_T1w")
    wf.connect(node, out, apply_xfm, "inputspec.input_image")

    node, out = strat_pool.get_data("T1w-brain-template")
    wf.connect(node, out, apply_xfm, "inputspec.reference")

    node, out = strat_pool.get_data("from-longitudinal_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, "inputspec.transform")

    outputs = {"space-template_desc-brain_T1w": (apply_xfm, "outputspec.output_image")}

    return (wf, outputs)


@nodeblock(
    name="warp_longitudinal_seg_to_T1w",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=[
        (
            "space-longitudinal_desc-brain_T1w",
            [
                "from-longitudinal_to-T1w_mode-image_desc-linear_xfm",
                "from-T1w_to-longitudinal_mode-image_desc-linear_xfm",
            ],
            "space-longitudinal_label-CSF_mask",
            "space-longitudinal_label-GM_mask",
            "space-longitudinal_label-WM_mask",
            "space-longitudinal_label-CSF_desc-preproc_mask",
            "space-longitudinal_label-GM_desc-preproc_mask",
            "space-longitudinal_label-WM_desc-preproc_mask",
            "space-longitudinal_label-CSF_probseg",
            "space-longitudinal_label-GM_probseg",
            "space-longitudinal_label-WM_probseg",
        ),
        "T1w-brain-template",
    ],
    outputs=[
        "label-CSF_mask",
        "label-GM_mask",
        "label-WM_mask",
        "label-CSF_desc-preproc_mask",
        "label-GM_desc-preproc_mask",
        "label-WM_desc-preproc_mask",
        "label-CSF_probseg",
        "label-GM_probseg",
        "label-WM_probseg",
    ],
)
def warp_longitudinal_seg_to_T1w(
    wf: pe.Workflow,
    cfg: Configuration,
    strat_pool: ResourcePool,
    pipe_num: int,
    opt: Optional[str] = None,
) -> NODEBLOCK_RETURN:
    """Transform anatomical images from longitudinal space template space."""
    if strat_pool.check_rpool("from-longitudinal_to-T1w_mode-image_desc-linear_xfm"):
        xfm_prov = strat_pool.get_cpac_provenance(
            "from-longitudinal_to-T1w_mode-image_desc-linear_xfm"
        )
        reg_tool = check_prov_for_regtool(xfm_prov)
        xfm: tuple[pe.Node, str] = strat_pool.get_data(
            "from-longitudinal_to-T1w_mode-image_desc-linear_xfm"
        )
    else:
        xfm_prov = strat_pool.get_cpac_provenance(
            "from-T1w_to-longitudinal_mode-image_desc-linear_xfm"
        )
        reg_tool = check_prov_for_regtool(xfm_prov)
        # create inverse xfm if we don't have it
        invt = pe.Node(interface=fsl.ConvertXFM(), name="convert_xfm")
        invt.inputs.invert_xfm = True
        wf.connect(
            *strat_pool.get_data("from-T1w_to-longitudinal_mode-image_desc-linear_xfm"),
            invt,
            "in_file",
        )
        xfm = (invt, "out_file")

    num_cpus = cfg.pipeline_setup["system_config"]["max_cores_per_participant"]

    num_ants_cores = cfg.pipeline_setup["system_config"]["num_ants_threads"]

    outputs = {}

    labels = [
        "CSF_mask",
        "CSF_desc-preproc_mask",
        "CSF_probseg",
        "GM_mask",
        "GM_desc-preproc_mask",
        "GM_probseg",
        "WM_mask",
        "WM_desc-preproc_mask",
        "WM_probseg",
    ]

    for label in labels:
        apply_xfm = apply_transform(
            f"warp_longitudinal_seg_to_T1w_{label}_{pipe_num}",
            reg_tool,
            time_series=False,
            num_cpus=num_cpus,
            num_ants_cores=num_ants_cores,
        )

        if reg_tool == "ants":
            apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
                "anatomical_registration"
            ]["registration"]["ANTs"]["interpolation"]
        elif reg_tool == "fsl":
            apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
                "anatomical_registration"
            ]["registration"]["FSL-FNIRT"]["interpolation"]

        node, out = strat_pool.get_data("space-longitudinal_desc-brain_T1w")
        wf.connect(node, out, apply_xfm, "inputspec.input_image")

        node, out = strat_pool.get_data("T1w-brain-template")
        wf.connect(node, out, apply_xfm, "inputspec.reference")

        wf.connect(*xfm, apply_xfm, "inputspec.transform")

        outputs[f"label-{label}"] = (apply_xfm, "outputspec.output_image")

    return (wf, outputs)


def anat_longitudinal_wf(
    subject_id: str, sub_list: list[dict], config: Configuration, dry_run: bool = False
) -> None:
    """
    Create and run longitudinal workflows for anatomical data.

    Parameters
    ----------
    subject_id
        the id of the subject
    sub_list
        a list of sessions for one subject
    config
        a Configuration object containing the information for the participant pipeline
    dry_run
        build graph without running?
    """
    nipype_config.update_config(
        {
            "execution": {
                "crashfile_format": "txt",
                "stop_on_first_crash": config[
                    "pipeline_setup", "system_config", "fail_fast"
                ],
            }
        }
    )
    config["subject_id"] = subject_id
    session_id_list: list[str] = []
    """List of lists for every strategy"""
    session_wfs = {}

    orig_pipe_name: str = config.pipeline_setup["pipeline_name"]

    strats_dct: dict[str, list[tuple[pe.Node, str] | str]] = {
        "desc-brain_T1w": [],
        "desc-head_T1w": [],
    }
    for session in sub_list:
        # Loop over the sessions to create the input for the longitudinal algorithm
        unique_id: str = session["unique_id"]
        session_id_list.append(unique_id)
        input_creds_path = check_creds_path(session.get("creds_path"), subject_id)

        workflow: pe.Workflow = initialize_nipype_wf(
            config,
            subject_id,
            unique_id,
            name="anat_longitudinal_pre-preproc",
        )
        rpool: ResourcePool
        workflow, rpool = initiate_rpool(workflow, config, session)
        pipeline_blocks = build_anat_preproc_stack(rpool, config)
        workflow = connect_pipeline(workflow, config, rpool, pipeline_blocks)

        session_wfs[unique_id] = rpool

        rpool.gather_pipes(workflow, config)

        for key in strats_dct.keys():
            strats_dct[key].append(cast(tuple[pe.Node, str], rpool.get_data(key)))
        if not dry_run:
            workflow.run()
            for key in strats_dct.keys():  # get the outputs from run-nodes
                for index, data in enumerate(list(strats_dct[key])):
                    if isinstance(data, tuple):
                        strats_dct[key][index] = workflow.get_output(*data)

    wf = initialize_nipype_wf(
        config,
        subject_id,
        name="template_node_brain",
    )

    config.pipeline_setup["pipeline_name"] = f"longitudinal_{orig_pipe_name}"

    num_sessions = len(strats_dct["desc-brain_T1w"])
    merge_brains = pe.Node(Merge(num_sessions), name="merge_brains")
    merge_skulls = pe.Node(Merge(num_sessions), name="merge_skulls")
    wf.add_nodes([merge_brains, merge_skulls])
    for i in list(range(0, num_sessions)):
        wf._connect_node_or_path_for_merge(
            merge_brains, strats_dct, "desc-brain_T1w", i
        )
        wf._connect_node_or_path_for_merge(merge_skulls, strats_dct, "desc-head_T1w", i)

    long_id = f"longitudinal_{subject_id}_strat-desc-brain_T1w"

    wf, rpool = initiate_rpool(wf, config, part_id=long_id)

    match config["longitudinal_template_generation", "using"]:
        case "C-PAC legacy":
            brain_output = "brain_template"
            head_output = "skull_template"

            # This node will generate the longitudinal template (the functions are
            # in longitudinal_preproc)
            # Later other algorithms could be added to calculate it, like the
            # multivariate template from ANTS
            # It would just require to change it here.

            # multiple variable names here for compatibility with other options later in this function
            brain_template_node = wholehead_template_node = template_node = (
                subject_specific_template(workflow_name="longitudinal_anat_template")
            )

            template_node.inputs.set(
                avg_method=config.longitudinal_template_generation["average_method"],
                dof=config.longitudinal_template_generation["dof"],
                interp=config.longitudinal_template_generation["legacy-specific"][
                    "interp"
                ],
                cost=config.longitudinal_template_generation["legacy-specific"]["cost"],
                convergence_threshold=config.longitudinal_template_generation[
                    "legacy-specific"
                ]["convergence_threshold"],
                thread_pool=config.longitudinal_template_generation["legacy-specific"][
                    "thread_pool"
                ],
                unique_id_list=list(session_wfs.keys()),
            )

            wf.connect(merge_brains, "out", brain_template_node, "input_brain_list")
            wf.connect(merge_skulls, "out", wholehead_template_node, "input_skull_list")

        case "mri_robust_template":
            brain_output = head_output = "mri_robust_template.out_file"
            brain_template_node = mri_robust_template(
                f"mri_robust_template_brain_{subject_id}", config, len(sub_list)
            )
            wholehead_template_node = mri_robust_template(
                f"mri_robust_template_head_{subject_id}", config, len(sub_list)
            )
            wf.connect(
                merge_brains, "out", brain_template_node, "mri_robust_template.in_files"
            )
            wf.connect(
                merge_brains,
                "out",
                wholehead_template_node,
                "mri_robust_template.in_files",
            )

        case _:
            msg = ": ".join(
                [
                    "Invalid 'using' value for longitudinal template generation",
                    str(config["longitudinal_template_generation", "using"]),
                ]
            )
            raise ValueError(msg)

    for suffix in ["", "-template"]:
        rpool.set_data(
            f"space-longitudinal_desc-brain_T1w{suffix}",
            brain_template_node,
            brain_output,
            {},
            "",
            brain_template_node.name,
        )

        for desc in ["head", "reorient"]:
            rpool.set_data(
                f"space-longitudinal_desc-{desc}_T1w{suffix}",
                wholehead_template_node,
                head_output,
                {},
                "",
                wholehead_template_node.name,
            )

    pipeline_blocks = [mask_longitudinal_T1w_brain]
    pipeline_blocks = build_T1w_registration_stack(
        rpool, config, pipeline_blocks, space="longitudinal"
    )
    pipeline_blocks = build_segmentation_stack(rpool, config, pipeline_blocks)

    rpool.gather_pipes(
        wf,
        config,
        add_excl=[
            "space-longitudinal_desc-brain_T1w",
            "space-longitudinal_desc-reorient_T1w",
            "space-longitudinal_desc-brain_mask",
        ],
    )
    wf = connect_pipeline(wf, config, rpool, pipeline_blocks)
    if not dry_run:
        wf.run()

    # now, just write out a copy of the above to each session
    config.pipeline_setup["pipeline_name"] = orig_pipe_name
    longitudinal_rpool = rpool
    cpr = cross_pool_resources(
        f"fsl_longitudinal_{subject_id}"
    )  # "fsl" for check_prov_for_regtool
    for i, session in enumerate(sub_list):
        unique_id = session["unique_id"]
        input_creds_path = check_creds_path(session.get("creds_path"), subject_id)

        ses_wf = initialize_nipype_wf(config, subject_id, unique_id)

        ses_wf, rpool = initiate_rpool(ses_wf, config, session)
        config.pipeline_setup["pipeline_name"] = f"longitudinal_{orig_pipe_name}"
        if "derivatives_dir" in session:
            ses_wf, rpool = ingress_output_dir(
                ses_wf,
                config,
                rpool,
                long_id,
                data_paths=session,
                part_id=subject_id,
                ses_id=unique_id,
                creds_path=input_creds_path,
            )

        select_sess = select_session_node(unique_id)

        match config["longitudinal_template_generation", "using"]:
            case "C-PAC legacy":
                assert isinstance(brain_template_node, pe.Node)
                for input_name, output_name in [
                    ("output_brains", "output_brain_list"),
                    ("warps", "warp_list"),
                ]:
                    cross_graph_connections(
                        wf,
                        ses_wf,
                        brain_template_node,
                        select_sess,
                        output_name,
                        input_name,
                        dry_run,
                    )

            case "mri_robust_template":
                assert isinstance(brain_template_node, pe.Workflow)
                assert isinstance(wholehead_template_node, pe.Workflow)
                index = i + 1
                head_select_sess = select_session_node(unique_id, "-wholehead")
                select_sess.set_input("session", f"space-longitudinal{index}")
                head_select_sess.set_input("session", f"space-longitudinal{index}")
                for input_name, output_name in [
                    ("output_brains", "mri_robust_template.mapmov"),
                    ("warps", "convert-to-FSL_.out_fsl"),
                ]:
                    cross_graph_connections(
                        wf,
                        ses_wf,
                        brain_template_node,
                        select_sess,
                        output_name,
                        input_name,
                        dry_run,
                    )
                    cross_graph_connections(
                        wf,
                        ses_wf,
                        wholehead_template_node,
                        head_select_sess,
                        output_name,
                        input_name,
                        dry_run,
                    )

                rpool.set_data(
                    "space-longitudinal_desc-head_T1w",
                    head_select_sess,
                    "brain_path",
                    {},
                    "",
                    head_select_sess.name,
                )
                rpool.set_data(
                    "from-T1w_to-longitudinal_mode-image_desc-linear_xfm",
                    head_select_sess,
                    "warp_path",
                    {},
                    "",
                    head_select_sess.name,
                )

        rpool.set_data(
            "space-longitudinal_desc-brain_T1w",
            select_sess,
            "brain_path",
            {},
            "",
            select_sess.name,
        )
        rpool.set_data(
            "from-T1w_to-longitudinal_mode-image_desc-linear_xfm",
            select_sess,
            "warp_path",
            {},
            "",
            select_sess.name,
        )

        config.pipeline_setup["pipeline_name"] = orig_pipe_name
        excl = ["space-template_desc-brain_T1w", "space-T1w_desc-brain_mask"]
        rpool.gather_pipes(ses_wf, config, add_excl=excl)
        cross_pool_keys = ["from-longitudinal_to-template_mode-image_xfm"]
        for key in cross_pool_keys:
            node, out = longitudinal_rpool.get_data(key)
            cross_graph_connections(wf, ses_wf, node, cpr, out, key, dry_run)
            rpool.set_data(key, cpr, key, {}, "", cpr.name)
        if not dry_run:
            ses_wf.run()

        pipeline_blocks = [
            warp_longitudinal_T1w_to_template,
            warp_longitudinal_seg_to_T1w,
        ]

        ses_wf = connect_pipeline(ses_wf, config, rpool, pipeline_blocks)

        rpool.gather_pipes(ses_wf, config)

        # this is going to run multiple times!
        # once for every strategy!
        if not dry_run:
            ses_wf.run()
