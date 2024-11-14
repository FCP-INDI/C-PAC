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
import os

from nipype.interfaces import fsl
from indi_aws import aws_utils

from CPAC.longitudinal.preproc import subject_specific_template
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.cpac_pipeline import (
    build_anat_preproc_stack,
    build_segmentation_stack,
    build_T1w_registration_stack,
    connect_pipeline,
    initialize_nipype_wf,
)
from CPAC.pipeline.engine import ingress_output_dir, initiate_rpool
from CPAC.pipeline.nodeblock import nodeblock
from CPAC.registration.registration import apply_transform
from CPAC.utils.configuration.configuration import Configuration
from CPAC.utils.interfaces.datasink import DataSink
from CPAC.utils.interfaces.function import Function
from CPAC.utils.utils import check_prov_for_regtool


@nodeblock(
    name="mask_T1w_longitudinal_template",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=["desc-brain_T1w"],
    outputs=["space-T1w_desc-brain_mask"],
)
def mask_T1w_longitudinal_template(wf, cfg, strat_pool, pipe_num, opt=None):
    brain_mask = pe.Node(
        interface=fsl.maths.MathsCommand(),
        name=f"longitudinal_anatomical_brain_mask_{pipe_num}",
    )
    brain_mask.inputs.args = "-bin"

    node, out = strat_pool.get_data("desc-brain_T1w")
    wf.connect(node, out, brain_mask, "in_file")

    outputs = {"space-T1w_desc-brain_mask": (brain_mask, "out_file")}

    return (wf, outputs)


def create_datasink(
    datasink_name,
    config,
    subject_id,
    session_id="",
    strat_name="",
    map_node_iterfield=None,
) -> pe.Node | pe.MapNode:
    """
    Parameters
    ----------
    datasink_name
    config
    subject_id
    session_id
    strat_name
    map_node_iterfield
    """
    encrypt_data = config.pipeline_setup["Amazon-AWS"]["s3_encryption"]

    # TODO Enforce value with schema validation
    # Extract credentials path for output if it exists
    try:
        # Get path to creds file
        creds_path = ""
        if config.pipeline_setup["Amazon-AWS"]["aws_output_bucket_credentials"]:
            creds_path = str(
                config.pipeline_setup["Amazon-AWS"]["aws_output_bucket_credentials"]
            )
            creds_path = os.path.abspath(creds_path)

        if (
            config.pipeline_setup["output_directory"]["path"]
            .lower()
            .startswith("s3://")
        ):
            # Test for s3 write access
            s3_write_access = aws_utils.test_bucket_access(
                creds_path, config.pipeline_setup["output_directory"]["path"]
            )

            if not s3_write_access:
                msg = "Not able to write to bucket!"
                raise Exception(msg)

    except Exception as e:
        if (
            config.pipeline_setup["output_directory"]["path"]
            .lower()
            .startswith("s3://")
        ):
            err_msg = (
                "There was an error processing credentials or "
                "accessing the S3 bucket. Check and try again.\n"
                "Error: %s" % e
            )
            raise Exception(err_msg)

    if map_node_iterfield is not None:
        ds = pe.MapNode(
            DataSink(infields=map_node_iterfield),
            name=f"sinker_{datasink_name}",
            iterfield=map_node_iterfield,
        )
    else:
        ds = pe.Node(DataSink(), name=f"sinker_{datasink_name}")

    ds.inputs.base_directory = config.pipeline_setup["output_directory"]["path"]
    ds.inputs.creds_path = creds_path
    ds.inputs.encrypt_bucket_keys = encrypt_data
    ds.inputs.container = os.path.join(
        "pipeline_%s_%s" % (config.pipeline_setup["pipeline_name"], strat_name),
        subject_id,
        session_id,
    )
    return ds


def connect_anat_preproc_inputs(
    strat, anat_preproc, strat_name, strat_nodes_list_list, workflow
):
    """
    Parameters
    ----------
    strat : Strategy
        the strategy object you want to fork
    anat_preproc : Workflow
        the anat_preproc workflow node to be connected and added to the resource pool
    strat_name : str
        name of the strategy
    strat_nodes_list_list : list
        a list of strat_nodes_list
    workflow : Workflow
        main longitudinal workflow

    Returns
    -------
    new_strat : Strategy
        the fork of strat with the resource pool updated
    strat_nodes_list_list : list
        a list of strat_nodes_list
    """
    new_strat = strat.fork()

    tmp_node, out_key = new_strat["anatomical"]
    workflow.connect(tmp_node, out_key, anat_preproc, "inputspec.anat")

    tmp_node, out_key = new_strat["template_cmass"]
    workflow.connect(tmp_node, out_key, anat_preproc, "inputspec.template_cmass")

    new_strat.append_name(anat_preproc.name)

    new_strat.update_resource_pool(
        {
            "anatomical_brain": (anat_preproc, "outputspec.brain"),
            "anatomical_skull_leaf": (anat_preproc, "outputspec.reorient"),
            "anatomical_brain_mask": (anat_preproc, "outputspec.brain_mask"),
        }
    )

    try:
        strat_nodes_list_list[strat_name].append(new_strat)
    except KeyError:
        strat_nodes_list_list[strat_name] = [new_strat]

    return new_strat, strat_nodes_list_list


def pick_map(file_list, index, file_type):
    if isinstance(file_list, list):
        if len(file_list) == 1:
            file_list = file_list[0]
        for file_name in file_list:
            if file_name.endswith(f"{file_type}_{index}.nii.gz"):
                return file_name
    return None


def select_session(session, output_brains, warps):
    brain_path = None
    warp_path = None
    for brain_path in output_brains:
        if f"{session}_" in brain_path:
            break
    for warp_path in warps:
        if f"{session}_" in warp_path:
            break
    return (brain_path, warp_path)


@nodeblock(
    name="mask_longitudinal_T1w_brain",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=["space-longitudinal_desc-brain_T1w"],
    outputs=["space-longitudinal_desc-brain_mask"],
)
def mask_longitudinal_T1w_brain(wf, cfg, strat_pool, pipe_num, opt=None):
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
    inputs=[
        (
            "space-longitudinal_desc-brain_T1w",
            "from-longitudinal_to-template_mode-image_xfm",
        )
    ],
    outputs=["space-template_desc-brain_T1w"],
)
def warp_longitudinal_T1w_to_template(wf, cfg, strat_pool, pipe_num, opt=None):
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

    node, out = strat_pool.get_data("T1w_brain_template")
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
            "from-longitudinal_to-T1w_mode-image_desc-linear_xfm",
            "space-longitudinal_label-CSF_mask",
            "space-longitudinal_label-GM_mask",
            "space-longitudinal_label-WM_mask",
            "space-longitudinal_label-CSF_desc-preproc_mask",
            "space-longitudinal_label-GM_desc-preproc_mask",
            "space-longitudinal_label-WM_desc-preproc_mask",
            "space-longitudinal_label-CSF_probseg",
            "space-longitudinal_label-GM_probseg",
            "space-longitudinal_label-WM_probseg",
        )
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
def warp_longitudinal_seg_to_T1w(wf, cfg, strat_pool, pipe_num, opt=None):
    xfm_prov = strat_pool.get_cpac_provenance(
        "from-longitudinal_to-T1w_mode-image_desc-linear_xfm"
    )
    reg_tool = check_prov_for_regtool(xfm_prov)

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

        node, out = strat_pool.get_data("T1w_brain_template")
        wf.connect(node, out, apply_xfm, "inputspec.reference")

        node, out = strat_pool.get_data("from-longitudinal_to-template_mode-image_xfm")
        wf.connect(node, out, apply_xfm, "inputspec.transform")

        outputs[f"label-{label}"] = (apply_xfm, "outputspec.output_image")

    return (wf, outputs)


def anat_longitudinal_wf(
    subject_id: str, sub_list: list[dict], config: Configuration
) -> None:
    """
    Create and run anatomical longitudinal workflow(s).

    Parameters
    ----------
    subject_id
        the id of the subject
    sub_list
        this is a list of sessions for one subject and each session if the same dictionary as the one given to
        prep_workflow
    config
        a configuration object containing the information of the pipeline config. (Same as for prep_workflow)
    """
    config["subject_id"] = subject_id
    session_id_list: list[list] = []
    """List of lists for every strategy"""
    session_wfs = {}

    cpac_dirs = []
    out_dir = config.pipeline_setup["output_directory"]["path"]

    orig_pipe_name = config.pipeline_setup["pipeline_name"]

    # Loop over the sessions to create the input for the longitudinal
    # algorithm
    for session in sub_list:
        unique_id = session["unique_id"]
        session_id_list.append(unique_id)

        try:
            creds_path = session["creds_path"]
            if creds_path and "none" not in creds_path.lower():
                if os.path.exists(creds_path):
                    input_creds_path = os.path.abspath(creds_path)
                else:
                    err_msg = (
                        'Credentials path: "%s" for subject "%s" '
                        'session "%s" was not found. Check this path '
                        "and try again." % (creds_path, subject_id, unique_id)
                    )
                    raise Exception(err_msg)
            else:
                input_creds_path = None
        except KeyError:
            input_creds_path = None

        workflow = initialize_nipype_wf(
            config,
            sub_list[0],
            # just grab the first one for the name
            name="anat_longitudinal_pre-preproc",
        )

        workflow, rpool = initiate_rpool(workflow, config, session)
        pipeline_blocks = build_anat_preproc_stack(rpool, config)
        workflow = connect_pipeline(workflow, config, rpool, pipeline_blocks)

        session_wfs[unique_id] = rpool

        rpool.gather_pipes(workflow, config)

        workflow.run()

        cpac_dir = os.path.join(
            out_dir, f"pipeline_{orig_pipe_name}", f"{subject_id}_{unique_id}"
        )
        cpac_dirs.append(os.path.join(cpac_dir, "anat"))

    # Now we have all the anat_preproc set up for every session
    # loop over the different anat preproc strategies
    strats_brain_dct = {}
    strats_head_dct = {}
    for cpac_dir in cpac_dirs:
        if os.path.isdir(cpac_dir):
            for filename in os.listdir(cpac_dir):
                if "T1w.nii" in filename:
                    for tag in filename.split("_"):
                        if "desc-" in tag and "brain" in tag:
                            if tag not in strats_brain_dct:
                                strats_brain_dct[tag] = []
                            strats_brain_dct[tag].append(
                                os.path.join(cpac_dir, filename)
                            )
                            if tag not in strats_head_dct:
                                strats_head_dct[tag] = []
                            head_file = filename.replace(tag, "desc-reorient")
                            strats_head_dct[tag].append(
                                os.path.join(cpac_dir, head_file)
                            )

    for strat in strats_brain_dct.keys():
        wf = initialize_nipype_wf(
            config,
            sub_list[0],
            # just grab the first one for the name
            name=f"template_node_{strat}",
        )

        config.pipeline_setup["pipeline_name"] = f"longitudinal_{orig_pipe_name}"

        template_node_name = f"longitudinal_anat_template_{strat}"

        # This node will generate the longitudinal template (the functions are
        # in longitudinal_preproc)
        # Later other algorithms could be added to calculate it, like the
        # multivariate template from ANTS
        # It would just require to change it here.
        template_node = subject_specific_template(workflow_name=template_node_name)

        template_node.inputs.set(
            avg_method=config.longitudinal_template_generation["average_method"],
            dof=config.longitudinal_template_generation["dof"],
            interp=config.longitudinal_template_generation["interp"],
            cost=config.longitudinal_template_generation["cost"],
            convergence_threshold=config.longitudinal_template_generation[
                "convergence_threshold"
            ],
            thread_pool=config.longitudinal_template_generation["thread_pool"],
            unique_id_list=list(session_wfs.keys()),
        )

        template_node.inputs.input_brain_list = strats_brain_dct[strat]
        template_node.inputs.input_skull_list = strats_head_dct[strat]

        long_id = f"longitudinal_{subject_id}_strat-{strat}"

        wf, rpool = initiate_rpool(wf, config, part_id=long_id)

        rpool.set_data(
            "space-longitudinal_desc-brain_T1w",
            template_node,
            "brain_template",
            {},
            "",
            template_node_name,
        )

        rpool.set_data(
            "space-longitudinal_desc-brain_T1w-template",
            template_node,
            "brain_template",
            {},
            "",
            template_node_name,
        )

        rpool.set_data(
            "space-longitudinal_desc-reorient_T1w",
            template_node,
            "skull_template",
            {},
            "",
            template_node_name,
        )

        rpool.set_data(
            "space-longitudinal_desc-reorient_T1w-template",
            template_node,
            "skull_template",
            {},
            "",
            template_node_name,
        )

        pipeline_blocks = [mask_longitudinal_T1w_brain]

        pipeline_blocks = build_T1w_registration_stack(rpool, config, pipeline_blocks)

        pipeline_blocks = build_segmentation_stack(rpool, config, pipeline_blocks)

        wf = connect_pipeline(wf, config, rpool, pipeline_blocks)

        excl = [
            "space-longitudinal_desc-brain_T1w",
            "space-longitudinal_desc-reorient_T1w",
            "space-longitudinal_desc-brain_mask",
        ]
        rpool.gather_pipes(wf, config, add_excl=excl)

        # this is going to run multiple times!
        # once for every strategy!
        wf.run()

        # now, just write out a copy of the above to each session
        config.pipeline_setup["pipeline_name"] = orig_pipe_name
        for session in sub_list:
            unique_id = session["unique_id"]

            try:
                creds_path = session["creds_path"]
                if creds_path and "none" not in creds_path.lower():
                    if os.path.exists(creds_path):
                        input_creds_path = os.path.abspath(creds_path)
                    else:
                        err_msg = (
                            'Credentials path: "%s" for subject "%s" '
                            'session "%s" was not found. Check this path '
                            "and try again." % (creds_path, subject_id, unique_id)
                        )
                        raise Exception(err_msg)
                else:
                    input_creds_path = None
            except KeyError:
                input_creds_path = None

            wf = initialize_nipype_wf(config, sub_list[0])

            wf, rpool = initiate_rpool(wf, config, session)

            config.pipeline_setup["pipeline_name"] = f"longitudinal_{orig_pipe_name}"
            rpool = ingress_output_dir(
                config, rpool, long_id, creds_path=input_creds_path
            )

            select_node_name = f"select_{unique_id}"
            select_sess = pe.Node(
                Function(
                    input_names=["session", "output_brains", "warps"],
                    output_names=["brain_path", "warp_path"],
                    function=select_session,
                ),
                name=select_node_name,
            )
            select_sess.inputs.session = unique_id

            wf.connect(template_node, "output_brain_list", select_sess, "output_brains")
            wf.connect(template_node, "warp_list", select_sess, "warps")

            rpool.set_data(
                "space-longitudinal_desc-brain_T1w",
                select_sess,
                "brain_path",
                {},
                "",
                select_node_name,
            )

            rpool.set_data(
                "from-T1w_to-longitudinal_mode-image_desc-linear_xfm",
                select_sess,
                "warp_path",
                {},
                "",
                select_node_name,
            )

            config.pipeline_setup["pipeline_name"] = orig_pipe_name
            excl = ["space-template_desc-brain_T1w", "space-T1w_desc-brain_mask"]

            rpool.gather_pipes(wf, config, add_excl=excl)
            wf.run()

    # begin single-session stuff again
    for session in sub_list:
        unique_id = session["unique_id"]

        try:
            creds_path = session["creds_path"]
            if creds_path and "none" not in creds_path.lower():
                if os.path.exists(creds_path):
                    input_creds_path = os.path.abspath(creds_path)
                else:
                    err_msg = (
                        'Credentials path: "%s" for subject "%s" '
                        'session "%s" was not found. Check this path '
                        "and try again." % (creds_path, subject_id, unique_id)
                    )
                    raise Exception(err_msg)
            else:
                input_creds_path = None
        except KeyError:
            input_creds_path = None

        wf = initialize_nipype_wf(config, sub_list[0])

        wf, rpool = initiate_rpool(wf, config, session)

        pipeline_blocks = [
            warp_longitudinal_T1w_to_template,
            warp_longitudinal_seg_to_T1w,
        ]

        wf = connect_pipeline(wf, config, rpool, pipeline_blocks)

        rpool.gather_pipes(wf, config)

        # this is going to run multiple times!
        # once for every strategy!
        wf.run()
