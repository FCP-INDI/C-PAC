# -*- coding: utf-8 -*-
# Copyright (C) 2020-2025  C-PAC Developers

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

import os

from nipype.interfaces import fsl

from CPAC.longitudinal.preproc import subject_specific_template
from CPAC.longitudinal.wf.utils import select_session
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.cpac_pipeline import (
    build_anat_preproc_stack,
    build_segmentation_stack,
    build_T1w_registration_stack,
    connect_pipeline,
    initialize_nipype_wf,
)
from CPAC.pipeline.engine import ingress_output_dir, initiate_rpool
from CPAC.pipeline.nodeblock import nodeblock, NODEBLOCK_RETURN
from CPAC.registration.registration import apply_transform
from CPAC.utils.configuration import Configuration
from CPAC.utils.interfaces.function import Function
from CPAC.utils.utils import check_prov_for_regtool


@nodeblock(
    name="mask_longitudinal_T1w_brain",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=["space-longitudinal_desc-preproc_T1w"],
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

    node, out = strat_pool.get_data("space-longitudinal_desc-preproc_T1w")
    wf.connect(node, out, brain_mask, "in_file")

    outputs = {"space-longitudinal_desc-brain_mask": (brain_mask, "out_file")}

    return wf, outputs


@nodeblock(
    name="warp_longitudinal_T1w_to_template",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=[
        (
            "space-longitudinal_desc-preproc_T1w",
            "from-longitudinal_to-template_mode-image_xfm",
            "T1w-brain-template",
        )
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

    node, out = strat_pool.get_data("space-longitudinal_desc-preproc_T1w")
    wf.connect(node, out, apply_xfm, "inputspec.input_image")

    node, out = strat_pool.get_data("T1w-brain-template")
    wf.connect(node, out, apply_xfm, "inputspec.reference")

    node, out = strat_pool.get_data("from-longitudinal_to-template_mode-image_xfm")
    wf.connect(node, out, apply_xfm, "inputspec.transform")

    outputs = {"space-template_desc-brain_T1w": (apply_xfm, "outputspec.output_image")}

    return wf, outputs


@nodeblock(
    name="warp_longitudinal_seg_to_T1w",
    config=["longitudinal_template_generation"],
    switch=["run"],
    inputs=[
        (
            "from-longitudinal_to-symtemplate_mode-image_xfm",
            "space-longitudinal_label-CSF_mask",
            "space-longitudinal_label-GM_mask",
            "space-longitudinal_label-WM_mask",
            "space-longitudinal_label-CSF_desc-preproc_mask",
            "space-longitudinal_label-GM_desc-preproc_mask",
            "space-longitudinal_label-WM_desc-preproc_mask",
            "space-longitudinal_label-CSF_probseg",
            "space-longitudinal_label-GM_probseg",
            "space-longitudinal_label-WM_probseg",
            "space-longitudinal_desc-preproc_T1w",
            "T1w-brain-template",
            "from-longitudinal_to-template_mode-image_xfm",
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
def warp_longitudinal_seg_to_T1w(
    wf, cfg, strat_pool, pipe_num, opt=None
) -> NODEBLOCK_RETURN:
    """Transform anatomical segmentation from longitudinal template to T1w space."""
    xfm_prov = strat_pool.get_cpac_provenance(
        "from-longitudinal_to-symtemplate_mode-image_xfm"
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

        node, out = strat_pool.get_data("space-longitudinal_desc-preproc_T1w")
        wf.connect(node, out, apply_xfm, "inputspec.input_image")

        node, out = strat_pool.get_data("T1w-brain-template")
        wf.connect(node, out, apply_xfm, "inputspec.reference")

        node, out = strat_pool.get_data("from-longitudinal_to-template_mode-image_xfm")
        wf.connect(node, out, apply_xfm, "inputspec.transform")

        outputs[f"label-{label}"] = (apply_xfm, "outputspec.output_image")

    return wf, outputs


def anat_longitudinal_wf(
    subject_id: str, sub_list: list[dict], config: Configuration
) -> None:
    """
    Create and run longitudinal workflows for anatomical data.

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
    # list of lists for every strategy
    session_id_list = []
    session_wfs = {}

    cpac_dirs = []
    out_dir = config.pipeline_setup["output_directory"]["path"]
    config.__setitem__("subject_id", subject_id)
    orig_pipe_name = config.pipeline_setup["pipeline_name"]

    # Loop over the sessions to create the input for the longitudinal
    # algorithm
    for i in range(len(sub_list)):
        unique_id = sub_list[i]["unique_id"]
        session_id_list.append(unique_id)

        try:
            creds_path = sub_list[i]["creds_path"]
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
            sub_list[i],
            # just grab the first one for the name
            name="anat_longitudinal_pre-preproc",
        )

        workflow, rpool = initiate_rpool(workflow, config, sub_list[i])
        pipeline_blocks = build_anat_preproc_stack(rpool, config)
        workflow = connect_pipeline(workflow, config, rpool, pipeline_blocks)

        session_wfs[unique_id] = rpool

        rpool.gather_pipes(workflow, config)

        workflow.run()

        cpac_dir = os.path.join(
            out_dir, f"pipeline_{orig_pipe_name}", subject_id, unique_id
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
                        if "desc-" in tag and "preproc" in tag:
                            if tag not in strats_brain_dct:
                                strats_brain_dct[tag] = []
                            strats_brain_dct[tag].append(
                                os.path.join(cpac_dir, filename)
                            )
                            if tag not in strats_head_dct:
                                strats_head_dct[tag] = []
                            head_file = filename.replace(tag, "desc-head")
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
        # in preproc)
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
            "space-longitudinal_desc-preproc_T1w",
            template_node,
            "brain_template",
            {},
            "",
            template_node_name,
        )

        rpool.set_data(
            "space-longitudinal_desc-preproc_T1w-template",
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
            "space-longitudinal_desc-preproc_T1w",
            "space-longitudinal_desc-reorient_T1w",
            "space-longitudinal_desc-brain_mask",
        ]
        rpool.gather_pipes(wf, config, add_excl=excl)

        # this is going to run multiple times!
        # once for every strategy!
        wf.run()

        # now, just write out a copy of the above to each session
        config.pipeline_setup["pipeline_name"] = orig_pipe_name
        for i in range(len(sub_list)):
            unique_id = sub_list[i]["unique_id"]
            session_id_list.append(unique_id)

            try:
                creds_path = sub_list[i]["creds_path"]
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
                sub_list[i],
            )

            wf, rpool = initiate_rpool(wf, config, sub_list[i])
            long_id = f"sub-longitudinal_ses-{subject_id}"
            config.pipeline_setup["pipeline_name"] = f"longitudinal_{orig_pipe_name}"
            long_id = long_id.replace("ses-sub-", "ses-sub")
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
                "space-longitudinal_desc-preproc_T1w",
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

            # rpool.gather_pipes(wf, config, add_excl=excl)

            # wf.run()

            # # begin single-session stuff again
            # for session in sub_list:
            #     unique_id = session["unique_id"]

            #     try:
            #         creds_path = session["creds_path"]
            #         if creds_path and "none" not in creds_path.lower():
            #             if os.path.exists(creds_path):
            #                 input_creds_path = os.path.abspath(creds_path)
            #             else:
            #                 err_msg = (
            #                     'Credentials path: "%s" for subject "%s" '
            #                     'session "%s" was not found. Check this path '
            #                     "and try again." % (creds_path, subject_id, unique_id)
            #                 )
            #                 raise Exception(err_msg)
            #         else:
            #             input_creds_path = None
            #     except KeyError:
            #         input_creds_path = None

            # nothing in rpool hmm

            # wf = initialize_nipype_wf(config, sub_list[0])

            # wf, rpool = initiate_rpool(wf, config, session)
            # print('first', rpool.get_resources())

            # print('second', rpool.get_resources())
            pipeline_blocks = [
                warp_longitudinal_T1w_to_template,
                warp_longitudinal_seg_to_T1w,
            ]

            wf = connect_pipeline(wf, config, rpool, pipeline_blocks)

            rpool.gather_pipes(wf, config)

            # this is going to run multiple times!
            # once for every strategy!
            wf.run()
