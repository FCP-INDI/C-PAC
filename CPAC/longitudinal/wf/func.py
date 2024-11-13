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

import os

import nipype.interfaces.io as nio

from CPAC.longitudinal.preproc import subject_specific_template
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.registration import (
    create_fsl_flirt_linear_reg,
    create_fsl_fnirt_nonlinear_reg,
    create_wf_calculate_ants_warp,
)
from CPAC.utils.datasource import (
    resolve_resolution,
)
from CPAC.utils.interfaces.function import Function
from CPAC.utils.strategy import Strategy
from CPAC.utils.utils import check_config_resources


# TODO check:
# 1 func alone works
# 2 anat + func works, pass anat strategy list?
def func_preproc_longitudinal_wf(subject_id, sub_list, config):
    """
    Parameters
    ----------
    subject_id : string
        the id of the subject
    sub_list : list of dict
        this is a list of sessions for one subject and each session if the same dictionary as the one given to
        prep_workflow
    config : configuration
        a configuration object containing the information of the pipeline config. (Same as for prep_workflow)

    Returns
    -------
    strat_list_ses_list : list of list
        a list of strategies; within each strategy, a list of sessions
    """
    datasink = pe.Node(nio.DataSink(), name="sinker")
    datasink.inputs.base_directory = config.pipeline_setup["working_directory"]["path"]

    session_id_list = []
    ses_list_strat_list = {}

    workflow_name = "func_preproc_longitudinal_" + str(subject_id)
    workflow = pe.Workflow(name=workflow_name)
    workflow.base_dir = config.pipeline_setup["working_directory"]["path"]
    workflow.config["execution"] = {
        "hash_method": "timestamp",
        "crashdump_dir": os.path.abspath(
            config.pipeline_setup["crash_directory"]["path"]
        ),
    }

    for sub_dict in sub_list:
        if "func" in sub_dict or "rest" in sub_dict:
            if "func" in sub_dict:
                sub_dict["func"]
            else:
                sub_dict["rest"]

            unique_id = sub_dict["unique_id"]
            session_id_list.append(unique_id)

            try:
                creds_path = sub_dict["creds_path"]
                if creds_path and "none" not in creds_path.lower():
                    if os.path.exists(creds_path):
                        input_creds_path = os.path.abspath(creds_path)
                    else:
                        err_msg = (
                            'Credentials path: "%s" for subject "%s" was not '
                            "found. Check this path and try again."
                            % (creds_path, subject_id)
                        )
                        raise Exception(err_msg)
                else:
                    input_creds_path = None
            except KeyError:
                input_creds_path = None

            strat = Strategy()
            strat_list = [strat]
            node_suffix = "_".join([subject_id, unique_id])

            # Functional Ingress Workflow
            # add optional flag
            workflow, diff, blip, fmap_rp_list = connect_func_ingress(
                workflow,
                strat_list,
                config,
                sub_dict,
                subject_id,
                input_creds_path,
                node_suffix,
            )

            # Functional Initial Prep Workflow
            workflow, strat_list = connect_func_init(
                workflow, strat_list, config, node_suffix
            )

            # Functional Image Preprocessing Workflow
            workflow, strat_list = connect_func_preproc(
                workflow, strat_list, config, node_suffix
            )

            # Distortion Correction
            workflow, strat_list = connect_distortion_correction(
                workflow, strat_list, config, diff, blip, fmap_rp_list, node_suffix
            )

            ses_list_strat_list[node_suffix] = strat_list

    # Here we have all the func_preproc set up for every session of the subject

    # TODO create a list of list ses_list_strat_list
    # a list of skullstripping strategies,
    # a list of sessions within each strategy list
    # TODO rename and reorganize dict
    # TODO update strat name
    strat_list_ses_list = {}
    strat_list_ses_list["func_default"] = []

    for sub_ses_id, strat_nodes_list in ses_list_strat_list.items():
        strat_list_ses_list["func_default"].append(strat_nodes_list[0])

    workflow.run()

    return strat_list_ses_list


def merge_func_preproc(working_directory):
    """
    Parameters
    ----------
    working_directory : string
        a path to the working directory

    Returns
    -------
    brain_list : list
        a list of func preprocessed brain
    skull_list : list
        a list of func preprocessed skull
    """
    brain_list = []
    skull_list = []

    for dirpath, dirnames, filenames in os.walk(working_directory):
        for f in filenames:
            if "func_get_preprocessed_median" in dirpath and ".nii.gz" in f:
                filepath = os.path.join(dirpath, f)
                brain_list.append(filepath)
            if "func_get_motion_correct_median" in dirpath and ".nii.gz" in f:
                filepath = os.path.join(dirpath, f)
                skull_list.append(filepath)

    brain_list.sort()
    skull_list.sort()

    return brain_list, skull_list


def register_func_longitudinal_template_to_standard(
    longitudinal_template_node, c, workflow, strat_init, strat_name
):
    (
        sub_mem_gb,
        num_cores_per_sub,
        num_ants_cores,
        num_omp_cores,
    ) = check_config_resources(c)

    strat_init_new = strat_init.fork()

    strat_init_new.update_resource_pool(
        {
            "functional_preprocessed_median": (
                longitudinal_template_node,
                "brain_template",
            ),
            "motion_correct_median": (longitudinal_template_node, "skull_template"),
        }
    )

    strat_list = [strat_init_new]

    new_strat_list = []

    regOption = c.anatomical_preproc["registration_workflow"]["registration"]["using"]

    if "FSL" in regOption:
        for num_strat, strat in enumerate(strat_list):
            flirt_reg_func_mni = create_fsl_flirt_linear_reg(
                "func_mni_flirt_register_%s_%d" % (strat_name, num_strat)
            )

            if c.functional_registration["2-func_registration_to_template"][
                "FNIRT_pipelines"
            ]["interpolation"] not in ["trilinear", "sinc", "spline"]:
                err_msg = 'The selected FSL interpolation method may be in the list of values: "trilinear", "sinc", "spline"'
                raise Exception(err_msg)

            # Input registration parameters
            flirt_reg_func_mni.inputs.inputspec.interp = c.functional_registration[
                "2-func_registration_to_template"
            ]["FNIRT_pipelines"]["interpolation"]

            node, out_file = strat["functional_preprocessed_median"]
            workflow.connect(
                node, out_file, flirt_reg_func_mni, "inputspec.input_brain"
            )

            # pass the reference files
            node, out_file = strat["template_brain_for_func_preproc"]
            workflow.connect(
                node, out_file, flirt_reg_func_mni, "inputspec.reference_brain"
            )

            if "ANTS" in regOption:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(flirt_reg_func_mni.name)

            strat.update_resource_pool(
                {
                    "registration_method": "FSL",
                    "func_longitudinal_to_mni_linear_xfm": (
                        flirt_reg_func_mni,
                        "outputspec.linear_xfm",
                    ),
                    "mni_to_func_longitudinal_linear_xfm": (
                        flirt_reg_func_mni,
                        "outputspec.invlinear_xfm",
                    ),
                    "func_longitudinal_template_to_standard": (
                        flirt_reg_func_mni,
                        "outputspec.output_brain",
                    ),
                }
            )

    strat_list += new_strat_list

    new_strat_list = []

    try:
        fsl_linear_reg_only = c.fsl_linear_reg_only
    except AttributeError:
        fsl_linear_reg_only = [0]

    if "FSL" in regOption and 0 in fsl_linear_reg_only:
        for num_strat, strat in enumerate(strat_list):
            if strat.get("registration_method") == "FSL":
                fnirt_reg_func_mni = create_fsl_fnirt_nonlinear_reg(
                    "func_mni_fnirt_register_%s_%d" % (strat_name, num_strat)
                )

                # brain input
                node, out_file = strat["functional_preprocessed_median"]
                workflow.connect(
                    node, out_file, fnirt_reg_func_mni, "inputspec.input_brain"
                )

                # brain reference
                node, out_file = strat["template_brain_for_func_preproc"]
                workflow.connect(
                    node, out_file, fnirt_reg_func_mni, "inputspec.reference_brain"
                )

                # skull input
                node, out_file = strat["motion_correct_median"]
                workflow.connect(
                    node, out_file, fnirt_reg_func_mni, "inputspec.input_skull"
                )

                # skull reference
                node, out_file = strat["template_skull_for_func_preproc"]
                workflow.connect(
                    node, out_file, fnirt_reg_func_mni, "inputspec.reference_skull"
                )

                node, out_file = strat["func_longitudinal_to_mni_linear_xfm"]
                workflow.connect(
                    node, out_file, fnirt_reg_func_mni, "inputspec.linear_aff"
                )

                node, out_file = strat["template_ref_mask"]
                workflow.connect(
                    node, out_file, fnirt_reg_func_mni, "inputspec.ref_mask"
                )

                # assign the FSL FNIRT config file specified in pipeline
                # config.yml
                fnirt_reg_func_mni.inputs.inputspec.fnirt_config = c.anatomical_preproc[
                    "registration_workflow"
                ]["registration"]["FSL-FNIRT"]["fnirt_config"]

                if 1 in fsl_linear_reg_only:
                    strat = strat.fork()
                    new_strat_list.append(strat)

                strat.append_name(fnirt_reg_func_mni.name)

                strat.update_resource_pool(
                    {
                        "func_longitudinal_to_mni_nonlinear_xfm": (
                            fnirt_reg_func_mni,
                            "outputspec.nonlinear_xfm",
                        ),
                        "func_longitudinal_template_to_standard": (
                            fnirt_reg_func_mni,
                            "outputspec.output_brain",
                        ),
                    },
                    override=True,
                )

    strat_list += new_strat_list

    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):
        # or run ANTS anatomical-to-MNI registration instead
        if "ANTS" in regOption and strat.get("registration_method") != "FSL":
            ants_reg_func_mni = create_wf_calculate_ants_warp(
                "func_mni_ants_register_%s_%d" % (strat_name, num_strat),
                num_threads=num_ants_cores,
                reg_ants_skull=c.anatomical_preproc["registration_workflow"][
                    "reg_with_skull"
                ],
            )

            if c.functional_registration["2-func_registration_to_template"][
                "ANTs_pipelines"
            ]["interpolation"] not in ["Linear", "BSpline", "LanczosWindowedSinc"]:
                err_msg = 'The selected ANTS interpolation method may be in the list of values: "Linear", "BSpline", "LanczosWindowedSinc"'
                raise Exception(err_msg)

            # Input registration parameters
            ants_reg_func_mni.inputs.inputspec.interp = c.functional_registration[
                "2-func_registration_to_template"
            ]["ANTs_pipelines"]["interpolation"]

            # calculating the transform with the skullstripped is
            # reported to be better, but it requires very high
            # quality skullstripping. If skullstripping is imprecise
            # registration with skull is preferred
            if c.anatomical_preproc["registration_workflow"]["reg_with_skull"]:
                # get the skull-stripped anatomical from resource pool
                node, out_file = strat["functional_preprocessed_median"]

                # pass the anatomical to the workflow
                workflow.connect(
                    node, out_file, ants_reg_func_mni, "inputspec.moving_brain"
                )

                # get the reorient skull-on anatomical from resource pool
                node, out_file = strat["motion_correct_median"]

                # pass the anatomical to the workflow
                workflow.connect(
                    node, out_file, ants_reg_func_mni, "inputspec.moving_skull"
                )

                # pass the reference file
                node, out_file = strat["template_brain_for_func_preproc"]
                workflow.connect(
                    node, out_file, ants_reg_func_mni, "inputspec.reference_brain"
                )

                # pass the reference file
                node, out_file = strat["template_skull_for_func_preproc"]
                workflow.connect(
                    node, out_file, ants_reg_func_mni, "inputspec.reference_skull"
                )

            else:
                node, out_file = strat["functional_preprocessed_median"]

                workflow.connect(
                    node, out_file, ants_reg_func_mni, "inputspec.moving_brain"
                )

                # pass the reference file
                node, out_file = strat["template_brain_for_func_preproc"]
                workflow.connect(
                    node, out_file, ants_reg_func_mni, "inputspec.reference_brain"
                )

            # pass the reference mask file
            node, out_file = strat["template_brain_mask_for_func_preproc"]
            workflow.connect(
                node, out_file, ants_reg_func_mni, "inputspec.reference_mask"
            )

            # pass the reference mask file
            node, out_file = strat["functional_brain_mask"]
            workflow.connect(node, out_file, ants_reg_func_mni, "inputspec.moving_mask")

            ants_reg_func_mni.inputs.inputspec.ants_para = c.anatomical_preproc[
                "registration_workflow"
            ]["registration"]["ANTs"]["T1_registration"]
            ants_reg_func_mni.inputs.inputspec.fixed_image_mask = None

            strat.append_name(ants_reg_func_mni.name)

            strat.update_resource_pool(
                {
                    "registration_method": "ANTS",
                    "ants_initial_xfm": (
                        ants_reg_func_mni,
                        "outputspec.ants_initial_xfm",
                    ),
                    "ants_rigid_xfm": (ants_reg_func_mni, "outputspec.ants_rigid_xfm"),
                    "ants_affine_xfm": (
                        ants_reg_func_mni,
                        "outputspec.ants_affine_xfm",
                    ),
                    "func_longitudinal_to_mni_nonlinear_xfm": (
                        ants_reg_func_mni,
                        "outputspec.warp_field",
                    ),
                    "mni_to_func_longitudinal_nonlinear_xfm": (
                        ants_reg_func_mni,
                        "outputspec.inverse_warp_field",
                    ),
                    "func_longitudinal_to_mni_ants_composite_xfm": (
                        ants_reg_func_mni,
                        "outputspec.composite_transform",
                    ),
                    "func_longitudinal_template_to_standard": (
                        ants_reg_func_mni,
                        "outputspec.normalized_output_brain",
                    ),
                }
            )

    strat_list += new_strat_list

    """
    # Func -> T1 Registration (Initial Linear Reg)
    workflow, strat_list, diff_complete = connect_func_to_anat_init_reg(workflow, strat_list, c)

    # Func -> T1 Registration (BBREG)
    workflow, strat_list = connect_func_to_anat_bbreg(workflow, strat_list, c, diff_complete)

    # Func -> T1/EPI Template
    workflow, strat_list = connect_func_to_template_reg(workflow, strat_list, c)
    """

    return workflow, strat_list


def func_longitudinal_template_wf(subject_id, strat_list, config):
    """
    Parameters
    ----------
    subject_id : string
        the id of the subject
    strat_list : list of list
        first level strategy, second level session
    config : configuration
        a configuration object containing the information of the pipeline config.

    Returns
    -------
        None
    """
    workflow_name = "func_longitudinal_template_" + str(subject_id)
    workflow = pe.Workflow(name=workflow_name)
    workflow.base_dir = config.pipeline_setup["working_directory"]["path"]
    workflow.config["execution"] = {
        "hash_method": "timestamp",
        "crashdump_dir": os.path.abspath(
            config.pipeline_setup["crash_directory"]["path"]
        ),
    }

    # strat_nodes_list = strat_list['func_default']
    strat_init = Strategy()

    templates_for_resampling = [
        (
            config.resolution_for_func_preproc,
            config.template_brain_only_for_func,
            "template_brain_for_func_preproc",
            "resolution_for_func_preproc",
        ),
        (
            config.resolution_for_func_preproc,
            config.template_skull_for_func,
            "template_skull_for_func_preproc",
            "resolution_for_func_preproc",
        ),
        (
            config.resolution_for_func_preproc,
            config.ref_mask_for_func,
            "template_ref_mask",
            "resolution_for_func_preproc",
        ),
        # TODO check float resolution
        (
            config.resolution_for_func_preproc,
            config.functional_registration["2-func_registration_to_template"][
                "target_template"
            ]["EPI_template"]["template_epi"],
            "template_epi",
            "resolution_for_func_preproc",
        ),
        (
            config.resolution_for_func_derivative,
            config.functional_registration["2-func_registration_to_template"][
                "target_template"
            ]["EPI_template"]["template_epi"],
            "template_epi_derivative",
            "resolution_for_func_derivative",
        ),
        (
            config.resolution_for_func_derivative,
            config.template_brain_only_for_func,
            "template_brain_for_func_derivative",
            "resolution_for_func_preproc",
        ),
        (
            config.resolution_for_func_derivative,
            config.template_skull_for_func,
            "template_skull_for_func_derivative",
            "resolution_for_func_preproc",
        ),
    ]

    for resolution, template, template_name, tag in templates_for_resampling:
        resampled_template = pe.Node(
            Function(
                input_names=["resolution", "template", "template_name", "tag"],
                output_names=["resampled_template"],
                function=resolve_resolution,
                as_module=True,
            ),
            name="resampled_" + template_name,
        )

        resampled_template.inputs.resolution = resolution
        resampled_template.inputs.template = template
        resampled_template.inputs.template_name = template_name
        resampled_template.inputs.tag = tag
        resampled_template.inputs.orientation = config["desired_orientation"]

        strat_init.update_resource_pool(
            {template_name: (resampled_template, "resampled_template")}
        )

    merge_func_preproc_node = pe.Node(
        Function(
            input_names=["working_directory"],
            output_names=["brain_list", "skull_list"],
            function=merge_func_preproc,
            as_module=True,
        ),
        name="merge_func_preproc",
    )

    merge_func_preproc_node.inputs.working_directory = config.pipeline_setup[
        "working_directory"
    ]["path"]

    template_node = subject_specific_template(
        workflow_name="subject_specific_func_template_" + subject_id
    )

    template_node.inputs.set(
        avg_method=config.longitudinal_template_average_method,
        dof=config.longitudinal_template_dof,
        interp=config.longitudinal_template_interp,
        cost=config.longitudinal_template_cost,
        convergence_threshold=config.longitudinal_template_convergence_threshold,
        thread_pool=config.longitudinal_template_thread_pool,
    )

    workflow.connect(
        merge_func_preproc_node, "brain_list", template_node, "input_brain_list"
    )

    workflow.connect(
        merge_func_preproc_node, "skull_list", template_node, "input_skull_list"
    )

    workflow, strat_list = register_func_longitudinal_template_to_standard(
        template_node, config, workflow, strat_init, "default"
    )

    workflow.run()
