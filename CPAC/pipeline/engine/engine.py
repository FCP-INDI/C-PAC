# Copyright (C) 2021-2024  C-PAC Developers

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
import ast
import copy
import hashlib
from itertools import chain
import json
import logging
import os
import re
from typing import Optional
import warnings
import pandas as pd


from nipype import config
from nipype.interfaces.utility import Rename

from CPAC.image_utils.spatial_smoothing import spatial_smoothing
from CPAC.image_utils.statistical_transforms import (
    fisher_z_score_standardize,
    z_score_standardize,
)
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.check_outputs import ExpectedOutputs
from CPAC.pipeline.nodeblock import NodeBlockFunction
from CPAC.pipeline.utils import MOVEMENT_FILTER_KEYS, name_fork, source_set
from CPAC.registration.registration import transform_derivative
from CPAC.resources.templates.lookup_table import lookup_identifier
from CPAC.utils.bids_utils import res_in_filename
from CPAC.utils.configuration import Configuration
from CPAC.utils.datasource import (
    create_anat_datasource,
    create_func_datasource,
    create_general_datasource,
    ingress_func_metadata,
    resolve_resolution,
)
from CPAC.utils.interfaces.datasink import DataSink
from CPAC.utils.interfaces.function import Function
from CPAC.utils.monitoring import (
    getLogger,
    LOGTAIL,
    WARNING_FREESURFER_OFF_WITH_DATA,
    WFLOGGER,
)
from CPAC.utils.outputs import Outputs
from CPAC.utils.utils import (
    check_prov_for_regtool,
    create_id_string,
    get_last_prov_entry,
    read_json,
    write_output_json,
)



class NodeBlock:
    def __init__(self, node_block_functions, debug=False):
        if not isinstance(node_block_functions, list):
            node_block_functions = [node_block_functions]

        self.node_blocks = {}

        for node_block_function in node_block_functions:  # <---- sets up the NodeBlock object in case you gave it a list of node blocks instead of a single one - for option forking.
            self.input_interface = []
            if isinstance(node_block_function, tuple):
                self.input_interface = node_block_function[1]
                node_block_function = node_block_function[0]
                if not isinstance(self.input_interface, list):
                    self.input_interface = [self.input_interface]

            if not isinstance(node_block_function, NodeBlockFunction):
                # If the object is a plain function `__name__` will be more useful than `str()`
                obj_str = (
                    node_block_function.__name__
                    if hasattr(node_block_function, "__name__")
                    else str(node_block_function)
                )
                msg = f'Object is not a nodeblock: "{obj_str}"'
                raise TypeError(msg)

            name = node_block_function.name
            self.name = name
            self.node_blocks[name] = {}

            if self.input_interface:
                for interface in self.input_interface:
                    for orig_input in node_block_function.inputs:
                        if isinstance(orig_input, tuple):
                            list_tup = list(orig_input)
                            if interface[0] in list_tup:
                                list_tup.remove(interface[0])
                                list_tup.append(interface[1])
                                node_block_function.inputs.remove(orig_input)
                                node_block_function.inputs.append(tuple(list_tup))
                        elif orig_input == interface[0]:
                            node_block_function.inputs.remove(interface[0])
                            node_block_function.inputs.append(interface[1])

            for key, val in node_block_function.legacy_nodeblock_dict().items():
                self.node_blocks[name][key] = val

            self.node_blocks[name]["block_function"] = node_block_function

            # TODO: fix/replace below
            self.outputs = {}
            for out in node_block_function.outputs:
                self.outputs[out] = None

            self.options = ["base"]
            if node_block_function.outputs is not None:
                self.options = node_block_function.outputs

            WFLOGGER.info("Connecting %s...", name)
            if debug:
                config.update_config({"logging": {"workflow_level": "DEBUG"}})
                logging.update_logging(config)
                WFLOGGER.debug(
                    '"inputs": %s\n\t "outputs": %s%s',
                    node_block_function.inputs,
                    list(self.outputs.keys()),
                    f'\n\t"options": {self.options}'
                    if self.options != ["base"]
                    else "",
                )
                config.update_config({"logging": {"workflow_level": "INFO"}})
                logging.update_logging(config)

    def get_name(self):
        return self.name

    def check_null(self, val):
        if isinstance(val, str):
            val = None if val.lower() == "none" else val
        return val

    def check_output(self, outputs, label, name):
        if label not in outputs:
            msg = (
                f'\n[!] Output name "{label}" in the block '
                "function does not match the outputs list "
                f'{outputs} in Node Block "{name}"\n'
            )
            raise NameError(msg)

    def grab_tiered_dct(self, cfg, key_list):
        cfg_dct = cfg.dict()
        for key in key_list:
            try:
                cfg_dct = cfg_dct.get(key, {})
            except KeyError as ke:
                msg = "[!] The config provided to the node block is not valid"
                raise KeyError(msg) from ke
        return cfg_dct

    def connect_block(self, wf, cfg, rpool):
        debug = cfg.pipeline_setup["Debugging"]["verbose"]
        all_opts = []
        for name, block_dct in self.node_blocks.items():
            opts = []
            config = self.check_null(block_dct["config"])
            option_key = self.check_null(block_dct["option_key"])
            option_val = self.check_null(block_dct["option_val"])
            if option_key and option_val:
                if not isinstance(option_key, list):
                    option_key = [option_key]
                if not isinstance(option_val, list):
                    option_val = [option_val]
                if config:
                    key_list = config + option_key
                else:
                    key_list = option_key
                if "USER-DEFINED" in option_val:
                    # load custom config data into each 'opt'
                    opts = self.grab_tiered_dct(cfg, key_list)
                else:
                    for option in option_val:
                        try:
                            if option in self.grab_tiered_dct(cfg, key_list):
                                # goes over the option_vals in the node block docstring, and checks if the user's pipeline config included it in the forking list
                                opts.append(option)
                        except AttributeError as err:
                            msg = f"{err}\nNode Block: {name}"
                            raise Exception(msg)

                if opts is None:
                    opts = [opts]

            elif option_key and not option_val:
                # enables multiple config forking entries
                if not isinstance(option_key[0], list):
                    msg = (
                        f"[!] The option_key field ({option_key}) "
                        f"for {name} exists but there is no "
                        "option_val.\n\nIf you are trying to "
                        "populate multiple option keys, the "
                        "option_val field must contain a list of "
                        "a list.\n"
                    )
                    raise ValueError(msg)
                for option_config in option_key:
                    # option_config is a list of pipe config levels down to the option
                    if config:
                        key_list = config + option_config
                    else:
                        key_list = option_config
                    option_val = option_config[-1]
                    if option_val in self.grab_tiered_dct(cfg, key_list[:-1]):
                        opts.append(option_val)
            else:  # AND, if there are multiple option-val's (in a list) in the docstring, it gets iterated below in 'for opt in option' etc. AND THAT'S WHEN YOU HAVE TO DELINEATE WITHIN THE NODE BLOCK CODE!!!
                opts = [None]
            all_opts += opts

        sidecar_additions = {
            "CpacConfigHash": hashlib.sha1(
                json.dumps(cfg.dict(), sort_keys=True).encode("utf-8")
            ).hexdigest(),
            "CpacConfig": cfg.dict(),
        }

        if cfg["pipeline_setup"]["output_directory"].get("user_defined"):
            sidecar_additions["UserDefined"] = cfg["pipeline_setup"][
                "output_directory"
            ]["user_defined"]

        for name, block_dct in self.node_blocks.items():
            # iterates over either the single node block in the sequence, or a list of node blocks within the list of node blocks, i.e. for option forking.
            switch = self.check_null(block_dct["switch"])
            config = self.check_null(block_dct["config"])
            option_key = self.check_null(block_dct["option_key"])
            option_val = self.check_null(block_dct["option_val"])
            inputs = self.check_null(block_dct["inputs"])
            outputs = self.check_null(block_dct["outputs"])

            block_function = block_dct["block_function"]

            opts = []
            if option_key and option_val:
                if not isinstance(option_key, list):
                    option_key = [option_key]
                if not isinstance(option_val, list):
                    option_val = [option_val]
                if config:
                    key_list = config + option_key
                else:
                    key_list = option_key
                if "USER-DEFINED" in option_val:
                    # load custom config data into each 'opt'
                    opts = self.grab_tiered_dct(cfg, key_list)
                else:
                    for option in option_val:
                        if option in self.grab_tiered_dct(cfg, key_list):
                            # goes over the option_vals in the node block docstring, and checks if the user's pipeline config included it in the forking list
                            opts.append(option)
            else:  # AND, if there are multiple option-val's (in a list) in the docstring, it gets iterated below in 'for opt in option' etc. AND THAT'S WHEN YOU HAVE TO DELINEATE WITHIN THE NODE BLOCK CODE!!!
                opts = [None]
                # THIS ALSO MEANS the multiple option-val's in docstring node blocks can be entered once in the entire node-block sequence, not in a list of multiples
            if not opts:
                # for node blocks where the options are split into different
                # block functions - opts will be empty for non-selected
                # options, and would waste the get_strats effort below
                continue

            if not switch:
                switch = [True]
            else:
                if config:
                    try:
                        key_list = config + switch
                    except TypeError as te:
                        msg = (
                            "\n\n[!] Developer info: Docstring error "
                            f"for {name}, make sure the 'config' or "
                            "'switch' fields are lists.\n\n"
                        )
                        raise TypeError(msg) from te
                    switch = self.grab_tiered_dct(cfg, key_list)
                elif isinstance(switch[0], list):
                    # we have multiple switches, which is designed to only work if
                    # config is set to "None"
                    switch_list = []
                    for key_list in switch:
                        val = self.grab_tiered_dct(cfg, key_list)
                        if isinstance(val, list):
                            # fork switches
                            if True in val:
                                switch_list.append(True)
                            if False in val:
                                switch_list.append(False)
                        else:
                            switch_list.append(val)
                    if False in switch_list:
                        switch = [False]
                    else:
                        switch = [True]
                else:
                    # if config is set to "None"
                    key_list = switch
                    switch = self.grab_tiered_dct(cfg, key_list)
                if not isinstance(switch, list):
                    switch = [switch]
            if True in switch:
                for (
                    pipe_idx,
                    strat_pool,  # strat_pool is a ResourcePool like {'desc-preproc_T1w': { 'json': info, 'data': (node, out) }, 'desc-brain_mask': etc.}
                ) in rpool.get_strats(inputs, debug).items():
                    # keep in mind rpool.get_strats(inputs) = {pipe_idx1: {'desc-preproc_T1w': etc.}, pipe_idx2: {..} }
                    fork = False in switch
                    for opt in opts:  # it's a dictionary of ResourcePools called strat_pools, except those sub-ResourcePools only have one level! no pipe_idx strat keys.
                        # remember, you can get 'data' or 'json' from strat_pool with member functions
                        # strat_pool has all of the JSON information of all the inputs!
                        # so when we set_data below for the TOP-LEVEL MAIN RPOOL (not the strat_pool), we can generate new merged JSON information for each output.
                        # particularly, our custom 'CpacProvenance' field.
                        node_name = name
                        pipe_x = rpool.get_pipe_number(pipe_idx)

                        replaced_inputs = []
                        for interface in self.input_interface:
                            if isinstance(interface[1], list):
                                for input_name in interface[1]:
                                    if strat_pool.check_rpool(input_name):
                                        break
                            else:
                                input_name = interface[1]
                            strat_pool.copy_resource(input_name, interface[0])
                            replaced_inputs.append(interface[0])
                        try:
                            wf, outs = block_function(wf, cfg, strat_pool, pipe_x, opt)
                        except IOError as e:  # duplicate node
                            WFLOGGER.warning(e)
                            continue

                        if not outs:
                            if block_function.__name__ == "freesurfer_postproc":
                                WFLOGGER.warning(WARNING_FREESURFER_OFF_WITH_DATA)
                                LOGTAIL["warnings"].append(
                                    WARNING_FREESURFER_OFF_WITH_DATA
                                )
                            continue

                        if opt and len(option_val) > 1:
                            node_name = f"{node_name}_{opt}"
                        elif opt and "USER-DEFINED" in option_val:
                            node_name = f'{node_name}_{opt["Name"]}'

                        if debug:
                            verbose_logger = getLogger("CPAC.engine")
                            verbose_logger.debug("\n=======================")
                            verbose_logger.debug("Node name: %s", node_name)
                            prov_dct = rpool.get_resource_strats_from_prov(
                                ast.literal_eval(pipe_idx)
                            )
                            for key, val in prov_dct.items():
                                verbose_logger.debug("-------------------")
                                verbose_logger.debug("Input - %s:", key)
                                sub_prov_dct = rpool.get_resource_strats_from_prov(val)
                                for sub_key, sub_val in sub_prov_dct.items():
                                    sub_sub_dct = rpool.get_resource_strats_from_prov(
                                        sub_val
                                    )
                                    verbose_logger.debug("  sub-input - %s:", sub_key)
                                    verbose_logger.debug("    prov = %s", sub_val)
                                    verbose_logger.debug(
                                        "    sub_sub_inputs = %s", sub_sub_dct.keys()
                                    )

                        for label, connection in outs.items():
                            self.check_output(outputs, label, name)
                            new_json_info = copy.deepcopy(strat_pool.get("json"))

                            # transfer over data-specific json info
                            # for example, if the input data json is _bold and the output is also _bold
                            data_type = label.split("_")[-1]
                            if data_type in new_json_info["subjson"]:
                                if (
                                    "SkullStripped"
                                    in new_json_info["subjson"][data_type]
                                ):
                                    new_json_info["SkullStripped"] = new_json_info[
                                        "subjson"
                                    ][data_type]["SkullStripped"]

                            # determine sources for the outputs, i.e. all input data into the node block
                            new_json_info["Sources"] = [
                                x
                                for x in strat_pool.get_entire_rpool()
                                if x != "json" and x not in replaced_inputs
                            ]

                            if isinstance(outputs, dict):
                                new_json_info.update(outputs[label])
                                if "Description" not in outputs[label]:
                                    # don't propagate old Description
                                    try:
                                        del new_json_info["Description"]
                                    except KeyError:
                                        pass
                                if "Template" in outputs[label]:
                                    template_key = outputs[label]["Template"]
                                    if template_key in new_json_info["Sources"]:
                                        # only if the pipeline config template key is entered as the 'Template' field
                                        # otherwise, skip this and take in the literal 'Template' string
                                        try:
                                            new_json_info["Template"] = new_json_info[
                                                "subjson"
                                            ][template_key]["Description"]
                                        except KeyError:
                                            pass
                                    try:
                                        new_json_info["Resolution"] = new_json_info[
                                            "subjson"
                                        ][template_key]["Resolution"]
                                    except KeyError:
                                        pass
                            else:
                                # don't propagate old Description
                                try:
                                    del new_json_info["Description"]
                                except KeyError:
                                    pass

                            if "Description" in new_json_info:
                                new_json_info["Description"] = " ".join(
                                    new_json_info["Description"].split()
                                )

                            for sidecar_key, sidecar_value in sidecar_additions.items():
                                if sidecar_key not in new_json_info:
                                    new_json_info[sidecar_key] = sidecar_value

                            try:
                                del new_json_info["subjson"]
                            except KeyError:
                                pass

                            if fork or len(opts) > 1 or len(all_opts) > 1:
                                if "CpacVariant" not in new_json_info:
                                    new_json_info["CpacVariant"] = {}
                                raw_label = rpool.get_raw_label(label)
                                if raw_label not in new_json_info["CpacVariant"]:
                                    new_json_info["CpacVariant"][raw_label] = []
                                new_json_info["CpacVariant"][raw_label].append(
                                    node_name
                                )

                            rpool.set_data(
                                label,
                                connection[0],
                                connection[1],
                                new_json_info,
                                pipe_idx,
                                node_name,
                                fork,
                            )

                            wf, post_labels = rpool.post_process(
                                wf,
                                label,
                                connection,
                                new_json_info,
                                pipe_idx,
                                pipe_x,
                                outs,
                            )

                            if rpool.func_reg:
                                for postlabel in post_labels:
                                    connection = (postlabel[1], postlabel[2])
                                    wf = rpool.derivative_xfm(
                                        wf,
                                        postlabel[0],
                                        connection,
                                        new_json_info,
                                        pipe_idx,
                                        pipe_x,
                                    )
        return wf


def wrap_block(node_blocks, interface, wf, cfg, strat_pool, pipe_num, opt):
    """Wrap a list of node block functions to use within other node blocks.

    Example usage:

        # This calls the 'bold_mask_afni' and 'bold_masking' node blocks to
        # skull-strip an EPI field map, without having to invoke the NodeBlock
        # connection system.

        # The interface dictionary tells wrap_block to set the EPI field map
        # in the parent node block's throw-away strat_pool as 'bold', so that
        # the 'bold_mask_afni' and 'bold_masking' node blocks will see that as
        # the 'bold' input.

        # It also tells wrap_block to set the 'desc-brain_bold' output of
        # the 'bold_masking' node block to 'opposite_pe_epi_brain' (what it
        # actually is) in the parent node block's strat_pool, which gets
        # returned.

        # Note 'bold' and 'desc-brain_bold' (all on the left side) are the
        # labels that 'bold_mask_afni' and 'bold_masking' understand/expect
        # through their interfaces and docstrings.

        # The right-hand side (the values of the 'interface' dictionary) are
        # what 'make sense' within the current parent node block - in this
        # case, the distortion correction node block dealing with field maps.

        interface = {'bold': (match_epi_fmaps_node, 'opposite_pe_epi'),
                     'desc-brain_bold': 'opposite_pe_epi_brain'}
        wf, strat_pool = wrap_block([bold_mask_afni, bold_masking],
                                    interface, wf, cfg, strat_pool,
                                    pipe_num, opt)

        ...further downstream in the parent node block:

        node, out = strat_pool.get_data('opposite_pe_epi_brain')

        # The above line will connect the output of the 'bold_masking' node
        # block (which is the skull-stripped version of 'opposite_pe_epi') to
        # the next node.

    """
    for block in node_blocks:
        # new_pool = copy.deepcopy(strat_pool)
        for in_resource, val in interface.items():
            if isinstance(val, tuple):
                strat_pool.set_data(
                    in_resource, val[0], val[1], {}, "", "", fork=True
                )  #
        if "sub_num" not in strat_pool.get_pool_info():
            strat_pool.set_pool_info({"sub_num": 0})
        sub_num = strat_pool.get_pool_info()["sub_num"]

        wf, outputs = block(wf, cfg, strat_pool, f"{pipe_num}-{sub_num}", opt)  #
        for out, val in outputs.items():
            if out in interface and isinstance(interface[out], str):
                strat_pool.set_data(
                    interface[out], outputs[out][0], outputs[out][1], {}, "", ""
                )
            else:
                strat_pool.set_data(out, outputs[out][0], outputs[out][1], {}, "", "")
        sub_num += 1
        strat_pool.set_pool_info({"sub_num": sub_num})

    return (wf, strat_pool)


def ingress_raw_anat_data(wf, rpool, cfg, data_paths, unique_id, part_id, ses_id):
    if "anat" not in data_paths[1]["ent__datatype"].values:
        WFLOGGER.warning("No anatomical data present.")
        return rpool

    # if "creds_path" not in data_paths:
    #     data_paths["creds_path"] = None

    anat_flow = create_anat_datasource(f"anat_T1w_gather_{part_id}_{ses_id}")

    anat = {}
    anat_data = data_paths[1].loc[data_paths[1]["ent__datatype"] == "anat"]
    if "T1w" in anat_data["ent__suffix"].values:
        anat["T1"] = anat_data["finfo__file_path"].values[0]

    # if isinstance(data_paths["anat"], str):
    #     anat["T1"] = data_paths["anat"]
    # elif "T1w" in data_paths["anat"]:
    #     anat["T1"] = data_paths["anat"]["T1w"]

    if "T1" in anat:
        anat_flow.inputs.inputnode.set(
            subject=part_id,
            anat=anat["T1"],
            creds_path=None,
            dl_dir=cfg.pipeline_setup["working_directory"]["path"],
            img_type="anat",
        )
        rpool.set_data("T1w", anat_flow, "outputspec.anat", {}, "", "anat_ingress")

    # if "T2w" in data_paths["anat"]:
    #     anat_flow_T2 = create_anat_datasource(f"anat_T2w_gather_{part_id}_{ses_id}")
    #     anat_flow_T2.inputs.inputnode.set(
    #         subject=part_id,
    #         anat=data_paths["anat"]["T2w"],
    #         creds_path=data_paths["creds_path"],
    #         dl_dir=cfg.pipeline_setup["working_directory"]["path"],
    #         img_type="anat",
    #     )
    #     rpool.set_data("T2w", anat_flow_T2, "outputspec.anat", {}, "", "anat_ingress")

    if cfg.surface_analysis["freesurfer"]["ingress_reconall"]:
        rpool = ingress_freesurfer(
            wf, rpool, cfg, data_paths, unique_id, part_id, ses_id
        )

    return rpool


def ingress_freesurfer(wf, rpool, cfg, data_paths, unique_id, part_id, ses_id):
    try:
        fs_path = os.path.join(cfg.pipeline_setup["freesurfer_dir"], part_id)
    except KeyError:
        WFLOGGER.warning("No FreeSurfer data present.")
        return rpool

    # fs_path = os.path.join(cfg.pipeline_setup['freesurfer_dir'], part_id)
    if not os.path.exists(fs_path):
        if "sub" in part_id:
            fs_path = os.path.join(
                cfg.pipeline_setup["freesurfer_dir"], part_id.replace("sub-", "")
            )
        else:
            fs_path = os.path.join(
                cfg.pipeline_setup["freesurfer_dir"], ("sub-" + part_id)
            )

        # patch for flo-specific data
        if not os.path.exists(fs_path):
            subj_ses = part_id + "-" + ses_id
            fs_path = os.path.join(cfg.pipeline_setup["freesurfer_dir"], subj_ses)
            if not os.path.exists(fs_path):
                WFLOGGER.info("No FreeSurfer data found for subject %s", part_id)
                return rpool

    # Check for double nested subj names
    if os.path.exists(os.path.join(fs_path, os.path.basename(fs_path))):
        fs_path = os.path.join(fs_path, part_id)

    fs_ingress = create_general_datasource("gather_freesurfer_dir")
    fs_ingress.inputs.inputnode.set(
        unique_id=unique_id,
        data=fs_path,
        creds_path=data_paths["creds_path"],
        dl_dir=cfg.pipeline_setup["working_directory"]["path"],
    )
    rpool.set_data(
        "freesurfer-subject-dir",
        fs_ingress,
        "outputspec.data",
        {},
        "",
        "freesurfer_config_ingress",
    )

    recon_outs = {
        "pipeline-fs_raw-average": "mri/rawavg.mgz",
        "pipeline-fs_subcortical-seg": "mri/aseg.mgz",
        "pipeline-fs_brainmask": "mri/brainmask.mgz",
        "pipeline-fs_wmparc": "mri/wmparc.mgz",
        "pipeline-fs_T1": "mri/T1.mgz",
        "pipeline-fs_hemi-L_desc-surface_curv": "surf/lh.curv",
        "pipeline-fs_hemi-R_desc-surface_curv": "surf/rh.curv",
        "pipeline-fs_hemi-L_desc-surfaceMesh_pial": "surf/lh.pial",
        "pipeline-fs_hemi-R_desc-surfaceMesh_pial": "surf/rh.pial",
        "pipeline-fs_hemi-L_desc-surfaceMesh_smoothwm": "surf/lh.smoothwm",
        "pipeline-fs_hemi-R_desc-surfaceMesh_smoothwm": "surf/rh.smoothwm",
        "pipeline-fs_hemi-L_desc-surfaceMesh_sphere": "surf/lh.sphere",
        "pipeline-fs_hemi-R_desc-surfaceMesh_sphere": "surf/rh.sphere",
        "pipeline-fs_hemi-L_desc-surfaceMap_sulc": "surf/lh.sulc",
        "pipeline-fs_hemi-R_desc-surfaceMap_sulc": "surf/rh.sulc",
        "pipeline-fs_hemi-L_desc-surfaceMap_thickness": "surf/lh.thickness",
        "pipeline-fs_hemi-R_desc-surfaceMap_thickness": "surf/rh.thickness",
        "pipeline-fs_hemi-L_desc-surfaceMap_volume": "surf/lh.volume",
        "pipeline-fs_hemi-R_desc-surfaceMap_volume": "surf/rh.volume",
        "pipeline-fs_hemi-L_desc-surfaceMesh_white": "surf/lh.white",
        "pipeline-fs_hemi-R_desc-surfaceMesh_white": "surf/rh.white",
        "pipeline-fs_xfm": "mri/transforms/talairach.lta",
    }

    for key, outfile in recon_outs.items():
        fullpath = os.path.join(fs_path, outfile)
        if os.path.exists(fullpath):
            fs_ingress = create_general_datasource(f"gather_fs_{key}_dir")
            fs_ingress.inputs.inputnode.set(
                unique_id=unique_id,
                data=fullpath,
                creds_path=data_paths["creds_path"],
                dl_dir=cfg.pipeline_setup["working_directory"]["path"],
            )
            rpool.set_data(
                key, fs_ingress, "outputspec.data", {}, "", f"fs_{key}_ingress"
            )
        else:
            warnings.warn(
                str(LookupError(f"\n[!] Path does not exist for {fullpath}.\n"))
            )

    return rpool


def ingress_raw_func_data(wf, rpool, cfg, data_paths, unique_id, part_id, ses_id):
    func_paths_dct = data_paths[1].loc[data_paths[1]["ent__datatype"] == "func"]

    func_wf = create_func_datasource(
        func_paths_dct, rpool, f"func_ingress_{part_id}_{ses_id}"
    )
    func_wf.inputs.inputnode.set(
        subject=part_id,
        creds_path=None,
        dl_dir=cfg.pipeline_setup["working_directory"]["path"],
    )
    func_wf.get_node("inputnode").iterables = ("scan", list(func_paths_dct.keys()))

    rpool.set_data("subject", func_wf, "outputspec.subject", {}, "", "func_ingress")
    rpool.set_data("bold", func_wf, "outputspec.rest", {}, "", "func_ingress")
    rpool.set_data("scan", func_wf, "outputspec.scan", {}, "", "func_ingress")
    rpool.set_data(
        "scan-params", func_wf, "outputspec.scan_params", {}, "", "scan_params_ingress"
    )

    # TODO: CHECK FOR PARAMETERS
    diff = None
    blip = None
    fmap_rp_list = None
    # wf, rpool, diff, blip, fmap_rp_list = ingress_func_metadata(
    #     wf, cfg, rpool, data_paths, part_id, None, ses_id
    # )

    # Memoize list of local functional scans
    # TODO: handle S3 files
    # Skip S3 files for now

    local_func_scans = (
        [file_path for file_path in func_paths_dct["finfo__file_path"].values]
        if not func_paths_dct.empty
        else []
    )

    # local_func_scans = [
    #     func_paths_dct[scan]["scan"]
    #     for scan in func_paths_dct.keys()
    #  #   if not func_paths_dct[scan]["scan"].startswith("s3://")
    # ]
    if local_func_scans:
        # pylint: disable=protected-access
        wf._local_func_scans = local_func_scans
        if cfg.pipeline_setup["Debugging"]["verbose"]:
            verbose_logger = getLogger("CPAC.engine")
            verbose_logger.debug("local_func_scans: %s", local_func_scans)
    del local_func_scans

    return (wf, rpool, diff, blip, fmap_rp_list)


def ingress_output_dir(
    wf, cfg, rpool, unique_id, data_paths, part_id, ses_id, creds_path=None
):
    dir_path = data_paths["derivatives_dir"]

    WFLOGGER.info("\nPulling outputs from %s.\n", dir_path)

    anat = os.path.join(dir_path, "anat")
    func = os.path.join(dir_path, "func")

    exts = [".nii", ".gz", ".mat", ".1D", ".txt", ".csv", ".rms", ".tsv"]

    outdir_anat = []
    outdir_func = []
    func_paths = {}
    func_dict = {}

    for subdir in [anat, func]:
        if os.path.isdir(subdir):
            for filename in os.listdir(subdir):
                for ext in exts:
                    if ext in filename:
                        if subdir == anat:
                            outdir_anat.append(os.path.join(subdir, filename))
                        else:
                            outdir_func.append(os.path.join(subdir, filename))

    # Add derivatives directory to rpool
    ingress = create_general_datasource("gather_derivatives_dir")
    ingress.inputs.inputnode.set(
        unique_id=unique_id,
        data=dir_path,
        creds_path=creds_path,
        dl_dir=cfg.pipeline_setup["working_directory"]["path"],
    )
    rpool.set_data(
        "derivatives-dir", ingress, "outputspec.data", {}, "", "outdir_config_ingress"
    )

    for subdir in [outdir_anat, outdir_func]:
        for filepath in subdir:
            filename = str(filepath)
            for ext in exts:
                filename = filename.split("/")[-1].replace(ext, "")

            data_label = filename.split(unique_id)[1].lstrip("_")

            if len(filename) == len(data_label):
                msg = (
                    "\n\n[!] Possibly wrong participant or "
                    "session in this directory?\n\n"
                    f"Filepath: {filepath}\n\n"
                )
                raise Exception(msg)

            bidstag = ""
            for tag in data_label.split("_"):
                for prefix in ["task-", "run-", "acq-", "rec"]:
                    if tag.startswith(prefix):
                        bidstag += f"{tag}_"
                        data_label = data_label.replace(f"{tag}_", "")
            data_label, json = strip_template(data_label, dir_path, filename)

            rpool, json_info, pipe_idx, node_name, data_label = json_outdir_ingress(
                rpool, filepath, exts, data_label, json
            )

            if (
                "template" in data_label
                and not json_info["Template"]
                == cfg.pipeline_setup["outdir_ingress"]["Template"]
            ):
                continue
            # Rename confounds to avoid confusion in nuisance regression
            if data_label.endswith("desc-confounds_timeseries"):
                data_label = "pipeline-ingress_desc-confounds_timeseries"

            if len(bidstag) > 1:
                # Remove tail symbol
                bidstag = bidstag[:-1]
                if bidstag.startswith("task-"):
                    bidstag = bidstag.replace("task-", "")

            # Rename bold mask for CPAC naming convention
            # and to avoid collision with anat brain mask
            if data_label.endswith("desc-brain_mask") and filepath in outdir_func:
                data_label = data_label.replace("brain_mask", "bold_mask")

            try:
                pipe_x = rpool.get_pipe_number(pipe_idx)
            except ValueError:
                pipe_x = len(rpool.pipe_list)
            if filepath in outdir_anat:
                ingress = create_general_datasource(
                    f"gather_anat_outdir_{data_label!s}_{pipe_x}"
                )
                ingress.inputs.inputnode.set(
                    unique_id=unique_id,
                    data=filepath,
                    creds_path=creds_path,
                    dl_dir=cfg.pipeline_setup["working_directory"]["path"],
                )
                rpool.set_data(
                    data_label,
                    ingress,
                    "outputspec.data",
                    json_info,
                    pipe_idx,
                    node_name,
                    f"outdir_{data_label}_ingress",
                    inject=True,
                )
            else:
                if data_label.endswith("desc-preproc_bold"):
                    func_key = data_label
                    func_dict[bidstag] = {}
                    func_dict[bidstag]["scan"] = str(filepath)
                    func_dict[bidstag]["scan_parameters"] = json_info
                    func_dict[bidstag]["pipe_idx"] = pipe_idx
                if data_label.endswith("desc-brain_mask"):
                    data_label = data_label.replace("brain_mask", "bold_mask")
                try:
                    func_paths[data_label].append(filepath)
                except:
                    func_paths[data_label] = []
                    func_paths[data_label].append(filepath)

    if func_dict:
        wf, rpool = func_outdir_ingress(
            wf,
            cfg,
            func_dict,
            rpool,
            unique_id,
            creds_path,
            part_id,
            func_key,
            func_paths,
        )

    if cfg.surface_analysis["freesurfer"]["ingress_reconall"]:
        rpool = ingress_freesurfer(
            wf, rpool, cfg, data_paths, unique_id, part_id, ses_id
        )
    return wf, rpool


def json_outdir_ingress(rpool, filepath, exts, data_label, json):
    desc_val = None
    for tag in data_label.split("_"):
        if "desc-" in tag:
            desc_val = tag
            break
    jsonpath = str(filepath)
    for ext in exts:
        jsonpath = jsonpath.replace(ext, "")
    jsonpath = f"{jsonpath}.json"

    if not os.path.exists(jsonpath):
        WFLOGGER.info(
            "\n\n[!] No JSON found for file %s.\nCreating %s..\n\n", filepath, jsonpath
        )
        json_info = {
            "Description": "This data was generated elsewhere and "
            "supplied by the user into this C-PAC run's "
            "output directory. This JSON file was "
            "automatically generated by C-PAC because a "
            "JSON file was not supplied with the data."
        }
        json_info = {**json_info, **json}
        write_output_json(json_info, jsonpath)
    else:
        json_info = read_json(jsonpath)
        json_info = {**json_info, **json}
    if "CpacProvenance" in json_info:
        if desc_val:
            # it's a C-PAC output, let's check for pipe_idx/strat integer
            # suffixes in the desc- entries.
            only_desc = str(desc_val)

            if only_desc[-1].isdigit():
                for idx in range(0, 3):
                    # let's stop at 3, please don't run >999 strategies okay?
                    if only_desc[-1].isdigit():
                        only_desc = only_desc[:-1]

                if only_desc[-1] == "-":
                    only_desc = only_desc.rstrip("-")
                else:
                    msg = (
                        "\n[!] Something went wrong with either "
                        "reading in the output directory or when "
                        "it was written out previously.\n\nGive "
                        "this to your friendly local C-PAC "
                        f"developer:\n\n{data_label!s}\n"
                    )
                    raise IOError(msg)

            # remove the integer at the end of the desc-* variant, we will
            # get the unique pipe_idx from the CpacProvenance below
            data_label = data_label.replace(desc_val, only_desc)

        # preserve cpac provenance/pipe_idx
        pipe_idx = rpool.generate_prov_string(json_info["CpacProvenance"])
        node_name = ""

    else:
        json_info["CpacProvenance"] = [f"{data_label}:Non-C-PAC Origin: {filepath}"]
        if "Description" not in json_info:
            json_info["Description"] = (
                "This data was generated elsewhere and "
                "supplied by the user into this C-PAC run's "
                "output directory. This JSON file was "
                "automatically generated by C-PAC because a "
                "JSON file was not supplied with the data."
            )
        pipe_idx = rpool.generate_prov_string(json_info["CpacProvenance"])
        node_name = f"{data_label}_ingress"

    return rpool, json_info, pipe_idx, node_name, data_label


def func_outdir_ingress(
    wf, cfg, func_dict, rpool, unique_id, creds_path, part_id, key, func_paths
):
    pipe_x = len(rpool.pipe_list)
    ingress = create_func_datasource(
        func_dict, rpool, f"gather_func_outdir_{key}_{pipe_x}"
    )
    ingress.inputs.inputnode.set(
        subject=unique_id,
        creds_path=creds_path,
        dl_dir=cfg.pipeline_setup["working_directory"]["path"],
    )
    rpool.set_data("subject", ingress, "outputspec.subject", {}, "", "func_ingress")
    ingress.get_node("inputnode").iterables = ("scan", list(func_dict.keys()))
    rpool.set_data(key, ingress, "outputspec.rest", {}, "", "func_ingress")

    rpool.set_data("scan", ingress, "outputspec.scan", {}, "", "func_ingress")
    rpool.set_data(
        "scan-params", ingress, "outputspec.scan_params", {}, "", "scan_params_ingress"
    )
    wf, rpool, diff, blip, fmap_rp_list = ingress_func_metadata(
        wf, cfg, rpool, func_dict, part_id, creds_path, key
    )

    # Have to do it this weird way to save the parsed BIDS tag & filepath
    mask_paths_key = (
        "desc-bold_mask"
        if "desc-bold_mask" in func_paths
        else "space-template_desc-bold_mask"
    )
    ts_paths_key = "pipeline-ingress_desc-confounds_timeseries"

    # Connect func data with approproate scan name
    iterables = pe.Node(
        Function(
            input_names=["scan", "mask_paths", "ts_paths"],
            output_names=["out_scan", "mask", "confounds"],
            function=set_iterables,
        ),
        name=f"set_iterables_{pipe_x}",
    )
    iterables.inputs.mask_paths = func_paths[mask_paths_key]
    iterables.inputs.ts_paths = func_paths[ts_paths_key]
    wf.connect(ingress, "outputspec.scan", iterables, "scan")

    for key in func_paths:
        if key in (mask_paths_key, ts_paths_key):
            ingress_func = create_general_datasource(f"ingress_func_data_{key}")
            ingress_func.inputs.inputnode.set(
                unique_id=unique_id,
                creds_path=creds_path,
                dl_dir=cfg.pipeline_setup["working_directory"]["path"],
            )
            wf.connect(iterables, "out_scan", ingress_func, "inputnode.scan")
            if key == mask_paths_key:
                wf.connect(iterables, "mask", ingress_func, "inputnode.data")
                rpool.set_data(
                    key, ingress_func, "inputnode.data", {}, "", f"outdir_{key}_ingress"
                )
            elif key == ts_paths_key:
                wf.connect(iterables, "confounds", ingress_func, "inputnode.data")
                rpool.set_data(
                    key, ingress_func, "inputnode.data", {}, "", f"outdir_{key}_ingress"
                )

    return wf, rpool


def set_iterables(scan, mask_paths=None, ts_paths=None):
    # match scan with filepath to get filepath
    mask_path = [path for path in mask_paths if scan in path]
    ts_path = [path for path in ts_paths if scan in path]

    return (scan, mask_path[0], ts_path[0])


def strip_template(data_label, dir_path, filename):
    json = {}
    # rename to template
    for prefix in ["space-", "from-", "to-"]:
        for bidstag in data_label.split("_"):
            if bidstag.startswith(prefix):
                template_key, template_val = bidstag.split("-")
                template_name, _template_desc = lookup_identifier(template_val)
                if template_name:
                    json["Template"] = template_val
                    data_label = data_label.replace(template_val, "template")
            elif bidstag.startswith("res-"):
                res_key, res_val = bidstag.split("-")
                json["Resolution"] = res_val
                data_label = data_label.replace(bidstag, "")
    if data_label.find("__"):
        data_label = data_label.replace("__", "_")
    return data_label, json


def ingress_pipeconfig_paths(cfg, rpool, unique_id, creds_path=None):
    # ingress config file paths
    # TODO: may want to change the resource keys for each to include one level up in the YAML as well

    import pandas as pd
    import pkg_resources as p

    template_csv = p.resource_filename("CPAC", "resources/cpac_templates.csv")
    template_df = pd.read_csv(template_csv, keep_default_na=False)

    for row in template_df.itertuples():
        key = row.Key
        val = row.Pipeline_Config_Entry
        val = cfg.get_nested(cfg, [x.lstrip() for x in val.split(",")])
        resolution = row.Intended_Resolution_Config_Entry
        desc = row.Description

        if not val:
            continue

        if resolution:
            res_keys = [x.lstrip() for x in resolution.split(",")]
            tag = res_keys[-1]
        json_info = {}

        if "$FSLDIR" in val:
            val = val.replace("$FSLDIR", cfg.pipeline_setup["system_config"]["FSLDIR"])
        if "$priors_path" in val:
            priors_path = (
                cfg.segmentation["tissue_segmentation"]["FSL-FAST"]["use_priors"][
                    "priors_path"
                ]
                or ""
            )
            if "$FSLDIR" in priors_path:
                priors_path = priors_path.replace(
                    "$FSLDIR", cfg.pipeline_setup["system_config"]["FSLDIR"]
                )
            val = val.replace("$priors_path", priors_path)
        if "${resolution_for_anat}" in val:
            val = val.replace(
                "${resolution_for_anat}",
                cfg.registration_workflows["anatomical_registration"][
                    "resolution_for_anat"
                ],
            )
        if "${func_resolution}" in val:
            val = val.replace(
                "${func_resolution}",
                cfg.registration_workflows["functional_registration"][
                    "func_registration_to_template"
                ]["output_resolution"][tag],
            )

        if desc:
            template_name, _template_desc = lookup_identifier(val)
            if template_name:
                desc = f"{template_name} - {desc}"
            json_info["Description"] = f"{desc} - {val}"
        if resolution:
            resolution = cfg.get_nested(cfg, res_keys)
            json_info["Resolution"] = resolution

            resampled_template = pe.Node(
                Function(
                    input_names=["resolution", "template", "template_name", "tag"],
                    output_names=["resampled_template"],
                    function=resolve_resolution,
                    as_module=True,
                ),
                name="resampled_" + key,
            )

            resampled_template.inputs.resolution = resolution
            resampled_template.inputs.template = val
            resampled_template.inputs.template_name = key
            resampled_template.inputs.tag = tag

            # the set_data below is set up a little differently, because we are
            # injecting and also over-writing already-existing entries
            # other alternative would have been to ingress into the
            # resampled_template node from the already existing entries, but we
            # didn't do that here
            rpool.set_data(
                key,
                resampled_template,
                "resampled_template",
                json_info,
                "",
                "template_resample",
            )  # pipe_idx (after the blank json {}) should be the previous strat that you want deleted! because you're not connecting this the regular way, you have to do it manually

        elif val:
            config_ingress = create_general_datasource(f"gather_{key}")
            config_ingress.inputs.inputnode.set(
                unique_id=unique_id,
                data=val,
                creds_path=creds_path,
                dl_dir=cfg.pipeline_setup["working_directory"]["path"],
            )
            rpool.set_data(
                key,
                config_ingress,
                "outputspec.data",
                json_info,
                "",
                f"{key}_config_ingress",
            )
    # templates, resampling from config
    """
    template_keys = [
        ("anat", ["network_centrality", "template_specification_file"]),
        ("anat", ["nuisance_corrections", "2-nuisance_regression",
                  "lateral_ventricles_mask"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "FSL-FAST", "use_priors",
          "CSF_path"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "FSL-FAST", "use_priors",
          "GM_path"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "FSL-FAST", "use_priors",
          "WM_path"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "Template_Based", "CSF"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "Template_Based", "GRAY"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "Template_Based", "WHITE"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T1w_ACPC_template"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T1w_brain_ACPC_template"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T2w_ACPC_template"]),
        ("anat", ["anatomical_preproc", "acpc_alignment", "T2w_brain_ACPC_template"])]

    def get_nested_attr(c, template_key):
        attr = getattr(c, template_key[0])
        keys = template_key[1:]

        def _get_nested(attr, keys):
            if len(keys) > 1:
                return (_get_nested(attr[keys[0]], keys[1:]))
            elif len(keys):
                return (attr[keys[0]])
            else:
                return (attr)

        return (_get_nested(attr, keys))

    def set_nested_attr(c, template_key, value):
        attr = getattr(c, template_key[0])
        keys = template_key[1:]

        def _set_nested(attr, keys):
            if len(keys) > 1:
                return (_set_nested(attr[keys[0]], keys[1:]))
            elif len(keys):
                attr[keys[0]] = value
            else:
                return (attr)

        return (_set_nested(attr, keys))

    for key_type, key in template_keys:
        attr = cfg.get_nested(cfg, key)
        if isinstance(attr, str) or attr == None:
            node = create_check_for_s3_node(
                key[-1],
                attr, key_type,
                data_paths['creds_path'],
                cfg.pipeline_setup['working_directory']['path'],
                map_node=False
            )
            cfg.set_nested(cfg, key, node)

    template_keys_in_list = [
        ("anat",
         ["segmentation", "tissue_segmentation", "ANTs_Prior_Based",
          "template_brain_list"]),
        ("anat",
         ["segmentation", "tissue_segmentation", "ANTs_Prior_Based",
          "template_segmentation_list"]),
    ]

    for key_type, key in template_keys_in_list:
        node = create_check_for_s3_node(
            key[-1],
            cfg.get_nested(cfg, key), key_type,
            data_paths['creds_path'],
            cfg.pipeline_setup['working_directory']['path'],
            map_node=True
        )
        cfg.set_nested(cfg, key, node)
    """

    return rpool

def ingress_all_data(wf, rpool, cfg, data_paths, unique_id, part_id, ses_id):




 #### One way to do it

    # for data in data_paths[1].iterrows():
    #     suffix = data[1]["ent__suffix"]
    #     datatype = data[1]["ent__datatype"]
    #     filepath = data[1]["finfo__file_path"]
    #     desc = data[1]["ent__desc"]

    #     data_flow = create_general_datasource(f"gather_{datatype}_{suffix}")
    #     data_flow.inputs.inputnode.set(
    #         unique_id=unique_id,
    #         data=filepath,
    #         creds_path=None,
    #         dl_dir=cfg.pipeline_setup["working_directory"]["path"],
    #     )
    #     rpool.set_data(
    #         f"{datatype}_{suffix}",
    #         data_flow,
    #         "outputspec.data",
    #         {},
    #         "",
    #         f"{datatype}_{suffix}_ingress",
    #     )

    return rpool

def initiate_rpool(wf, cfg, data_paths=None, part_id=None):
    """
    Initialize a new ResourcePool.

    data_paths format:
      {'anat': {
            'T1w': '{T1w path}',
            'T2w': '{T2w path}'
        },
       'creds_path': {None OR path to credentials CSV},
       'func': {
           '{scan ID}':
               {
                   'scan': '{path to BOLD}',
                   'scan_parameters': {scan parameter dictionary}
               }
       },
       'site_id': 'site-ID',
       'subject_id': 'sub-01',
       'unique_id': 'ses-1',
       'derivatives_dir': '{derivatives_dir path}'}
    """
    # TODO: refactor further, integrate with the ingress_data functionality
    # TODO: used for BIDS-Derivatives (below), and possible refactoring of
    # TODO: the raw data config to use 'T1w' label instead of 'anat' etc.

    if data_paths:
        part_id = data_paths[0][0]
        ses_id = data_paths[0][1]
        unique_id = f"{part_id}_{ses_id}"

    elif part_id:
        unique_id = part_id
        creds_path = None
    from .resource_pool import ResourcePool
    rpool = ResourcePool(name=unique_id, cfg=cfg)

    # if data_paths:
    #     rpool = ingress_all_data(
    #         wf, rpool, cfg, data_paths, unique_id, part_id, ses_id
    #     )
    rpool.build_rpool(data_paths)

    # grab any file paths from the pipeline config YAML
    # creds_path = None
    # rpool = ingress_pipeconfig_paths(cfg, rpool, unique_id, creds_path)

    # output files with 4 different scans
    resource_description = {
        "suffix": "T1w",
        "desc": "preproc",
        "space": "MNI152NLin6ASym"
    }
    resource_content = rpool.get_resource(resource_description)
    #print(dir(rpool.get_resource("T1w")[0]))
    #rpool.write_to_disk(cfg.pipeline_setup["working_directory"]["path"])
    #print(rpool.get_resource("T1w"))

    # Ensure the directory exists
    os.makedirs('/code/output', exist_ok=True)

    # Now, safely open the file. It will be created if it does not exist.
    with open('/code/output/output.txt', 'w') as file:
      
        # Write the content to the file
        file.write(str(resource_content))
    import sys
    sys.exit()


    return (wf, rpool)


def run_node_blocks(blocks, data_paths, cfg=None):
    import os

    from CPAC.pipeline import nipype_pipeline_engine as pe
    from CPAC.pipeline.engine import NodeBlock

    if not cfg:
        cfg = {
            "pipeline_setup": {
                "working_directory": {"path": os.getcwd()},
                "log_directory": {"path": os.getcwd()},
            }
        }

    # TODO: WE HAVE TO PARSE OVER UNIQUE ID'S!!!
    _, rpool = initiate_rpool(cfg, data_paths)

    wf = pe.Workflow(name="node_blocks")
    wf.base_dir = cfg.pipeline_setup["working_directory"]["path"]
    wf.config["execution"] = {
        "hash_method": "timestamp",
        "crashdump_dir": cfg.pipeline_setup["log_directory"]["path"],
    }

    run_blocks = []
    if rpool.check_rpool("desc-preproc_T1w"):
        WFLOGGER.info("Preprocessed T1w found, skipping anatomical preprocessing.")
    else:
        run_blocks += blocks[0]
    if rpool.check_rpool("desc-preproc_bold"):
        WFLOGGER.info("Preprocessed BOLD found, skipping functional preprocessing.")
    else:
        run_blocks += blocks[1]

    for block in run_blocks:
        wf = NodeBlock(
            block, debug=cfg["pipeline_setup", "Debugging", "verbose"]
        ).connect_block(wf, cfg, rpool)
    rpool.gather_pipes(wf, cfg)

    wf.run()


class NodeData:
    r"""Attribute access for ResourcePool.get_data outputs.

    Class to hold outputs of CPAC.pipeline.engine.ResourcePool().get_data(), so one can
    do ``node_data = strat_pool.node_data(resource)`` and have ``node_data.node`` and
    ``node_data.out`` instead of doing ``node, out = strat_pool.get_data(resource)``
    and needing two variables (``node`` and ``out``) to store that information.

    Also includes ``variant`` attribute providing the resource's self-keyed value
    within its ``CpacVariant`` dictionary.

    Examples
    --------
    >>> rp = ResourcePool()
    >>> rp.node_data(None)
    NotImplemented (NotImplemented)

    >>> rp.set_data('test',
    ...             pe.Node(Function(input_names=[]), 'test'),
    ...             'b', [], 0, 'test')
    >>> rp.node_data('test')
    test (b)
    >>> rp.node_data('test').out
    'b'

    >>> try:
    ...     rp.node_data('b')
    ... except LookupError as lookup_error:
    ...     print(str(lookup_error).strip().split('\n')[0].strip())
    [!] C-PAC says: None of the listed resources are in the resource pool:
    """

    # pylint: disable=too-few-public-methods
    def __init__(self, strat_pool=None, resource=None, **kwargs):
        self.node = NotImplemented
        self.out = NotImplemented
        if strat_pool is not None and resource is not None:
            self.node, self.out = strat_pool.get_data(resource, **kwargs)

    def __repr__(self):  # noqa: D105
        return f'{getattr(self.node, "name", str(self.node))} ({self.out})'
