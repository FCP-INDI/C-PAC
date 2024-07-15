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
"""C-PAC pipeline engine."""

import ast
import hashlib
import json
import os
from typing import Any, Optional

from nipype import config, logging  # type: ignore [import-untyped]

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.engine.nodeblock import NODEBLOCK_INPUTS, NodeBlockFunction
from CPAC.pipeline.engine.resource import ResourcePool
from CPAC.utils.configuration.configuration import Configuration
from CPAC.utils.monitoring import (
    getLogger,
    LOGTAIL,
    WARNING_FREESURFER_OFF_WITH_DATA,
    WFLOGGER,
)

PIPELINE_BLOCKS = list["NodeBlockFunction | PIPELINE_BLOCKS"]


class NodeBlock:
    def __init__(
        self,
        node_block_functions: NodeBlockFunction | PIPELINE_BLOCKS,
        debug: bool = False,
    ) -> None:
        """Create a ``NodeBlock`` from a list of py:class:`~CPAC.pipeline.engine.nodeblock.NodeBlockFunction`s."""
        if not isinstance(node_block_functions, list):
            node_block_functions = [node_block_functions]

        self.node_blocks: dict[str, Any] = {}

        for node_block_function in node_block_functions:  # <---- sets up the NodeBlock object in case you gave it a list of node blocks instead of a single one - for option forking.
            self.input_interface = []
            if isinstance(node_block_function, tuple):
                self.input_interface = node_block_function[1]
                node_block_function = node_block_function[0]  # noqa: PLW2901
                if not isinstance(self.input_interface, list):
                    self.input_interface = [self.input_interface]

            if not isinstance(node_block_function, NodeBlockFunction):
                # If the object is a plain function `__name__` will be more useful than `str()`
                obj_str = (
                    node_block_function.__name__  # type: ignore [attr-defined]
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
            self.outputs: dict[str, Optional[str]] = {}
            for out in node_block_function.outputs:
                self.outputs[out] = None

            self.options: list[str] | dict[str, Any] = ["base"]
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

    def connect_block(self, wf: pe.Workflow, cfg: Configuration, rpool: ResourcePool):
        debug = bool(cfg.pipeline_setup["Debugging"]["verbose"])  # type: ignore [attr-defined]
        all_opts: list[str] = []

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
            inputs: NODEBLOCK_INPUTS = self.check_null(block_dct["inputs"])
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
                # THIS ALSO MEANS the multiple option-val's in docstring node blocks can be entered once in the entire node-block sequence, not in a list of multiples
            if not opts:
                # for node blocks where the options are split into different
                # block functions - opts will be empty for non-selected
                # options, and would waste the get_strats effort below
                continue
            all_opts += opts

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
                                ast.literal_eval(str(pipe_idx))
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
                            new_json_info = strat_pool.json

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


def run_node_blocks(blocks, data_paths, cfg=None):
    from CPAC.pipeline.engine import NodeBlock
    from CPAC.pipeline.engine.resource import ResourcePool

    if not cfg:
        cfg = {
            "pipeline_setup": {
                "working_directory": {"path": os.getcwd()},
                "log_directory": {"path": os.getcwd()},
            }
        }

    # TODO: WE HAVE TO PARSE OVER UNIQUE ID'S!!!

    wf = pe.Workflow(name="node_blocks")
    rpool = ResourcePool(wf=wf, cfg=cfg, data_paths=data_paths)
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
