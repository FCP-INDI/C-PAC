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

import os

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.monitoring import WFLOGGER


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
    from CPAC.pipeline.engine.nodeblock import NodeBlock
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
        wf = rpool.connect_block(
            wf, NodeBlock(block, debug=cfg["pipeline_setup", "Debugging", "verbose"])
        )
    rpool.gather_pipes(wf, cfg)

    wf.run()
