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
