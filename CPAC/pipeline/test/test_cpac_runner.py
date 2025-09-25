# Copyright (C) 2021-2025  C-PAC Developers

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
"""Run C-PAC in a container."""

import os
from pathlib import Path

import pytest

from CPAC.pipeline.cpac_pipeline import load_cpac_pipe_config
from CPAC.pipeline.cpac_runner import run_T1w_longitudinal
from CPAC.pipeline.utils import get_shell
from CPAC.resources.configs import CONFIGS_PATH
from CPAC.utils.bids_utils import create_cpac_data_config


def test_shell() -> None:
    """Test that ``get_shell`` returns a path to an executable BASH."""
    shell: str = get_shell()
    assert shell.lower().endswith("bash"), "Default shell isn't BASH?"
    assert Path(shell).exists(), "No default shell found."
    assert os.access(shell, os.X_OK), "Default shell not executable."


@pytest.mark.skip(reason="not a pytest test")
def test_run_T1w_longitudinal(bids_dir, cfg, test_dir, part_id):
    sub_data_list = create_cpac_data_config(
        bids_dir, participant_labels=[part_id], skip_bids_validator=True
    )
    cfg = load_cpac_pipe_config(cfg)

    cfg.pipeline_setup["output_directory"]["path"] = os.path.join(test_dir, "out")
    cfg.pipeline_setup["working_directory"]["path"] = os.path.join(test_dir, "work")

    run_T1w_longitudinal(sub_data_list, cfg)


if __name__ == "__main__":
    bids_dir = "/Users/steven.giavasis/data/neurodata_hnu"
    test_dir = "/test_dir"
    part_id = "0025427"
    cfg = str(CONFIGS_PATH / "pipeline_config_default.yml")
    test_run_T1w_longitudinal(bids_dir, cfg, test_dir, part_id)
