# Copyright (C) 2022-2025  C-PAC Developers

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
"""Tests for surface configuration."""

from pathlib import Path
from typing import cast

import pytest
import yaml

from CPAC.pipeline.cpac_pipeline import run_workflow
from CPAC.resources.configs import CONFIGS_PATH
from CPAC.utils.configuration import Configuration


@pytest.mark.skip(reason="timing out for unrelated reasons")
@pytest.mark.timeout(60)
def test_duplicate_freesurfer(tmp_path: Path) -> None:
    """The pipeline should build fast if freesurfer is not self-duplicating."""
    config = Configuration(yaml.safe_load("FROM: abcd-options"))
    with (CONFIGS_PATH / "data_config_S3-BIDS-ABIDE.yml").open("r") as data_config:
        sub_dict = yaml.safe_load(data_config)[0]
    for directory in ["output", "working", "log", "crash_log"]:
        directory_key = ["pipeline_setup", f"{directory}_directory", "path"]
        item = cast(str, config[directory_key])
        config[directory_key] = str(tmp_path / item.lstrip("/"))
    run_workflow(sub_dict, config, False, test_config=True)
