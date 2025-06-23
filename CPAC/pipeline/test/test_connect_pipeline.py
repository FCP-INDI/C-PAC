#!/usr/bin/env python3
# Copyright (C) 2025  C-PAC Developers

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
"""Test pipeline connections."""

from logging import INFO
from pathlib import Path
from typing import Callable

import pytest
import yaml

from CPAC.pipeline.cpac_runner import run
from CPAC.utils.configuration import preconfig_yaml
from CPAC.utils.monitoring import log_nodes_cb


@pytest.mark.parametrize("preconfig", ["abcd-options"])
def test_config(
    caplog: pytest.LogCaptureFixture, preconfig: str, tmp_path: Path
) -> None:
    """Run 'test_config' analysis level."""
    caplog.set_level(INFO)
    data_config_file = tmp_path / "data_config.yaml"
    with data_config_file.open("w") as _f:
        yaml.dump(
            [
                {
                    "anat": "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/KKI/sub-1019436/ses-1/anat/sub-1019436_ses-1_run-1_T1w.nii.gz",
                    "func": {
                        "rest_acq-1_run-1": {
                            "scan": "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/KKI/sub-1019436/ses-1/func/sub-1019436_ses-1_task-rest_acq-1_run-1_bold.nii.gz",
                            "scan_parameters": "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/KKI/task-rest_acq-1_bold.json",
                        }
                    },
                    "site": "KKI",
                    "subject_id": "1019436",
                    "unique_id": "1",
                }
            ],
            _f,
        )
    plugin = "MultiProc"
    plugin_args: dict[str, int | bool | Callable] = {
        "n_procs": 2,
        "memory_gb": 10,
        "raise_insufficient": True,
        "status_callback": log_nodes_cb,
    }
    tracking = False
    exitcode = run(
        str(data_config_file),
        preconfig_yaml(preconfig),
        plugin=plugin,
        plugin_args=plugin_args,
        tracking=tracking,
        test_config=True,
    )
    if exitcode != 0:
        records = list(caplog.records)
        msg: str
        msg = str(records[-1])
        if hasattr(records[-1], "exc_info"):
            exc_info = records[-1].exc_info
            if (
                exc_info
                and exc_info[0]
                and exc_info[1]
                and hasattr(exc_info[1], "args")
            ):
                msg = exc_info[1].args[0]
                raise exc_info[0](exc_info[1])
        raise AssertionError(msg)
