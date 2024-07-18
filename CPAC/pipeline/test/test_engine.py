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
"""Tests for C-PAC pipeline engine."""

from pathlib import Path

import pytest

from CPAC.pipeline.cpac_pipeline import (
    build_anat_preproc_stack,
    build_workflow,
)
from CPAC.pipeline.engine import ResourcePool
from CPAC.utils.bids_utils import create_cpac_data_config
from CPAC.utils.configuration import Configuration, Preconfiguration


def _set_up_test(
    bids_examples: Path, preconfig: str, tmp_path: Path
) -> tuple[Configuration, dict]:
    """Set up `cfg` and `sub_data` for engine tests."""
    bids_dir = str(bids_examples / "ds051")
    sub_data = create_cpac_data_config(bids_dir, skip_bids_validator=True)[0]
    cfg = Preconfiguration(preconfig)
    cfg.pipeline_setup["output_directory"]["path"] = str(tmp_path / "out")
    cfg.pipeline_setup["working_directory"]["path"] = str(tmp_path / "work")
    cfg.pipeline_setup["log_directory"]["path"] = str(tmp_path / "logs")
    return cfg, sub_data


@pytest.mark.parametrize("preconfig", ["default"])
def test_ingress_func_raw_data(
    bids_examples: Path, preconfig: str, tmp_path: Path
) -> None:
    """Test :py:method:`ResourcePool.ingress_raw_func_data` ."""
    cfg, sub_data_dct = _set_up_test(bids_examples, preconfig, tmp_path)
    rpool = ResourcePool(cfg=cfg, data_paths=sub_data_dct)
    rpool.gather_pipes(rpool.wf, cfg, all_types=True)


@pytest.mark.parametrize("preconfig", ["default"])
def test_ingress_anat_raw_data(
    bids_examples: Path, preconfig: str, tmp_path: Path
) -> None:
    """Test :py:method:`ResourcePool.ingress_raw_anat_data` ."""
    cfg, sub_data_dct = _set_up_test(bids_examples, preconfig, tmp_path)
    rpool = ResourcePool(
        cfg=cfg,
        data_paths=sub_data_dct,
    )
    rpool.ingress_raw_anat_data()
    rpool.gather_pipes(rpool.wf, cfg, all_types=True)


@pytest.mark.parametrize("preconfig", ["default"])
def test_ingress_pipeconfig_data(
    bids_examples: Path, preconfig: str, tmp_path: Path
) -> None:
    """Test :py:method:`ResourcePool.ingress_pipeconfig_paths` ."""
    cfg, sub_data_dct = _set_up_test(bids_examples, preconfig, tmp_path)
    rpool = ResourcePool(
        cfg=cfg,
        data_paths=sub_data_dct,
    )
    rpool.gather_pipes(rpool.wf, cfg, all_types=True)


@pytest.mark.parametrize("preconfig", ["anat-only"])
def test_build_anat_preproc_stack(
    bids_examples: Path, preconfig: str, tmp_path: Path
) -> None:
    """Test :py:func:`~build_anat_preproc_stack` ."""
    cfg, sub_data_dct = _set_up_test(bids_examples, preconfig, tmp_path)

    rpool = ResourcePool(cfg=cfg, data_paths=sub_data_dct)
    pipeline_blocks = build_anat_preproc_stack(rpool, cfg)
    wf = rpool.connect_pipeline(rpool.wf, cfg, pipeline_blocks)
    rpool.gather_pipes(wf, cfg)


@pytest.mark.parametrize("preconfig", ["default"])
def test_build_workflow(bids_examples: Path, preconfig: str, tmp_path: Path) -> None:
    """Test :py:func:`~build_workflow` ."""
    cfg, sub_data_dct = _set_up_test(bids_examples, preconfig, tmp_path)
    rpool = ResourcePool(cfg=cfg, data_paths=sub_data_dct)
    wf = build_workflow(sub_data_dct["subject_id"], sub_data_dct, cfg)
    rpool.gather_pipes(wf, cfg)
