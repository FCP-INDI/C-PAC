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
"""Unit tests for the C-PAC pipeline engine."""

from argparse import Namespace
import os
from pathlib import Path
from typing import cast

from _pytest.logging import LogCaptureFixture
import pytest

from CPAC.pipeline.cpac_pipeline import (
    build_anat_preproc_stack,
    build_workflow,
    connect_pipeline,
    initialize_nipype_wf,
    load_cpac_pipe_config,
)
from CPAC.pipeline.engine import (
    ingress_pipeconfig_paths,
    ingress_raw_anat_data,
    ingress_raw_func_data,
    initiate_rpool,
    ResourcePool,
)
from CPAC.utils.bids_utils import create_cpac_data_config


@pytest.mark.skip(reason="not a pytest test")
def test_ingress_func_raw_data(pipe_config, bids_dir, test_dir):
    sub_data_dct = create_cpac_data_config(bids_dir, skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup["output_directory"]["path"] = os.path.join(test_dir, "out")
    cfg.pipeline_setup["working_directory"]["path"] = os.path.join(test_dir, "work")

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    part_id = sub_data_dct["subject_id"]
    ses_id = sub_data_dct["unique_id"]

    unique_id = f"{part_id}_{ses_id}"

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    if "func" in sub_data_dct:
        wf, rpool, diff, blip, fmap_rp_list = ingress_raw_func_data(
            wf, rpool, cfg, sub_data_dct, unique_id, part_id, ses_id
        )

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()


@pytest.mark.skip(reason="not a pytest test")
def test_ingress_anat_raw_data(pipe_config, bids_dir, test_dir):
    sub_data_dct = create_cpac_data_config(bids_dir, skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup["output_directory"]["path"] = os.path.join(test_dir, "out")
    cfg.pipeline_setup["working_directory"]["path"] = os.path.join(test_dir, "work")

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    part_id = sub_data_dct["subject_id"]
    ses_id = sub_data_dct["unique_id"]

    unique_id = f"{part_id}_{ses_id}"

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    rpool = ingress_raw_anat_data(
        wf, rpool, cfg, sub_data_dct, unique_id, part_id, ses_id
    )

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()


@pytest.mark.skip(reason="not a pytest test")
def test_ingress_pipeconfig_data(pipe_config, bids_dir, test_dir):
    sub_data_dct = create_cpac_data_config(bids_dir, skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup["output_directory"]["path"] = os.path.join(test_dir, "out")
    cfg.pipeline_setup["working_directory"]["path"] = os.path.join(test_dir, "work")
    cfg.pipeline_setup["log_directory"]["path"] = os.path.join(test_dir, "logs")

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    part_id = sub_data_dct["subject_id"]
    ses_id = sub_data_dct["unique_id"]

    unique_id = f"{part_id}_{ses_id}"

    rpool = ResourcePool(name=unique_id, cfg=cfg)

    wf, rpool = ingress_pipeconfig_paths(wf, cfg, rpool, sub_data_dct, unique_id)

    rpool.gather_pipes(wf, cfg, all=True)

    wf.run()


@pytest.mark.skip(reason="not a pytest test")
def test_build_anat_preproc_stack(pipe_config, bids_dir, test_dir):
    sub_data_dct = create_cpac_data_config(bids_dir, skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup["output_directory"]["path"] = os.path.join(test_dir, "out")
    cfg.pipeline_setup["working_directory"]["path"] = os.path.join(test_dir, "work")
    cfg.pipeline_setup["log_directory"]["path"] = os.path.join(test_dir, "logs")

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    wf, rpool = initiate_rpool(wf, cfg, sub_data_dct)

    pipeline_blocks = build_anat_preproc_stack(rpool, cfg)
    wf = connect_pipeline(wf, cfg, rpool, pipeline_blocks)

    rpool.gather_pipes(wf, cfg)

    wf.run()


@pytest.mark.skip(reason="not a pytest test")
def test_build_workflow(pipe_config, bids_dir, test_dir):
    sub_data_dct = create_cpac_data_config(bids_dir, skip_bids_validator=True)[0]
    cfg = load_cpac_pipe_config(pipe_config)

    cfg.pipeline_setup["output_directory"]["path"] = os.path.join(test_dir, "out")
    cfg.pipeline_setup["working_directory"]["path"] = os.path.join(test_dir, "work")
    cfg.pipeline_setup["log_directory"]["path"] = os.path.join(test_dir, "logs")

    wf = initialize_nipype_wf(cfg, sub_data_dct)

    wf, rpool = initiate_rpool(wf, cfg, sub_data_dct)

    wf, _, _ = build_workflow(sub_data_dct["subject_id"], sub_data_dct, cfg)

    rpool.gather_pipes(wf, cfg)

    wf.run()


def test_missing_resource(
    bids_examples: Path, caplog: LogCaptureFixture, tmp_path: Path
) -> None:
    """Test the error message thrown when a resource is missing."""
    from datetime import datetime

    import yaml

    from CPAC.pipeline.cpac_runner import run
    from CPAC.utils.bids_utils import sub_list_filter_by_labels
    from CPAC.utils.configuration import Preconfiguration, set_subject
    from CPAC.utils.configuration.yaml_template import create_yaml_from_template

    st = datetime.now().strftime("%Y-%m-%dT%H-%M-%SZ")
    namespace = Namespace(
        bids_dir=str(bids_examples / "ds113b"),
        output_dir=str(tmp_path / "output"),
        analysis_level="test_config",
        participant_label="sub-01",
    )
    c = Preconfiguration("anat-only")
    c["pipeline_setup", "output_directory", "path"] = namespace.output_dir
    c["pipeline_setup", "log_directory", "path"] = str(tmp_path / "logs")
    c["pipeline_setup", "working_directory", "path"] = str(tmp_path / "work")
    c["pipeline_setup", "system_config", "maximum_memory_per_participant"] = 1.0
    c["pipeline_setup", "system_config", "max_cores_per_participant"] = 1
    c["pipeline_setup", "system_config", "num_participants_at_once"] = 1
    c["pipeline_setup", "system_config", "num_ants_threads"] = 1
    c["pipeline_setup", "working_directory", "remove_working_dir"] = True
    sub_list = create_cpac_data_config(
        namespace.bids_dir,
        namespace.participant_label,
        None,
        True,
        only_one_anat=False,
    )
    sub_list = sub_list_filter_by_labels(list(sub_list), {"T1w": None, "bold": None})
    for i, sub in enumerate(sub_list):
        if isinstance(sub.get("anat"), dict):
            for anat_key in sub["anat"]:
                if isinstance(sub["anat"][anat_key], list) and len(
                    sub["anat"][anat_key]
                ):
                    sub_list[i]["anat"][anat_key] = sub["anat"][anat_key][0]
        if isinstance(sub.get("anat"), list) and len(sub["anat"]):
            sub_list[i]["anat"] = sub["anat"][0]
    data_config_file = f"cpac_data_config_{st}.yml"
    sublogdirs = [set_subject(sub, c)[2] for sub in sub_list]
    # write out the data configuration file
    data_config_file = os.path.join(sublogdirs[0], data_config_file)
    with open(data_config_file, "w", encoding="utf-8") as _f:
        noalias_dumper = yaml.dumper.SafeDumper
        noalias_dumper.ignore_aliases = lambda self, data: True
        yaml.dump(sub_list, _f, default_flow_style=False, Dumper=noalias_dumper)

    # update and write out pipeline config file
    pipeline_config_file = os.path.join(sublogdirs[0], f"cpac_pipeline_config_{st}.yml")
    with open(pipeline_config_file, "w", encoding="utf-8") as _f:
        _f.write(create_yaml_from_template(c))
    minimized_config = f"{pipeline_config_file[:-4]}_min.yml"
    with open(minimized_config, "w", encoding="utf-8") as _f:
        _f.write(create_yaml_from_template(c, import_from="blank"))
    for config_file in (data_config_file, pipeline_config_file, minimized_config):
        os.chmod(config_file, 0o444)  # Make config files readonly

    if len(sublogdirs) > 1:
        # If more than one run is included in the given data config
        # file, an identical copy of the data and pipeline config
        # will be included in the log directory for each run
        for sublogdir in sublogdirs[1:]:
            for config_file in (
                data_config_file,
                pipeline_config_file,
                minimized_config,
            ):
                try:
                    os.link(config_file, config_file.replace(sublogdirs[0], sublogdir))
                except FileExistsError:
                    pass

    run(
        data_config_file,
        pipeline_config_file,
        plugin="Linear",
        plugin_args={
            "n_procs": int(
                cast(
                    int | str,
                    c["pipeline_setup", "system_config", "max_cores_per_participant"],
                )
            ),
            "memory_gb": int(
                cast(
                    int | str,
                    c[
                        "pipeline_setup",
                        "system_config",
                        "maximum_memory_per_participant",
                    ],
                )
            ),
            "raise_insufficient": c[
                "pipeline_setup", "system_config", "raise_insufficient"
            ],
        },
        tracking=False,
        test_config=namespace.analysis_level == "test_config",
    )

    assert "can be output from" in caplog.text


# bids_dir = "/Users/steven.giavasis/data/HBN-SI_dataset/rawdata"
# test_dir = "/test_dir"

# cfg = "/Users/hecheng.jin/GitHub/DevBranch/CPAC/resources/configs/pipeline_config_monkey-ABCD.yml"

# test_ingress_func_raw_data(cfg, bids_dir, test_dir)
# test_ingress_anat_raw_data(cfg, bids_dir, test_dir)
# test_ingress_pipeconfig_data(cfg, bids_dir, test_dir)
# test_build_anat_preproc_stack(cfg, bids_dir, test_dir)
if __name__ == "__main__":
    cfg = "/Users/hecheng.jin/GitHub/pipeline_config_monkey-ABCDlocal.yml"
    bids_dir = "/Users/hecheng.jin/Monkey/monkey_data_oxford/site-ucdavis"
    test_dir = "/Users/hecheng.jin/GitHub/Test/T2preproc"
    test_build_workflow(cfg, bids_dir, test_dir)
