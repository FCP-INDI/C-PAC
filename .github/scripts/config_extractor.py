import json
import os
import pathlib as pl
import sys
from typing import Literal

import yaml

from utils import cd, filesafe, multi_get


def _download_cpac_repo(cpac_dir: pl.Path, checkout_sha: str) -> None:
    """Downloads C-PAC configs from github and extracts them to the specified directory"""

    print(f"Check out C-PAC ({checkout_sha}) from github...")
    print("-------------------------------------------")
    os.system(f"git clone https://github.com/FCP-INDI/C-PAC.git {cpac_dir}")
    with cd(cpac_dir):
        if os.system(f'git checkout "{checkout_sha}"') != 0:
            print(f"Could not checkout {checkout_sha}")
            exit(1)
    print("-------------------------------------------")


def fetch_and_expand_cpac_configs(
    cpac_dir: pl.Path,
    output_dir: pl.Path,
    checkout_sha: str,
    config_names_ids: dict[str, str],
) -> None:
    """
    Fetches C-PAC configs from github, fully expands them (FROM: parent),
    and then saves them to the specified directory.
    """
    if not (cpac_dir / "CPAC").exists():
        cpac_dir.mkdir(parents=True, exist_ok=True)
        _download_cpac_repo(cpac_dir=cpac_dir, checkout_sha=checkout_sha)

    output_dir.mkdir(parents=True, exist_ok=True)

    cpac_module_path = str(cpac_dir.absolute())

    if cpac_module_path not in sys.path:
        sys.path.append(cpac_module_path)

    from CPAC.utils.configuration.configuration import Preconfiguration  # noqa
    from CPAC.utils.configuration.yaml_template import create_yaml_from_template  # noqa

    for config_name, config_id in config_names_ids.items():
        conf = Preconfiguration(config_id)
        config_yaml_string = create_yaml_from_template(conf.dict(), "blank")

        with open(output_dir / (filesafe(config_name) + ".yml"), "w", encoding="utf-8") as handle:
            handle.write(config_yaml_string)


def check_cpac_config(config: dict) -> tuple[Literal[True], None] | tuple[Literal[False], Exception]:
    """Checks if the specified file is a valid C-PAC config file"""
    from CPAC.utils.configuration.configuration import Configuration  # noqa

    try:
        Configuration(config)
    except Exception as e:
        return False, e
    return True, None


def get_cpac_config_ids() -> list[str]:
    from CPAC.pipeline import ALL_PIPELINE_CONFIGS
    return ALL_PIPELINE_CONFIGS

def fetch_and_expand_all_cpac_configs(
    cpac_dir: pl.Path,
    output_dir: pl.Path,
    checkout_sha: str,
):
    config_names_ids = {i: i for i in get_cpac_config_ids()}
    fetch_and_expand_cpac_configs(
        cpac_dir=cpac_dir,
        output_dir=output_dir,
        checkout_sha=checkout_sha,
        config_names_ids=config_names_ids,
    )

def _normalize_index_union(indices) -> list[list[str]]:
    if not indices:
        return []
    if isinstance(indices[0], list):
        return indices
    return [indices]

def normalize_index_union(indices) -> list[list[str]]:
    re = _normalize_index_union(indices)
    assert all(isinstance(item, list) for item in re)
    assert all(isinstance(i, str) for item in re for i in item)
    return re


if __name__ == "__main__":    
    import sys

    sys.path.append("cpac_source")

    CPAC_DIR = pl.Path("cpac_source")
    CONFIG_DIR = pl.Path("configs")
    CPAC_SHA = "dc41bf4f94da07dd78aeaf2fb894e11999f34748"


    with open("nodeblock_index.json") as f:
        nbs = json.load(f)

    configs = {}
    
    for config_path in CONFIG_DIR.glob("*.yml"):
        with open(config_path, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle)

        configs[config_path.stem] = config

    def _any_true_in_config(config, multi_index_union):  
        for path in paths:
            if multi_get(config, path):
                return True
        return False

    for nb in nbs:
        nb_configs = normalize_index_union(nb["decorator_args"].get("config"))
        nb_switchs = normalize_index_union(nb["decorator_args"].get("switch"))

        # multiply

        if not nb_configs and not nb_switchs:
            continue

        paths: list[list[str]] =  []
        if not nb_configs:
            paths = nb_switchs
        else:
            for nb_config in nb_configs:
                for nb_switch in nb_switchs:
                    paths.append(nb_config + nb_switch)

        assert all(isinstance(item, list) for item in paths), paths
        assert all(isinstance(i, str) for item in paths for i in item), paths

        configs_with_this_enabled = []
        for config_name, config in configs.items():
            if _any_true_in_config(config, paths):
                configs_with_this_enabled.append(config_name)

        nb['workflows'] = configs_with_this_enabled

    with open("nodeblock_index.json", "w", encoding="utf8") as handle:
        json.dump(nbs, handle, indent=2)
    




            
                


        

