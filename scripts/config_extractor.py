import json
import os
import pathlib as pl
import sys
from typing import Literal
import yaml
import base64
import hashlib
import re
from contextlib import contextmanager
from typing import Any, Generator, Sequence


@contextmanager
def cd(path: str | os.PathLike[str]) -> Generator[None, None, None]:
    """Context manager for changing the working directory"""
    old_wd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_wd)


def filesafe(s: str, replacement: str = "-") -> str:
    """
    Converts a string to a file safe string.
    Removes all non-alphanumeric characters and
    replaces them with the replacement string.
    """
    return re.sub(r"[^\w\d-]", replacement, s).lower()


def print_warning(msg: str) -> None:
    """Prints a colored warning message to the console"""
    print(f"\033[93mWARNING: {msg}\033[0m")


def multi_get(obj: dict, index: Sequence) -> Any | None:  # noqa: ANN401
    """
    Gets a value from a nested dictionary.
    Returns None if the path does not exist.
    """
    for i in index:
        if not isinstance(obj, dict) or i not in obj:
            return None
        obj = obj[i]
    return obj


def multi_set(obj: dict, index: Sequence, value: Any) -> bool:  # noqa: ANN401
    """
    Sets a value in a nested dictionary.
    Returns True if the path exists or was able to be created
    and the value was set.
    """
    for idx, i in enumerate(index):
        if not isinstance(obj, dict):
            return False

        if idx == len(index) - 1:
            obj[i] = value
            return True

        if i not in obj:
            obj[i] = {}

        obj = obj[i]
    assert False


def multi_del(obj: dict, index: Sequence) -> Any | None:  # noqa: ANN401
    """
    Deletes a value from a nested dictionary.
    Returns the value if the path exists and
    the value was deleted.
    """
    for idx, i in enumerate(index):
        if not isinstance(obj, dict):
            return None

        if idx == len(index) - 1:
            if i in obj:
                val = obj[i]
                del obj[i]
                return val
            return None

        if i not in obj:
            return None

        obj = obj[i]
    assert False


def aslist(obj: Any) -> list:  # noqa: ANN401
    """
    Converts an object to a list. If the object is
    already a list, it is returned as is.
    """
    if isinstance(obj, list):
        return obj
    return [obj]


def b64_urlsafe_hash(s: str) -> str:
    """
    Hashes a string and returns a base64 urlsafe encoded version of the hash.
    """
    return base64.urlsafe_b64encode(hashlib.sha1(s.encode()).digest()).decode().replace("=", "")


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

    sys.path.append(".")

    CPAC_DIR = pl.Path(".")
    CONFIG_DIR = pl.Path(".")
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
    




            
                


        

