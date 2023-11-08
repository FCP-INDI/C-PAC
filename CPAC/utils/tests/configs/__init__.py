"""Configs for testing"""
from pathlib import Path
from pkg_resources import resource_filename
import yaml

_TEST_CONFIGS_PATH = Path(resource_filename("CPAC", "utils/tests/configs"))
with open(_TEST_CONFIGS_PATH / "neurostars_23786.yml", "r", encoding="utf-8"
          ) as _f:
    # A loaded YAML file to test https://tinyurl.com/neurostars23786
    NEUROSTARS_23786 = _f.read()
with open(_TEST_CONFIGS_PATH / "neurostars_24035.yml", "r", encoding="utf-8"
          ) as _f:
    # A loaded YAML file to test https://tinyurl.com/neurostars24035
    NEUROSTARS_24035 = _f.read()
# A loaded YAML file to test https://tinyurl.com/cmicnlslack420349
SLACK_420349 = {
    "filepath": yaml.dump({
        "FROM": str(_TEST_CONFIGS_PATH / "neurostars_23786.yml"),
        "pipeline_setup": {"pipeline_name": "slack_420349_filepath"}}),
    "preconfig": yaml.dump({"FROM": "preproc",
        "pipeline_setup": {"pipeline_name": "slack_420349_preconfig"}})}
__all__ = ["NEUROSTARS_23786", "NEUROSTARS_24035", "SLACK_420349"]
