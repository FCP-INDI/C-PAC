"""Configs for testing."""

from importlib import resources
from importlib.resources.abc import Traversable

import yaml

_TEST_CONFIGS_PATH: Traversable = resources.files("CPAC").joinpath(
    "utils/tests/configs"
)
with (_TEST_CONFIGS_PATH / "neurostars_23786.yml").open("r", encoding="utf-8") as _f:
    # A loaded YAML file to test https://tinyurl.com/neurostars23786
    NEUROSTARS_23786 = _f.read()
with (_TEST_CONFIGS_PATH / "neurostars_24035.yml").open("r", encoding="utf-8") as _f:
    # A loaded YAML file to test https://tinyurl.com/neurostars24035
    NEUROSTARS_24035 = _f.read()
# A loaded YAML file to test https://tinyurl.com/cmicnlslack420349
SLACK_420349 = {
    "filepath": yaml.dump(
        {
            "FROM": str(_TEST_CONFIGS_PATH / "neurostars_23786.yml"),
            "pipeline_setup": {"pipeline_name": "slack_420349_filepath"},
        }
    ),
    "preconfig": yaml.dump(
        {
            "FROM": "preproc",
            "pipeline_setup": {"pipeline_name": "slack_420349_preconfig"},
        }
    ),
}
__all__ = ["NEUROSTARS_23786", "NEUROSTARS_24035", "SLACK_420349"]
