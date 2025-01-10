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
"""Tests for CLI utilities."""

from logging import INFO
import os
from pathlib import Path
import sys

import click
import pytest
import semver

CPAC_DIR = str(Path(__file__).parent.parent.parent)
sys.path.append(CPAC_DIR)
DATA_DIR = os.path.join(CPAC_DIR, "dev", "circleci_data")

from CPAC.__main__ import utils as CPAC_main_utils  # noqa: E402


def _click_backport(command, key):
    """Switch back to underscores for older versions of click."""
    return _resolve_alias(command, key.replace("-", "_").replace("opt_out", "opt-out"))


def _compare_version(left: str, right: str) -> int:
    """Handle required verbosity after ``semver.compare`` deprecation."""
    return semver.Version.compare(
        semver.Version.parse(left), semver.Version.parse(right)
    )


try:
    _BACKPORT_CLICK = _compare_version(click.__version__, "7.0.0") < 0
except ValueError:
    try:
        _BACKPORT_CLICK = _compare_version(f"{click.__version__}.0", "7.0.0") < 0
    except ValueError:
        _BACKPORT_CLICK = True


def _resolve_alias(command, key):
    """Resolve alias if possible."""
    return command.resolve_alias(key) if hasattr(command, "resolve_alias") else key


@pytest.mark.parametrize("multiword_connector", ["-", "_"])
def test_build_data_config(caplog, cli_runner, multiword_connector):
    """
    Test CLI ``utils data-config new-settings-template``...

    ...and ``utils data_config new_settings_template``.
    """
    caplog.set_level(INFO)
    if multiword_connector == "-" and _BACKPORT_CLICK:
        return
    os.chdir(DATA_DIR)
    test_yaml = os.path.join(DATA_DIR, "data_settings.yml")
    _delete_test_yaml(test_yaml)
    if multiword_connector == "_":
        data_config = CPAC_main_utils.commands[
            _click_backport(CPAC_main_utils, "data-config")
        ]
        result = cli_runner.invoke(
            data_config.commands[_click_backport(data_config, "new-settings-template")]
        )
    else:
        result = cli_runner.invoke(
            CPAC_main_utils.commands["data-config"].commands["new-settings-template"]
        )

    assert result.exit_code == 0
    assert "\n".join(caplog.messages).startswith(
        "\nGenerated a default data_settings YAML file for editing"
    )
    assert os.path.exists(test_yaml)
    _delete_test_yaml(test_yaml)


def test_new_settings_template(bids_examples: Path, caplog, cli_runner):
    """Test CLI ``utils new-settings-template``."""
    caplog.set_level(INFO)
    os.chdir(CPAC_DIR)
    assert bids_examples.exists()

    result = cli_runner.invoke(
        CPAC_main_utils.commands[
            _click_backport(CPAC_main_utils, "data-config")
        ].commands["build"],
        [os.path.join(DATA_DIR, "data_settings_bids_examples_ds051_default_BIDS.yml")],
    )

    participant_yaml = os.path.join(DATA_DIR, "data_config_ds051.yml")
    group_yaml = os.path.join(DATA_DIR, "group_analysis_participants_ds051.txt")

    assert result.exit_code == 0
    assert "\n".join(caplog.messages).startswith(
        "\nGenerating data configuration file.."
    )
    assert os.path.exists(participant_yaml)
    assert os.path.exists(group_yaml)
    _delete_test_yaml(participant_yaml)
    _delete_test_yaml(group_yaml)


def test_repickle(cli_runner):  # noqa
    fn = "python_2_pickle.pkl"
    pickle_path = os.path.join(DATA_DIR, fn)
    backups = [_Backup(pickle_path), _Backup(f"{pickle_path}z")]

    result = cli_runner.invoke(CPAC_main_utils.commands["repickle"], [DATA_DIR])

    assert result.exit_code == 0
    assert (
        f"Converted pickle {fn} from a Python 2 pickle to a Python 3 "
        "pickle." in result.output
    )

    result = cli_runner.invoke(CPAC_main_utils.commands["repickle"], [DATA_DIR])
    assert result.exit_code == 0
    assert f"Pickle {fn} is a Python 3 pickle." in result.output

    [backup.restore() for backup in backups]


class _Backup:
    def __init__(self, filepath):
        self.path = filepath
        with open(self.path, "rb") as r:
            self.data = r.read()

    def restore(self):
        with open(self.path, "wb") as w:
            w.write(self.data)


def _delete_test_yaml(test_yaml):
    if os.path.exists(test_yaml):
        os.remove(test_yaml)


def _test_repickle(pickle_path, gzipped=False):
    # pylint: disable=import-outside-toplevel,unused-import
    backup = _Backup(pickle_path)
    if gzipped:
        import gzip  # noqa: F401
    backup.restore()
