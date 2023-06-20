import os
import sys

from pathlib import Path

import click
import semver

CPAC_DIR = str(Path(__file__).parent.parent.parent)
sys.path.append(CPAC_DIR)
DATA_DIR = os.path.join(CPAC_DIR, 'dev', 'circleci_data')

from CPAC.__main__ import utils as CPAC_main_utils \
    # noqa: E402  # pylint: disable=wrong-import-position

try:
    _BACKPORT_CLICK = semver.compare(click.__version__, '8.0.0') < 0
except ValueError:
    try:
        _BACKPORT_CLICK = semver.compare(f'{click.__version__}.0', '8.0.0') < 0
    except ValueError:
        _BACKPORT_CLICK = True

def _click_backport(key):
    """Switch back to underscores for older versions of click"""
    return key.replace('-', '_') if _BACKPORT_CLICK else key


def test_build_data_config(cli_runner):
    os.chdir(DATA_DIR)
    test_yaml = os.path.join(DATA_DIR, "data_settings.yml")
    _delete_test_yaml(test_yaml)

    result = cli_runner.invoke(
        CPAC_main_utils.commands[_click_backport('data-config')].commands[
            _click_backport('new-settings-template')
        ]
    )

    assert result.exit_code == 0
    assert result.output.startswith(
        "\nGenerated a default data_settings YAML file for editing"
    )
    assert os.path.exists(test_yaml)
    _delete_test_yaml(test_yaml)


def test_new_settings_template(cli_runner):
    os.chdir(CPAC_DIR)

    example_dir = os.path.join(CPAC_DIR, 'bids-examples')
    if not os.path.exists(example_dir):
        from git import Repo
        Repo.clone_from(
            "https://github.com/bids-standard/bids-examples.git",
            example_dir
        )

    result = cli_runner.invoke(
        CPAC_main_utils.commands[_click_backport('data-config')
                                 ].commands['build'],
        [os.path.join(
            DATA_DIR,
            "data_settings_bids_examples_ds051_default_BIDS.yml"
        )]
    )

    participant_yaml = os.path.join(DATA_DIR, "data_config_ds051.yml")
    group_yaml = os.path.join(DATA_DIR,
                              "group_analysis_participants_ds051.txt")

    assert result.exit_code == 0
    assert result.output.startswith("\nGenerating data configuration file..")
    assert os.path.exists(participant_yaml)
    assert os.path.exists(group_yaml)
    _delete_test_yaml(participant_yaml)
    _delete_test_yaml(group_yaml)


def test_repickle(cli_runner):
    fn = 'python_2_pickle.pkl'
    pickle_path = os.path.join(DATA_DIR, fn)
    backups = [_Backup(pickle_path), _Backup(f'{pickle_path}z')]

    result = cli_runner.invoke(
        CPAC_main_utils.commands['repickle'],
        [DATA_DIR]
    )

    assert result.exit_code == 0
    assert (
        f'Converted pickle {fn} from a Python 2 pickle to a Python 3 '
        'pickle.' in result.output
    )

    result = cli_runner.invoke(
        CPAC_main_utils.commands['repickle'],
        [DATA_DIR]
    )
    assert result.exit_code == 0
    assert f'Pickle {fn} is a Python 3 pickle.' in result.output

    [backup.restore() for backup in backups]


class _Backup():
    def __init__(self, filepath):
        self.path = filepath
        with open(self.path, 'rb') as r:
            self.data = r.read()

    def restore(self):
        with open(self.path, 'wb') as w:
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
