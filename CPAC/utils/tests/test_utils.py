"""Tests of CPAC utility functions."""

import multiprocessing
from unittest import mock

from _pytest.logging import LogCaptureFixture
import pytest

from CPAC.func_preproc import get_motion_ref
from CPAC.pipeline.engine.nodeblock import NodeBlockFunction
from CPAC.utils.configuration import Configuration
from CPAC.utils.monitoring.custom_logging import log_subprocess
from CPAC.utils.tests import old_functions
from CPAC.utils.utils import (
    check_config_resources,
    check_system_deps,
    ScanParameters,
)

SCAN_PARAMS = {
    "BIDS": {
        "params": {
            "RepetitionTime": 2.0,
            "ScanOptions": "FS",
            "SliceAcquisitionOrder": "Interleaved Ascending",
        },
        "expected_TR": 2.0,
    },
    "C-PAC": {
        "params": {
            "tr": 2.5,
            "acquisition": "seq+z",
            "reference": "24",
            "first_TR": 1,
            "last_TR": "",
        },
        "expected_TR": 2.5,
    },
    "nested": {
        "params": {
            "TR": {"scan": 3},
            "first_TR": {"scan": 0},
            "last_TR": {"scan": 450},
        },
        "expected_TR": 3,
    },
}


def _installation_check(command: str, flag: str) -> None:
    """Test that command is installed by running specified version or help flag.

    Parameters
    ----------
    command : str
        Command to run.

    flag : str
        Version or help flag to run.

    Raises
    ------
    AssertionError
        If command is not installed.

    Returns
    -------
    None
    """
    _, exit_code = log_subprocess([command, flag])
    assert exit_code == 0


def test_check_config_resources():
    """Test check_config_resources function."""
    with (
        mock.patch.object(multiprocessing, "cpu_count", return_value=2),
        pytest.raises(SystemError) as system_error,
    ):
        check_config_resources(
            Configuration(
                {"pipeline_setup": {"system_config": {"max_cores_per_participant": 10}}}
            )
        )
    error_string = str(system_error.value)
    assert "threads running in parallel (10)" in error_string
    assert "threads available (2)" in error_string


@pytest.mark.parametrize("scan_params", ["BIDS", "C-PAC", "nested"])
@pytest.mark.parametrize("convert_to", [int, float, str])
def test_fetch_and_convert(
    caplog: LogCaptureFixture, scan_params: str, convert_to: type
) -> None:
    """Test functionality to fetch and convert scan parameters."""
    params = ScanParameters(SCAN_PARAMS[scan_params]["params"], "subject", "scan")
    TR = params.fetch_and_convert(
        keys=["TR", "RepetitionTime"],
        convert_to=convert_to,
    )
    if TR and "RepetitionTime" in params.params:
        old_TR = convert_to(
            old_functions.check(
                params.params, params.subject, params.scan, "RepetitionTime", False
            )
        )
        assert TR == old_TR
    try:
        old_TR = convert_to(
            old_functions.try_fetch_parameter(
                params.params, params.subject, params.scan, ["TR", "RepetitionTime"]
            )
        )
    except TypeError:
        old_TR = None
    assert (
        (TR == convert_to(SCAN_PARAMS[scan_params]["expected_TR"]))
        and isinstance(TR, convert_to)
        and TR == old_TR
    )
    if scan_params == "C-PAC":
        assert "Using case-insenitive match: 'TR' ≅ 'tr'." in caplog.text
    else:
        assert "Using case-insenitive match: 'TR' ≅ 'tr'." not in caplog.text
    not_TR = params.fetch_and_convert(
        keys=["NotTR", "NotRepetitionTime"],
        convert_to=convert_to,
    )
    assert not_TR is None
    if "first_TR" in params.params:
        first_tr = params.fetch_and_convert(["first_TR"], int, 1, False)
        old_first_tr = old_functions.check(
            params.params, params.subject, params.scan, "first_TR", False
        )
        if old_first_tr:
            old_first_tr = old_functions.check2(old_first_tr)
        assert first_tr == old_first_tr
    if "last_TR" in params.params:
        last_tr = params.fetch_and_convert(["last_TR"], int, "", False)
        old_last_tr = old_functions.check(
            params.params, params.subject, params.scan, "last_TR", False
        )
        if old_last_tr:
            old_last_tr = old_functions.check2(old_last_tr)
        assert last_tr == old_last_tr


@pytest.mark.parametrize("executable", ["Xvfb"])
def test_executable(executable):
    """Make sure executable is installed."""
    _installation_check(executable, "-help")


def test_NodeBlock_option_SSOT():  # pylint: disable=invalid-name
    """Test using NodeBlock dictionaries for SSOT for options."""
    assert isinstance(get_motion_ref, NodeBlockFunction)
    with pytest.raises(ValueError) as value_error:
        get_motion_ref(None, None, None, None, opt="chaos")
    error_message = str(value_error.value).rstrip()
    for opt in get_motion_ref.option_val:
        assert f"'{opt}'" in error_message
    assert error_message.endswith("Tool input: 'chaos'")


def test_system_deps():
    """Test system dependencies.

    Raises an exception if dependencies are not met.
    """
    check_system_deps(*([True] * 4))
