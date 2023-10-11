"""Tests of CPAC utility functions"""
import multiprocessing
from unittest import mock
import pytest
from CPAC.func_preproc import get_motion_ref
from CPAC.pipeline.nodeblock import NodeBlockFunction
from CPAC.utils.configuration import Configuration
from CPAC.utils.monitoring.custom_logging import log_subprocess
from CPAC.utils.utils import check_config_resources, check_system_deps, \
                             try_fetch_parameter

scan_params_bids = {
    'RepetitionTime': 2.0,
    'ScanOptions': 'FS',
    'SliceAcquisitionOrder': 'Interleaved Ascending',
}
scan_params_cpac = {
    'tr': 2.5,
    'acquisition': 'seq+z',
    'reference': '24',
    'first_tr': '',
    'last_tr': '',
}


def _installation_check(command: str, flag: str) -> None:
    """Test that command is installed by running specified version or
    help flag.

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
    with mock.patch.object(multiprocessing, 'cpu_count', return_value=2), \
         pytest.raises(SystemError) as system_error:
        check_config_resources(Configuration({'pipeline_setup': {
            'system_config': {'max_cores_per_participant': 10}}}))
    error_string = str(system_error.value)
    assert 'threads running in parallel (10)' in error_string
    assert 'threads available (2)' in error_string


def test_function():
    TR = try_fetch_parameter(scan_params_bids, '0001', 'scan',
                             ['TR', 'RepetitionTime'])
    assert TR == 2.0

    TR = try_fetch_parameter(scan_params_cpac, '0001', 'scan',
                             ['TR', 'RepetitionTime'])
    assert TR == 2.5


@pytest.mark.parametrize("executable", ["Xvfb"])
def test_executable(executable):
    """Make sure executable is installed"""
    _installation_check(executable, "-help")


def test_NodeBlock_option_SSOT():  # pylint: disable=invalid-name
    '''Test using NodeBlock dictionaries for SSOT for options'''
    assert isinstance(get_motion_ref, NodeBlockFunction)
    with pytest.raises(ValueError) as value_error:
        get_motion_ref(None, None, None, None, opt='chaos')
    error_message = str(value_error.value).rstrip()
    for opt in get_motion_ref.option_val:
        assert f"'{opt}'" in error_message
    assert error_message.endswith('Tool input: \'chaos\'')


def test_system_deps():
    """Test system dependencies.
    Raises an exception if dependencies are not met.
    """
    check_system_deps(*([True] * 4))
