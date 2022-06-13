"""Tests of CPAC utility functions"""
import pytest
from CPAC.func_preproc.func_preproc import get_motion_ref
from CPAC.utils.configuration import Configuration
from CPAC.utils.docs import grab_docstring_dct
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


def test_NodeBlock_option_SSOT():  # pylint: disable=invalid-name
    '''Test using NodeBlock dictionaries for SSOT for options'''
    nodebock_opts = grab_docstring_dct(get_motion_ref).get('option_val')
    with pytest.raises(ValueError) as value_error:
        get_motion_ref(None, None, None, None, opt='chaos')
    error_message = str(value_error.value).rstrip()
    for opt in nodebock_opts:
        assert f"'{opt}'" in error_message
    assert error_message.endswith('Tool input: \'chaos\'')


def test_system_deps():
    """Test system dependencies.
    Raises an exception if dependencies are not met.
    """
    check_system_deps(*([True] * 4))

