"""Tests of CPAC utility functions"""
import pytest
from CPAC.func_preproc.func_preproc import get_motion_ref
from CPAC.utils.docs import grab_docstring_dct
from CPAC.utils.utils import check_system_deps, try_fetch_parameter

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

