"""Tests for C-PAC utilities"""
from CPAC.utils.utils import check_system_deps, try_fetch_parameter

scan_params_bids = {
    "RepetitionTime": 2.0,
    "ScanOptions": "FS",
    "SliceAcquisitionOrder": "Interleaved Ascending",
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


def test_system_deps():
    """Test system dependencies.
    Raises an exception if dependencies are not met.
    """
    check_system_deps(*([True] * 4))
