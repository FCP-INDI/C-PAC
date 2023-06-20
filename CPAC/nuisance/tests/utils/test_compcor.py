import pytest

from CPAC.nuisance.utils import compcor


def test_TR_string_to_float():
    assert abs(compcor.TR_string_to_float('123') - 123.) < 1e-08
    assert abs(compcor.TR_string_to_float('123s') - 123.) < 1e-08
    assert abs(compcor.TR_string_to_float('123ms') - 123. / 1000) < 1e-08

    assert abs(compcor.TR_string_to_float('1.23') - 1.23) < 1e-08
    assert abs(compcor.TR_string_to_float('1.23s') - 1.23) < 1e-08
    assert abs(compcor.TR_string_to_float('1.23ms') - 1.23 / 1000) < 1e-08

    with pytest.raises(Exception):
        compcor.TR_string_to_float(None)

    with pytest.raises(Exception):
        compcor.TR_string_to_float(['123'])

    with pytest.raises(Exception):
        compcor.TR_string_to_float(123)

    with pytest.raises(Exception):
        compcor.TR_string_to_float('ms')
