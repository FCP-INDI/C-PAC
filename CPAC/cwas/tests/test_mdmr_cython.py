import os
import pytest


@pytest.mark.skip(reason='possibly deprecated')
def test_mdmr():

    from CPAC.cwas.mdmr import mdmr
    from CPAC.cwas.cwas import calc_cwas
    import numpy as np

    X = np.genfromtxt(os.path.join(os.path.dirname(__file__), 'X.csv'), delimiter=',')
    Y = np.genfromtxt(os.path.join(os.path.dirname(__file__), 'Y.csv'), delimiter=',')

    X = X.reshape((X.shape[0], X.shape[1], 1))

    F_value, p_value = calc_cwas(
        X, Y, np.array([0, 1, 2], dtype=int), 1000, [0])

    assert np.isclose(p_value.mean(), 1.0, rtol=0.1)
