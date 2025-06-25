# Copyright (C) 2018-2025  C-PAC Developers

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
"""MDMR Cython tests."""

from importlib.resources import as_file, files

import numpy as np
import pytest

from CPAC.cwas.cwas import calc_cwas


@pytest.mark.skip(reason="possibly deprecated")
def test_mdmr() -> None:
    with as_file(files("CPAC").joinpath("cwas/tests")) as _f:
        X = np.genfromtxt(_f / "X.csv", delimiter=",")
        Y = np.genfromtxt(_f / "Y.csv", delimiter=",")

    X = X.reshape((X.shape[0], X.shape[1], 1))

    _F_value, p_value = calc_cwas(X, Y, np.array([0, 1, 2], dtype=int), 1000, [0])
    assert np.isclose(p_value.mean(), 1.0, rtol=0.1)
