# Copyright (C) 2022 - 2024  C-PAC Developers

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
"""Tests for bandpass filters."""

from importlib.abc import Traversable
from importlib.resources import files
from pathlib import Path

from numpy.typing import NDArray
import pytest

from CPAC.nuisance.bandpass import read_1D

RAW_ONE_D: Traversable = files("CPAC").joinpath("nuisance/tests/regressors.1D")


@pytest.mark.parametrize("start_line", list(range(6)))
def test_read_1D(start_line: int, tmp_path: Path) -> None:
    """Test the correct number of rows are read when reading a 1D file."""
    regressor: Path = tmp_path / f"regressor_startAtL{start_line}.1D"
    # create a regressor.1D file with (5 - ``start_line``) lines of header
    with (
        RAW_ONE_D.open("r", encoding="utf-8") as _raw,
        regressor.open("w", encoding="utf-8") as _test_file,
    ):
        for line in _raw.readlines()[start_line:]:
            _test_file.write(line)
    header: list[str]
    data: NDArray
    header, data = read_1D(regressor)
    # should get the same array no matter how many lines of header
    assert data.shape == (10, 29)
    # all header lines should be captured
    assert len(header) == 5 - start_line
