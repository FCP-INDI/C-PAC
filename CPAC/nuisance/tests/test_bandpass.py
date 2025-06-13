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

import numpy as np
from numpy.typing import NDArray
import pytest
from scipy.fft import fft

from CPAC.nuisance.bandpass import ideal_bandpass, read_1D

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


@pytest.mark.parametrize("sample_period", [1.0, 1000.0])
@pytest.mark.parametrize(
    "lowcut, highcut, in_freq, out_freq",
    [
        (0.005, 0.05, 0.01, 0.2),
        (0.01, 0.1, 0.02, 0.15),
        (0.02, 0.08, 0.04, 0.12),
    ],
)
def test_ideal_bandpass_with_various_cutoffs(lowcut, highcut, in_freq, out_freq, sample_period):
    """Test the ideal bandpass filter with various cutoff frequencies."""
    t = np.arange(512) * sample_period
    signal = np.sin(2 * np.pi * in_freq * t) + np.sin(2 * np.pi * out_freq * t)

    filtered = ideal_bandpass(signal, sample_period, (lowcut, highcut))

    freqs = np.fft.fftfreq(len(signal), d=sample_period)
    orig_fft = np.abs(fft(signal))
    filt_fft = np.abs(fft(filtered))

    idx_in = np.argmin(np.abs(freqs - in_freq))
    idx_out = np.argmin(np.abs(freqs - out_freq))

    assert filt_fft[idx_in] > 0.5 * orig_fft[idx_in]
    assert filt_fft[idx_out] < 0.1 * orig_fft[idx_out]
