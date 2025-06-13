# Copyright (C) 2022 - 2025  C-PAC Developers

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
from os import getenv
from pathlib import Path

from networkx import DiGraph
import numpy as np
from numpy.typing import NDArray
import pytest
import nibabel as nib
from scipy.fft import fft

from CPAC.nuisance.bandpass import ideal_bandpass, read_1D
from CPAC.nuisance.nuisance import filtering_bold_and_regressors
from CPAC.nuisance.utils.utils import load_censor_tsv
from CPAC.pipeline.engine import ResourcePool
from CPAC.pipeline.nipype_pipeline_engine import Workflow
from CPAC.pipeline.test.test_engine import _download
from CPAC.utils.configuration import Preconfiguration
from CPAC.utils.tests.osf import download_file

RAW_ONE_D: Traversable = files("CPAC").joinpath("nuisance/tests/regressors.1D")


class TestResourcePool(ResourcePool):
    """ResourcePool with OSF download function."""

    def osf(self, resource: str, file: str, destination: Path, index: int) -> None:
        """Download a file from the Open Science Framework."""
        _download(self, resource, download_file, file, destination, index)


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


@pytest.mark.parametrize(
    "lowcut, highcut, in_freq, out_freq",
    [
        (0.005, 0.05, 0.01, 0.2),
        (0.01, 0.1, 0.02, 0.15),
        (0.02, 0.08, 0.04, 0.12),
    ],
)
def test_ideal_bandpass_with_various_cutoffs(lowcut, highcut, in_freq, out_freq):
    """Test the ideal bandpass filter with various cutoff frequencies."""
    sample_period = 1.0
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


@pytest.mark.skipif(
    not getenv("OSF_DATA"),
    reason="OSF API key not set in OSF_DATA environment variable",
)
def test_frequency_filter(tmp_path: Path) -> None:
    """Test that the bandpass filter works as expected."""
    cfg = Preconfiguration("benchmark-FNIRT")
    rpool = TestResourcePool(cfg)
    wf = Workflow("bandpass_filtering", base_dir=str(tmp_path))
    index = 0
    for resource, file in {
        "realigned_file": "residuals.nii.gz",
        "regressor_file": "regressors.1D",
    }.items():
        rpool.osf(resource, file, tmp_path, index)
        index += 1

    filt = filtering_bold_and_regressors(
        cfg["nuisance_corrections", "2-nuisance_regression", "Regressors"][0]
    )
    residuals = rpool.node_data("realigned_file")
    regressors = rpool.node_data("regressor_file")
    wf.connect(
        [
            (residuals.node, filt, [(residuals.out, "inputspec.functional_file_path")]),
            (
                regressors.node,
                filt,
                [(regressors.out, "inputspec.regressors_file_path")],
            ),
        ]
    )
    res: DiGraph = wf.run()
    out_node = next(iter(res.nodes))
    output = out_node.run()
    trs = nib.load(output.outputs.bandpassed_file).header["dim"][4]  # type: ignore[reportPrivateImportUsage]
    array = load_censor_tsv(output.outputs.regressor_file, trs)
    assert not all(
        [array.min() == 0, array.max() == 0, array.sum() == 0]
    ), "Bandpass filter filtered all signals."
