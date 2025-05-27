# Copyright (C) 2019-2025  C-PAC Developers

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
"""Bandpass filtering utilities."""

import os
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
import nibabel as nib
from scipy.fftpack import fft, ifft


def ideal_bandpass(data, sample_period, bandpass_freqs):
    """
    Apply ideal bandpass filtering to a 1D time series data using FFT. Derived from YAN Chao-Gan 120504 based on REST.

    Parameters
    ----------
    data : NDArray
        1D time series data to be filtered.
    sample_period : float
        Length of sampling period in seconds.
    bandpass_freqs : tuple
        Tuple containing the bandpass frequencies (LowCutoff, HighCutoff).

    Returns
    -------
    NDArray
        Filtered time series data.

    """
    sample_freq = 1.0 / sample_period
    sample_length = data.shape[0]
    nyquist_freq = sample_freq / 2.0

    # Length of zero-padded data for efficient FFT
    N = int(2 ** np.ceil(np.log2(len(data))))
    data_p = np.zeros(N)
    data_p[:sample_length] = data

    LowCutoff, HighCutoff = bandpass_freqs

    if LowCutoff is None:  # No lower cutoff (low-pass filter)
        low_cutoff_i = 0
    elif LowCutoff > nyquist_freq:
        # Cutoff beyond fs/2 (all-stop filter)
        low_cutoff_i = int(N / 2)
    else:
        low_cutoff_i = np.ceil(LowCutoff * N * sample_period).astype("int")

    if HighCutoff > nyquist_freq or HighCutoff is None:
        # Cutoff beyond fs/2 or unspecified (become a highpass filter)
        high_cutoff_i = int(N / 2)
    else:
        high_cutoff_i = np.fix(HighCutoff * N * sample_period).astype("int")

    freq_mask = np.zeros_like(data_p, dtype="bool")
    freq_mask[low_cutoff_i : high_cutoff_i + 1] = True
    freq_mask[N - high_cutoff_i : N + 1 - low_cutoff_i] = True

    f_data = fft(data_p)
    f_data[~freq_mask] = 0.0
    return np.real_if_close(ifft(f_data)[:sample_length])


def read_1D(one_D: Path | str) -> tuple[list[str], NDArray]:
    """Parse a header from a 1D file, returing that header and a Numpy Array."""
    header = []
    with open(one_D, "r") as _f:
        # Each leading line that doesn't start with a number goes into the header
        for line in _f.readlines():
            try:
                float(line.split()[0])
                break
            except ValueError:
                header.append(line)

    regressor = np.loadtxt(one_D, skiprows=len(header))
    return header, regressor


def bandpass_voxels(realigned_file, regressor_file, bandpass_freqs, sample_period=None):
    """Perform ideal bandpass filtering on each voxel time-series.

    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    bandpass_freqs : tuple
        Tuple containing the bandpass frequencies. (LowCutoff_HighPass HighCutoff_LowPass)
    sample_period : float, optional
        Length of sampling period in seconds.  If not specified,
        this value is read from the nifti file provided.

    Returns
    -------
    bandpassed_file : string
        Path of filtered output (nifti file).

    """
    nii = nib.load(realigned_file)
    data = nii.get_fdata().astype("float64")
    mask = (data != 0).sum(-1) != 0
    Y = data[mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))

    if not sample_period:
        hdr = nii.header
        sample_period = float(hdr.get_zooms()[3])
        # Sketchy check to convert TRs in millisecond units
        if sample_period > 20.0:  # noqa: PLR2004
            sample_period /= 1000.0

    Y_bp = np.zeros_like(Y)
    for j in range(Y.shape[1]):
        Y_bp[:, j] = ideal_bandpass(Yc[:, j], sample_period, bandpass_freqs)

    data[mask] = Y_bp.T
    img = nib.Nifti1Image(data, header=nii.header, affine=nii.affine)
    bandpassed_file = os.path.join(os.getcwd(), "bandpassed_demeaned_filtered.nii.gz")
    img.to_filename(bandpassed_file)

    regressor_bandpassed_file = None

    if regressor_file is not None:
        if regressor_file.endswith(".nii.gz") or regressor_file.endswith(".nii"):
            nii = nib.load(regressor_file)
            data = nii.get_fdata().astype("float64")
            mask = (data != 0).sum(-1) != 0
            Y = data[mask].T
            Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
            Y_bp = np.zeros_like(Y)
            for j in range(Y.shape[1]):
                Y_bp[:, j] = ideal_bandpass(Yc[:, j], sample_period, bandpass_freqs)
            data[mask] = Y_bp.T

            img = nib.Nifti1Image(data, header=nii.header, affine=nii.affine)
            regressor_bandpassed_file = os.path.join(
                os.getcwd(), "regressor_bandpassed_demeaned_filtered.nii.gz"
            )
            img.to_filename(regressor_bandpassed_file)

        else:
            header: list[str]
            regressor: NDArray
            header, regressor = read_1D(regressor_file)
            Yc = regressor - np.tile(regressor.mean(0), (regressor.shape[0], 1))
            Y_bp = np.zeros_like(Yc)

            # Modify to allow just 1 regressor column
            shape = (
                regressor.shape[0] if len(regressor.shape) < 1 else regressor.shape[1]
            )
            for j in range(shape):
                Y_bp[:, j] = ideal_bandpass(Yc[:, j], sample_period, bandpass_freqs)

            regressor_bandpassed_file = os.path.join(
                os.getcwd(), "regressor_bandpassed_demeaned_filtered.1D"
            )
            with open(regressor_bandpassed_file, "w") as ofd:
                # write out the header information
                for line in header:
                    ofd.write(line)

                nuisance_regressors = np.array(Y_bp)
                np.savetxt(ofd, nuisance_regressors, fmt="%.18f", delimiter="\t")
    return bandpassed_file, regressor_bandpassed_file


def afni_1dBandpass(in_file, highpass, lowpass, tr=1):
    """
    Perform AFNI 1dBandpass.

    Parameters
    ----------
    in_file : string
        Path of an input 1D file
    highpass : float
        LowCutoff/HighPass
    lowpass : float
        HighCutoff/LowPass

    Returns
    -------
    out_file : string
        Path of an output 1D file
    """
    import os

    basename = os.path.basename(in_file)
    filename, file_extension = os.path.splitext(basename)
    out_file = os.path.join(os.getcwd(), filename + "_bp" + file_extension)

    cmd = "1dBandpass -dt %f %f %f %s > %s" % (tr, highpass, lowpass, in_file, out_file)
    os.system(cmd)

    return out_file
