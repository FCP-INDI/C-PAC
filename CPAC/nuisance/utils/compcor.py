# Copyright (C) 2019-2024  C-PAC Developers

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
import os

import numpy as np
import nibabel as nib
from nibabel.filebasedimages import ImageFileError
from nipype import logging
from scipy import signal
from scipy.linalg import svd

from CPAC.utils import safe_shape

iflogger = logging.getLogger("nipype.interface")


def calc_compcor_components(data_filename, num_components, mask_filename):
    if num_components < 1:
        msg = f"Improper value for num_components ({num_components}), should be >= 1."
        raise ValueError(msg)

    try:
        image_data = nib.load(data_filename).get_fdata().astype(np.float64)
    except (ImageFileError, MemoryError, OSError, TypeError, ValueError) as e:
        msg = f"Unable to load data from {data_filename}"
        raise ImageFileError(msg) from e

    try:
        binary_mask = nib.load(mask_filename).get_fdata().astype(np.int16)
    except (ImageFileError, MemoryError, OSError, TypeError, ValueError) as e:
        msg = f"Unable to load data from {mask_filename}"
        raise ImageFileError(msg) from e

    if not safe_shape(image_data, binary_mask):
        msg = (
            f"The data in {data_filename} and {mask_filename} do not have a"
            " consistent shape"
        )
        raise ValueError(msg)

    # make sure that the values in binary_mask are binary
    binary_mask[binary_mask > 0] = 1
    binary_mask[binary_mask != 1] = 0

    # reduce the image data to only the voxels in the binary mask
    image_data = image_data[binary_mask == 1, :]

    # filter out any voxels whose variance equals 0
    iflogger.info("Removing zero variance components")
    image_data = image_data[image_data.std(1) != 0, :]

    if image_data.shape.count(0):
        err = (
            "\n\n[!] No wm or csf signals left after removing those "
            "with zero variance.\n\n"
        )
        raise Exception(err)

    iflogger.info("Detrending and centering data")
    Y = signal.detrend(image_data, axis=1, type="linear").T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yc = Yc / np.tile(np.array(Yc.std(0)).reshape(1, Yc.shape[1]), (Yc.shape[0], 1))
    iflogger.info("Calculating SVD decomposition of Y*Y'")
    U, S, Vh = np.linalg.svd(Yc, full_matrices=False)

    # write out the resulting regressor file
    regressor_file = os.path.join(os.getcwd(), "compcor_regressors.1D")
    np.savetxt(regressor_file, U[:, :num_components], delimiter="\t", fmt="%16g")

    return regressor_file


# cosine_filter adapted from nipype 'https://github.com/nipy/nipype/blob/d353f0d879826031334b09d33e9443b8c9b3e7fe/nipype/algorithms/confounds.py'
def cosine_filter(
    input_image_path,
    timestep,
    period_cut=128,
    remove_mean=True,
    axis=-1,
    failure_mode="error",
):
    """
    input_image_path: string
            Bold image to be filtered.
    timestep: float
            'Repetition time (TR) of series (in sec) - derived from image header if unspecified'
    period_cut: float
            Minimum period (in sec) for DCT high-pass filter, nipype default value: 128.

    """
    from CPAC.nuisance.utils.compcor import _cosine_drift, _full_rank

    input_img = nib.load(input_image_path)
    input_data = input_img.get_fdata()

    datashape = input_data.shape
    timepoints = datashape[axis]
    if datashape[0] == 0 and failure_mode != "error":
        return input_data, np.array([])

    input_data = input_data.reshape((-1, timepoints))

    frametimes = timestep * np.arange(timepoints)
    X = _full_rank(_cosine_drift(period_cut, frametimes))[0]
    X[:, :-1] if X.shape[1] > 1 else np.array([])

    betas = np.linalg.lstsq(X, input_data.T)[0]

    if not remove_mean:
        X = X[:, :-1]
        betas = betas[:-1]

    residuals = input_data - X.dot(betas).T

    output_data = residuals.reshape(datashape)

    hdr = input_img.header
    output_img = nib.Nifti1Image(output_data, header=hdr, affine=input_img.affine)

    file_name = input_image_path[input_image_path.rindex("/") + 1 :]

    cosfiltered_img = os.path.join(os.getcwd(), file_name)

    output_img.to_filename(cosfiltered_img)

    return cosfiltered_img


# _cosine_drift and _full_rank copied from nipype 'https://github.com/nipy/nipype/blob/d353f0d879826031334b09d33e9443b8c9b3e7fe/nipype/algorithms/confounds.py'
def _cosine_drift(period_cut, frametimes):
    """
     Create a cosine drift matrix with periods greater or equals to period_cut.

    Parameters.
    ----------
    period_cut : float
         Cut period of the low-pass filter (in sec)
    frametimes : array of shape(nscans)
         The sampling times (in sec)

    Returns
    -------
    cdrift :  array of shape(n_scans, n_drifts)
             cosin drifts plus a constant regressor at cdrift[:,0]
    Ref: http://en.wikipedia.org/wiki/Discrete_cosine_transform DCT-II
    """
    len_tim = len(frametimes)
    n_times = np.arange(len_tim)
    hfcut = 1.0 / period_cut  #  input parameter is the period

    # frametimes.max() should be (len_tim-1)*dt
    dt = frametimes[1] - frametimes[0]
    # hfcut = 1/(2*dt) yields len_time
    # If series is too short, return constant regressor
    order = max(int(np.floor(2 * len_tim * hfcut * dt)), 1)
    cdrift = np.zeros((len_tim, order))
    nfct = np.sqrt(2.0 / len_tim)

    for k in range(1, order):
        cdrift[:, k - 1] = nfct * np.cos((np.pi / len_tim) * (n_times + 0.5) * k)

    cdrift[:, order - 1] = 1.0  #  or 1./sqrt(len_tim) to normalize
    return cdrift


def _full_rank(X, cmax=1e15):
    """
    This function possibly adds a scalar matrix to X
    to guarantee that the condition number is smaller than a given threshold.

    Parameters.
    ----------
    X : array of shape(nrows, ncols)
    cmax=1.e-15, float tolerance for condition number

    Returns
    -------
    X : array of shape(nrows, ncols) after regularization
    cmax=1.e-15, float tolerance for condition number
    """
    U, s, V = fallback_svd(X, full_matrices=False)
    smax, smin = s.max(), s.min()
    c = smax / smin
    if c < cmax:
        return X, c
    iflogger.warning("Matrix is singular at working precision, regularizing...")
    lda = (smax - cmax * smin) / (cmax - 1)
    s = s + lda
    X = np.dot(U, np.dot(np.diag(s), V))
    return X, cmax


def fallback_svd(a, full_matrices=True, compute_uv=True):
    try:
        return np.linalg.svd(a, full_matrices=full_matrices, compute_uv=compute_uv)
    except np.linalg.LinAlgError:
        pass

    return svd(
        a, full_matrices=full_matrices, compute_uv=compute_uv, lapack_driver="gesvd"
    )


def TR_string_to_float(tr):
    """
    Convert TR string to seconds (float). Suffixes 's' or 'ms' to indicate
    seconds or milliseconds.

    Parameters
    ----------
    tr : TR string representation. May use suffixes 's' or 'ms' to indicate
    seconds or milliseconds.

    Returns
    -------
    tr in seconds (float)
    """
    if not isinstance(tr, str):
        msg = f"Improper type for TR_string_to_float ({tr})."
        raise TypeError(msg)

    tr_str = tr.replace(" ", "")

    try:
        if tr_str.endswith("ms"):
            tr_numeric = float(tr_str[:-2]) * 0.001
        elif tr.endswith("s"):
            tr_numeric = float(tr_str[:-1])
        else:
            tr_numeric = float(tr_str)
    except Exception as exc:
        msg = f'Can not convert TR string to float: "{tr}".'
        raise ValueError(msg) from exc

    return tr_numeric
