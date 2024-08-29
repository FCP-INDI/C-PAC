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
from scipy import signal
from scipy.linalg import svd

from CPAC.utils import safe_shape
from CPAC.utils.monitoring import IFLOGGER


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
    IFLOGGER.info("Removing zero variance components")
    image_data = image_data[image_data.std(1) != 0, :]

    if image_data.shape.count(0):
        err = (
            "\n\n[!] No wm or csf signals left after removing those "
            "with zero variance.\n\n"
        )
        raise Exception(err)

    IFLOGGER.info("Detrending and centering data")
    Y = signal.detrend(image_data, axis=1, type="linear").T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yc = Yc / np.tile(np.array(Yc.std(0)).reshape(1, Yc.shape[1]), (Yc.shape[0], 1))
    IFLOGGER.info("Calculating SVD decomposition of Y*Y'")
    U, S, Vh = np.linalg.svd(Yc, full_matrices=False)

    # write out the resulting regressor file
    regressor_file = os.path.join(os.getcwd(), "compcor_regressors.1D")
    np.savetxt(regressor_file, U[:, :num_components], delimiter="\t", fmt="%16g")

    return regressor_file


def cosine_filter(
    input_image_path,
    timestep,
    period_cut=128,
    remove_mean=True,
    axis=-1,
    failure_mode="error",
):
    """
    The function applies a cosine filter to the input BOLD image using the discrete cosine transform (DCT) method.
    
    Adapted from nipype implementation. https://github.com/nipy/nipype/blob/d353f0d/nipype/algorithms/confounds.py#L1086-L1107
    It removes the low-frequency drift from the voxel time series. The filtered image is saved to disk.


    Parameters
    ----------
    input_image_path : str
        Path to the BOLD image to be filtered.
    timestep : float
        Repetition time (TR) of the series (in seconds). Derived from image header if unspecified.
    period_cut : float, optional
        Minimum period (in seconds) for the DCT high-pass filter. Default value is 128.
    remove_mean : bool, optional
        Whether to remove the mean from the voxel time series before filtering. Default is True.
    axis : int, optional
        The axis along which to apply the filter. Default is -1 (last axis).
    failure_mode : {'error', 'ignore'}, optional
        Specifies how to handle failure modes. If set to 'error', the function raises an error.
        If set to 'ignore', it returns the input data unchanged in case of failure. Default is 'error'.

    Returns
    -------
    cosfiltered_img : str
        Path to the filtered BOLD image.

    """
    # STATEMENT OF CHANGES:
    #     This function is derived from sources licensed under the Apache-2.0 terms,
    #     and this function has been changed.

    # CHANGES:
    #     * Refactored to take and return filepaths instead of loaded data
    #     * Removed caluclation and return of `non_constant_regressors`
    #     * Modified docstring to reflect local changes
    #     * Updated style to match C-PAC codebase
    #     * Updated to use generator and iterate over voxel time series to optimize memory usage.

    # ORIGINAL WORK'S ATTRIBUTION NOTICE:
    #    Copyright (c) 2009-2016, Nipype developers

    #    Licensed under the Apache License, Version 2.0 (the "License");
    #    you may not use this file except in compliance with the License.
    #    You may obtain a copy of the License at

    #        http://www.apache.org/licenses/LICENSE-2.0

    #    Unless required by applicable law or agreed to in writing, software
    #    distributed under the License is distributed on an "AS IS" BASIS,
    #    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    #    See the License for the specific language governing permissions and
    #    limitations under the License.

    #    Prior to release 0.12, Nipype was licensed under a BSD license.

    # Modifications copyright (C) 2019 - 2024  C-PAC Developers

    def voxel_generator():
        for i in range(datashape[0]):
            for j in range(datashape[1]):
                for k in range(datashape[2]):
                    yield input_data[i, j, k, :]
    
    from nipype.algorithms.confounds import _cosine_drift, _full_rank
    
    input_img = nib.load(input_image_path)
    input_data = input_img.get_fdata()
    datashape = input_data.shape
    timepoints = datashape[axis]
    if datashape[0] == 0 and failure_mode != "error":
        return input_data, np.array([])

    frametimes = timestep * np.arange(timepoints)
    X = _full_rank(_cosine_drift(period_cut, frametimes))[0]

    # Reshape the input data to bring the time dimension to the last axis if it's not already
    if axis != -1:
        reshaped_data = np.moveaxis(input_data, axis, -1)
    else:
        reshaped_data = input_data

    reshaped_output_data = np.zeros_like(reshaped_data)

    voxel_gen = voxel_generator()

    for i in range(reshaped_data.shape[0]):
        print(f"calculating {i+1} of {reshaped_data.shape[0]} row of voxels")
        for j in range(reshaped_data.shape[1]):
            for k in range(reshaped_data.shape[2]):
                voxel_time_series = next(voxel_gen)
                betas = np.linalg.lstsq(X, voxel_time_series.T, rcond=None)[0]

                if not remove_mean:
                    X = X[:, :-1]
                    betas = betas[:-1]

                residuals = voxel_time_series - X.dot(betas)
                reshaped_output_data[i, j, k, :] = residuals

    # Move the time dimension back to its original position if it was reshaped
    if axis != -1:
        output_data = np.moveaxis(reshaped_output_data, -1, axis)
    else:
        output_data = reshaped_output_data

    hdr = input_img.header
    output_img = nib.Nifti1Image(output_data, header=hdr, affine=input_img.affine)
    file_name = input_image_path[input_image_path.rindex("/") + 1 :]

    cosfiltered_img = os.path.join(os.getcwd(), file_name)

    output_img.to_filename(cosfiltered_img)

    return cosfiltered_img


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
