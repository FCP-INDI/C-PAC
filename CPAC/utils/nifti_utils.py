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
"""Utlities for NIfTI images."""

import os
from typing import Literal, overload

import numpy as np
import nibabel as nib
from nipype.interfaces.afni import utils as afni_utils

from CPAC.pipeline import nipype_pipeline_engine as pe


def nifti_image_input(
    image: str | nib.nifti1.Nifti1Image,
) -> nib.nifti1.Nifti1Image:
    """Test if an input is a path or a nifti.image and the image loaded through nibabel.

    Parameters
    ----------
    image : str or nibabel.nifti1.Nifti1Image
        path to the nifti file or the image already loaded through nibabel

    Returns
    -------
    img : nibabel.nifti1.Nifti1Image
        load and return the nifti image if image is a path, otherwise simply
        return image
    """
    if isinstance(image, nib.nifti1.Nifti1Image):
        img = image
    elif isinstance(image, str):
        if not os.path.exists(image):
            msg = f"{image} does not exist."
            raise FileNotFoundError(msg)
        img = nib.load(image)
    else:
        msg = "Image can be either a string or a nifti1.Nifti1Image"
        raise TypeError(msg)
    return img


def more_zeros_than_ones(image):
    """
    Return True if there are more zeros than other values in a given nifti image.

    Parameters
    ----------
    image : str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be inverted or
        the image already loaded through nibabel

    Returns
    -------
    more_zeros : boolean
    """
    img = nifti_image_input(image)
    data = img.get_fdata()
    nb_zeros = len(np.where(data == 0)[0])
    size = data.size
    return nb_zeros > size - nb_zeros


def inverse_nifti_values(image):
    """
    Replace zeros by ones and non-zero values by 1.

    Parameters
    ----------
    image : str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be inverted or
        the image already loaded through nibabel

    Returns
    -------
    output : Nibabel Nifti1Image
    """
    img = nifti_image_input(image)
    data = img.get_fdata()
    zeros = np.where(data)
    out_data = np.ones(data.shape)
    out_data[zeros] = 0

    return nib.nifti1.Nifti1Image(out_data, img.affine)


@overload
def orientation_node(
    name: str,
    orientation: Literal["RPI", "LPI", "RAI", "LAI", "RAS", "LAS", "RPS", "LPS"],
    node_type: type[pe.MapNode],
) -> pe.MapNode: ...
@overload
def orientation_node(
    name: str,
    orientation: Literal["RPI", "LPI", "RAI", "LAI", "RAS", "LAS", "RPS", "LPS"],
    node_type: type[pe.Node],
) -> pe.Node: ...
def orientation_node(
    name: str,
    orientation: Literal["RPI", "LPI", "RAI", "LAI", "RAS", "LAS", "RPS", "LPS"],
    node_type: type[pe.Node | pe.MapNode] = pe.Node,
) -> pe.Node | pe.MapNode:
    """Return a node configured to resample an input with AFNI 3dresample."""
    kwargs = {
        "interface": afni_utils.Resample(
            orientation=orientation,
            outputtype="NIFTI_GZ",
        ),
        "name": name,
        "mem_gb": 0,
        "mem_x": (0.0115, "in_file", "t"),
    }
    if node_type == pe.MapNode:
        kwargs["iterfield"] = ["in_file", "out_file"]
    return node_type(**kwargs)
