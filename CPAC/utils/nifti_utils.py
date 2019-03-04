import os
import numpy as np
import six

import nibabel as nib


def more_zeros_than_ones(image):
    """
    Return True is there is more zeros than other values in a given nifti image.
    Parameters
    ----------
    image: str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be inverted or
        the image already loaded through nibabel

    Returns
    -------
    more_zeros: boolean
    """
    if isinstance(image, nib.nifti1.Nifti1Image):
        img = image
    elif isinstance(image, six.string_types):
        if not os.path.exists(image):
            raise ValueError(str(image) + " does not exist.")
        else:
            img = nib.load(image)
    else:
        raise TypeError("Image can be either a string or a nifti1.Nifti1Image")

    data = img.get_data()
    nb_zeros = len(np.where(data == 0)[0])
    size = data.size
    more_zeros = nb_zeros > size - nb_zeros
    return more_zeros


def inverse_nifti_values(image):
    """
    Replace zeros by ones and non-zero values by 1
    Parameters
    ----------
    image: str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be inverted or
        the image already loaded through nibabel

    Returns
    -------
    output: Nibabel Nifti1Image
    """
    if isinstance(image, nib.nifti1.Nifti1Image):
        img = image
    elif isinstance(image, six.string_types):
        if not os.path.exists(image):
            raise ValueError(str(image) + " does not exist.")
        else:
            img = nib.load(image)
    else:
        raise TypeError("Image can be either a string or a nifti1.Nifti1Image")

    data = img.get_data()
    zeros = np.where(data)
    out_data = np.ones(data.shape)
    out_data[zeros] = 0

    return nib.nifti1.Nifti1Image(out_data, img.affine)
