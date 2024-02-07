import os

import numpy as np
import six
import nibabel as nib


def nifti_image_input(image):
    """
    Test if an input is a path or a nifti.image and the image loaded through
    nibabel.

    Parameters.
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
    elif isinstance(image, six.string_types):
        if not os.path.exists(image):
            raise ValueError(str(image) + " does not exist.")
        else:
            img = nib.load(image)
    else:
        msg = "Image can be either a string or a nifti1.Nifti1Image"
        raise TypeError(msg)
    return img


def more_zeros_than_ones(image):
    """
    Return True if there are more zeros than other values in a given nifti image.

    Parameters.
    ----------
    image : str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be inverted or
        the image already loaded through nibabel

    Returns
    -------
    more_zeros : boolean
    """
    if isinstance(image, nib.nifti1.Nifti1Image):
        img = image
    elif isinstance(image, six.string_types):
        if not os.path.exists(image):
            raise ValueError(str(image) + " does not exist.")
        else:
            img = nib.load(image)
    else:
        msg = "Image can be either a string or a nifti1.Nifti1Image"
        raise TypeError(msg)

    data = img.get_fdata()
    nb_zeros = len(np.where(data == 0)[0])
    size = data.size
    return nb_zeros > size - nb_zeros


def inverse_nifti_values(image):
    """
    Replace zeros by ones and non-zero values by 1.

    Parameters.
    ----------
    image : str or nibabel.nifti1.Nifti1Image
        path to the nifti file to be inverted or
        the image already loaded through nibabel

    Returns
    -------
    output : Nibabel Nifti1Image
    """
    if isinstance(image, nib.nifti1.Nifti1Image):
        img = image
    elif isinstance(image, six.string_types):
        if not os.path.exists(image):
            raise ValueError(str(image) + " does not exist.")
        else:
            img = nib.load(image)
    else:
        msg = "Image can be either a string or a nifti1.Nifti1Image"
        raise TypeError(msg)

    data = img.get_fdata()
    zeros = np.where(data)
    out_data = np.ones(data.shape)
    out_data[zeros] = 0

    return nib.nifti1.Nifti1Image(out_data, img.affine)
