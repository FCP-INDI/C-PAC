import os
import numpy as np

import nibabel as nib


def more_zeros_than_ones(image_path=None, nifti=None):
    """
    Return True is there is more zeros than other values in a given nifti image.
    Parameters
    ----------
    image_path: str
        path to the nifti file to be checked
    nifti: nibabel.nifti1.Nifti1Image
        the image already loaded through nibabel

    Returns
    -------
    more_zeros: boolean
    """
    if nifti is not None:
        img = nifti
    elif image_path is not None:
        if not os.path.exists(image_path):
            raise ValueError(str(image_path) + " does not exist.")
        else:
            img = nib.load(image_path)
    else:
        raise ValueError("Empty parameters")

    data = img.get_data()
    nb_zeros = len(np.where(data == 0)[0])
    size = data.size
    more_zeros = nb_zeros > size - nb_zeros
    return more_zeros


def inverse_nifti_values(image_path=None, nifti=None):
    """
    Replace zeros by ones and non-zero values by 1
    Parameters
    ----------
    image_path: str
        path to the nifti file to be inverted
    nifti: nibabel.nifti1.Nifti1Image
        the image already loaded through nibabel

    Returns
    -------
    output: Nibabel Nifti1Image
    """
    if nifti is not None:
        img = nifti
    elif image_path is not None:
        if not os.path.exists(image_path):
            raise ValueError(str(image_path) + " does not exist.")
        else:
            img = nib.load(image_path)
    else:
        raise ValueError("Empty parameters")

    data = img.get_data()
    zeros = np.where(data)
    out_data = np.ones(data.shape)
    out_data[zeros] = 0

    return nib.nifti1.Nifti1Image(out_data, img.affine)
