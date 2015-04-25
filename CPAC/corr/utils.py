

def compute_corr(in_file, mask_file):

    """
    Computes the Network Correlation Matrix

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)

    Returns
    -------

    out_file : nifti file
        ReHo map of the input EPI image

    """

    import nibabel as nb
    import numpy as np
    import os
    import sys

    out_file = None

    res_fname = (in_file)
    res_mask_fname = (mask_file)
    CUTNUMBER = 10


    out_file = corr_file

    return out_file


