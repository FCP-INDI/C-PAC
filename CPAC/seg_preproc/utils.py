# Import packages
import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import scipy.ndimage as nd
import numpy as np
import nibabel as nb

def check_if_file_is_empty(in_file):
    """
    Raise exception if regressor fie is empty.

    Parameters
    ----------

    in_file : nii file (string)
        regressor file

    Returns
    -------

    in_file : string
        return same file

    """
    import nibabel as nb
    import numpy as np
    nii = nb.load(in_file)
    data = nii.get_data()
    if data.size == 0 or np.all(data==0) or np.all(data==np.nan):
        raise ValueError('File {0} is empty. Use a lower threshold or turn '
                         'off regressors.'.format(in_file))
    return in_file


def pick_wm_0(probability_maps):

    """
    Returns the csf probability map from the list of segmented probability maps

    Parameters
    ----------

    probability_maps : list (string)
        List of Probability Maps

    Returns
    -------

    file : string
        Path to segment_prob_0.nii.gz is returned

    """

    if isinstance(probability_maps, list):
        if len(probability_maps) == 1:
            probability_maps = probability_maps[0]
        for filename in probability_maps:
            if filename.endswith("prob_0.nii.gz"):
                return filename
    return None


def pick_wm_1(probability_maps):

    """
    Returns the gray matter probability map from the list of segmented probability maps

    Parameters
    ----------

    probability_maps : list (string)
        List of Probability Maps

    Returns
    -------

    file : string
        Path to segment_prob_1.nii.gz is returned

    """

    if isinstance(probability_maps, list):
        if len(probability_maps) == 1:
            probability_maps = probability_maps[0]
        for filename in probability_maps:
            if filename.endswith("prob_1.nii.gz"):
                return filename
    return None


def pick_wm_2(probability_maps):

    """
    Returns the white matter probability map from the list of segmented probability maps

    Parameters
    ----------

    probability_maps : list (string)
        List of Probability Maps

    Returns
    -------

    file : string
        Path to segment_prob_2.nii.gz is returned

    """

    if isinstance(probability_maps, list):
        if len(probability_maps) == 1:
            probability_maps = probability_maps[0]
        for filename in probability_maps:
            if filename.endswith("prob_2.nii.gz"):
                return filename
    return None


def erosion(roi_mask, erosion_prop):
    roi_mask_img = nb.load(roi_mask)
    roi_mask_data = roi_mask_img.get_fdata()
    orig_vol = np.sum(roi_mask_data > 0)
    
    while np.sum(roi_mask_data > 0) / orig_vol > erosion_prop :
        roi_mask_data = nd.binary_erosion(roi_mask_data, iterations=1)
    
    hdr = roi_mask_img.get_header()
    output_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                 affine=roi_mask_img.get_affine())
    out_file = os.path.join(os.getcwd(), 'segment_tissue_mask.nii.gz')

    output_img.to_filename(out_file)

    return out_file



