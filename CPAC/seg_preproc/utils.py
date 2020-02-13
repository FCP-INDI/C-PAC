# Import packages
import os
import sys
import re
import subprocess
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


def pick_wm_prob_0(probability_maps):

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


def pick_wm_prob_1(probability_maps):

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


def pick_wm_prob_2(probability_maps):

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


def pick_wm_class_0(tissue_class_files):

    """
    Returns the csf tissu class file from the list of segmented tissue class files

    Parameters
    ----------

    tissue_class_files : list (string)
        List of tissue class files

    Returns
    -------

    file : string
        Path to segment_seg_0.nii.gz is returned

    """

    if isinstance(tissue_class_files, list):
        if len(tissue_class_files) == 1:
            tissue_class_files = tissue_class_files[0]
        for filename in tissue_class_files:
            if filename.endswith("seg_0.nii.gz"):
                return filename
    return None


def pick_wm_class_1(tissue_class_files):

    """
    Returns the gray matter tissue class file from the list of segmented tissue class files

    Parameters
    ----------

    tissue_class_files : list (string)
        List of tissue class files

    Returns
    -------

    file : string
        Path to segment_seg_1.nii.gz is returned

    """

    if isinstance(tissue_class_files, list):
        if len(tissue_class_files) == 1:
            tissue_class_files = tissue_class_files[0]
        for filename in tissue_class_files:
            if filename.endswith("seg_1.nii.gz"):
                return filename
    return None


def pick_wm_class_2(tissue_class_files):

    """
    Returns the white matter tissue class file from the list of segmented tissue class files

    Parameters
    ----------

    tissue_class_files : list (string)
        List of tissue class files

    Returns
    -------

    file : string
        Path to segment_seg_2.nii.gz is returned

    """

    if isinstance(tissue_class_files, list):
        if len(tissue_class_files) == 1:
            tissue_class_files = tissue_class_files[0]
        for filename in tissue_class_files:
            if filename.endswith("seg_2.nii.gz"):
                return filename
    return None

# This functionality is adapted from poldracklab/niworkflows:
# https://github.com/poldracklab/niworkflows/blob/master/niworkflows/interfaces/utils.py
# https://fmriprep.readthedocs.io/
# https://poldracklab.stanford.edu/
# We are temporarily maintaining our own copy for more granular control.

def mask_erosion(roi_mask = None, skullstrip_mask = None, mask_erosion_mm = None, mask_erosion_prop = None):

    """
    Returns eroded segment mask and skull-stripped brain mask 

    Parameters
    ----------

    roi_mask : string 
        Path to binarized segment mask

    skullstrip_mask : string 
        Path to skull-stripped brain mask 

    mask_erosion_prop : float
        Proportion of erosion skull-stripped brain mask

    Returns
    -------

    output_roi_mask : string
        Path to eroded segment mask

    eroded_skullstrip_mask : string
        Path to eroded skull-stripped brain mask 

    """
    skullstrip_mask_img =  nb.load(skullstrip_mask)
    skullstrip_mask_data = skullstrip_mask_img.get_fdata()

    roi_mask_img = nb.load(roi_mask)
    roi_mask_data = roi_mask_img.get_fdata()
    erode_in = (mask_erosion_mm is not None and mask_erosion_mm > 0 or
                mask_erosion_prop is not None and mask_erosion_prop < 1 and mask_erosion_prop > 0)
    if erode_in:
        if mask_erosion_mm:
            iter_n = max(int(mask_erosion_mm / max(skullstrip_mask_img.header.get_zooms())), 1)
            skullstrip_mask_data = nd.binary_erosion(skullstrip_mask_data, iterations=iter_n)
        else :
            orig_vol = np.sum(skullstrip_mask_data > 0)
            while np.sum(skullstrip_mask_data > 0) / (orig_vol*1.0) > mask_erosion_prop :
                skullstrip_mask_data = nd.binary_erosion(skullstrip_mask_data, iterations=1)
        
        roi_mask_data[~skullstrip_mask_data] = 0
        
    hdr = roi_mask_img.get_header()
    output_roi_mask_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                 affine=roi_mask_img.get_affine())
    output_roi_mask = os.path.join(os.getcwd(), 'segment_tissue_eroded_mask.nii.gz')
    output_roi_mask_img.to_filename(output_roi_mask)

    hdr = skullstrip_mask_img.get_header()
    output_skullstrip_mask_img = nb.Nifti1Image(skullstrip_mask_data, header=hdr,
                                 affine=skullstrip_mask_img.get_affine())
    eroded_skullstrip_mask = os.path.join(os.getcwd(), 'eroded_skullstrip_mask.nii.gz')

    output_skullstrip_mask_img.to_filename(eroded_skullstrip_mask)

    return output_roi_mask, eroded_skullstrip_mask


# This functionality is adapted from poldracklab/niworkflows:
# https://github.com/poldracklab/niworkflows/blob/master/niworkflows/interfaces/utils.py
# https://fmriprep.readthedocs.io/
# https://poldracklab.stanford.edu/
# We are temporarily maintaining our own copy for more granular control.

def erosion(roi_mask = None, erosion_mm = None, erosion_prop = None):


    """
    Returns eroded tissue segment mask 

    Parameters
    ----------

    roi_mask : string 
        Path to binarized segment (ROI) mask

    erosion_prop : float
        Proportion of erosion segment mask

    Returns
    -------

    eroded_roi_mask : string
        Path to eroded segment mask

    """

    roi_mask_img = nb.load(roi_mask)
    roi_mask_data = roi_mask_img.get_fdata()
    orig_vol = np.sum(roi_mask_data > 0)

    erode_out = (erosion_mm is not None and erosion_mm > 0 or
                 erosion_prop is not None and erosion_prop < 1 and erosion_prop > 0)
    if erode_out:
        if erosion_mm:
            iter_n = max(int(erosion_mm / max(roi_mask_img.header.get_zooms())), 1)
            iter_n = int(erosion_mm / max(roi_mask_img.header.get_zooms()))
            roi_mask_data = nd.binary_erosion(roi_mask_data, iterations=iter_n)
        else:
            while np.sum(roi_mask_data > 0) / (orig_vol*1.0) > erosion_prop :
                roi_mask_data = nd.binary_erosion(roi_mask_data, iterations=1)

    hdr = roi_mask_img.get_header()
    output_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                 affine=roi_mask_img.get_affine())
    eroded_roi_mask = os.path.join(os.getcwd(), 'segment_tissue_mask.nii.gz')

    output_img.to_filename(eroded_roi_mask)

    return eroded_roi_mask