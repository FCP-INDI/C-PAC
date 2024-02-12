# Import packages
import os
from CPAC.pipeline import nipype_pipeline_engine as pe
import scipy.ndimage as nd
import numpy as np
import nibabel as nb


def check_if_file_is_empty(in_file):
    """Raise exception if regressor fie is empty.

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
    data = nii.get_fdata()
    if data.size == 0 or np.all(data == 0) or np.all(data == np.nan):
        raise ValueError('File {0} is empty. Use a lower threshold or turn '
                         'off regressors.'.format(in_file))
    return in_file


def _erode(roi_mask, erosion_mm, erosion_prop):
    """Function to perform in-common erosion steps

    Parameters
    ----------
    mask : str
        Path to mask

    erosion_mm : int, float or None
        Kernel width in mm

    erosion_prop : float or None
        Target volume ratio, 0 < erosion_prop < 1

    Returns
    -------
    mask_img : nibabel image
        Original mask image

    erode : bool
        erosion_bool or erosion_mm

    mask_data : numpy array
        eroded mask data
    """
    mask_img = nb.load(roi_mask)
    mask_data = mask_img.get_fdata()
    orig_vol = np.sum(mask_data > 0)
    erode = ((erosion_mm is not None and erosion_mm > 0) or
             (erosion_prop is not None and 0 < erosion_prop < 1))
    if erode:
        if erosion_mm:
            iter_n = max(int(erosion_mm / max(mask_img.header.get_zooms())), 1)
            mask_data = nd.binary_erosion(mask_data, iterations=iter_n)
        else:
            while np.sum(mask_data > 0) / (orig_vol * 1.0) > erosion_prop:
                mask_data = nd.binary_erosion(mask_data, iterations=1)
    return mask_img, erode, mask_data


def pick_wm_prob_0(probability_maps):
    """Returns the csf probability map from the list of segmented
    probability maps

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
    """Returns the gray matter probability map from the list of segmented probability maps

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
    """Returns the white matter probability map from the list of segmented probability maps

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
    """Returns the csf tissue class file from the list of segmented tissue class files

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
    """Returns the gray matter tissue class file from the list of segmented tissue class files

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
    """Returns the white matter tissue class file from the list of segmented tissue class files

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


def mask_erosion(roi_mask=None, skullstrip_mask=None, mask_erosion_mm=None,
                 mask_erosion_prop=None):
    """Returns eroded segment mask and skull-stripped brain mask

    # This functionality is adapted from poldracklab/niworkflows:
    #   https://github.com/poldracklab/niworkflows/blob/master/niworkflows/interfaces/utils.py
    #   https://fmriprep.readthedocs.io/
    #   https://poldracklab.stanford.edu/
    # We are temporarily maintaining our own copy for more granular control.

    Parameters
    ----------

    roi_mask : string
        Path to binarized segment mask

    skullstrip_mask : string
        Path to skull-stripped brain mask

    mask_erosion_prop : float
        Target volume ratio for skull-stripped brain mask

    Returns
    -------

    output_roi_mask : string
        Path to eroded segment mask

    eroded_skullstrip_mask : string
        Path to eroded skull-stripped brain mask
    """
    roi_mask_img = nb.load(roi_mask)
    roi_mask_data = roi_mask_img.get_fdata()
    skullstrip_mask_img, erode_in, skullstrip_mask_data = _erode(
        skullstrip_mask, mask_erosion_mm, mask_erosion_prop)

    if erode_in:
        # pylint: disable=invalid-unary-operand-type
        roi_mask_data[~skullstrip_mask_data] = 0

    hdr = roi_mask_img.header
    output_roi_mask_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                         affine=roi_mask_img.affine)
    output_roi_mask = os.path.join(os.getcwd(),
                                   'segment_tissue_eroded_mask.nii.gz')
    output_roi_mask_img.to_filename(output_roi_mask)

    hdr = skullstrip_mask_img.header
    output_skullstrip_mask_img = nb.Nifti1Image(
        skullstrip_mask_data, header=hdr,
        affine=skullstrip_mask_img.affine)
    eroded_skullstrip_mask = os.path.join(os.getcwd(),
                                          'eroded_skullstrip_mask.nii.gz')

    output_skullstrip_mask_img.to_filename(eroded_skullstrip_mask)

    return output_roi_mask, eroded_skullstrip_mask


# This functionality is adapted from poldracklab/niworkflows:
# https://github.com/poldracklab/niworkflows/blob/master/niworkflows/interfaces/utils.py
# https://fmriprep.readthedocs.io/
# https://poldracklab.stanford.edu/
# We are temporarily maintaining our own copy for more granular control.
def erosion(roi_mask=None, erosion_mm=None, erosion_prop=None):
    """Returns eroded tissue segment mask

    Parameters
    ----------
    roi_mask : string
        Path to binarized segment (ROI) mask

    erosion_prop : float
        Target volume ratio for erosion segment mask

    Returns
    -------
    eroded_roi_mask : string
        Path to eroded segment mask
    """
    roi_mask_img, _, roi_mask_data = _erode(roi_mask, erosion_mm, erosion_prop)

    hdr = roi_mask_img.header
    output_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                affine=roi_mask_img.affine)
    eroded_roi_mask = os.path.join(os.getcwd(), 'segment_tissue_mask.nii.gz')

    output_img.to_filename(eroded_roi_mask)

    return eroded_roi_mask


def hardcoded_antsJointLabelFusion(anatomical_brain, anatomical_brain_mask,
                                   template_brain_list,
                                   template_segmentation_list):
    """run antsJointLabelFusion.sh

    Parameters
    ----------

    anatomical_brain : string (nifti file)
        Target image to be labeled.

    anatomical_brain_mask : string (nifti file)
        Target mask image

    template_brain_list : list
        Atlas to be warped to target image.

    template_segmentation_list : list
        Labels corresponding to atlas.

    Returns
    -------

    multiatlas_Intensity : string (nifti file)

    multiatlas_Labels : string (nifti file)
    """
    import os
    import subprocess

    cmd = ["${ANTSPATH}${ANTSPATH:+/}antsJointLabelFusion.sh"]
    cmd.append(
        " -d 3 -o ants_multiatlas_ -t {0} -x {1} -y b -c 0".format(
            anatomical_brain, anatomical_brain_mask))

    if (not len(template_brain_list) == len(template_segmentation_list)):
        err_msg = '\n\n[!] C-PAC says: Please check ANTs Prior-based ' \
                  'Segmentation setting. For performing ANTs Prior-based ' \
                  'segmentation method the number of specified ' \
                  'segmentations should be identical to the number of atlas ' \
                  'image sets.\n\n'
        raise Exception(err_msg)
    else:
        for index in range(len(template_brain_list)):
            cmd.append(
                " -g {0} -l {1}".format(
                    template_brain_list[index],
                    template_segmentation_list[index]))

    # write out the actual command-line entry for testing/validation later
    command_file = os.path.join(os.getcwd(), 'command.txt')
    with open(command_file, 'wt') as f:
        f.write(' '.join(cmd))

    str = ""
    bash_cmd = str.join(cmd)

    try:
        retcode = subprocess.check_output(bash_cmd, shell=True) \
            # noqa: F841  # pylint: disable=unused-variable
    except Exception as e:  # pylint: disable=broad-except,invalid-name
        # pylint: disable=raise-missing-from
        raise Exception('[!] antsJointLabel segmentation method did not '
                        'complete successfully.\n\nError '
                        'details:\n{0}\n{1}\n'.format(
                            e,
                            getattr(e, 'output', '')))

    multiatlas_Intensity = None
    multiatlas_Labels = None

    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:
        if "Intensity" in f:
            multiatlas_Intensity = os.getcwd() + "/" + f
        if "Labels" in f:
            multiatlas_Labels = os.getcwd() + "/" + f

    if not multiatlas_Labels:
        raise Exception("\n\n[!] No multiatlas labels file found. "
                        "antsJointLabelFusion may not have completed "
                        "successfully.\n\n")

    return multiatlas_Intensity, multiatlas_Labels


def pick_tissue_from_labels_file(multiatlas_Labels, csf_label=[4,14,15,24,43],
                                 gm_label=[3,42], wm_label=[2,41]):
    """Pick tissue mask from multiatlas labels file
    based off of FreeSurferColorLUT https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
    or user provided label value

    Parameters
    ----------
    multiatlas_Labels : string (nifti file)

    csf_label : list
        a list of integer label values corresponding to CSF in multiatlas file

    gm_label : list
        a list of integer label value corresponding to Gray Matter in multiatlas file

    wm_label : list
        a list of integer label value corresponding to White Matter in multiatlas file

    Returns
    -------
    csf_mask : string (nifti file)

    gm_mask : string (nifti file)

    wm_mask : string (nifti file)
    """
    # pylint: disable=import-outside-toplevel,redefined-outer-name,reimported
    import os
    import nibabel as nb
    import numpy as np

    img = nb.load(multiatlas_Labels)
    data = img.get_fdata()

    # pick tissue mask from multiatlas labels file
    # based off of FreeSurferColorLUT or user provided label values
    # hard-coded csf/gm/wm label values are based off of FreeSurferColorLUT

    # FreeSurfer Ventricle Labels:
    # Left-Lateral-Ventricle 4, 3rd-Ventricle 14, 4th-Ventricle 15, Right-Lateral-Ventricle 43

    csf = np.zeros(np.size(data))
    csf[np.where(np.in1d(data, np.array(csf_label)))] = 1
    csf = csf.reshape(data.shape)

    gm = np.zeros(np.size(data))
    gm[np.where(np.in1d(data, np.array(gm_label)))] = 1
    gm = gm.reshape(data.shape)

    wm = np.zeros(np.size(data))
    wm[np.where(np.in1d(data, np.array(wm_label)))] = 1
    wm = wm.reshape(data.shape)

    save_img_csf = nb.Nifti1Image(csf, header=img.header, affine=img.affine)
    save_img_gm = nb.Nifti1Image(gm, header=img.header, affine=img.affine)
    save_img_wm = nb.Nifti1Image(wm, header=img.header, affine=img.affine)

    save_img_csf.to_filename('csf_mask.nii.gz')
    save_img_gm.to_filename('gm_mask.nii.gz')
    save_img_wm.to_filename('wm_mask.nii.gz')

    csf_mask = os.path.join(os.getcwd(), 'csf_mask.nii.gz')
    gm_mask = os.path.join(os.getcwd(), 'gm_mask.nii.gz')
    wm_mask = os.path.join(os.getcwd(), 'wm_mask.nii.gz')

    return csf_mask, gm_mask, wm_mask
