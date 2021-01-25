# Import packages
import os
import nipype.pipeline.engine as pe
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
    data = nii.get_data()
    if data.size == 0 or np.all(data == 0) or np.all(data == np.nan):
        raise ValueError('File {0} is empty. Use a lower threshold or turn '
                         'off regressors.'.format(in_file))
    return in_file


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
    """  # noqa
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
    """  # noqa
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
    """  # noqa
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
    """  # noqa
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
    """  # noqa
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
    """  # noqa
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
        Proportion of erosion skull-stripped brain mask

    Returns
    -------

    output_roi_mask : string
        Path to eroded segment mask

    eroded_skullstrip_mask : string
        Path to eroded skull-stripped brain mask
    """
    skullstrip_mask_img = nb.load(skullstrip_mask)
    skullstrip_mask_data = skullstrip_mask_img.get_fdata()

    roi_mask_img = nb.load(roi_mask)
    roi_mask_data = roi_mask_img.get_fdata()
    erode_in = (mask_erosion_mm is not None and mask_erosion_mm > 0 or
                mask_erosion_prop is not None and mask_erosion_prop < 1 and
                mask_erosion_prop > 0)
    if erode_in:
        if mask_erosion_mm:
            iter_n = max(
                int(mask_erosion_mm / max(
                    skullstrip_mask_img.header.get_zooms()
                )), 1)
            skullstrip_mask_data = nd.binary_erosion(
                skullstrip_mask_data, iterations=iter_n)
        else:
            orig_vol = np.sum(skullstrip_mask_data > 0)
            while (
                np.sum(skullstrip_mask_data > 0) / (orig_vol*1.0) >
                mask_erosion_prop
            ):
                skullstrip_mask_data = nd.binary_erosion(
                    skullstrip_mask_data, iterations=1)

        roi_mask_data[~skullstrip_mask_data] = 0

    hdr = roi_mask_img.get_header()
    output_roi_mask_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                         affine=roi_mask_img.get_affine())
    output_roi_mask = os.path.join(os.getcwd(),
                                   'segment_tissue_eroded_mask.nii.gz')
    output_roi_mask_img.to_filename(output_roi_mask)

    hdr = skullstrip_mask_img.get_header()
    output_skullstrip_mask_img = nb.Nifti1Image(
        skullstrip_mask_data, header=hdr,
        affine=skullstrip_mask_img.get_affine())
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
                 erosion_prop is not None and erosion_prop < 1 and
                 erosion_prop > 0)
    if erode_out:
        if erosion_mm:
            iter_n = max(
                int(erosion_mm / max(roi_mask_img.header.get_zooms())),
                1
            )
            iter_n = int(erosion_mm / max(roi_mask_img.header.get_zooms()))
            roi_mask_data = nd.binary_erosion(roi_mask_data, iterations=iter_n)
        else:
            while np.sum(roi_mask_data > 0) / (orig_vol*1.0) > erosion_prop:
                roi_mask_data = nd.binary_erosion(roi_mask_data, iterations=1)

    hdr = roi_mask_img.get_header()
    output_img = nb.Nifti1Image(roi_mask_data, header=hdr,
                                affine=roi_mask_img.get_affine())
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

    anatomical_brain_mask: string (nifti file)
        Target mask image

    template_brain_list: list
        Atlas to be warped to target image.

    template_segmentation_list: list
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
        retcode = subprocess.check_output(bash_cmd, shell=True)  # noqa F841
    except Exception as e:
        raise Exception('[!] antsJointLabel segmentation method did not '
                        'complete successfully.\n\nError '
                        'details:\n{0}\n{1}\n'.format(e, e.output))

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


def pick_tissue_from_labels_file(multiatlas_Labels, csf_label=24,
                                 left_gm_label=3, left_wm_label=2,
                                 right_gm_label=42, right_wm_label=41,
                                 include_ventricles=False):
    """Pick tissue mask from multiatlas labels file
    based off of FreeSurferColorLUT https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
    or user provided label value

    Parameters
    ----------

    multiatlas_Labels : string (nifti file)

    csf_label: integer
        label value corresponding to CSF in multiatlas file

    left_gm_label: integer
        label value corresponding to Left Gray Matter in multiatlas file

    left_wm_label: integer
        label value corresponding to Left White Matter in multiatlas file

    right_gm_label: integer
        label value corresponding to Right Gray Matter in multiatlas file

    right_wm_label: integer
        label value corresponding to Right White Matter in multiatlas file

    include_ventricles: boolean
        whether include labels of ventricles in CSF or not

    Returns
    -------

    csf_mask : string (nifti file)

    gm_mask : string (nifti file)

    wm_mask : string (nifti file)
    """  # noqa
    import os
    import nibabel as nb
    import numpy as np

    img = nb.load(multiatlas_Labels)
    data = img.get_data()

    # pick tissue mask from multiatlas labels file
    # based off of FreeSurferColorLUT or user provided label values
    # hard-coded csf/gm/wm label values are based off of FreeSurferColorLUT

    # FreeSurfer Ventricle Labels:
    # Left-Lateral-Ventricle 4, 3rd-Ventricle 14, 4th-Ventricle 15, Right-Lateral-Ventricle 43

    csf = data.copy()
    if include_ventricles:
        csf[np.logical_and(np.logical_and(np.logical_and(np.logical_and(csf != csf_label, csf != 4), csf != 14), csf != 15), csf != 43)] = 0
        csf[np.logical_or(np.logical_or(np.logical_or(np.logical_or(csf == csf_label, csf == 4), csf == 14), csf == 15), csf == 43)] = 1
    else:
        csf[csf != csf_label] = 0
        csf[csf == csf_label] = 1

    gm = data.copy()

    gm[np.logical_and(gm != right_gm_label, gm != left_gm_label)] = 0
    gm[np.logical_or(gm == right_gm_label, gm == left_gm_label)] = 1

    wm = data.copy()
    wm[np.logical_and(wm != right_wm_label, wm != left_wm_label)] = 0
    wm[np.logical_or(wm == right_wm_label, wm == left_wm_label)] = 1

    save_img_csf = nb.Nifti1Image(csf, header=img.get_header(), affine=img.get_affine())
    save_img_gm = nb.Nifti1Image(gm, header=img.get_header(), affine=img.get_affine())
    save_img_wm = nb.Nifti1Image(wm, header=img.get_header(), affine=img.get_affine())

    save_img_csf.to_filename('csf_mask.nii.gz')
    save_img_gm.to_filename('gm_mask.nii.gz')
    save_img_wm.to_filename('wm_mask.nii.gz')

    csf_mask = os.path.join(os.getcwd(), 'csf_mask.nii.gz')
    gm_mask = os.path.join(os.getcwd(), 'gm_mask.nii.gz')
    wm_mask = os.path.join(os.getcwd(), 'wm_mask.nii.gz')

    return csf_mask, gm_mask, wm_mask
