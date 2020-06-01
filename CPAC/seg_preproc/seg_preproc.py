import numpy as np

from nipype.interfaces.utility import Function
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype.interfaces.ants as ants
from nipype.interfaces.ants import WarpImageMultiTransform
from CPAC.seg_preproc.utils import (
    check_if_file_is_empty,
    pick_wm_prob_0,
    pick_wm_prob_1,
    pick_wm_prob_2,
    pick_wm_class_0,
    pick_wm_class_1,
    pick_wm_class_2,
    erosion,
    mask_erosion,
    hardcoded_antsJointLabelFusion,
    pick_tissue_from_labels_file)

import nipype.pipeline.engine as pe
import scipy.ndimage as nd
import numpy as np
from CPAC.registration.utils import check_transforms, generate_inverse_transform_flags


def create_seg_preproc(use_ants,
                        use_priors,
                        use_threshold,
                        csf_use_erosion=False,
                        wm_use_erosion=False,
                        gm_use_erosion=False,
                        wf_name='seg_preproc'):

    """Segment the subject's anatomical brain into cerebral spinal fluids,
    white matter and gray matter and refine them using template-space tissue
    priors, if selected to do so.

    Parameters
    ----------
    use_ants: boolean
        Whether we are using ANTs or FSL-FNIRT for registration purposes.
    wf_name : string
        name of the workflow

    Returns
    -------
    seg_preproc : workflow
        Workflow Object for Segmentation Workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/seg_preproc/seg_preproc.py>`_

    Workflow Inputs: ::

        csf_threshold.csf_threshold : list (float)
            Threshold of Cerebral Spinal Fluid probabilities

        wm_threshold.wm_threshold : list (float)
            Threshold of White Matter probabilities

        gm_threshold.gm_threshold : list (float)
            Threshold of Gray Matter probabilities

        inputspec.brain : string (existing nifti file)
            Anatomical image(without skull)

        inputspec.standard2highres_mat : string (existing affine transformation .mat file)
            File for transformation from mni space to anatomical space

        inputspec.PRIOR_CSF : string (existing nifti file)
            FSL Standard CSF Tissue prior image , binarized with threshold of 0.4

        inputspec.PRIOR_GRAY : string (existing nifti file)
            FSL Standard GRAY Matter Tissue prior image , binarized with threshold of 0.66

        inputspec.PRIOR_WHITE : string (existing nifti file)
            FSL Standard White Matter Tissue prior image , binarized with threshold of 0.2

    Workflow Outputs: ::

        outputspec.csf_mni2t1 : string (nifti file)
            outputs CSF prior template(in MNI space) registered to anatomical space

        outputspec.gm_mni2t1 : string (nifti file)
            outputs gray matter prior template registered to anatomical space

        outputspec.gm_mask : string (nifti file)
            outputs image after masking gm_combo with gm prior in t1 space

        outputspec.wm_mni2t1 : string (nifti file)
            outputs White Matter prior template(in MNI space) registered to anatomical space

        outputspec.wm_mask : string (nifti file)
            outputs image after masking wm_combo with white matter(wm) prior in t1 space

        outputspec.probability_maps : string (nifti file)
            outputs individual probability maps (output from brain segmentation using FAST)

        outputspec.mixeltype : string (nifti file)
            outputs mixeltype volume file _mixeltype (output from brain segmentation using FAST)

        outputspec.partial_volume_map : string (nifti file)
            outputs partial volume file _pveseg (output from brain segmentation using FAST)

        outputspec.partial_volume_files : string (nifti file)
            outputs partial volume estimate files _pve_ (output from brain segmentation using FAST)


    Order of commands:

    - Segment the Anatomical brain. For details see `fast <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FAST>`_::

        fast
        -t 1
        -g
        -p
        -o segment
        mprage_brain.nii.gz

    - Register CSF template in template space to t1 space. For details see `flirt <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT>`_::

        flirt
        -in PRIOR_CSF
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres_inv.mat
        -out csf_mni2t1

    - Find overlap between csf probability map and csf_mni2t1. For details see  `fslmaths <http://www.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro/index.htm>`_::

        fslmaths
        segment_prob_0.nii.gz
        -mas csf_mni2t1.nii.gz
        csf_combo.nii.gz

    - Threshold and binarize CSF probability map ::

        fslmaths
        csf_combo.nii.gz
        -thr 0.4 (threshold value can be changeable by user)
        -bin csf_bin.nii.gz

    - Generate CSF csf_mask, by applying csf prior in t1 space to binarized csf probability map ::

        fslmaths
        csf_bin.nii.gz
        -mas csf_mni2t1
        csf_mask

    - Register WM template in template space to t1 space ::

        flirt
        -in PRIOR_WM
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres.mat
        -out wm_mni2t1

    - Find overlap between wm probability mask and wm_mni2t1 ::

        fslmaths
        segment_prob_2.nii.gz
        -mas wm_mni2t1.nii.gz
        wm_combo.nii.gz

    - Threshold and binarize WM probability map ::

        fslmaths
        wm_combo.nii.gz
        -thr 0.4 (threshold value can be changeable by user)
        -bin wm_bin.nii.gz

    - Generate WM csf_mask, by applying wm_prior in t1 space to binarized wm map ::

        fslmaths
        wm_bin.nii.gz
        -mas wm_mni2t1
        wm_mask

    - Register GM template in template space to t1 space ::

        flirt
        -in PRIOR_GM
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres.mat
        -out gm_mni2t1

    - Find overlap between gm probability map and gm_mni2t1 ::

        fslmaths
        segment_prob_1.nii.gz
        -mas gm_mni2t1.nii.gz
        gm_combo.nii.gz

    - Threshold and binarize GM probability map ::

        fslmaths
        gm_combo.nii.gz
        -thr 0.4 (threshold value can be changeable by user)
        -bin gm_bin.nii.gz

    - Generate GM csf_mask, by applying gm prior in t1 space to thresholded binarized gm probability map ::

        fslmaths
        gm_bin.nii.gz
        -mas gm_mni2t1
        gm_mask


    Examples
    --------
    >>> import CPAC.seg_preproc as seg_wflow
    >>> seg = seg_wflow.create_seg_preproc()
    >>> seg.inputs.inputspec.standard2highres_mat = '/home/data/Projects/C-PAC/working_directory/s1001/reg_preproc/standard2highres.mat'
    >>> seg.inputs.inputspec.PRIOR_CSF = '/home/data/Projects/C-PAC/tissuepriors/2mm/avg152T1_csf_bin.nii.gz'
    >>> seg.inputs.inputspec.PRIOR_WHITE = '/home/data/Projects/C-PAC/tissuepriors/2mm/avg152T1_white_bin.nii.gz'
    >>> seg.inputs.inputspec.PRIOR_GRAY = '/home/data/Projects/C-PAC/tissuepriors/2mm/avg152T1_gray_bin.nii.gz'
    >>> seg.inputs.inputspec.brain = '/home/data/Projects/C-PAC/working_directory/s1001/anat_preproc/mprage_brain.nii.gz'
    >>> seg_preproc.run() # doctest: +SKIP


    High Level Graph:

    .. image:: ../images/seg_preproc.dot.png
        :width: 1100
        :height: 480

    Detailed Graph:

    .. image:: ../images/seg_preproc_detailed.dot.png
        :width: 1100
        :height: 480
    """

    preproc = pe.Workflow(name = wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['brain',
                                                       'brain_mask',
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig',
                                                       'PRIOR_CSF',
                                                       'PRIOR_GRAY',
                                                       'PRIOR_WHITE']),
                        name='inputspec')

    inputnode_csf_threshold = pe.Node(util.IdentityInterface(
                                    fields=['csf_threshold']),
                             name='csf_threshold')

    inputnode_wm_threshold = pe.Node(util.IdentityInterface(
                                    fields=['wm_threshold']),
                             name='wm_threshold')

    inputnode_gm_threshold = pe.Node(util.IdentityInterface(
                                    fields=['gm_threshold']),
                             name='gm_threshold')

    inputnode_csf_erosion_prop = pe.Node(util.IdentityInterface(
                                    fields=['csf_erosion_prop']),
                             name='csf_erosion_prop')

    inputnode_wm_erosion_prop = pe.Node(util.IdentityInterface(
                                    fields=['wm_erosion_prop']),
                             name='wm_erosion_prop')

    inputnode_gm_erosion_prop = pe.Node(util.IdentityInterface(
                                    fields=['gm_erosion_prop']),
                             name='gm_erosion_prop')

    inputnode_csf_mask_erosion_mm = pe.Node(util.IdentityInterface(
                                    fields=['csf_mask_erosion_mm']),
                             name='csf_mask_erosion_mm')

    inputnode_wm_mask_erosion_mm = pe.Node(util.IdentityInterface(
                                    fields=['wm_mask_erosion_mm']),
                             name='wm_mask_erosion_mm')

    inputnode_gm_mask_erosion_mm = pe.Node(util.IdentityInterface(
                                    fields=['gm_mask_erosion_mm']),
                             name='gm_mask_erosion_mm')

    inputnode_csf_erosion_mm = pe.Node(util.IdentityInterface(
                                    fields=['csf_erosion_mm']),
                             name='csf_erosion_mm')

    inputnode_wm_erosion_mm = pe.Node(util.IdentityInterface(
                                    fields=['wm_erosion_mm']),
                             name='wm_erosion_mm')

    inputnode_gm_erosion_mm = pe.Node(util.IdentityInterface(
                                    fields=['gm_erosion_mm']),
                             name='gm_erosion_mm')

    outputNode = pe.Node(util.IdentityInterface(fields=['csf_mask',
                                                        'gm_mask',
                                                        'wm_mask',
                                                        'csf_probability_map',
                                                        'probability_maps',
                                                        'tissue_class_files',
                                                        'mixeltype',
                                                        'partial_volume_map',
                                                        'partial_volume_files']),
                        name='outputspec')

    segment = pe.Node(interface=fsl.FAST(), name='segment', mem_gb=1.5)
    segment.inputs.img_type = 1
    segment.inputs.segments = True
    segment.inputs.probability_maps = True
    segment.inputs.out_basename = 'segment'

    check_wm = pe.Node(name='check_wm', interface=Function(function=check_if_file_is_empty, input_names=['in_file'], output_names=['out_file']))
    check_gm = pe.Node(name='check_gm', interface=Function(function=check_if_file_is_empty, input_names=['in_file'], output_names=['out_file']))
    check_csf = pe.Node(name='check_csf', interface=Function(function=check_if_file_is_empty, input_names=['in_file'], output_names=['out_file']))

    preproc.connect(inputNode, 'brain', segment, 'in_files')

    preproc.connect(segment, 'mixeltype',
                    outputNode, 'mixeltype')
    preproc.connect(segment, 'partial_volume_files',
                    outputNode, 'partial_volume_files')
    preproc.connect(segment, 'partial_volume_map',
                    outputNode, 'partial_volume_map')
    preproc.connect(segment, 'tissue_class_files',
                    outputNode, 'tissue_class_files')
    preproc.connect(segment, 'probability_maps',
                    outputNode, 'probability_maps')

    process_csf = process_segment_map('CSF', use_ants, use_priors, use_threshold, use_erosion=csf_use_erosion)


    if use_ants:
        preproc.connect(inputNode, 'standard2highres_init',
                        process_csf, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        process_csf, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    process_csf, 'inputspec.brain')
    preproc.connect(inputNode, 'brain_mask',
                    process_csf, 'inputspec.brain_mask')
    preproc.connect(inputnode_csf_threshold, 'csf_threshold',
                    process_csf, 'inputspec.threshold')
    preproc.connect(inputnode_csf_erosion_prop, 'csf_erosion_prop',
                    process_csf, 'inputspec.erosion_prop')
    preproc.connect(inputnode_csf_mask_erosion_mm, 'csf_mask_erosion_mm',
                    process_csf, 'inputspec.mask_erosion_mm')
    preproc.connect(inputnode_csf_erosion_mm, 'csf_erosion_mm',
                    process_csf, 'inputspec.erosion_mm')
    preproc.connect(inputNode, 'PRIOR_CSF',
                    process_csf, 'inputspec.tissue_prior')
    preproc.connect(segment, ('tissue_class_files', pick_wm_class_0),
                    process_csf, 'inputspec.tissue_class_file')
    preproc.connect(segment, ('probability_maps', pick_wm_prob_0),
                    process_csf, 'inputspec.probability_tissue_map')
    preproc.connect(inputNode, 'standard2highres_mat',
                    process_csf, 'inputspec.standard2highres_mat')
    preproc.connect(process_csf, 'outputspec.segment_mask',
                    outputNode, 'csf_mask')
    preproc.connect(process_csf, 'outputspec.probability_tissue_map',
                    outputNode, 'csf_probability_map')

    process_wm = process_segment_map('WM', use_ants, use_priors, use_threshold, use_erosion=wm_use_erosion)


    if use_ants:
        preproc.connect(inputNode, 'standard2highres_init',
                        process_wm, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        process_wm, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    process_wm, 'inputspec.brain')
    preproc.connect(inputNode, 'brain_mask',
                    process_wm, 'inputspec.brain_mask')
    preproc.connect(inputnode_wm_threshold, 'wm_threshold',
                    process_wm, 'inputspec.threshold')
    preproc.connect(inputnode_wm_erosion_prop, 'wm_erosion_prop',
                    process_wm, 'inputspec.erosion_prop')
    preproc.connect(inputnode_wm_mask_erosion_mm, 'wm_mask_erosion_mm',
                    process_wm, 'inputspec.mask_erosion_mm')
    preproc.connect(inputnode_wm_erosion_mm, 'wm_erosion_mm',
                    process_wm, 'inputspec.erosion_mm')
    preproc.connect(inputNode, 'PRIOR_WHITE',
                    process_wm, 'inputspec.tissue_prior')
    preproc.connect(segment, ('tissue_class_files', pick_wm_class_2),
                    process_wm, 'inputspec.tissue_class_file')
    preproc.connect(segment, ('probability_maps', pick_wm_prob_2),
                    process_wm, 'inputspec.probability_tissue_map')
    preproc.connect(inputNode, 'standard2highres_mat',
                    process_wm, 'inputspec.standard2highres_mat')
    preproc.connect(process_wm, 'outputspec.segment_mask',
                    outputNode, 'wm_mask')

    process_gm = process_segment_map('GM', use_ants, use_priors, use_threshold, use_erosion=gm_use_erosion)

    if use_ants:
        preproc.connect(inputNode, 'standard2highres_init',
                        process_gm, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        process_gm, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    process_gm, 'inputspec.brain')
    preproc.connect(inputNode, 'brain_mask',
                    process_gm, 'inputspec.brain_mask')
    preproc.connect(inputnode_gm_threshold, 'gm_threshold',
                    process_gm, 'inputspec.threshold')
    preproc.connect(inputnode_gm_erosion_prop, 'gm_erosion_prop',
                    process_gm, 'inputspec.erosion_prop')
    preproc.connect(inputnode_gm_mask_erosion_mm, 'gm_mask_erosion_mm',
                    process_gm, 'inputspec.mask_erosion_mm')
    preproc.connect(inputnode_gm_erosion_mm, 'gm_erosion_mm',
                    process_gm, 'inputspec.erosion_mm')
    preproc.connect(inputNode, 'PRIOR_GRAY',
                    process_gm, 'inputspec.tissue_prior')
    preproc.connect(segment, ('tissue_class_files', pick_wm_class_1),
                    process_gm, 'inputspec.tissue_class_file')
    preproc.connect(segment, ('probability_maps', pick_wm_prob_1),
                    process_gm, 'inputspec.probability_tissue_map')
    preproc.connect(inputNode, 'standard2highres_mat',
                    process_gm, 'inputspec.standard2highres_mat')
    preproc.connect(process_gm, 'outputspec.segment_mask',
                    outputNode, 'gm_mask')

    return preproc


def process_segment_map(wf_name,
                        use_ants,
                        use_priors,
                        use_threshold,
                        use_erosion):
    """This is a sub workflow used inside segmentation workflow to process
    probability maps obtained in segmentation. Steps include overlapping
    of the prior tissue with probability maps, thresholding and binarizing
    it and creating a mask that is used in further analysis.

    Parameters
    ----------
    wf_name : string
        Workflow Name
    use_priors: boolean
        Whether or not to use template-space tissue priors to further refine
        the resulting segmentation tissue masks.
    use_threshold: String
        Choose threshold to further refine
        the resulting segmentation tissue masks.
    use_erosion: boolean
        Whether or not to erode the resulting segmentation tissue masks.
    use_ants : boolean
        Whether or not to use ANTs or FSL for transform application.

    Returns
    -------
    preproc : workflow
        Workflow Object for process_segment_map Workflow


    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/seg_preproc/seg_preproc.py>`_


    Workflow Inputs::

        inputspec.brain : string (existing nifti file)
            Anatomical image(without skull)

        inputspec.standard2highres_mat : string (existing affine transformation .mat file)
            path to transformation matrix from mni space to anatomical space

        inputspec.threshold : float
            threshold value

        inputspec.tissue_prior : string (existing nifti file)
            path to FSL Standard Tissue prior image

        inputspec.probability_tissue_map : string (nifti file)
            tissue Probability map obtained from fsl FAST

    Workflow Outputs::

        outputspec.segment_mni2t1 : string (nifti file)
            path to output CSF prior template(in MNI space) registered to anatomical space

        outputspec.segment_combo : string (nifti file)
            path to output image containing overlap between csf probability map and segment_mni2t1

        outputspec.segment_thresh : string (nifti file)
            path to output image after Thresholding segment_combo

        outputspec.segment_bin : string (nifti file)
            path to output image after binarizing segment_thresh

        outputspec.segment_erosion : string (nifti file)
            path to output image after eroding segment_bin

        outputspec.segment_mask : string (nifti file)
            path to output image after masking segment_combo with its tissue prior in t1 space


    Order of commands:

    - Register tissue prior in MNI space to t1 space.

    - Threshold segment probability map

    - Binarize threshed segment probability map

    - Erose binarized segment mask

    - Generate segment mask, by applying tissue prior in t1 space to thresholded binarized segment probability map


    High Level Graph:

    .. image:: ../images/process_segment_map.dot.png
        :width: 1100
        :height: 480

    Detailed Graph:

    .. image:: ../images/process_segment_map_detailed.dot.png
        :width: 1100
        :height: 480

    """

    import nipype.interfaces.utility as util

    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['tissue_prior',
                                                       'threshold',
                                                       'erosion_prop',
                                                       'mask_erosion_mm',
                                                       'erosion_mm',
                                                       'brain',
                                                       'brain_mask',
                                                       'tissue_class_file',
                                                       'probability_tissue_map',
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['segment_mask',
                                                        'probability_tissue_map']),
                        name='outputspec')

    def form_threshold_string(threshold):
        return '-thr %f ' % (threshold)

    def form_mask_erosion_prop(erosion_prop):
        return erosion_prop**3

    if use_ants:
        collect_linear_transforms = pe.Node(util.Merge(3),
                                            name='{0}_collect_linear_transforms'.format(wf_name))
        preproc.connect(inputNode, 'standard2highres_init', collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig', collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat', collect_linear_transforms, 'in3')

        # check transform list to exclude Nonetype (missing) init/rig/affine
        check_transform = pe.Node(util.Function(input_names=['transform_list'], 
                                                output_names=['checked_transform_list', 'list_length'],
                                                function=check_transforms), name='{0}_check_transforms'.format(wf_name))
        
        preproc.connect(collect_linear_transforms, 'out', check_transform, 'transform_list')

        # generate inverse transform flags, which depends on the number of transforms
        inverse_transform_flags = pe.Node(util.Function(input_names=['transform_list'], 
                                                        output_names=['inverse_transform_flags'],
                                                        function=generate_inverse_transform_flags), 
                                                        name='{0}_inverse_transform_flags'.format(wf_name))

        preproc.connect(check_transform, 'checked_transform_list', inverse_transform_flags, 'transform_list')

        # mni to t1
        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='{0}_prior_mni_to_t1'.format(wf_name))

        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'

        preproc.connect(inverse_transform_flags, 'inverse_transform_flags', tissueprior_mni_to_t1, 'invert_transform_flags')
        preproc.connect(inputNode, 'tissue_prior', tissueprior_mni_to_t1, 'input_image')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference_image')
        preproc.connect(check_transform, 'checked_transform_list', tissueprior_mni_to_t1, 'transforms')

        if 'FSL-FAST Thresholding' in use_threshold:
            input_1, value_1 = (inputNode, 'tissue_class_file')
        else:
            input_1, value_1 = (inputNode, 'probability_tissue_map')

        if use_priors:
            overlap_segmentmap_with_prior = pe.Node(interface=fsl.MultiImageMaths(),
                                                    name='overlap_%s_map_with_prior' % (wf_name))
            overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '

            preproc.connect(input_1, value_1, overlap_segmentmap_with_prior, 'in_file')

            preproc.connect(tissueprior_mni_to_t1, 'output_image', overlap_segmentmap_with_prior, 'operand_files')

            input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')


        if 'Customized Thresholding' in use_threshold:
            segmentmap_threshold = pe.Node(interface=fsl.ImageMaths(),
                                                name='threshold_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, ('threshold', form_threshold_string), segmentmap_threshold, 'op_string')

            preproc.connect(input_1, value_1, segmentmap_threshold, 'in_file')

            input_1, value_1 = (segmentmap_threshold, 'out_file')


        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(),
                                                name='binarize_%s' % (wf_name))
        binarize_threshold_segmentmap.inputs.op_string = '-bin '

        preproc.connect(input_1, value_1, binarize_threshold_segmentmap, 'in_file')

        input_1, value_1 = (binarize_threshold_segmentmap, 'out_file')

        preproc.connect(input_1, value_1, outputNode, 'probability_tissue_map')

        ero_imports = ['import scipy.ndimage as nd' , 'import numpy as np', 'import nibabel as nb', 'import os']

        if use_erosion:
            # mask erosion
            eroded_mask = pe.Node(util.Function(input_names = ['roi_mask', 'skullstrip_mask', 'mask_erosion_mm', 'mask_erosion_prop'],
                                                output_names = ['output_roi_mask', 'eroded_skullstrip_mask'],
                                                function = mask_erosion,
                                                imports = ero_imports),
                                                name='erode_skullstrip_mask_%s' % (wf_name))
            preproc.connect(inputNode, ('erosion_prop', form_mask_erosion_prop), eroded_mask, 'mask_erosion_prop')
            preproc.connect(inputNode, 'mask_erosion_mm', eroded_mask, 'mask_erosion_mm')
            preproc.connect(inputNode, 'brain_mask', eroded_mask, 'skullstrip_mask')
            preproc.connect(input_1, value_1, eroded_mask, 'roi_mask')

            input_1, value_1 = (eroded_mask, 'output_roi_mask')

            # erosion
            erosion_segmentmap = pe.Node(util.Function(input_names = ['roi_mask', 'erosion_mm', 'erosion_prop'],
                                                output_names = ['eroded_roi_mask'],
                                                function = erosion,
                                                imports = ero_imports),
                                                name='erosion_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, 'erosion_prop', erosion_segmentmap, 'erosion_prop')
            preproc.connect(inputNode, 'erosion_mm', erosion_segmentmap, 'erosion_mm')
            preproc.connect(input_1, value_1, erosion_segmentmap, 'roi_mask')
            input_1, value_1 = (erosion_segmentmap, 'eroded_roi_mask')

        preproc.connect (input_1, value_1, outputNode, 'segment_mask')


    else:
        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='{0}_prior_mni_to_t1'.format(wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'

        # mni to t1
        preproc.connect(inputNode, 'tissue_prior', tissueprior_mni_to_t1, 'in_file')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference')

        preproc.connect(inputNode, 'standard2highres_mat', tissueprior_mni_to_t1, 'in_matrix_file')

        if 'FSL-FAST Thresholding' in use_threshold:
            input_1, value_1 = (inputNode, 'tissue_class_file')
        else:
            input_1, value_1 = (inputNode, 'probability_tissue_map')

        if use_priors:
            overlap_segmentmap_with_prior = pe.Node(interface=fsl.MultiImageMaths(),
                                                    name='overlap_%s_map_with_prior' % (wf_name))
            overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '

            preproc.connect(input_1, value_1, overlap_segmentmap_with_prior, 'in_file')

            preproc.connect(tissueprior_mni_to_t1, 'out_file', overlap_segmentmap_with_prior, 'operand_files')

            input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')


        if 'Customized Thresholding' in use_threshold:
            segmentmap_threshold = pe.Node(interface=fsl.ImageMaths(),
                                                name='threshold_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, ('threshold', form_threshold_string), segmentmap_threshold, 'op_string')

            preproc.connect(input_1, value_1, segmentmap_threshold, 'in_file')

            input_1, value_1 = (segmentmap_threshold, 'out_file')


        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(),
                                                name='binarize_%s' % (wf_name))
        binarize_threshold_segmentmap.inputs.op_string = '-bin '

        preproc.connect(input_1, value_1, binarize_threshold_segmentmap, 'in_file')

        input_1, value_1 = (binarize_threshold_segmentmap, 'out_file')

        ero_imports = ['import scipy.ndimage as nd' , 'import numpy as np', 'import nibabel as nb', 'import os']

        if use_erosion:
            # mask erosion
            eroded_mask = pe.Node(util.Function(input_names = ['roi_mask', 'skullstrip_mask', 'mask_erosion_mm', 'mask_erosion_prop'],
                                                output_names = ['output_roi_mask', 'eroded_skullstrip_mask'],
                                                function = mask_erosion,
                                                imports = ero_imports),
                                                name='erode_skullstrip_mask_%s' % (wf_name))
            preproc.connect(inputNode, ('erosion_prop', form_mask_erosion_prop), eroded_mask, 'mask_erosion_prop')
            preproc.connect(inputNode, 'mask_erosion_mm', eroded_mask, 'mask_erosion_mm')
            preproc.connect(inputNode, 'brain_mask', eroded_mask, 'skullstrip_mask')
            preproc.connect(input_1, value_1, eroded_mask, 'roi_mask')

            input_1, value_1 = (eroded_mask, 'output_roi_mask')

            # erosion
            erosion_segmentmap = pe.Node(util.Function(input_names = ['roi_mask', 'erosion_mm', 'erosion_prop'],
                                                output_names = ['eroded_roi_mask'],
                                                function = erosion,
                                                imports = ero_imports),
                                                name='erosion_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, 'erosion_prop', erosion_segmentmap, 'erosion_prop')
            preproc.connect(inputNode, 'erosion_mm', erosion_segmentmap, 'erosion_mm')
            preproc.connect(input_1, value_1, erosion_segmentmap, 'roi_mask')
            input_1, value_1 = (erosion_segmentmap, 'eroded_roi_mask')

        preproc.connect (input_1, value_1, outputNode, 'segment_mask')

    return preproc

def create_seg_preproc_template_based(use_ants,
                                    wf_name='seg_preproc_templated_based'):

    """Generate the subject's cerebral spinal fluids,
    white matter and gray matter mask based on provided template, if selected to do so.

    Parameters
    ----------
    use_ants: boolean
        Whether we are using ANTs or FSL-FNIRT for registration purposes.
    wf_name : string
        name of the workflow

    Returns
    -------
    seg_preproc_templated_based : workflow
        Workflow Object for Segmentation Workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/seg_preproc/seg_preproc.py>`_

    Workflow Inputs: ::

        inputspec.brain : string (existing nifti file)
            Anatomical image(without skull)
                Note: Mean EPI will replace anatomical image, if anatomical data doesn't exist.

        inputspec.standard2highres_mat : string (existing affine transformation .mat file)
            File for transformation from mni space to anatomical space

        inputspec.CSF_template : string (existing nifti file)
            CSF tissue mask on template space

        inputspec.GRAY_template : string (existing nifti file)
            GRAY Matter CSF tissue mask on template space

        inputspec.WHITE_template : string (existing nifti file)
            White Matter tissue mask on template space

    Workflow Outputs: ::

        outputspec.csf_mni2t1 : string (nifti file)
            outputs CSF prior template(in MNI space) registered to anatomical space

        outputspec.gm_mni2t1 : string (nifti file)
            outputs gray matter prior template registered to anatomical space

        outputspec.wm_mni2t1 : string (nifti file)
            outputs White Matter prior template(in MNI space) registered to anatomical space


    Order of commands:

    - Register CSF template in template space to t1(or mean EPI) space. For details see `flirt <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT>`_::

        flirt
        -in CSF_template
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres_inv.mat
        -out csf_mni2t1

        -bin csf_bin.nii.gz


    - Register WM template in template space to t1 space ::

        flirt
        -in WHITE_template
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres.mat
        -out wm_mni2t1


    - Register GM template in template space to t1 space ::

        flirt
        -in GRAY_template
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres.mat
        -out gm_mni2t1

    """

    preproc = pe.Workflow(name = wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['brain',
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig',
                                                       'CSF_template',
                                                       'WHITE_template',
                                                       'GRAY_template']),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=['csf_mask',
                                                        'gm_mask',
                                                        'wm_mask']),
                        name='outputspec')


    csf_template2t1 = tissue_mask_template_to_t1('CSF', use_ants)

    if use_ants:
        preproc.connect(inputNode, 'standard2highres_init',
                        csf_template2t1, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        csf_template2t1, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    csf_template2t1, 'inputspec.brain')
    preproc.connect(inputNode, 'CSF_template',
                    csf_template2t1, 'inputspec.tissue_mask_template')
    preproc.connect(inputNode, 'standard2highres_mat',
                    csf_template2t1, 'inputspec.standard2highres_mat')
    preproc.connect(csf_template2t1, 'outputspec.segment_mask_temp2t1',
                    outputNode, 'csf_mask')


    wm_template2t1 = tissue_mask_template_to_t1('WM', use_ants)

    if use_ants:
        preproc.connect(inputNode, 'standard2highres_init',
                        wm_template2t1, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        wm_template2t1, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    wm_template2t1, 'inputspec.brain')
    preproc.connect(inputNode, 'WHITE_template',
                    wm_template2t1, 'inputspec.tissue_mask_template')
    preproc.connect(inputNode, 'standard2highres_mat',
                    wm_template2t1, 'inputspec.standard2highres_mat')
    preproc.connect(wm_template2t1, 'outputspec.segment_mask_temp2t1',
                    outputNode, 'wm_mask')


    gm_template2t1 = tissue_mask_template_to_t1('GM', use_ants)

    if use_ants:
        preproc.connect(inputNode, 'standard2highres_init',
                        gm_template2t1, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        gm_template2t1, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    gm_template2t1, 'inputspec.brain')
    preproc.connect(inputNode, 'GRAY_template',
                    gm_template2t1, 'inputspec.tissue_mask_template')
    preproc.connect(inputNode, 'standard2highres_mat',
                    gm_template2t1, 'inputspec.standard2highres_mat')
    preproc.connect(gm_template2t1, 'outputspec.segment_mask_temp2t1',
                    outputNode, 'gm_mask')

    return preproc


def tissue_mask_template_to_t1(wf_name,
                                use_ants):

    import nipype.interfaces.utility as util

    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['brain',
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig',
                                                       'tissue_mask_template']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['segment_mask_temp2t1']),
                        name='outputspec')

    if use_ants:
        collect_linear_transforms = pe.Node(util.Merge(3),
                                            name='{0}_collect_linear_transforms'.format(wf_name))

        preproc.connect(inputNode, 'standard2highres_init', collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig', collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat', collect_linear_transforms, 'in3')

        # check transform list to exclude Nonetype (missing) init/rig/affine
        check_transform = pe.Node(util.Function(input_names=['transform_list'], 
                                                output_names=['checked_transform_list', 'list_length'],
                                                function=check_transforms), name='{0}_check_transforms'.format(wf_name))
        
        preproc.connect(collect_linear_transforms, 'out', check_transform, 'transform_list')

        # generate inverse transform flags, which depends on the number of transforms
        inverse_transform_flags = pe.Node(util.Function(input_names=['transform_list'], 
                                                        output_names=['inverse_transform_flags'],
                                                        function=generate_inverse_transform_flags), 
                                                        name='{0}_inverse_transform_flags'.format(wf_name))

        preproc.connect(check_transform, 'checked_transform_list', inverse_transform_flags, 'transform_list')

        # mni to t1
        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='{0}_mni_to_t1'.format(wf_name))

        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'

        preproc.connect(inverse_transform_flags, 'inverse_transform_flags', tissueprior_mni_to_t1, 'invert_transform_flags')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference_image')
        preproc.connect(check_transform, 'checked_transform_list', tissueprior_mni_to_t1, 'transforms')
        preproc.connect(inputNode, 'tissue_mask_template', tissueprior_mni_to_t1, 'input_image')

        preproc.connect (tissueprior_mni_to_t1, 'output_image', outputNode, 'segment_mask_temp2t1')

    else:
        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='{0}_mni_to_t1'.format(wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'

        # mni to t1
        preproc.connect(inputNode, 'tissue_mask_template', tissueprior_mni_to_t1, 'in_file')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference')
        preproc.connect(inputNode, 'standard2highres_mat', tissueprior_mni_to_t1, 'in_matrix_file')

        preproc.connect (tissueprior_mni_to_t1, 'out_file', outputNode, 'segment_mask_temp2t1')

    return preproc


def create_seg_preproc_antsJointLabel_method(wf_name='seg_preproc_templated_based'):

    """Generate the subject's cerebral spinal fluids,
    white matter and gray matter mask based on provided template, if selected to do so.

    Parameters
    ----------
    wf_name : string
        name of the workflow

    Returns
    -------
    seg_preproc_templated_based : workflow
        Workflow Object for Segmentation Workflow

    Notes
    -----

    Workflow Inputs: ::

        inputspec.brain : string (existing nifti file)
            Anatomical image(without skull)
               
        inputspec.template_brain : string (existing nifti file)
            Template anatomical image(without skull)

        inputspec.template_segmentation : string (existing nifti file)
            Template segmentation image(without skull)

    Workflow Outputs: ::

        outputspec.csf_mask : string (nifti file)
            outputs CSF mask

        outputspec.gm_mask : string (nifti file)
            outputs gray matter mask

        outputspec.wm_mask : string (nifti file)
            outputs White Matter mask


    """

    preproc = pe.Workflow(name = wf_name)
    inputNode = pe.Node(util.IdentityInterface(fields=['anatomical_brain',
                                                       'anatomical_brain_mask',
                                                       'template_brain_list',
                                                       'template_segmentation_list',
                                                       'csf_label', 
                                                       'left_gm_label', 
                                                       'left_wm_label', 
                                                       'right_gm_label', 
                                                       'right_wm_label']),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=['csf_mask',
                                                        'gm_mask',
                                                        'wm_mask']),
                        name='outputspec')


    seg_preproc_antsJointLabel = pe.Node(util.Function(input_names=['anatomical_brain', 'anatomical_brain_mask', 'template_brain_list', 'template_segmentation_list'], 
                                                        output_names=['multiatlas_Intensity', 'multiatlas_Labels'],
                                                        function=hardcoded_antsJointLabelFusion), 
                                                        name='{0}_antsJointLabel'.format(wf_name))

    preproc.connect(inputNode, 'anatomical_brain',
                    seg_preproc_antsJointLabel, 'anatomical_brain')
    preproc.connect(inputNode, 'anatomical_brain_mask',
                    seg_preproc_antsJointLabel, 'anatomical_brain_mask')
    preproc.connect(inputNode, 'template_brain_list',
                    seg_preproc_antsJointLabel, 'template_brain_list')
    preproc.connect(inputNode, 'template_segmentation_list',
                    seg_preproc_antsJointLabel, 'template_segmentation_list')                 
    
    pick_tissue = pe.Node(util.Function(input_names=['multiatlas_Labels', 'csf_label', 'left_gm_label', 'left_wm_label', 'right_gm_label', 'right_wm_label'], 
                                        output_names=['csf_mask', 'gm_mask', 'wm_mask'],
                                        function=pick_tissue_from_labels_file), 
                                        name='{0}_tissue_mask'.format(wf_name))
    
    preproc.connect(seg_preproc_antsJointLabel, 'multiatlas_Labels',
                    pick_tissue, 'multiatlas_Labels') 
    preproc.connect(inputNode, 'csf_label',
                    pick_tissue, 'csf_label') 
    preproc.connect(inputNode, 'left_gm_label',
                    pick_tissue, 'left_gm_label') 
    preproc.connect(inputNode, 'left_wm_label',
                    pick_tissue, 'left_wm_label')                    
    preproc.connect(inputNode, 'right_gm_label',
                    pick_tissue, 'right_gm_label') 
    preproc.connect(inputNode, 'right_wm_label',
                    pick_tissue, 'right_wm_label') 

    preproc.connect(pick_tissue, 'csf_mask',
                    outputNode, 'csf_mask')
    preproc.connect(pick_tissue, 'gm_mask',
                    outputNode, 'gm_mask')   
    preproc.connect(pick_tissue, 'wm_mask',
                    outputNode, 'wm_mask')

    return preproc
