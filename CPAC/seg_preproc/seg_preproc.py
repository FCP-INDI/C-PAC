from nipype.interfaces.utility import Function
import nipype.algorithms.rapidart as ra
from nipype.interfaces import afni, ants, freesurfer, fsl, utility as util
from nipype.interfaces.ants import WarpImageMultiTransform
from CPAC.anat_preproc.utils import mri_convert
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
from nipype.interfaces import freesurfer
from CPAC.anat_preproc.utils import mri_convert
from CPAC.registration.utils import (
    check_transforms,
    generate_inverse_transform_flags)
from CPAC.pipeline.schema import valid_options


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
    use_priors: boolean

    use_threshold: list

    csf_use_erosion: boolean

    wm_use_erosion: boolean

    gm_use_erosion: boolean

    wf_name: string
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

    .. exec::
        from CPAC.seg_preproc import create_seg_preproc
        wf = create_seg_preproc(False, False, ['FSL-FAST Thresholding'])
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/seg_preproc.dot'
        )

    High Level Graph:

    .. image:: ../../images/generated/seg_preproc.png
        :width: 1100
        :height: 100

    Detailed Graph:

    .. image:: ../../images/generated/seg_preproc_detailed.png
        :width: 1100
        :height: 480
    """  # noqa
    preproc = pe.Workflow(name=wf_name)
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

    outputNode = pe.Node(
        util.IdentityInterface(fields=['csf_mask',
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

    check_wm = pe.Node(name='check_wm', interface=Function(
        function=check_if_file_is_empty, input_names=['in_file'],
        output_names=['out_file']))
    check_gm = pe.Node(name='check_gm', interface=Function(
        function=check_if_file_is_empty, input_names=['in_file'],
        output_names=['out_file']))
    check_csf = pe.Node(name='check_csf', interface=Function(
        function=check_if_file_is_empty, input_names=['in_file'],
        output_names=['out_file']))

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

    process_csf = process_segment_map('CSF', use_ants, use_priors,
                                      use_threshold,
                                      use_erosion=csf_use_erosion)

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

    process_wm = process_segment_map('WM', use_ants, use_priors,
                                     use_threshold, use_erosion=wm_use_erosion)

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

    process_gm = process_segment_map('GM', use_ants, use_priors,
                                     use_threshold, use_erosion=gm_use_erosion)

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
    use_threshold: list
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

    .. exec::
        from CPAC.seg_preproc import process_segment_map
        wf = process_segment_map('segment_map_wf',
                                False,
                                False,
                                ['FSL-FAST Thresholding'],
                                False)
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/process_segment_map.dot'
        )

    High Level Graph:

    .. image:: ../../images/generated/process_segment_map.png
        :width: 1100
        :height: 480

    Detailed Graph:

    .. image:: ../../images/generated/process_segment_map_detailed.png
        :width: 1100
        :height: 480

    """  # noqa
    import nipype.interfaces.utility as util

    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['tissue_prior',
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

    outputNode = pe.Node(
        util.IdentityInterface(fields=['segment_mask',
                                       'probability_tissue_map']),
        name='outputspec')

    def form_threshold_string(threshold):
        return '-thr %f ' % (threshold)

    def form_mask_erosion_prop(erosion_prop):
        return erosion_prop**3

    if use_ants:
        collect_linear_transforms = pe.Node(
            util.Merge(3), name='{0}_collect_linear_transforms'.format(
                wf_name))
        preproc.connect(inputNode, 'standard2highres_init',
                        collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig',
                        collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat',
                        collect_linear_transforms, 'in3')

        # check transform list to exclude Nonetype (missing) init/rig/affine
        check_transform = pe.Node(
            util.Function(input_names=['transform_list'],
                          output_names=['checked_transform_list',
                                        'list_length'],
                          function=check_transforms),
            name='{0}_check_transforms'.format(wf_name))

        preproc.connect(collect_linear_transforms, 'out',
                        check_transform, 'transform_list')

        # generate inverse transform flags, which depends on the
        # number of transforms
        inverse_transform_flags = pe.Node(
            util.Function(input_names=['transform_list'],
                          output_names=['inverse_transform_flags'],
                          function=generate_inverse_transform_flags),
            name='{0}_inverse_transform_flags'.format(wf_name))

        preproc.connect(check_transform, 'checked_transform_list',
                        inverse_transform_flags, 'transform_list')

        # mni to t1
        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='{0}_prior_mni_to_t1'.format(
                                            wf_name))

        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'

        preproc.connect(inverse_transform_flags, 'inverse_transform_flags',
                        tissueprior_mni_to_t1, 'invert_transform_flags')
        preproc.connect(inputNode, 'tissue_prior',
                        tissueprior_mni_to_t1, 'input_image')
        preproc.connect(inputNode, 'brain',
                        tissueprior_mni_to_t1, 'reference_image')
        preproc.connect(check_transform, 'checked_transform_list',
                        tissueprior_mni_to_t1, 'transforms')

        if 'FSL-FAST' in use_threshold:
            input_1, value_1 = (inputNode, 'tissue_class_file')
        else:
            input_1, value_1 = (inputNode, 'probability_tissue_map')

        if use_priors:
            overlap_segmentmap_with_prior = pe.Node(
                interface=fsl.MultiImageMaths(),
                name='overlap_%s_map_with_prior' % (wf_name))
            overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '

            preproc.connect(input_1, value_1,
                            overlap_segmentmap_with_prior, 'in_file')

            preproc.connect(tissueprior_mni_to_t1, 'output_image',
                            overlap_segmentmap_with_prior, 'operand_files')

            input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')

        if 'Custom' in use_threshold:
            segmentmap_threshold = pe.Node(
                interface=fsl.ImageMaths(),
                name='threshold_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, ('threshold', form_threshold_string),
                            segmentmap_threshold, 'op_string')

            preproc.connect(input_1, value_1, segmentmap_threshold, 'in_file')

            input_1, value_1 = (segmentmap_threshold, 'out_file')

        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(),
                                                name='binarize_%s' % (wf_name))
        binarize_threshold_segmentmap.inputs.op_string = '-bin '

        preproc.connect(input_1, value_1,
                        binarize_threshold_segmentmap, 'in_file')

        input_1, value_1 = (binarize_threshold_segmentmap, 'out_file')

        preproc.connect(input_1, value_1, outputNode, 'probability_tissue_map')

        ero_imports = ['import scipy.ndimage as nd', 'import numpy as np',
                       'import nibabel as nb', 'import os']

        if use_erosion:
            # mask erosion
            eroded_mask = pe.Node(
                util.Function(input_names=['roi_mask', 'skullstrip_mask',
                                           'mask_erosion_mm',
                                           'mask_erosion_prop'],
                              output_names=['output_roi_mask',
                                            'eroded_skullstrip_mask'],
                              function=mask_erosion,
                              imports=ero_imports),
                name='erode_skullstrip_mask_%s' % (wf_name))
            preproc.connect(
                inputNode, ('erosion_prop', form_mask_erosion_prop),
                eroded_mask, 'mask_erosion_prop')
            preproc.connect(inputNode, 'mask_erosion_mm',
                            eroded_mask, 'mask_erosion_mm')
            preproc.connect(inputNode, 'brain_mask',
                            eroded_mask, 'skullstrip_mask')
            preproc.connect(input_1, value_1, eroded_mask, 'roi_mask')

            input_1, value_1 = (eroded_mask, 'output_roi_mask')

            # erosion
            erosion_segmentmap = pe.Node(
                util.Function(input_names=['roi_mask', 'erosion_mm',
                                           'erosion_prop'],
                              output_names=['eroded_roi_mask'],
                              function=erosion,
                              imports=ero_imports),
                name='erosion_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, 'erosion_prop',
                            erosion_segmentmap, 'erosion_prop')
            preproc.connect(inputNode, 'erosion_mm',
                            erosion_segmentmap, 'erosion_mm')
            preproc.connect(input_1, value_1, erosion_segmentmap, 'roi_mask')
            input_1, value_1 = (erosion_segmentmap, 'eroded_roi_mask')

        preproc.connect(input_1, value_1, outputNode, 'segment_mask')

    else:
        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='{0}_prior_mni_to_t1'.format(
                                            wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'

        # mni to t1
        preproc.connect(inputNode, 'tissue_prior',
                        tissueprior_mni_to_t1, 'in_file')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference')

        preproc.connect(inputNode, 'standard2highres_mat',
                        tissueprior_mni_to_t1, 'in_matrix_file')

        if 'FSL-FAST' in use_threshold:
            input_1, value_1 = (inputNode, 'tissue_class_file')
        else:
            input_1, value_1 = (inputNode, 'probability_tissue_map')

        if use_priors:
            overlap_segmentmap_with_prior = pe.Node(
                interface=fsl.MultiImageMaths(),
                name='overlap_%s_map_with_prior' % (wf_name))
            overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '

            preproc.connect(input_1, value_1,
                            overlap_segmentmap_with_prior, 'in_file')

            preproc.connect(tissueprior_mni_to_t1, 'out_file',
                            overlap_segmentmap_with_prior, 'operand_files')

            input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')

        if 'Custom' in use_threshold:
            segmentmap_threshold = pe.Node(interface=fsl.ImageMaths(),
                                           name='threshold_segmentmap_%s' % (
                                               wf_name))
            preproc.connect(inputNode, ('threshold', form_threshold_string),
                            segmentmap_threshold, 'op_string')

            preproc.connect(input_1, value_1, segmentmap_threshold, 'in_file')

            input_1, value_1 = (segmentmap_threshold, 'out_file')

        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(),
                                                name='binarize_%s' % (wf_name))
        binarize_threshold_segmentmap.inputs.op_string = '-bin '

        preproc.connect(input_1, value_1,
                        binarize_threshold_segmentmap, 'in_file')

        input_1, value_1 = (binarize_threshold_segmentmap, 'out_file')

        ero_imports = ['import scipy.ndimage as nd', 'import numpy as np',
                       'import nibabel as nb', 'import os']

        if use_erosion:
            # mask erosion
            eroded_mask = pe.Node(
                util.Function(input_names=['roi_mask', 'skullstrip_mask',
                                           'mask_erosion_mm',
                                           'mask_erosion_prop'],
                              output_names=['output_roi_mask',
                                            'eroded_skullstrip_mask'],
                              function=mask_erosion,
                              imports=ero_imports),
                name='erode_skullstrip_mask_%s' % (wf_name))
            preproc.connect(
                inputNode, ('erosion_prop', form_mask_erosion_prop),
                eroded_mask, 'mask_erosion_prop')
            preproc.connect(inputNode, 'mask_erosion_mm',
                            eroded_mask, 'mask_erosion_mm')
            preproc.connect(inputNode, 'brain_mask',
                            eroded_mask, 'skullstrip_mask')
            preproc.connect(input_1, value_1, eroded_mask, 'roi_mask')

            input_1, value_1 = (eroded_mask, 'output_roi_mask')

            # erosion
            erosion_segmentmap = pe.Node(
                util.Function(input_names=['roi_mask', 'erosion_mm',
                                           'erosion_prop'],
                              output_names=['eroded_roi_mask'],
                              function=erosion,
                              imports=ero_imports),
                name='erosion_segmentmap_%s' % (wf_name))
            preproc.connect(inputNode, 'erosion_prop',
                            erosion_segmentmap, 'erosion_prop')
            preproc.connect(inputNode, 'erosion_mm',
                            erosion_segmentmap, 'erosion_mm')
            preproc.connect(input_1, value_1, erosion_segmentmap, 'roi_mask')
            input_1, value_1 = (erosion_segmentmap, 'eroded_roi_mask')

        preproc.connect(input_1, value_1, outputNode, 'segment_mask')

    return preproc


def create_seg_preproc_template_based(
    use_ants, wf_name='seg_preproc_templated_based'
):
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

    """  # noqa
    preproc = pe.Workflow(name=wf_name)
    inputNode = pe.Node(
        util.IdentityInterface(fields=['brain',
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


def tissue_mask_template_to_t1(wf_name, use_ants):

    import nipype.interfaces.utility as util

    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(
        util.IdentityInterface(fields=['brain',
                                       'standard2highres_init',
                                       'standard2highres_mat',
                                       'standard2highres_rig',
                                       'tissue_mask_template']),
        name='inputspec')

    outputNode = pe.Node(
        util.IdentityInterface(fields=['segment_mask_temp2t1']),
        name='outputspec')

    if use_ants:
        collect_linear_transforms = pe.Node(
            util.Merge(3),
            name='{0}_collect_linear_transforms'.format(wf_name))

        preproc.connect(inputNode, 'standard2highres_init',
                        collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig',
                        collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat',
                        collect_linear_transforms, 'in3')

        # check transform list to exclude Nonetype (missing) init/rig/affine
        check_transform = pe.Node(
            util.Function(input_names=['transform_list'],
                          output_names=['checked_transform_list',
                                        'list_length'],
                          function=check_transforms),
            name='{0}_check_transforms'.format(wf_name))

        preproc.connect(collect_linear_transforms, 'out',
                        check_transform, 'transform_list')

        # generate inverse transform flags, which depends on the
        # number of transforms
        inverse_transform_flags = pe.Node(
            util.Function(input_names=['transform_list'],
                          output_names=['inverse_transform_flags'],
                          function=generate_inverse_transform_flags),
            name='{0}_inverse_transform_flags'.format(wf_name))

        preproc.connect(check_transform, 'checked_transform_list',
                        inverse_transform_flags, 'transform_list')

        # mni to t1
        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='{0}_mni_to_t1'.format(wf_name))

        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'

        preproc.connect(inverse_transform_flags, 'inverse_transform_flags',
                        tissueprior_mni_to_t1, 'invert_transform_flags')
        preproc.connect(inputNode, 'brain',
                        tissueprior_mni_to_t1, 'reference_image')
        preproc.connect(check_transform, 'checked_transform_list',
                        tissueprior_mni_to_t1, 'transforms')
        preproc.connect(inputNode, 'tissue_mask_template',
                        tissueprior_mni_to_t1, 'input_image')

        preproc.connect(tissueprior_mni_to_t1, 'output_image',
                        outputNode, 'segment_mask_temp2t1')

    else:
        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='{0}_mni_to_t1'.format(wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'

        # mni to t1
        preproc.connect(inputNode, 'tissue_mask_template',
                        tissueprior_mni_to_t1, 'in_file')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference')
        preproc.connect(inputNode, 'standard2highres_mat',
                        tissueprior_mni_to_t1, 'in_matrix_file')

        preproc.connect(tissueprior_mni_to_t1, 'out_file',
                        outputNode, 'segment_mask_temp2t1')

    return preproc


def create_seg_preproc_antsJointLabel_method(wf_name='seg_preproc_templated_based'):
    """
    Generate the subject's cerebral spinal fluids,
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
                                                       'template_segmentation'
                                                       '_list',
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

    seg_preproc_antsJointLabel = pe.Node(
        util.Function(input_names=['anatomical_brain',
                                   'anatomical_brain_mask',
                                   'template_brain_list',
                                   'template_segmentation_list'],
                      output_names=['multiatlas_Intensity',
                                    'multiatlas_Labels'],
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

    pick_tissue = pe.Node(util.Function(input_names=['multiatlas_Labels',
                                                     'csf_label',
                                                     'left_gm_label',
                                                     'left_wm_label',
                                                     'right_gm_label',
                                                     'right_wm_label'],
                                        output_names=['csf_mask', 'gm_mask',
                                                      'wm_mask'],
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


def create_seg_preproc_freesurfer(config=None,
                                  wf_name='seg_preproc_freesurfer'):
    """
    Generate the subject's segmentations based on freesurfer.

    Parameters
    ----------
    wf_name : string
        name of the workflow

    Returns
    -------
    seg_preproc_freesurfer : workflow
        workflow object for segmentation workflow

    Notes
    -----

    Workflow Inputs: ::

        inputspec.subject_dir : string (existing nifti file)
            FreeSurfer autorecon1 dir

    Workflow Outputs: ::

        outputspec.wm_mask : string (nifti file)
            outputs White Matter mask
    """  # noqa
    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['subject_dir']),
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['wm_mask',
                                                        'gm_mask',
                                                        'csf_mask',
                                                        'subject_id']),
                         name='outputspec')

    reconall2 = pe.Node(interface=freesurfer.ReconAll(),
                        name='anat_autorecon2')

    reconall2.inputs.directive = 'autorecon2'
    reconall2.inputs.openmp = config.num_omp_threads  # TODO: update nested

    if config.autorecon2_args is not None:  # TODO: update nested
        reconall2.inputs.args = config.autorecon2_args  # TODO: update nested

    preproc.connect(inputnode, 'subject_dir',
                    reconall2, 'subjects_dir')

    preproc.connect(reconall2, 'subject_id',
                    outputnode, 'subject_id')

    # register FS segmentations (aseg.mgz) to native space
    fs_aseg_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                                name='fs_aseg_to_native')

    fs_aseg_to_native.inputs.reg_header = True
    fs_aseg_to_native.inputs.interp = 'nearest'

    preproc.connect(reconall2, 'aseg',
                    fs_aseg_to_native, 'source_file')

    preproc.connect(reconall2, 'rawavg',
                    fs_aseg_to_native, 'target_file')

    preproc.connect(inputnode, 'subject_dir',
                    fs_aseg_to_native, 'subjects_dir')

    # convert registered FS segmentations from .mgz to .nii.gz
    fs_aseg_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=mri_convert),
                               name='fs_aseg_to_nifti')

    fs_aseg_to_nifti.inputs.args = '-rt nearest'

    preproc.connect(fs_aseg_to_native, 'transformed_file',
                    fs_aseg_to_nifti, 'in_file')

    pick_tissue = pe.Node(util.Function(input_names=['multiatlas_Labels'],
                                        output_names=['csf_mask', 'gm_mask',
                                                      'wm_mask'],
                                        function=pick_tissue_from_labels_file),
                          name=f'{wf_name}_tissue_mask')

    pick_tissue.inputs.include_ventricles = True

    preproc.connect(fs_aseg_to_nifti, 'out_file',
                    pick_tissue, 'multiatlas_Labels')

    preproc.connect(pick_tissue, 'wm_mask',
                    outputnode, 'wm_mask')

    preproc.connect(pick_tissue, 'gm_mask',
                    outputnode, 'gm_mask')

    preproc.connect(pick_tissue, 'csf_mask',
                    outputnode, 'csf_mask')

    return preproc


def connect_anat_segmentation(workflow, strat_list, c, strat_name=None):
    """Segmentation Preprocessing Workflow
    
    Parameters
    ----------
        workflow : workflow
            main preprocessing workflow
        strat_list : list
            list of strategies
        c : configuration
            pipeline configuration
        strat_name : str
            name of strategy

    Returns
    -------
        workflow : workflow
            updated main preprocessing workflow
        strat_list : list
            list of updated strategies
    """  # noqa
    from CPAC.seg_preproc.seg_preproc import (
        create_seg_preproc,
        create_seg_preproc_template_based
    )

    new_strat_list = []

    if True in c.anatomical_preproc['segmentation_workflow']['run']:

        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            seg_preproc = None

            if not any(
                o in c.anatomical_preproc['segmentation_workflow'][
                    '1-segmentation'
                ]['using'] for o in valid_options['segmentation']['using']
            ):
                err = '\n\n[!] C-PAC says: Your segmentation thresholding ' \
                      'options setting does not include any of {0}.\n\n' \
                      'Options you provided:\nanatomical_preproc' \
                      '[\'segmentation_workflow\'][\'1-segmentation\']' \
                      '[\'using\']: {1}\n\n'.format(
                          str(valid_options['segmentation']['using']),
                          str(c.anatomical_preproc['segmentation_workflow'][
                              '1-segmentation'
                          ]['using'])
                      )
                raise Exception(err)

            if strat.get('registration_method') == 'FSL':
                use_ants = False
            elif strat.get('registration_method') == 'ANTS':
                use_ants = True

            if strat_name is not None:
                seg_preproc_wf_name = f'seg_preproc_{strat_name}_{num_strat}'
            else:
                seg_preproc_wf_name = f'seg_preproc_{num_strat}'

            use_threshold = []
            if c.anatomical_preproc['segmentation_workflow'][
                '3-custom_thresholding'
            ]['run'] is True:
                if 'FSL-FAST' in c.anatomical_preproc[
                    'segmentation_workflow'
                ]['1-segmentation']['using']:
                    use_threshold.append('FSL-FAST')
                if c.anatomical_preproc[
                    'segmentation_workflow'
                ]['3-custom_thresholding']['run'] is True:
                    use_threshold.append('Custom')

            seg_preproc = create_seg_preproc(
                use_ants=use_ants,
                use_priors=c.anatomical_preproc['segmentation_workflow'][
                    '2-use_priors']['run'],
                use_threshold=use_threshold,
                csf_use_erosion=c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_csf']['run'],
                wm_use_erosion=c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_wm']['run'],
                gm_use_erosion=c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_gm']['run'],
                wf_name=seg_preproc_wf_name)

            seg_preproc.inputs.csf_threshold.csf_threshold = \
                c.anatomical_preproc['segmentation_workflow'][
                    '3-custom_thresholding']['CSF_threshold_value']
            seg_preproc.inputs.wm_threshold.wm_threshold = \
                c.anatomical_preproc['segmentation_workflow'][
                    '3-custom_thresholding']['WM_threshold_value']
            seg_preproc.inputs.gm_threshold.gm_threshold = \
                c.anatomical_preproc['segmentation_workflow'][
                    '3-custom_thresholding']['GM_threshold_value']

            seg_preproc.inputs.csf_erosion_prop.csf_erosion_prop = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_csf']['csf_erosion_prop']
            seg_preproc.inputs.wm_erosion_prop.wm_erosion_prop = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_wm']['wm_erosion_prop']
            seg_preproc.inputs.gm_erosion_prop.gm_erosion_prop = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_gm']['gm_erosion_prop']

            seg_preproc.inputs.csf_mask_erosion_mm.csf_mask_erosion_mm = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_csf']['csf_mask_erosion_mm']
            seg_preproc.inputs.wm_mask_erosion_mm.wm_mask_erosion_mm = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_wm']['wm_mask_erosion_mm']
            seg_preproc.inputs.gm_mask_erosion_mm.gm_mask_erosion_mm = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_gm']['gm_mask_erosion_mm']

            seg_preproc.inputs.csf_erosion_mm.csf_erosion_mm = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_csf']['csf_erosion_mm']
            seg_preproc.inputs.wm_erosion_mm.wm_erosion_mm = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_wm']['wm_erosion_mm']
            seg_preproc.inputs.gm_erosion_mm.gm_erosion_mm = \
                c.anatomical_preproc['segmentation_workflow'][
                    '4-erosion']['erode_gm']['gm_erosion_mm']

            node, out_file = strat['anatomical_brain_mask']
            workflow.connect(node, out_file,
                             seg_preproc, 'inputspec.brain_mask')

            if seg_preproc is None:
                continue

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             seg_preproc, 'inputspec.brain')

            if strat.get('registration_method') == 'FSL':
                node, out_file = strat['mni_to_anatomical_linear_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_mat')

            elif strat.get('registration_method') == 'ANTS':
                node, out_file = strat['ants_initial_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_init')

                node, out_file = strat['ants_rigid_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_rig')

                node, out_file = strat['ants_affine_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc,
                                 'inputspec.standard2highres_mat')

            priors_c = c.anatomical_preproc['segmentation_workflow'][
                '2-use_priors']

            workflow.connect(priors_c['CSF_path'], 'local_path',
                             seg_preproc, 'inputspec.PRIOR_CSF')

            workflow.connect(priors_c['GM_path'], 'local_path',
                             seg_preproc, 'inputspec.PRIOR_GRAY')

            workflow.connect(priors_c['WM_path'], 'local_path',
                             seg_preproc, 'inputspec.PRIOR_WHITE')

            if False in c.anatomical_preproc['segmentation_workflow']['run']:
                strat = strat.fork()
                new_strat_list.append(strat)

            if c.anatomical_preproc['segmentation_workflow'][
                '4-erosion'
            ]['erode_anatomical_brain_mask']['run'] is True:
                ero_imports = ['import scipy.ndimage as nd',
                               'import numpy as np', 'import nibabel as nb',
                               'import os']
                eroded_mask = pe.Node(
                    util.Function(input_names=['roi_mask', 'skullstrip_mask',
                                               'mask_erosion_mm',
                                               'mask_erosion_prop'],
                                  output_names=['output_roi_mask',
                                                'eroded_skullstrip_mask'],
                                  function=mask_erosion,
                                  imports=ero_imports),
                    name='erode_skullstrip_brain_mask')
                eroded_mask.inputs.mask_erosion_mm = c.anatomical_preproc[
                    'segmentation_workflow'
                ]['4-erosion']['erode_anatomical_brain_mask'][
                    'brain_mask_erosion_mm']

                node, out_file = strat['anatomical_brain_mask']
                workflow.connect(node, out_file,
                                 eroded_mask, 'skullstrip_mask')

                workflow.connect(
                    seg_preproc, 'outputspec.csf_probability_map',
                    eroded_mask, 'roi_mask')

                strat.update_resource_pool({
                    'anatomical_eroded_brain_mask': (
                        eroded_mask, 'eroded_skullstrip_mask')
                })

            strat.append_name(seg_preproc.name)
            strat.update_resource_pool({
                'anatomical_gm_mask': (seg_preproc, 'outputspec.gm_mask'),
                'anatomical_csf_mask': (seg_preproc, 'outputspec.csf_mask'),
                'anatomical_wm_mask': (seg_preproc, 'outputspec.wm_mask'),
                'seg_probability_maps': (seg_preproc,
                                         'outputspec.probability_maps'),
                'seg_mixeltype': (seg_preproc, 'outputspec.mixeltype'),
                'seg_partial_volume_map': (seg_preproc,
                                           'outputspec.partial_volume_map'),
                'seg_partial_volume_files': (seg_preproc,
                                             'outputspec.partial_volume_files')
            })

    strat_list += new_strat_list

    ANTs_Prior_c = c.anatomical_preproc['segmentation_workflow'][
        '1-segmentation']['ANTs_Prior_Based']

    if True in ANTs_Prior_c['run']:

        for num_strat, strat in enumerate(strat_list):

            seg_preproc_ants_prior_based = \
                create_seg_preproc_antsJointLabel_method(
                    wf_name='seg_preproc_ants_prior_{0}'.format(num_strat))

            if seg_preproc_ants_prior_based is None:
                continue

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             seg_preproc_ants_prior_based,
                             'inputspec.anatomical_brain')

            node, out_file = strat['anatomical_brain_mask']
            workflow.connect(node, out_file,
                             seg_preproc_ants_prior_based,
                             'inputspec.anatomical_brain_mask')

            workflow.connect(
                ANTs_Prior_c['template_brain_list'], 'local_path',
                seg_preproc_ants_prior_based, 'inputspec.template_brain_list')
            workflow.connect(
                ANTs_Prior_c['template_segmentation_list'], 'local_path',
                seg_preproc_ants_prior_based,
                'inputspec.template_segmentation_list')
            seg_preproc_ants_prior_based.inputs.inputspec.csf_label = \
                ANTs_Prior_c['CSF_label']
            seg_preproc_ants_prior_based.inputs.inputspec.left_gm_label = \
                ANTs_Prior_c['left_GM_label']
            seg_preproc_ants_prior_based.inputs.inputspec.right_gm_label = \
                ANTs_Prior_c['right_GM_label']
            seg_preproc_ants_prior_based.inputs.inputspec.left_wm_label = \
                ANTs_Prior_c['left_WM_label']
            seg_preproc_ants_prior_based.inputs.inputspec.right_wm_label = \
                ANTs_Prior_c['right_WM_label']

            if False in ANTs_Prior_c['run']:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(seg_preproc_ants_prior_based.name)
            strat.update_resource_pool({
                'anatomical_gm_mask': (
                    seg_preproc_ants_prior_based, 'outputspec.gm_mask'),
                'anatomical_csf_mask': (
                    seg_preproc_ants_prior_based, 'outputspec.csf_mask'),
                'anatomical_wm_mask': (
                    seg_preproc_ants_prior_based, 'outputspec.wm_mask')
            })

    strat_list += new_strat_list

    template_based_c = c.anatomical_preproc['segmentation_workflow'][
        '1-segmentation']['Template_Based']

    if (
        True in template_based_c['run'] and
        'T1 Template' in template_based_c['template_for_segmentation']
    ):
        for num_strat, strat in enumerate(strat_list):

            nodes = strat.get_nodes_names()

            if strat.get('registration_method') == 'FSL':
                use_ants = False
            elif strat.get('registration_method') == 'ANTS':
                use_ants = True

            seg_preproc_template_based = create_seg_preproc_template_based(
                use_ants=use_ants,
                wf_name='seg_preproc_t1_template_{0}'.format(num_strat))

            if seg_preproc_template_based is None:
                continue

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file,
                             seg_preproc_template_based, 'inputspec.brain')

            if strat.get('registration_method') == 'FSL':
                node, out_file = strat['mni_to_anatomical_linear_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_mat')

            elif strat.get('registration_method') == 'ANTS':
                node, out_file = strat['ants_initial_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_init')

                node, out_file = strat['ants_rigid_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_rig')

                node, out_file = strat['ants_affine_xfm']
                workflow.connect(node, out_file,
                                 seg_preproc_template_based,
                                 'inputspec.standard2highres_mat')

            workflow.connect(template_based_c['CSF'], 'local_path',
                             seg_preproc_template_based,
                             'inputspec.CSF_template')

            workflow.connect(template_based_c['GRAY'], 'local_path',
                             seg_preproc_template_based,
                             'inputspec.GRAY_template')

            workflow.connect(template_based_c['WHITE'], 'local_path',
                             seg_preproc_template_based,
                             'inputspec.WHITE_template')

            if False in template_based_c['run']:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(seg_preproc_template_based.name)
            strat.update_resource_pool({
                'anatomical_gm_mask': (
                    seg_preproc_template_based, 'outputspec.gm_mask'),
                'anatomical_csf_mask': (
                    seg_preproc_template_based, 'outputspec.csf_mask'),
                'anatomical_wm_mask': (
                    seg_preproc_template_based, 'outputspec.wm_mask')
            })

    strat_list += new_strat_list

    return workflow, strat_list
