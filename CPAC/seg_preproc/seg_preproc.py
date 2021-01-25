from nipype.interfaces.utility import Function
import nipype.algorithms.rapidart as ra
from nipype.interfaces import afni, ants, freesurfer, fsl, utility as util
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
from nipype.interfaces import freesurfer

from CPAC.utils.utils import check_prov_for_regtool

from CPAC.anat_preproc.utils import mri_convert
from CPAC.registration.utils import (
    check_transforms,
    generate_inverse_transform_flags)
from CPAC.registration.registration import apply_transform
from CPAC.pipeline.schema import valid_options


def process_segment_map(wf_name, use_priors, use_custom_threshold, reg_tool):
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
                                       'template_to_T1_xfm']),
        name='inputspec')

    outputNode = pe.Node(
        util.IdentityInterface(fields=['segment_mask',
                                       'probability_tissue_map']),
        name='outputspec')

    # FSL-FAST
    #  'tissue_class_files' output is a list of individual binary tissue masks
    #      triggered by 'segments' boolean input (-g or --segments)
    #  'probability_maps' output is a list of individual probability maps
    #      triggered by 'probability_maps' boolean input (-p)

    def form_threshold_string(threshold):
        return '-thr %f ' % (threshold)

    def form_mask_erosion_prop(erosion_prop):
        return erosion_prop ** 3

    if not use_custom_threshold:
        # already binary tissue mask
        input_1, value_1 = (inputNode, 'tissue_class_file')
    else:
        # probability map
        input_1, value_1 = (inputNode, 'probability_tissue_map')

    if use_priors:
        apply_xfm = apply_transform(f'seg_tissue_priors_template_to_T1',
                                    reg_tool=reg_tool)
        apply_xfm.inputs.inputspec.interpolation = "NearestNeighbor"

        preproc.connect(inputNode, 'tissue_prior', apply_xfm,
                        'inputspec.input_image')
        preproc.connect(inputNode, 'brain', apply_xfm,
                        'inputspec.reference')
        preproc.connect(inputNode, 'template_to_T1_xfm', apply_xfm,
                        'inputspec.transform')

        overlap_segmentmap_with_prior = pe.Node(
            interface=fsl.MultiImageMaths(),
            name='overlap_%s_map_with_prior' % (wf_name))
        overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '

        preproc.connect(input_1, value_1,
                        overlap_segmentmap_with_prior, 'in_file')

        preproc.connect(apply_xfm, 'outputspec.output_image',
                        overlap_segmentmap_with_prior, 'operand_files')

        input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')

    if use_custom_threshold:
        segmentmap_threshold = pe.Node(
            interface=fsl.ImageMaths(),
            name='threshold_segmentmap_%s' % (wf_name))
        preproc.connect(inputNode, ('threshold', form_threshold_string),
                        segmentmap_threshold, 'op_string')

        preproc.connect(input_1, value_1, segmentmap_threshold, 'in_file')

        # these are the probability maps, not the binary tissue masks
        input_1, value_1 = (segmentmap_threshold, 'out_file')

        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(),
                                                name='binarize_%s' % (
                                                wf_name))
        binarize_threshold_segmentmap.inputs.op_string = '-bin '

        preproc.connect(input_1, value_1,
                        binarize_threshold_segmentmap, 'in_file')

        input_1, value_1 = (binarize_threshold_segmentmap, 'out_file')

    # regardless of input, they are binary tissue masks now
    preproc.connect(input_1, value_1, outputNode, 'segment_mask')

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
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1,
                        'reference')
        preproc.connect(inputNode, 'standard2highres_mat',
                        tissueprior_mni_to_t1, 'in_matrix_file')

        preproc.connect(tissueprior_mni_to_t1, 'out_file',
                        outputNode, 'segment_mask_temp2t1')

    return preproc


def create_seg_preproc_antsJointLabel_method(
        wf_name='seg_preproc_templated_based'):
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

    preproc = pe.Workflow(name=wf_name)

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


def tissue_seg_fsl_fast(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "tissue_seg_fsl_fast",
     "config": ["segmentation"],
     "switch": ["run"],
     "option_key": ["tissue_segmentation", "using"],
     "option_val": "FSL-FAST",
     "inputs": ["desc-brain_T1w",
                "space-T1w_desc-brain_mask",
                "from-template_to-T1w_mode-image_desc-linear_xfm",
                "CSF_path",
                "GM_path",
                "WM_path"],
     "outputs": ["label-CSF_mask",
                 "label-GM_mask",
                 "label-WM_mask",
                 "label-CSF_probseg",
                 "label-GM_probseg",
                 "label-WM_probseg"]}
    '''

    # FSL-FAST
    #  'tissue_class_files' output is a list of individual binary tissue masks
    #      triggered by 'segments' boolean input (-g or --segments)
    #  'probability_maps' output is a list of individual probability maps
    #      triggered by 'probability_maps' boolean input (-p)

    segment = pe.Node(interface=fsl.FAST(), name=f'segment_{pipe_num}',
                      mem_gb=1.5)
    segment.inputs.img_type = 1
    segment.inputs.segments = True
    segment.inputs.probability_maps = True
    segment.inputs.out_basename = 'segment'

    check_wm = pe.Node(name='check_wm',
                       interface=Function(function=check_if_file_is_empty,
                                          input_names=['in_file'],
                                          output_names=['out_file']))
    check_gm = pe.Node(name='check_gm',
                       interface=Function(function=check_if_file_is_empty,
                                          input_names=['in_file'],
                                          output_names=['out_file']))
    check_csf = pe.Node(name='check_csf', interface=Function(
        function=check_if_file_is_empty, input_names=['in_file'],
        output_names=['out_file']))

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, segment, 'in_files')

    use_custom_threshold = cfg['segmentation']['tissue_segmentation'][
                               'FSL-FAST']['thresholding'] == 'Custom'

    use_priors = cfg['segmentation']['tissue_segmentation']['FSL-FAST'][
        'use_priors']['run']

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-template_to-T1w_mode-image_desc-linear_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    process_csf = process_segment_map(f'CSF_{pipe_num}', use_priors,
                                      use_custom_threshold, reg_tool)
    process_csf.inputs.inputspec.threshold = cfg['segmentation'][
        'tissue_segmentation']['FSL-FAST']['thresholding']['Custom'][
        'CSF_threshold_value']

    if use_priors:
        node, out = strat_pool.get_data('CSF_path')
        wf.connect(node, out, process_csf, 'inputspec.tissue_prior')

    process_gm = process_segment_map(f'GM_{pipe_num}', use_priors,
                                     use_custom_threshold, reg_tool)
    process_gm.inputs.inputspec.threshold = cfg['segmentation'][
        'tissue_segmentation']['FSL-FAST']['thresholding']['Custom'][
        'GM_threshold_value']

    if use_priors:
        node, out = strat_pool.get_data('GM_path')
        wf.connect(node, out, process_gm, 'inputspec.tissue_prior')

    process_wm = process_segment_map(f'WM_{pipe_num}', use_priors,
                                     use_custom_threshold, reg_tool)
    process_wm.inputs.inputspec.threshold = cfg['segmentation'][
        'tissue_segmentation']['FSL-FAST']['thresholding']['Custom'][
        'WM_threshold_value']

    if use_priors:
        node, out = strat_pool.get_data('WM_path')
        wf.connect(node, out, process_wm, 'inputspec.tissue_prior')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, process_csf, 'inputspec.brain')
    wf.connect(node, out, process_gm, 'inputspec.brain')
    wf.connect(node, out, process_wm, 'inputspec.brain')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, process_csf, 'inputspec.brain_mask')
    wf.connect(node, out, process_gm, 'inputspec.brain_mask')
    wf.connect(node, out, process_wm, 'inputspec.brain_mask')

    node, out = strat_pool.get_data(
        'from-template_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out, process_csf, 'inputspec.template_to_T1_xfm')
    wf.connect(node, out, process_gm, 'inputspec.template_to_T1_xfm')
    wf.connect(node, out, process_wm, 'inputspec.template_to_T1_xfm')

    wf.connect(segment, ('tissue_class_files', pick_wm_class_0),
               process_csf, 'inputspec.tissue_class_file')
    wf.connect(segment, ('probability_maps', pick_wm_prob_0),
               process_csf, 'inputspec.probability_tissue_map')

    wf.connect(segment, ('tissue_class_files', pick_wm_class_1),
               process_gm, 'inputspec.tissue_class_file')
    wf.connect(segment, ('probability_maps', pick_wm_prob_1),
               process_gm, 'inputspec.probability_tissue_map')

    wf.connect(segment, ('tissue_class_files', pick_wm_class_2),
               process_wm, 'inputspec.tissue_class_file')
    wf.connect(segment, ('probability_maps', pick_wm_prob_2),
               process_wm, 'inputspec.probability_tissue_map')

    outputs = {
        'label-CSF_probseg': (segment, ('probability_maps', pick_wm_prob_0)),
        'label-GM_probseg': (segment, ('probability_maps', pick_wm_prob_1)),
        'label-WM_probseg': (segment, ('probability_maps', pick_wm_prob_2)),
        'label-CSF_mask': (process_csf, 'outputspec.segment_mask'),
        'label-GM_mask': (process_gm, 'outputspec.segment_mask'),
        'label-WM_mask': (process_wm, 'outputspec.segment_mask')
    }

    return (wf, outputs)


def tissue_seg_T1_template_based(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "tissue_seg_T1_template_based",
     "config": ["segmentation"],
     "switch": ["run"],
     "option_key": ["tissue_segmentation", "using"],
     "option_val": "Template_Based",
     "inputs": ["desc-brain_T1w",
                "from-template_to-T1w_mode-image_desc-linear_xfm"],
     "outputs": ["label-CSF_mask",
                 "label-GM_mask",
                 "label-WM_mask"]}
    '''

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-template_to-T1w_mode-image_desc-linear_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)
    use_ants = reg_tool == 'ants'

    csf_template2t1 = tissue_mask_template_to_t1('CSF', use_ants)
    csf_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['CSF']

    gm_template2t1 = tissue_mask_template_to_t1('GM', use_ants)
    gm_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['GRAY']

    wm_template2t1 = tissue_mask_template_to_t1('WM', use_ants)
    wm_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['WHITE']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, csf_template2t1, 'inputspec.brain')
    wf.connect(node, out, gm_template2t1, 'inputspec.brain')
    wf.connect(node, out, wm_template2t1, 'inputspec.brain')

    node, out = \
        strat_pool.get_data('from-template_to-T1w_mode-image_desc-linear_xfm')
    wf.connect(node, out,
               csf_template2t1, 'inputspec.standard2highres_mat')
    wf.connect(node, out,
               wm_template2t1, 'inputspec.standard2highres_mat')
    wf.connect(node, out,
               gm_template2t1, 'inputspec.standard2highres_mat')

    outputs = {
        'label-CSF_mask': (
        csf_template2t1, 'outputspec.segment_mask_temp2t1'),
        'label-GM_mask': (gm_template2t1, 'outputspec.segment_mask_temp2t1'),
        'label-WM_mask': (wm_template2t1, 'outputspec.segment_mask_temp2t1')
    }

    return (wf, outputs)


def tissue_seg_EPI_template_based(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "tissue_seg_EPI_template_based",
     "config": ["segmentation"],
     "switch": ["run"],
     "option_key": ["tissue_segmentation", "using"],
     "option_val": "Template_Based",
     "inputs": ["desc-mean_bold",
                "from-template_to-bold_mode-image_desc-linear_xfm"],
     "outputs": ["space-bold_label-CSF_mask",
                 "space-bold_label-GM_mask",
                 "space-bold_label-WM_mask"]}
    '''

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-template_to-T1w_mode-image_desc-linear_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)
    use_ants = reg_tool == 'ants'

    csf_template2t1 = tissue_mask_template_to_t1('CSF', use_ants)
    csf_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['CSF']

    gm_template2t1 = tissue_mask_template_to_t1('GM', use_ants)
    gm_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['GRAY']

    wm_template2t1 = tissue_mask_template_to_t1('WM', use_ants)
    wm_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['WHITE']

    node, out = strat_pool.get_data('desc-mean_bold')
    wf.connect(node, out, csf_template2t1, 'inputspec.brain')
    wf.connect(node, out, gm_template2t1, 'inputspec.brain')
    wf.connect(node, out, wm_template2t1, 'inputspec.brain')

    node, out = \
        strat_pool.get_data(
            'from-template_to-bold_mode-image_desc-linear_xfm')
    wf.connect(node, out,
               csf_template2t1, 'inputspec.standard2highres_mat')
    wf.connect(node, out,
               wm_template2t1, 'inputspec.standard2highres_mat')
    wf.connect(node, out,
               gm_template2t1, 'inputspec.standard2highres_mat')

    outputs = {
        'space-bold_label-CSF_mask': (csf_template2t1,
                                      'outputspec.segment_mask_temp2t1'),
        'space-bold_label-GM_mask': (gm_template2t1,
                                     'outputspec.segment_mask_temp2t1'),
        'space-bold_label-WM_mask': (wm_template2t1,
                                     'outputspec.segment_mask_temp2t1')
    }

    return (wf, outputs)


def tissue_seg_ants_prior(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "tissue_seg_ants_prior",
     "config": ["segmentation"],
     "switch": ["run"],
     "option_key": ["tissue_segmentation", "using"],
     "option_val": "ANTs_Prior_Based",
     "inputs": ["desc-brain_T1w",
                "space-T1w_desc-brain_mask"],
     "outputs": ["label-CSF_mask",
                 "label-GM_mask",
                 "label-WM_mask"]}
    '''

    seg_preproc_ants_prior_based = \
        create_seg_preproc_antsJointLabel_method(wf_name=f'seg_preproc_'
                                                         f'ants_prior_'
                                                         f'{pipe_num}')

    seg_preproc_ants_prior_based.inputs.inputspec.template_brain_list = \
        cfg['segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
            'template_brain_list']
    seg_preproc_ants_prior_based.inputs.inputspec.template_segmentation_list = \
        cfg['segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
            'template_segmentation_list']

    seg_preproc_ants_prior_based.inputs.inputspec.csf_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'ANTs_prior_seg_CSF_label']

    seg_preproc_ants_prior_based.inputs.inputspec.left_gm_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'ANTs_prior_seg_left_GM_label']
    seg_preproc_ants_prior_based.inputs.inputspec.right_gm_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'ANTs_prior_seg_right_GM_label']

    seg_preproc_ants_prior_based.inputs.inputspec.left_wm_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'ANTs_prior_seg_left_WM_label']
    seg_preproc_ants_prior_based.inputs.inputspec.right_wm_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'ANTs_prior_seg_right_WM_label']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out,
               seg_preproc_ants_prior_based, 'inputspec.anatomical_brain')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, seg_preproc_ants_prior_based,
               'inputspec.anatomical_brain_mask')

    outputs = {
        'label-CSF_mask': (
        seg_preproc_ants_prior_based, 'outputspec.csf_mask'),
        'label-GM_mask': (seg_preproc_ants_prior_based, 'outputspec.gm_mask'),
        'label-WM_mask': (seg_preproc_ants_prior_based, 'outputspec.wm_mask')
    }

    return (wf, outputs)

