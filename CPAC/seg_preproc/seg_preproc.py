from CPAC.pipeline.nodeblock import nodeblock
from nipype.interfaces.utility import Function
from nipype.interfaces import ants, freesurfer, fsl, utility as util

from CPAC.anat_preproc.utils import freesurfer_hemispheres, mri_convert
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.registration.registration import apply_transform
from CPAC.registration.utils import (
    check_transforms,
    generate_inverse_transform_flags)
from CPAC.seg_preproc.utils import (
    check_if_file_is_empty,
    pick_wm_prob_0,
    pick_wm_prob_1,
    pick_wm_prob_2,
    pick_wm_class_0,
    pick_wm_class_1,
    pick_wm_class_2,
    hardcoded_antsJointLabelFusion)
from CPAC.utils.interfaces.function.seg_preproc import \
    pick_tissue_from_labels_file_interface
from CPAC.utils.utils import check_prov_for_regtool


def process_segment_map(wf_name, use_priors, use_custom_threshold, reg_tool):
    """This is a sub workflow used inside segmentation workflow to process
    probability maps obtained in segmentation. Steps include overlapping
    of the prior tissue with probability maps, thresholding and binarizing
    it and creating a mask that is used in further analysis.

    Parameters
    ----------
    wf_name : string
        Workflow Name
    use_priors : boolean
        Whether or not to use template-space tissue priors to further refine
        the resulting segmentation tissue masks.
    use_threshold : list
        Choose threshold to further refine
        the resulting segmentation tissue masks.
    use_erosion : boolean
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

    """
    # pylint: disable=import-outside-toplevel,redefined-outer-name,reimported
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
            name='overlap_%s_map_with_prior' % (wf_name),
            mem_gb=1.775,
            mem_x=(5022839943792975 / 2417851639229258349412352, 'in_file'))
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
                                                       'gm_label',
                                                       'wm_label']),
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

    pick_tissue = pe.Node(pick_tissue_from_labels_file_interface(),
                          name='{0}_tissue_mask'.format(wf_name))

    preproc.connect(seg_preproc_antsJointLabel, 'multiatlas_Labels',
                    pick_tissue, 'multiatlas_Labels')
    preproc.connect(inputNode, 'csf_label',
                    pick_tissue, 'csf_label')
    preproc.connect(inputNode, 'gm_label',
                    pick_tissue, 'gm_label')
    preproc.connect(inputNode, 'wm_label',
                    pick_tissue, 'wm_label')

    preproc.connect(pick_tissue, 'csf_mask',
                    outputNode, 'csf_mask')
    preproc.connect(pick_tissue, 'gm_mask',
                    outputNode, 'gm_mask')
    preproc.connect(pick_tissue, 'wm_mask',
                    outputNode, 'wm_mask')

    return preproc


@nodeblock(
    name="tissue_seg_fsl_fast",
    config=["segmentation"],
    switch=["run"],
    option_key=["tissue_segmentation", "using"],
    option_val="FSL-FAST",
    inputs=[
        (
            ["desc-brain_T1w", "space-longitudinal_desc-brain_T1w"],
            ["space-T1w_desc-brain_mask", "space-longitudinal_desc-brain_mask"],
            [
                "from-template_to-T1w_mode-image_desc-linear_xfm",
                "from-template_to-longitudinal_mode-image_desc-linear_xfm",
            ],
        ),
        "CSF-path",
        "GM-path",
        "WM-path",
    ],
    outputs=[
        "label-CSF_mask",
        "label-GM_mask",
        "label-WM_mask",
        "label-CSF_desc-preproc_mask",
        "label-GM_desc-preproc_mask",
        "label-WM_desc-preproc_mask",
        "label-CSF_probseg",
        "label-GM_probseg",
        "label-WM_probseg",
        "label-CSF_pveseg",
        "label-GM_pveseg",
        "label-WM_pveseg",
        "space-longitudinal_label-CSF_mask",
        "space-longitudinal_label-GM_mask",
        "space-longitudinal_label-WM_mask",
        "space-longitudinal_label-CSF_desc-preproc_mask",
        "space-longitudinal_label-GM_desc-preproc_mask",
        "space-longitudinal_label-WM_desc-preproc_mask",
        "space-longitudinal_label-CSF_probseg",
        "space-longitudinal_label-GM_probseg",
        "space-longitudinal_label-WM_probseg",
    ],
)
def tissue_seg_fsl_fast(wf, cfg, strat_pool, pipe_num, opt=None):
    # FSL-FAST
    #  'tissue_class_files' output is a list of individual binary tissue masks
    #      triggered by 'segments' boolean input (-g or --segments)
    #  'probability_maps' output is a list of individual probability maps
    #      triggered by 'probability_maps' boolean input (-p)

    segment = pe.Node(interface=fsl.FAST(),
                      name=f'segment_{pipe_num}',
                      mem_gb=3.48,
                      mem_x=(3444233104315183 / 19342813113834066795298816,
                             'in_files'))
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

    connect, resource = \
        strat_pool.get_data(["desc-brain_T1w",
                             "space-longitudinal_desc-brain_T1w"],
                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, segment, 'in_files')


    use_custom_threshold = cfg['segmentation']['tissue_segmentation'][
                               'FSL-FAST']['thresholding'][
                               'use'] == 'Custom'

    use_priors = cfg['segmentation']['tissue_segmentation'][
        'FSL-FAST']['use_priors']['run']

    long = ''
    if 'space-longitudinal' in resource:
        long = 'space-longitudinal_'

    if use_priors:
        xfm = 'from-template_to-T1w_mode-image_desc-linear_xfm'
        if 'space-longitudinal' in resource:
            xfm = 'from-template_to-longitudinal_mode-image_desc-linear_xfm'
        xfm_prov = strat_pool.get_cpac_provenance(xfm)
        reg_tool = check_prov_for_regtool(xfm_prov)
    else:
        xfm_prov = None
        reg_tool = None
        xfm = None


    process_csf = process_segment_map(f'CSF_{pipe_num}', use_priors,
                                      use_custom_threshold, reg_tool)
    process_csf.inputs.inputspec.threshold = cfg['segmentation'][
        'tissue_segmentation']['FSL-FAST']['thresholding']['Custom'][
        'CSF_threshold_value']

    get_pve_csf = pe.Node(interface=fsl.maths.MathsCommand(),
                          name=f'get_pve_csf_{pipe_num}')
    get_pve_csf.inputs.args = '-thr 0.5 -uthr 1.5 -bin'
    wf.connect(segment, 'partial_volume_map', get_pve_csf, 'in_file')

    get_pve_gm = pe.Node(interface=fsl.maths.MathsCommand(),
                          name=f'get_pve_gm_{pipe_num}')
    get_pve_gm.inputs.args = '-thr 1.5 -uthr 2.5 -bin'
    wf.connect(segment, 'partial_volume_map', get_pve_gm, 'in_file')

    get_pve_wm = pe.Node(interface=fsl.maths.MathsCommand(),
                          name=f'get_pve_wm_{pipe_num}')
    get_pve_wm.inputs.args = '-thr 2.5 -uthr 3.5 -bin'
    wf.connect(segment, 'partial_volume_map', get_pve_wm, 'in_file')

    if use_priors:
        node, out = strat_pool.get_data('CSF-path')
        wf.connect(node, out, process_csf, 'inputspec.tissue_prior')

    process_gm = process_segment_map(f'GM_{pipe_num}', use_priors,
                                     use_custom_threshold, reg_tool)
    process_gm.inputs.inputspec.threshold = cfg['segmentation'][
        'tissue_segmentation']['FSL-FAST']['thresholding']['Custom'][
        'GM_threshold_value']

    if use_priors:
        node, out = strat_pool.get_data('GM-path')
        wf.connect(node, out, process_gm, 'inputspec.tissue_prior')

    process_wm = process_segment_map(f'WM_{pipe_num}', use_priors,
                                     use_custom_threshold, reg_tool)
    process_wm.inputs.inputspec.threshold = cfg['segmentation'][
        'tissue_segmentation']['FSL-FAST']['thresholding']['Custom'][
        'WM_threshold_value']

    if use_priors:
        node, out = strat_pool.get_data('WM-path')
        wf.connect(node, out, process_wm, 'inputspec.tissue_prior')

    node, out = strat_pool.get_data(["desc-brain_T1w",
                                     "space-longitudinal_desc-brain_T1w"])
    wf.connect(node, out, process_csf, 'inputspec.brain')
    wf.connect(node, out, process_gm, 'inputspec.brain')
    wf.connect(node, out, process_wm, 'inputspec.brain')

    node, out = strat_pool.get_data(["space-T1w_desc-brain_mask",
                                     "space-longitudinal_desc-brain_mask"])
    wf.connect(node, out, process_csf, 'inputspec.brain_mask')
    wf.connect(node, out, process_gm, 'inputspec.brain_mask')
    wf.connect(node, out, process_wm, 'inputspec.brain_mask')

    if use_priors:
        node, out = strat_pool.get_data(xfm)
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

    get_csf = pe.Node(util.Function(input_names=['probability_maps'],
                                    output_names=['filename'],
                                    function=pick_wm_prob_0),
                      name=f'get_csf_{pipe_num}')

    wf.connect(segment, 'probability_maps', get_csf, 'probability_maps')

    outputs = {
        f'{long}label-CSF_probseg': (get_csf, 'filename'),
        f'{long}label-GM_probseg':
            (segment, ('probability_maps', pick_wm_prob_1)),
        f'{long}label-WM_probseg':
            (segment, ('probability_maps', pick_wm_prob_2)),
        f'{long}label-CSF_mask':
            (segment, ('tissue_class_files', pick_wm_class_0)),
        f'{long}label-GM_mask':
            (segment, ('tissue_class_files', pick_wm_class_1)),
        f'{long}label-WM_mask':
            (segment, ('tissue_class_files', pick_wm_class_2)),
        f'{long}label-CSF_desc-preproc_mask':
            (process_csf, 'outputspec.segment_mask'),
        f'{long}label-GM_desc-preproc_mask':
            (process_gm, 'outputspec.segment_mask'),
        f'{long}label-WM_desc-preproc_mask':
            (process_wm, 'outputspec.segment_mask'),
        f'{long}label-CSF_pveseg': (get_pve_csf, 'out_file'),
        f'{long}label-GM_pveseg': (get_pve_gm, 'out_file'),
        f'{long}label-WM_pveseg': (get_pve_wm, 'out_file'),
    }

    return (wf, outputs)


@nodeblock(
    name="tissue_seg_T1_template_based",
    config=["segmentation"],
    switch=["run"],
    option_key=["tissue_segmentation", "using"],
    option_val="Template_Based",
    inputs=[("desc-brain_T1w", "from-template_to-T1w_mode-image_desc-linear_xfm")],
    outputs=["label-CSF_mask", "label-GM_mask", "label-WM_mask"],
)
def tissue_seg_T1_template_based(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-template_to-T1w_mode-image_desc-linear_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)
    use_ants = reg_tool == 'ants'

    csf_template2t1 = tissue_mask_template_to_t1(f'CSF_{pipe_num}', 
                                                 use_ants)
    csf_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['CSF']

    gm_template2t1 = tissue_mask_template_to_t1(f'GM_{pipe_num}',
                                                use_ants)
    gm_template2t1.inputs.inputspec.tissue_mask_template = cfg[
        'segmentation']['tissue_segmentation']['Template_Based']['GRAY']

    wm_template2t1 = tissue_mask_template_to_t1(f'WM_{pipe_num}',
                                                use_ants)
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


@nodeblock(
    name="tissue_seg_EPI_template_based",
    config=["segmentation"],
    switch=["run"],
    option_key=["tissue_segmentation", "using"],
    option_val="Template_Based",
    inputs=[("desc-mean_bold", "from-EPItemplate_to-bold_mode-image_desc-linear_xfm")],
    outputs=[
        "space-bold_label-CSF_mask",
        "space-bold_label-GM_mask",
        "space-bold_label-WM_mask",
    ],
)
def tissue_seg_EPI_template_based(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-EPItemplate_to-bold_mode-image_desc-linear_xfm')
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
            'from-EPItemplate_to-bold_mode-image_desc-linear_xfm')
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


@nodeblock(
    name="tissue_seg_ants_prior",
    config=["segmentation"],
    switch=["run"],
    option_key=["tissue_segmentation", "using"],
    option_val="ANTs_Prior_Based",
    inputs=[
        (
            "desc-brain_T1w",
            ["space-T1w_desc-brain_mask", "space-T1w_desc-acpcbrain_mask"],
        )
    ],
    outputs=["label-CSF_mask", "label-GM_mask", "label-WM_mask"],
)
def tissue_seg_ants_prior(wf, cfg, strat_pool, pipe_num, opt=None):

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
        'CSF_label']

    seg_preproc_ants_prior_based.inputs.inputspec.gm_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'GM_label']

    seg_preproc_ants_prior_based.inputs.inputspec.wm_label = cfg[
        'segmentation']['tissue_segmentation']['ANTs_Prior_Based'][
        'WM_label']

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out,
               seg_preproc_ants_prior_based, 'inputspec.anatomical_brain')

    node, out = strat_pool.get_data(['space-T1w_desc-brain_mask',
                                     'space-T1w_desc-acpcbrain_mask'])
    wf.connect(node, out, seg_preproc_ants_prior_based,
               'inputspec.anatomical_brain_mask')

    outputs = {
        'label-CSF_mask': (
        seg_preproc_ants_prior_based, 'outputspec.csf_mask'),
        'label-GM_mask': (seg_preproc_ants_prior_based, 'outputspec.gm_mask'),
        'label-WM_mask': (seg_preproc_ants_prior_based, 'outputspec.wm_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="tissue_seg_freesurfer",
    config=["segmentation"],
    switch=["run"],
    option_key=["tissue_segmentation", "using"],
    option_val="FreeSurfer",
    inputs=[
        "freesurfer-subject-dir",
        "pipeline-fs_raw-average",
        "pipeline-fs_subcortical-seg",
    ],
    outputs=[
        "pipeline-fs_hemi-L_desc-surface_curv",
        "pipeline-fs_hemi-R_desc-surface_curv",
        "pipeline-fs_hemi-L_desc-surfaceMesh_pial",
        "pipeline-fs_hemi-R_desc-surfaceMesh_pial",
        "pipeline-fs_hemi-L_desc-surfaceMesh_smoothwm",
        "pipeline-fs_hemi-R_desc-surfaceMesh_smoothwm",
        "pipeline-fs_hemi-L_desc-surfaceMesh_sphere",
        "pipeline-fs_hemi-R_desc-surfaceMesh_sphere",
        "pipeline-fs_hemi-L_desc-surfaceMap_sulc",
        "pipeline-fs_hemi-R_desc-surfaceMap_sulc",
        "pipeline-fs_hemi-L_desc-surfaceMap_thickness",
        "pipeline-fs_hemi-R_desc-surfaceMap_thickness",
        "pipeline-fs_hemi-L_desc-surfaceMap_volume",
        "pipeline-fs_hemi-R_desc-surfaceMap_volume",
        "pipeline-fs_hemi-L_desc-surfaceMesh_white",
        "pipeline-fs_hemi-R_desc-surfaceMesh_white",
        "label-CSF_mask",
        "label-GM_mask",
        "label-WM_mask",
    ],
)
def tissue_seg_freesurfer(wf, cfg, strat_pool, pipe_num, opt=None):

    node, out = strat_pool.get_data('freesurfer-subject-dir')

    fs_aseg_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                                name=f'fs_aseg_to_native_{pipe_num}')
    fs_aseg_to_native.inputs.reg_header = True
    fs_aseg_to_native.inputs.interp = 'nearest'

    wf.connect(node, out, fs_aseg_to_native, 'subjects_dir')

    node, out = strat_pool.get_data('pipeline-fs_subcortical-seg')
    wf.connect(node, out, fs_aseg_to_native, 'source_file')

    node, out = strat_pool.get_data('pipeline-fs_raw-average')
    wf.connect(node, out, fs_aseg_to_native, 'target_file')

    fs_aseg_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=mri_convert),
                               name=f'fs_aseg_to_nifti_{pipe_num}')
    fs_aseg_to_nifti.inputs.args = '-rt nearest'

    wf.connect(fs_aseg_to_native, 'transformed_file',
               fs_aseg_to_nifti, 'in_file')

    pick_tissue = pe.Node(pick_tissue_from_labels_file_interface(),
                          name=f'select_fs_tissue_{pipe_num}')

    pick_tissue.inputs.csf_label = cfg['segmentation'][
        'tissue_segmentation']['FreeSurfer']['CSF_label']
    pick_tissue.inputs.gm_label = cfg['segmentation'][
        'tissue_segmentation']['FreeSurfer']['GM_label']
    pick_tissue.inputs.wm_label = cfg['segmentation'][
        'tissue_segmentation']['FreeSurfer']['WM_label']

    wf.connect(fs_aseg_to_nifti, 'out_file', pick_tissue, 'multiatlas_Labels')

    erode_tissues = {}
    if cfg['segmentation']['tissue_segmentation']['FreeSurfer']['erode'] > 0:
        for tissue in ['csf', 'wm', 'gm']:
            erode_tissues[tissue] = pe.Node(
                interface=freesurfer.model.Binarize(),
                name=f'erode_{tissue}_{pipe_num}')
            erode_tissues[tissue].inputs.match = [1]
            erode_tissues[tissue].inputs.erode = cfg['segmentation'][
                'tissue_segmentation']['FreeSurfer']['erode']
            wf.connect(pick_tissue, f'{tissue}_mask', erode_tissues[tissue],
                       'in_file')

    if erode_tissues:
        outputs = {
            'label-CSF_mask': (erode_tissues['csf'], 'binary_file'), 
            'label-WM_mask': (erode_tissues['wm'], 'binary_file'),
            'label-GM_mask': (erode_tissues['gm'], 'binary_file')
            }

    else:
        outputs = {
            'label-CSF_mask': (pick_tissue, 'csf_mask'),
            'label-WM_mask': (pick_tissue, 'wm_mask'),
            'label-GM_mask': (pick_tissue, 'gm_mask')
        }

    return (wf, outputs)
