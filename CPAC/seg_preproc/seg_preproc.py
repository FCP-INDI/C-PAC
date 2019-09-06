# Import packages

import os
import sys
import commands
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
    pick_wm_0,
    pick_wm_1,
    pick_wm_2,
    erosion,
    mask_erosion)
    
import nipype.pipeline.engine as pe
import scipy.ndimage as nd
import numpy as np

def create_seg_preproc(use_ants,
                        use_priors,
                        use_threshold,
                        use_erosion,
                        erosion_prop,
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

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                        # 'csf_mni2t1',
                                                        # 'csf_combo',
                                                        # 'csf_thresh',
                                                        # 'csf_bin',
                                                        'csf_mask',
                                                        # 'gm_mni2t1',
                                                        # 'gm_combo',
                                                        # 'gm_thresh',
                                                        # 'gm_bin',
                                                        'gm_mask',
                                                        # 'wm_mni2t1',
                                                        # 'wm_combo',
                                                        # 'wm_thresh',
                                                        # 'wm_bin',
                                                        'wm_mask',
                                                        'probability_maps',
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

    preproc.connect(segment, 'probability_maps',
                    outputNode, 'probability_maps')
    preproc.connect(segment, 'mixeltype',
                    outputNode, 'mixeltype')
    preproc.connect(segment, 'partial_volume_files',
                    outputNode, 'partial_volume_files')
    preproc.connect(segment, 'partial_volume_map',
                    outputNode, 'partial_volume_map')

    process_csf = process_segment_map('CSF', use_ants, use_priors, use_threshold, use_erosion, erosion_prop)


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

    preproc.connect(inputNode, 'PRIOR_CSF',
                    process_csf, 'inputspec.tissue_prior')

    preproc.connect(segment, ('probability_maps', pick_wm_0),
                    process_csf, 'inputspec.probability_tissue_map')

    preproc.connect(inputNode, 'standard2highres_mat',
                    process_csf, 'inputspec.standard2highres_mat')

    # preproc.connect(process_csf, 'outputspec.tissueprior_mni2t1',
    #                 outputNode, 'csf_mni2t1')
    # preproc.connect(process_csf, 'outputspec.segment_combo',
    #                 outputNode, 'csf_combo')
    # preproc.connect(process_csf, 'outputspec.segment_thresh',
    #                 outputNode, 'csf_thresh')
    # preproc.connect(process_csf, 'outputspec.segment_bin',
    #                 outputNode, 'csf_bin')
    preproc.connect(process_csf, 'outputspec.segment_mask',
                    outputNode, 'csf_mask')

    process_wm = process_segment_map('WM', use_ants, use_priors, use_threshold, use_erosion, erosion_prop)


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
    preproc.connect(inputNode, 'PRIOR_WHITE',
                    process_wm, 'inputspec.tissue_prior')
    preproc.connect(segment, ('probability_maps', pick_wm_2),
                    process_wm, 'inputspec.probability_tissue_map')

    preproc.connect(inputNode, 'standard2highres_mat',
                    process_wm, 'inputspec.standard2highres_mat')
    # preproc.connect(process_wm, 'outputspec.tissueprior_mni2t1',
    #                 outputNode, 'wm_mni2t1')
    # preproc.connect(process_wm, 'outputspec.segment_combo',
    #                 outputNode, 'wm_combo')
    # preproc.connect(process_wm, 'outputspec.segment_thresh',
    #                 outputNode, 'wm_thresh')
    # preproc.connect(process_wm, 'outputspec.segment_bin',
    #                 outputNode, 'wm_bin')
    preproc.connect(process_wm, 'outputspec.segment_mask',
                    outputNode, 'wm_mask')

    process_gm = process_segment_map('GM', use_ants, use_priors, use_threshold, use_erosion, erosion_prop)


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
    preproc.connect(inputNode, 'PRIOR_GRAY',
                    process_gm, 'inputspec.tissue_prior')
    preproc.connect(segment, ('probability_maps', pick_wm_1),
                    process_gm, 'inputspec.probability_tissue_map')
    preproc.connect(inputNode, 'standard2highres_mat',
                    process_gm, 'inputspec.standard2highres_mat')
    # preproc.connect(process_gm, 'outputspec.tissueprior_mni2t1',
    #                 outputNode, 'gm_mni2t1')
    # preproc.connect(process_gm, 'outputspec.segment_combo',
    #                 outputNode, 'gm_combo')
    # preproc.connect(process_gm, 'outputspec.segment_thresh',
    #                 outputNode, 'gm_thresh')                
    # preproc.connect(process_gm, 'outputspec.segment_bin',
    #                 outputNode, 'gm_bin')
    preproc.connect(process_gm, 'outputspec.segment_mask',
                    outputNode, 'gm_mask')

    return preproc


def process_segment_map(wf_name,
                        use_ants,
                        use_priors,
                        use_threshold,
                        use_erosion,
                        erosion_prop):
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
    use_threshold: boolean
        Whether or not to use threshold to further refine
        the resulting segmentation tissue masks.
    thresh: float, default value: 0.95
        Thresh values, if use threshold.
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
                                                       'brain',
                                                       'brain_mask',
                                                       'probability_tissue_map',
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=[
                                                        # 'tissueprior_mni2t1',
                                                        # 'segment_combo',
                                                        # 'segment_thresh',
                                                        # 'segment_bin',
                                                        # 'segment_erosion',
                                                        'segment_mask']),
                        name='outputspec')
    
    def form_threshold_string(threshold):
        return '-thr %f ' % (threshold)


    if use_ants:
        collect_linear_transforms = pe.Node(util.Merge(3),
                                            name='{0}_collect_linear_transforms'.format(wf_name))

        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='{0}_prior_mni_to_t1'.format(wf_name))
        tissueprior_mni_to_t1.inputs.invert_transform_flags = [True, True, True]
        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'

        # mni to t1
        preproc.connect(inputNode, 'tissue_prior', tissueprior_mni_to_t1, 'input_image')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference_image')

        preproc.connect(inputNode, 'standard2highres_init', collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig', collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat', collect_linear_transforms, 'in3')

        preproc.connect(collect_linear_transforms, 'out', tissueprior_mni_to_t1, 'transforms')
    
        input_1, value_1 = (inputNode, 'probability_tissue_map')

        if use_priors:
            overlap_segmentmap_with_prior = pe.Node(interface=fsl.MultiImageMaths(), 
                                                    name='overlap_%s_map_with_prior' % (wf_name))
            overlap_segmentmap_with_prior.inputs.op_string = '-mas %s ' 

            preproc.connect(input_1, value_1, overlap_segmentmap_with_prior, 'in_file')
            
            preproc.connect(tissueprior_mni_to_t1, 'output_image', overlap_segmentmap_with_prior, 'operand_files')
            
            input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')


        if use_threshold:
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


        # # create segment mask
        # segment_mask = pe.Node(interface=fsl.MultiImageMaths(),
        #                     name='{0}_mask'.format(wf_name))
        # segment_mask.inputs.op_string = ' -mas %s '
        
        ero_imports = ['import scipy.ndimage as nd' , 'import numpy as np', 'import nibabel as nb', 'import os']

        if use_erosion:
            # mask erosion 
            eroded_mask = pe.Node(util.Function(input_names = ['roi_mask', 'skullstrip_mask', 'mask_erosion_prop'], 
                                                output_names = ['output_roi_mask', 'eroded_skullstrip_mask'], 
                                                function = mask_erosion,
                                                imports = ero_imports),                                    
                                                name='erode_skullstrip_mask_%s' % (wf_name))
            eroded_mask.inputs.mask_erosion_prop =  erosion_prop**3 
            preproc.connect(inputNode, 'brain_mask', eroded_mask, 'skullstrip_mask')
            preproc.connect(input_1, value_1, eroded_mask, 'roi_mask')
            
            input_1, value_1 = (eroded_mask, 'output_roi_mask')

            # erosion 
            erosion_segmentmap = pe.Node(util.Function(input_names = ['roi_mask', 'erosion_prop'], 
                                                output_names = ['eroded_roi_mask'], 
                                                function = erosion,
                                                imports = ero_imports),                                    
                                                name='erosion_segmentmap_%s' % (wf_name))

            erosion_segmentmap.inputs.erosion_prop =  erosion_prop   
            preproc.connect(input_1, value_1, erosion_segmentmap, 'roi_mask')
            input_1, value_1 = (erosion_segmentmap, 'eroded_roi_mask')

        #connect to output nodes
        # preproc.connect(tissueprior_mni_to_t1, 'output_image', outputNode, 'tissueprior_mni2t1')
        
        # preproc.connect(overlap_segmentmap_with_prior, 'out_file', outputNode, 'segment_combo')

        # preproc.connect(segmentmap_threshold, 'out_file', outputNode, 'segment_thresh')
    
        # preproc.connect(binarize_threshold_segmentmap, 'out_file', outputNode, 'segment_bin')
                
        # preproc.connect(erosion_segmentmap, 'out_file', outputNode, 'segment_erosion')  
        
        preproc.connect (input_1, value_1, outputNode, 'segment_mask')


    else:
        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='{0}_prior_mni_to_t1'.format(wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'
        
        # mni to t1
        preproc.connect(inputNode, 'tissue_prior', tissueprior_mni_to_t1, 'input_image')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference_image')

        preproc.connect(inputNode, 'standard2highres_init', collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig', collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat', collect_linear_transforms, 'in3')

        preproc.connect(collect_linear_transforms, 'out', tissueprior_mni_to_t1, 'transforms')
            
        input_1, value_1 = (inputNode, 'probability_tissue_map')

        if use_priors:
            overlap_segmentmap_with_prior = pe.Node(interface=fsl.ImageMaths(), 
                                                    name='overlap_%s_map_with_prior' % (wf_name))
            overlap_segmentmap_with_prior.inputs.op_string = '-mas %s ' 

            preproc.connect(input_1, value_1, overlap_segmentmap_with_prior, 'in_file')
            
            preproc.connect(tissueprior_mni_to_t1, 'output_image', overlap_segmentmap_with_prior, 'operand_files')
            
            input_1, value_1 = (overlap_segmentmap_with_prior, 'out_file')


        if use_threshold:
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


        # # create segment mask
        # segment_mask = pe.Node(interface=fsl.MultiImageMaths(),
        #                     name='{0}_mask'.format(wf_name))
        # segment_mask.inputs.op_string = ' -mas %s '
        
        ero_imports = ['import scipy.ndimage as nd' , 'import numpy as np', 'import nibabel as nb', 'import os']

        if use_erosion:
            # mask erosion 
            eroded_mask = pe.Node(util.Function(input_names = ['roi_mask', 'skullstrip_mask', 'mask_erosion_prop'], 
                                                output_names = ['output_roi_mask', 'eroded_skullstrip_mask'], 
                                                function = mask_erosion,
                                                imports = ero_imports),                                    
                                                name='erode_skullstrip_mask_%s' % (wf_name))
            eroded_mask.inputs.mask_erosion_prop =  erosion_prop**3 
            preproc.connect(inputNode, 'brain_mask', eroded_mask, 'skullstrip_mask')
            preproc.connect(input_1, value_1, eroded_mask, 'roi_mask')
            
            input_1, value_1 = (eroded_mask, 'output_roi_mask')

            # erosion 
            erosion_segmentmap = pe.Node(util.Function(input_names = ['roi_mask', 'erosion_prop'], 
                                                output_names = ['eroded_roi_mask'], 
                                                function = erosion,
                                                imports = ero_imports),                                    
                                                name='erosion_segmentmap_%s' % (wf_name))

            erosion_segmentmap.inputs.erosion_prop =  erosion_prop   
            preproc.connect(input_1, value_1, erosion_segmentmap, 'roi_mask')
            input_1, value_1 = (erosion_segmentmap, 'eroded_roi_mask')

        #connect to output nodes
        # preproc.connect(tissueprior_mni_to_t1, 'output_image', outputNode, 'tissueprior_mni2t1')
        
        # preproc.connect(overlap_segmentmap_with_prior, 'out_file', outputNode, 'segment_combo')
    
        # preproc.connect(binarize_threshold_segmentmap, 'out_file', outputNode, 'segment_bin')
                
        # preproc.connect(erosion_segmentmap, 'out_file', outputNode, 'segment_erosion')  
        
        preproc.connect (input_1, value_1, outputNode, 'segment_mask')

    return preproc


