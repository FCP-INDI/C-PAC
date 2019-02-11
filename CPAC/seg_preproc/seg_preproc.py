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
from CPAC.seg_preproc.utils import * 
import nipype.pipeline.engine as pe


def create_seg_preproc(use_ants, wf_name ='seg_preproc'):


    """
    Segment the subject's anatomical brain into cerebral spinal fluids, white matter and gray matter
    and binarize them.

    Parameters
    ----------

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
    
    - Register CSF template in MNI space to t1 space. For details see `flirt <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FLIRT>`_::
    
        flirt
        -in PRIOR_CSF
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres_inv.mat
        -out csf_mni2t1

    - Threshold and binarize CSF probability map ::

        fslmaths
        csf_combo.nii.gz
        -thr 0.4
        -bin csf_bin.nii.gz

    - Generate CSF csf_mask, by applying csf prior in t1 space to thresholded binarized csf probability map ::

        fslmaths
        csf_bin.nii.gz
        -mas csf_mni2t1
        csf_mask


    - Register WM template in MNI space to t1 space ::
        
        flirt
        -in PRIOR_WM
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres.mat
        -out wm_mni2t1

    - Threshold and binarize WM probability map ::

        fslmaths
        wm_combo.nii.gz
        -thr 0.4
        -bin wm_bin.nii.gz

    - Generate WM csf_mask, by applying wm_prior in t1 space to thresholded binarized wm probability map ::

        fslmaths
        wm_bin.nii.gz
        -mas wm_mni2t1
        wm_mask
 
    - Register GM template in MNI space to t1 space ::
    
        flirt
        -in PRIOR_GM
        -ref mprage_brain.nii.gz
        -applyxfm
        -init standard2highres.mat
        -out gm_mni2t1

    - Threshold and binarize GM probability map ::

        fslmaths
        gm_combo.nii.gz
        -thr 0.4
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
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig',
                                                       'PRIOR_CSF',
                                                       'PRIOR_GRAY',
                                                       'PRIOR_WHITE']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['csf_mni2t1',
                                                        'csf_mask',
                                                        'gm_mni2t1',
                                                        'gm_mask',
                                                        'wm_mni2t1',
                                                        'probability_maps',
                                                        'mixeltype',
                                                        'partial_volume_map',
                                                        'partial_volume_files',
                                                        'wm_mask']),
                        name='outputspec')

    segment = pe.Node(interface=fsl.FAST(), name='segment', mem_gb=1.5)
    segment.inputs.img_type = 1
    segment.inputs.segments = True
    segment.inputs.probability_maps = True
    segment.inputs.out_basename = 'segment'

    check_wm = pe.Node(name='check_wm', interface=Function(function=check_if_file_is_empty, input_names=['in_file'], output_names=['out_file']))
    check_gm = pe.Node(name='check_gm', interface=Function(function=check_if_file_is_empty, input_names=['in_file'], output_names=['out_file']))
    check_csf = pe.Node(name='check_csf', interface=Function(function=check_if_file_is_empty, input_names=['in_file'], output_names=['out_file']))

    #connections
    preproc.connect(inputNode, 'brain',
                    segment, 'in_files')

    preproc.connect(segment, 'probability_maps',
                    outputNode, 'probability_maps')
    preproc.connect(segment, 'mixeltype',
                    outputNode, 'mixeltype')
    preproc.connect(segment, 'partial_volume_files',
                    outputNode, 'partial_volume_files')
    preproc.connect(segment, 'partial_volume_map',
                    outputNode, 'partial_volume_map')

    ##get binarize thresholded csf mask
    process_csf = process_segment_map('CSF', use_ants)

    if use_ants == True:
        preproc.connect(inputNode, 'standard2highres_init',
                        process_csf, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        process_csf, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    process_csf, 'inputspec.brain',)
    preproc.connect(inputNode, 'PRIOR_CSF',
                    process_csf, 'inputspec.tissue_prior')

    #tissue_class_files = binary segmented volume file one val for each class
    preproc.connect(segment, ('tissue_class_files', pick_wm_0),
                    process_csf, 'inputspec.probability_map')
    

    preproc.connect(inputNode, 'standard2highres_mat',
                    process_csf, 'inputspec.standard2highres_mat')


    preproc.connect(process_csf, 'outputspec.segment_mask',
                    outputNode, 'csf_mask')

    #get binarize thresholded wm mask
    process_wm = process_segment_map('WM', use_ants)

    if use_ants == True:
        preproc.connect(inputNode, 'standard2highres_init',
                        process_wm, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        process_wm, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    process_wm, 'inputspec.brain',)
    preproc.connect(inputNode, 'PRIOR_WHITE',
                    process_wm, 'inputspec.tissue_prior')
    preproc.connect(segment, ('tissue_class_files', pick_wm_2),
                    process_wm, 'inputspec.probability_map')

    preproc.connect(inputNode, 'standard2highres_mat',
                    process_wm, 'inputspec.standard2highres_mat')

    preproc.connect(process_wm, 'outputspec.tissueprior_mni2t1',
                    outputNode, 'wm_mni2t1')
    preproc.connect(process_wm, 'outputspec.segment_mask',
                    outputNode, 'wm_mask')

    # get binarize thresholded gm mask
    process_gm = process_segment_map('GM', use_ants)

    if use_ants == True:
        preproc.connect(inputNode, 'standard2highres_init',
                        process_gm, 'inputspec.standard2highres_init')
        preproc.connect(inputNode, 'standard2highres_rig',
                        process_gm, 'inputspec.standard2highres_rig')

    preproc.connect(inputNode, 'brain',
                    process_gm, 'inputspec.brain',)
    preproc.connect(inputNode, 'PRIOR_GRAY', 
                    process_gm, 'inputspec.tissue_prior')
    preproc.connect(segment, ('tissue_class_files', pick_wm_1),
                    process_gm, 'inputspec.probability_map')
    preproc.connect(inputNode, 'standard2highres_mat',
                    process_gm, 'inputspec.standard2highres_mat')
    
    preproc.connect(process_gm, 'outputspec.tissueprior_mni2t1',
                    outputNode, 'gm_mni2t1')
    preproc.connect(process_gm, 'outputspec.segment_mask',
                    outputNode, 'gm_mask')

    return preproc


def process_segment_map(wf_name, use_ants):

    """
    This is a sub workflow used inside segmentation workflow to process 
    probability maps obtained in segmentation. Steps include overlapping 
    of the prior tissue with probability maps, thresholding and binarizing 
    it and creating a mask that is used in further analysis.


    Parameters
    ----------
    wf_name : string
        Workflow Name

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
    
        inputspec.tissue_prior : string (existing nifti file)
            path to FSL Standard Tissue prior image 
            
        inputspec.probability_map : string (nifti file)
            tissue Probability map obtained from fsl FAST
        
    Workflow Outputs::

        outputspec.segment_mni2t1 : string (nifti file)
            path to output CSF prior template(in MNI space) registered to anatomical space
    
        outputspec.segment_mask : string (nifti file)
            path to output image after masking segment_combo with its tissue prior in t1 space
        
        
    Order of commands:
 
    - Register tissue prior in MNI space to t1 space. 
    
    - Threshold and binarize segment probability map 
    
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
                                                       'brain',
                                                       'probability_map',
                                                       'standard2highres_init',
                                                       'standard2highres_mat',
                                                       'standard2highres_rig']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['tissueprior_mni2t1',
                                                        'segment_mask']),
                        name='outputspec')

    def form_threshold_string(threshold):
        return '-thr %f -bin ' % (threshold)

    if use_ants == True:

        collect_linear_transforms = pe.Node(util.Merge(3), name='%s_collect_linear_transforms' % (wf_name))

        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='%s_prior_mni_to_t1' % (wf_name))
        tissueprior_mni_to_t1.inputs.invert_transform_flags = [True, True, True]
        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'

        segment_mask = pe.Node(interface=fsl.MultiImageMaths(),
                               name='%s_mask' % (wf_name))
        segment_mask.inputs.op_string = '-mas %s '

        # mni to t1
        preproc.connect(inputNode, 'tissue_prior', tissueprior_mni_to_t1, 'input_image')
        preproc.connect(inputNode, 'brain', tissueprior_mni_to_t1, 'reference_image')

        preproc.connect(inputNode, 'standard2highres_init', collect_linear_transforms, 'in1')
        preproc.connect(inputNode, 'standard2highres_rig', collect_linear_transforms, 'in2')
        preproc.connect(inputNode, 'standard2highres_mat', collect_linear_transforms, 'in3')

        preproc.connect(collect_linear_transforms, 'out', tissueprior_mni_to_t1, 'transforms')

        #create segment mask
        preproc.connect(inputNode, 'probability_map',
                        segment_mask, 'in_file')
        preproc.connect(tissueprior_mni_to_t1, 'output_image', 
                        segment_mask, 'operand_files')


        #connect to output nodes
        preproc.connect(tissueprior_mni_to_t1, 'output_image', 
                        outputNode, 'tissueprior_mni2t1')
        preproc.connect(segment_mask, 'out_file', outputNode, 'segment_mask')

    else:

        tissueprior_mni_to_t1 = pe.Node(interface=fsl.FLIRT(),
                                        name='%s_prior_mni_to_t1' % (wf_name))
        tissueprior_mni_to_t1.inputs.apply_xfm = True
        tissueprior_mni_to_t1.inputs.interp = 'nearestneighbour'

        segment_mask = pe.Node(interface=fsl.MultiImageMaths(),
                               name='%s_mask' % (wf_name))
        segment_mask.inputs.op_string = ' -bin -mas %s'

        # mni to t1
        preproc.connect(inputNode, 'tissue_prior',
                        tissueprior_mni_to_t1, 'in_file')
        preproc.connect(inputNode, 'brain',
                        tissueprior_mni_to_t1, 'reference')
        preproc.connect(inputNode, 'standard2highres_mat',
                        tissueprior_mni_to_t1, 'in_matrix_file')

        # create segment mask
        preproc.connect(inputNode, 'probability_map',
                        segment_mask, 'in_file')
        preproc.connect(tissueprior_mni_to_t1, 'out_file',
                        segment_mask, 'operand_files')

        # connect to output nodes
        preproc.connect(tissueprior_mni_to_t1, 'out_file',
                        outputNode, 'tissueprior_mni2t1')
        preproc.connect(segment_mask, 'out_file',
                        outputNode, 'segment_mask')

    return preproc
