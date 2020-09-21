# -*- coding: utf-8 -*-
import os
from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces import freesurfer
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.anat_preproc.ants import init_brain_extraction_wf
from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string, mri_convert
from CPAC.utils.datasource import create_check_for_s3_node
from CPAC.unet.function import predict_volumes

def patch_cmass_output(lst, index=0):
    """
    Parameters
    ----------
    lst : list of tuples
        output of afni.CenterMass()
    index : int
        index in the list of tuples

    Returns
    -------
        tuple
            one set of center of mass coordinates
    """
    if len(lst) <= index:
        raise IndexError("lst index out of range")
    return lst[index]

def create_anat_preproc(method='afni', already_skullstripped=False, config=None, wf_name='anat_preproc', sub_dir=None):
    """
    The main purpose of this workflow is to process T1 scans. Raw mprage file is deobliqued, reoriented
    into RPI and skullstripped. Also, a whole brain only mask is generated from the skull stripped image
    for later use in registration.

    Parameters
    ----------
    method : string or None
        Skull-stripping method, None if no skullstripping required
    
    already_skullstripped : boolean

    config : configuration
        Pipeline configuration object

    wf_name : string
        Anatomical preprocessing workflow name

    sub_dir : path
        Subject-specific working directory

    Returns
    -------
    anat_preproc : workflow
        Anatomical preprocessing workflow

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/anat_preproc/anat_preproc.py>`_

    Workflow Inputs::
        inputspec.anat : string
            User input anatomical (T1) Image, in any of the 8 orientations

    Workflow Outputs::

        outputspec.refit : string
            Path to deobliqued anatomical image

        outputspec.reorient : string
            Path to RPI oriented anatomical image

        outputspec.skullstrip : string
            Path to skull stripped RPI oriented mprage file with normalized intensities.

        outputspec.brain : string
            Path to skull stripped RPI brain image with original intensity values and not normalized or scaled.

    Order of commands:
    - Deobliqing the scans. ::
        3drefit -deoblique mprage.nii.gz

    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior  (RPI) orientation ::
        3dresample -orient RPI
                   -prefix mprage_RPI.nii.gz
                   -inset mprage.nii.gz

    - Skull-Stripping the image ::
        Using AFNI ::
            3dSkullStrip -input mprage_RPI.nii.gz
                         -o_ply mprage_RPI_3dT.nii.gz
        or using BET ::
            bet mprage_RPI.nii.gz

    - The skull-stripping step modifies the intensity values. To get back the original intensity values, we do an element wise product of RPI data with step function of skull-stripped data ::
        3dcalc -a mprage_RPI.nii.gz
               -b mprage_RPI_3dT.nii.gz
               -expr 'a*step(b)'
               -prefix mprage_RPI_3dc.nii.gz

    High Level Workflow Graph:
    .. image:: ../images/anatpreproc_graph.dot.png
       :width: 500

    Detailed Workflow Graph:
    .. image:: ../images/anatpreproc_graph_detailed.dot.png
       :width: 500

    Examples
    --------
    >>> from CPAC.anat_preproc import create_anat_preproc
    >>> preproc = create_anat_preproc()
    >>> preproc.inputs.inputspec.anat = 'sub1/anat/mprage.nii.gz'
    >>> preproc.run() #doctest: +SKIP
    """

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat', 
                                                       'brain_mask',
                                                       'template_brain_only_for_anat',
                                                       'template_skull_for_anat',
                                                       'template_cmass']), 
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['refit',
                                                        'reorient',
                                                        'skullstrip',
                                                        'brain',
                                                        'brain_mask',
                                                        'freesurfer_subject_dir',
                                                        'center_of_mass']),
                        name='outputspec')

    anat_deoblique = pe.Node(interface=afni.Refit(),
                             name='anat_deoblique')
    anat_deoblique.inputs.deoblique = True
    preproc.connect(inputnode, 'anat', anat_deoblique, 'in_file')

    preproc.connect(anat_deoblique, 'out_file', outputnode, 'refit')    
    # Disable non_local_means_filtering and n4_bias_field_correction when run niworkflows-ants
    if method == 'niworkflows-ants':
        config.non_local_means_filtering = False 
        config.n4_bias_field_correction = False
        
    if config.non_local_means_filtering and config.n4_bias_field_correction:
        denoise = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise')
        preproc.connect(anat_deoblique, 'out_file', denoise, 'input_image')
        n4 = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, shrink_factor=2, copy_header=True),
            name='anat_n4')
        preproc.connect(denoise, 'output_image', n4, 'input_image')
    elif config.non_local_means_filtering and not config.n4_bias_field_correction:
        denoise = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise')
        preproc.connect(anat_deoblique, 'out_file', denoise, 'input_image')
    elif not config.non_local_means_filtering and config.n4_bias_field_correction:
        n4 = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, shrink_factor=2, copy_header=True),
            name='anat_n4')
        preproc.connect(anat_deoblique, 'out_file', n4, 'input_image')

    # Anatomical reorientation
    anat_reorient = pe.Node(interface=afni.Resample(),
                            name='anat_reorient')
    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'
    
    if config.n4_bias_field_correction:
        preproc.connect(n4, 'output_image', anat_reorient, 'in_file')
    elif config.non_local_means_filtering and not config.n4_bias_field_correction:
        preproc.connect(denoise, 'output_image', anat_reorient, 'in_file')
    else:
        preproc.connect(anat_deoblique, 'out_file', anat_reorient, 'in_file')

    preproc.connect(anat_reorient, 'out_file', outputnode, 'reorient')

    if already_skullstripped:

        anat_skullstrip = pe.Node(interface=util.IdentityInterface(fields=['out_file']),
                                    name='anat_skullstrip')

        preproc.connect(anat_reorient, 'out_file',
                        anat_skullstrip, 'out_file')

        preproc.connect(anat_skullstrip, 'out_file',
                        outputnode, 'skullstrip')
        
        # binarize skullstripped brain to get brain mask
        brain_mask = pe.Node(interface=fsl.maths.MathsCommand(), name='brain_mask')
        brain_mask.inputs.args = '-bin'

        preproc.connect(anat_skullstrip, 'out_file',
                        brain_mask, 'in_file')

        preproc.connect(brain_mask, 'out_file',
                        outputnode, 'brain_mask')

        preproc.connect(anat_skullstrip, 'out_file',
                        outputnode, 'brain')

    else:

        if method == 'afni':
            # Skull-stripping using AFNI 3dSkullStrip
            inputnode_afni = pe.Node(
                util.IdentityInterface(fields=['mask_vol',
                                               'shrink_factor',
                                               'var_shrink_fac',
                                               'shrink_fac_bot_lim',
                                               'avoid_vent',
                                               'niter',
                                               'pushout',
                                               'touchup',
                                               'fill_hole',
                                               'avoid_eyes',
                                               'use_edge',
                                               'exp_frac',
                                               'smooth_final',
                                               'push_to_edge',
                                               'use_skull',
                                               'perc_int',
                                               'max_inter_iter',
                                               'blur_fwhm',
                                               'fac',
                                               'monkey']),
                name='AFNI_options')

            skullstrip_args = pe.Node(util.Function(input_names=['spat_norm',
                                                                 'spat_norm_dxyz',
                                                                 'mask_vol',
                                                                 'shrink_fac',
                                                                 'var_shrink_fac',
                                                                 'shrink_fac_bot_lim',
                                                                 'avoid_vent',
                                                                 'niter',
                                                                 'pushout',
                                                                 'touchup',
                                                                 'fill_hole',
                                                                 'avoid_eyes',
                                                                 'use_edge',
                                                                 'exp_frac',
                                                                 'smooth_final',
                                                                 'push_to_edge',
                                                                 'use_skull',
                                                                 'perc_int',
                                                                 'max_inter_iter',
                                                                 'blur_fwhm',
                                                                 'fac',
                                                                 'monkey'],
                                                    output_names=['expr'],
                                                    function=create_3dskullstrip_arg_string),
                                      name='anat_skullstrip_args')

            preproc.connect([
                (inputnode_afni, skullstrip_args, [
                    ('mask_vol', 'mask_vol'),
                    ('shrink_factor', 'shrink_fac'),
                    ('var_shrink_fac', 'var_shrink_fac'),
                    ('shrink_fac_bot_lim', 'shrink_fac_bot_lim'),
                    ('avoid_vent', 'avoid_vent'),
                    ('niter', 'niter'),
                    ('pushout', 'pushout'),
                    ('touchup', 'touchup'),
                    ('fill_hole', 'fill_hole'),
                    ('avoid_eyes', 'avoid_eyes'),
                    ('use_edge', 'use_edge'),
                    ('exp_frac', 'exp_frac'),
                    ('smooth_final', 'smooth_final'),
                    ('push_to_edge', 'push_to_edge'),
                    ('use_skull', 'use_skull'),
                    ('perc_int', 'perc_int'),
                    ('max_inter_iter', 'max_inter_iter'),
                    ('blur_fwhm', 'blur_fwhm'),
                    ('fac', 'fac'),
                    ('monkey','monkey')
                ])
            ])

            anat_skullstrip = pe.Node(interface=afni.SkullStrip(),
                                      name='anat_skullstrip')

            anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip, 'in_file')

            preproc.connect(skullstrip_args, 'expr',
                            anat_skullstrip, 'args')

            # Generate anatomical brain mask

            anat_brain_mask = pe.Node(interface=afni.Calc(),
                                            name='anat_brain_mask')

            anat_brain_mask.inputs.expr = 'step(a)'
            anat_brain_mask.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_skullstrip, 'out_file',
                            anat_brain_mask, 'in_file_a')

            # Apply skull-stripping step mask to original volume
            anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                            name='anat_skullstrip_orig_vol')

            anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
            anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_a')

            preproc.connect(anat_brain_mask, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_b')

            preproc.connect(anat_brain_mask, 'out_file',
                            outputnode, 'brain_mask')

            preproc.connect(anat_skullstrip_orig_vol, 'out_file',
                            outputnode, 'brain')

        elif method == 'fsl':
            # Skull-stripping using FSL BET
            inputnode_bet = pe.Node(
                util.IdentityInterface(fields=['frac',
                                               'mask_boolean',
                                               'mesh_boolean',
                                               'outline',
                                               'padding',
                                               'radius',
                                               'reduce_bias',
                                               'remove_eyes',
                                               'robust',
                                               'skull',
                                               'surfaces',
                                               'threshold',
                                               'vertical_gradient']),
                name='BET_options')

            anat_skullstrip = pe.Node(
                interface=fsl.BET(), name='anat_skullstrip')
            anat_skullstrip.inputs.output_type = 'NIFTI_GZ'
            
            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip, 'in_file')

            preproc.connect([
                (inputnode_bet, anat_skullstrip, [
                    ('frac', 'frac'),
                    ('mask_boolean', 'mask'),
                    ('mesh_boolean', 'mesh'),
                    ('outline', 'outline'),
                    ('padding', 'padding'),
                    ('radius', 'radius'),
                    ('reduce_bias', 'reduce_bias'),
                    ('remove_eyes', 'remove_eyes'),
                    ('robust', 'robust'),
                    ('skull', 'skull'),
                    ('surfaces', 'surfaces'),
                    ('threshold', 'threshold'),
                    ('vertical_gradient', 'vertical_gradient'),
                ])
            ])

            preproc.connect(anat_skullstrip, 'out_file',
                            outputnode, 'skullstrip')

            # Apply skull-stripping step mask to original volume
            anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                            name='anat_skullstrip_orig_vol')

            anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
            anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_a')

            preproc.connect(anat_skullstrip, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_b')

            preproc.connect(anat_skullstrip, 'mask_file',
                            outputnode, 'brain_mask')

            preproc.connect(anat_skullstrip_orig_vol, 'out_file',
                            outputnode, 'brain')

        elif method == 'niworkflows-ants': 
            # Skull-stripping using niworkflows-ants  
            anat_skullstrip_ants = init_brain_extraction_wf(tpl_target_path=config.niworkflows_ants_template_path,
                                                            tpl_mask_path=config.niworkflows_ants_mask_path,
                                                            tpl_regmask_path=config.niworkflows_ants_regmask_path,
                                                            name='anat_skullstrip_ants')

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_ants, 'inputnode.in_files')

            preproc.connect(anat_skullstrip_ants, 'copy_xform.out_file',
                            outputnode, 'skullstrip')

            preproc.connect(anat_skullstrip_ants, 'copy_xform.out_file',
                            outputnode, 'brain')

            preproc.connect(anat_skullstrip_ants, 'atropos_wf.copy_xform.out_mask',
                            outputnode, 'brain_mask')
        
        elif method == 'freesurfer':

            reconall = pe.Node(interface=freesurfer.ReconAll(),
                            name='anat_freesurfer')
            reconall.inputs.directive = 'autorecon1'
            
            num_strat = len(config.skullstrip_option)-1 # the number of strategy that FreeSurfer will be

            freesurfer_subject_dir = os.path.join(sub_dir, f'anat_preproc_freesurfer_{num_strat}', 'anat_freesurfer')
            
            # create the node dir
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
            if not os.path.exists(os.path.join(sub_dir, f'anat_preproc_freesurfer_{num_strat}')):
                os.mkdir(os.path.join(sub_dir, f'anat_preproc_freesurfer_{num_strat}'))
            if not os.path.exists(freesurfer_subject_dir):
                os.mkdir(freesurfer_subject_dir)
            
            reconall.inputs.subjects_dir = freesurfer_subject_dir
            reconall.inputs.openmp = config.num_omp_threads
            preproc.connect(anat_reorient, 'out_file',
                            reconall, 'T1_files')
            preproc.connect(reconall, 'subjects_dir',
                            outputnode, 'freesurfer_subject_dir')

            # register FS brain mask to native space
            fs_brain_mask_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                            name='fs_brain_mask_to_native')
            fs_brain_mask_to_native.inputs.reg_header = True
            preproc.connect(reconall, 'brainmask',
                            fs_brain_mask_to_native, 'source_file')
            preproc.connect(reconall, 'rawavg',
                            fs_brain_mask_to_native, 'target_file')
              
            # convert brain mask file from .mgz to .nii.gz
            fs_brain_mask_to_nifti = pe.Node(util.Function(input_names=['in_file'], 
                                              output_names=['out_file'],
                                              function=mri_convert),                        
                                    name='fs_brain_mask_to_nifti')
            preproc.connect(fs_brain_mask_to_native, 'transformed_file',
                            fs_brain_mask_to_nifti, 'in_file')
            preproc.connect(fs_brain_mask_to_nifti, 'out_file',
                            outputnode, 'brain_mask')

            # apply mask
            fs_brain = pe.Node(interface=fsl.ApplyMask(),
                            name='anat_skullstrip')
            preproc.connect(anat_reorient, 'out_file',
                            fs_brain, 'in_file')
            preproc.connect(fs_brain_mask_to_nifti, 'out_file',
                            fs_brain, 'mask_file')
            preproc.connect(fs_brain, 'out_file',
                            outputnode, 'brain')

        elif method == 'mask':

            brain_mask_deoblique = pe.Node(interface=afni.Refit(),
                                    name='brain_mask_deoblique')
            brain_mask_deoblique.inputs.deoblique = True
            preproc.connect(inputnode, 'brain_mask',
                            brain_mask_deoblique, 'in_file')

            brain_mask_reorient = pe.Node(interface=afni.Resample(),
                                    name='brain_mask_reorient')
            brain_mask_reorient.inputs.orientation = 'RPI'
            brain_mask_reorient.inputs.outputtype = 'NIFTI_GZ'
            preproc.connect(brain_mask_deoblique, 'out_file',
                            brain_mask_reorient, 'in_file')


            anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                            name='anat_skullstrip_orig_vol')
            anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
            anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(anat_reorient, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_a')

            preproc.connect(brain_mask_reorient, 'out_file',
                            anat_skullstrip_orig_vol, 'in_file_b')

            preproc.connect(brain_mask_reorient, 'out_file',
                            outputnode, 'brain_mask')

            preproc.connect(anat_skullstrip_orig_vol, 'out_file',
                            outputnode, 'brain')

        elif method == 'unet':
            """
            UNet
            options (following numbers are default):
            input_slice: 3
            conv_block: 5
            kernel_root: 16
            rescale_dim: 256
            """
            # TODO: add options to pipeline_config
            unet_check_for_s3 = create_check_for_s3_node('unet', config.unet_model)
            unet_mask = pe.Node(util.Function(input_names=['model_path', 'cimg_in'], 
                                              output_names=['out_path'],
                                              function=predict_volumes),                        
                                name='unet_mask')
            
            preproc.connect(unet_check_for_s3, 'local_path', unet_mask, 'model_path')
            preproc.connect(anat_reorient, 'out_file', unet_mask, 'cimg_in')

            """
            Revised mask with ANTs
            """
            # fslmaths <whole head> -mul <mask> brain.nii.gz
            unet_masked_brain = pe.Node(interface=fsl.MultiImageMaths(), name='unet_masked_brain')
            unet_masked_brain.inputs.op_string = "-mul %s"
            preproc.connect(anat_reorient, 'out_file', unet_masked_brain, 'in_file')
            preproc.connect(unet_mask, 'out_path', unet_masked_brain, 'operand_files')

            # flirt -v -dof 6 -in brain.nii.gz -ref NMT_SS_0.5mm.nii.gz -o brain_rot2atl -omat brain_rot2atl.mat -interp sinc
            # TODO: antsRegistration -z 0 -d 3 -r [NMT_SS_0.5mm.nii.gz,brain.nii.gz,0] -o [transform,brain_rot2atl.nii.gz,brain_inv_rot2atl.nii.gz] -t Rigid[0.1] -m MI[NMT_SS_0.5mm.nii.gz,brain.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-08,10] -s 3.0x2.0x1.0x0.0 -f 8x4x2x1 -u 1 -t Affine[0.1] -m MI[NMT_SS_0.5mm.nii.gz,brain.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-08,10] -s 3.0x2.0x1.0x0.0 -f 8x4x2x1 -u 1
            native_brain_to_template_brain = pe.Node(interface=fsl.FLIRT(), name='native_brain_to_template_brain')
            native_brain_to_template_brain.inputs.dof = 6
            native_brain_to_template_brain.inputs.interp = 'sinc'
            preproc.connect(unet_masked_brain, 'out_file', native_brain_to_template_brain, 'in_file')
            preproc.connect(inputnode, 'template_brain_only_for_anat', native_brain_to_template_brain, 'reference')
            
            # flirt -in head.nii.gz -ref NMT_0.5mm.nii.gz -o head_rot2atl -applyxfm -init brain_rot2atl.mat
            # TODO: antsApplyTransforms -d 3 -i head.nii.gz -r NMT_0.5mm.nii.gz -n Linear -o head_rot2atl.nii.gz -v -t transform1Rigid.mat -t transform2Affine.mat -t transform0DerivedInitialMovingTranslation.mat 
            native_head_to_template_head = pe.Node(interface=fsl.FLIRT(), name='native_head_to_template_head')
            native_head_to_template_head.inputs.apply_xfm = True
            preproc.connect(anat_reorient, 'out_file', native_head_to_template_head, 'in_file')
            preproc.connect(native_brain_to_template_brain, 'out_matrix_file', native_head_to_template_head, 'in_matrix_file')
            preproc.connect(inputnode, 'template_skull_for_anat', native_head_to_template_head, 'reference')

            # fslmaths NMT_SS_0.5mm.nii.gz -bin templateMask.nii.gz
            template_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(), name='template_brain_mask')
            template_brain_mask.inputs.args = '-bin'
            preproc.connect(inputnode, 'template_brain_only_for_anat', template_brain_mask, 'in_file')

            # ANTS 3 -m  CC[head_rot2atl.nii.gz,NMT_0.5mm.nii.gz,1,5] -t SyN[0.25] -r Gauss[3,0] -o atl2T1rot -i 60x50x20 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
            ants_template_head_to_template = pe.Node(interface=ants.Registration(), name='template_head_to_template')
            ants_template_head_to_template.inputs.metric = ['CC']
            ants_template_head_to_template.inputs.metric_weight = [1,5]
            ants_template_head_to_template.inputs.transforms = ['SyN']
            ants_template_head_to_template.inputs.transform_parameters = [(0.25,)]
            ants_template_head_to_template.inputs.interpolation = 'NearestNeighbor'
            ants_template_head_to_template.inputs.number_of_iterations = [[60,50,20]] 
            ants_template_head_to_template.inputs.smoothing_sigmas = [[0.6,0.2,0.0]]
            ants_template_head_to_template.inputs.shrink_factors = [[4,2,1]] 
            ants_template_head_to_template.inputs.convergence_threshold = [1.e-8]
            preproc.connect(native_head_to_template_head, 'out_file', ants_template_head_to_template, 'fixed_image')
            preproc.connect(inputnode, 'template_skull_for_anat', ants_template_head_to_template, 'moving_image')
            # antsApplyTransforms -d 3 -i templateMask.nii.gz -t atl2T1rotWarp.nii.gz atl2T1rotAffine.txt -r brain_rot2atl.nii.gz -o brain_rot2atl_mask.nii.gz
            template_head_transform_to_template = pe.Node(interface=ants.ApplyTransforms(), name='template_head_transform_to_template')
            template_head_transform_to_template.inputs.dimension = 3
            preproc.connect(template_brain_mask, 'out_file', template_head_transform_to_template, 'input_image')
            preproc.connect(native_brain_to_template_brain, 'out_file', template_head_transform_to_template, 'reference_image')
            preproc.connect(ants_template_head_to_template, 'forward_transforms', template_head_transform_to_template, 'transforms')

            # TODO: replace convert_xfm and flirt with: 
            # antsApplyTransforms -d 3 -i brain_rot2atl_mask.nii.gz -r brain.nii.gz -n linear -o brain_mask.nii.gz -t [transform0DerivedInitialMovingTranslation.mat,1] -t [transform2Affine.mat,1] -t [transform1Rigid.mat,1] 
            # convert_xfm -omat brain_rot2native.mat -inverse brain_rot2atl.mat 
            invt = pe.Node(interface=fsl.ConvertXFM(), name='convert_xfm')
            invt.inputs.invert_xfm = True
            preproc.connect(native_brain_to_template_brain, 'out_matrix_file', invt, 'in_file')

            # flirt -in brain_rot2atl_mask.nii.gz -ref brain.nii.gz -o brain_mask.nii.gz -applyxfm -init brain_rot2native.mat
            template_brain_to_native_brain = pe.Node(interface=fsl.FLIRT(), name='template_brain_to_native_brain')
            template_brain_to_native_brain.inputs.apply_xfm = True
            preproc.connect(template_head_transform_to_template, 'output_image', template_brain_to_native_brain, 'in_file')
            preproc.connect(unet_masked_brain, 'out_file', template_brain_to_native_brain, 'reference')
            preproc.connect(invt, 'out_file', template_brain_to_native_brain, 'in_matrix_file')

            # fslmaths brain_mask.nii.gz -thr .5 -bin brain_mask_thr.nii.gz
            refined_mask = pe.Node(interface=fsl.Threshold(), name='refined_mask')
            refined_mask.inputs.thresh = 0.5
            refined_mask.inputs.args = '-bin'
            preproc.connect(template_brain_to_native_brain, 'out_file', refined_mask, 'in_file')

            # get a new brain with mask
            refined_brain = pe.Node(interface=fsl.MultiImageMaths(), name='refined_brain')
            refined_brain.inputs.op_string = "-mul %s"
            preproc.connect(anat_reorient, 'out_file', refined_brain, 'in_file')
            preproc.connect(refined_mask, 'out_file', refined_brain, 'operand_files')
            
            preproc.connect(refined_mask, 'out_file', outputnode, 'brain_mask')
            preproc.connect(refined_brain, 'out_file', outputnode, 'brain')

    return preproc

# TODO multi-thread doens't work - debug!
def reconstruct_surface(num_omp_threads):

    """
    Parameters
    ----------
    num_omp_threads : int
        number of openmp threads

    Returns
    -------
    surface_reconstruction : Workflow
        autorecon3 surface reconstruction workflow
    
    Notes
    -----
    wm_seg is not a necessary input for autorecon3 but make a fake connection 
    so that nipype can run autorecon2 and 3 sequentially
    """
    
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    from nipype.interfaces import freesurfer

    surface_reconstruction = pe.Workflow(name='surface_reconstruction')

    inputnode = pe.Node(util.IdentityInterface(fields=['subject_dir',
                                                        'wm_seg']), 
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['curv',
                                                        'pial',
                                                        'smoothwm',
                                                        'sphere',
                                                        'sulc',
                                                        'thickness',
                                                        'volume',
                                                        'white']),
                        name='outputspec')

    reconall3 = pe.Node(interface=freesurfer.ReconAll(),
                        name='anat_autorecon3')

    reconall3.inputs.directive = 'autorecon3'
    reconall3.inputs.openmp = num_omp_threads
    
    surface_reconstruction.connect(inputnode, 'subject_dir',
                                    reconall3, 'subjects_dir')

    surface_reconstruction.connect(reconall3, 'curv',
                                    outputnode, 'curv')

    surface_reconstruction.connect(reconall3, 'pial',
                                    outputnode, 'pial')

    surface_reconstruction.connect(reconall3, 'smoothwm',
                                    outputnode, 'smoothwm')

    surface_reconstruction.connect(reconall3, 'sphere',
                                    outputnode, 'sphere')

    surface_reconstruction.connect(reconall3, 'sulc',
                                    outputnode, 'sulc')

    surface_reconstruction.connect(reconall3, 'thickness',
                                    outputnode, 'thickness')

    surface_reconstruction.connect(reconall3, 'volume',
                                    outputnode, 'volume')

    surface_reconstruction.connect(reconall3, 'white',
                                    outputnode, 'white')

    return surface_reconstruction                           