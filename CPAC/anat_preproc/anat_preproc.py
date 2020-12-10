# -*- coding: utf-8 -*-
import os
from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces.fsl import utils as fsl_utils
from nipype.interfaces.fsl import preprocess as fsl_preproc
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.anat_preproc.ants import init_brain_extraction_wf
from nipype.interfaces import freesurfer
from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string, \
    fsl_aff_to_rigid, \
    mri_convert, \
    wb_command, \
    fslmaths_command
from CPAC.utils.datasource import create_check_for_s3_node
from CPAC.unet.function import predict_volumes

from CPAC.seg_preproc.utils import pick_tissue_from_labels_file


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


def acpc_alignment(skullstrip_tool='afni', config=None, acpc_target='whole-head', wf_name='acpc_align'):
               
    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_leaf', 
                                                       'brain_mask',
                                                       'template_brain_only_for_anat',
                                                       'template_skull_for_anat',
                                                       'template_brain_for_acpc',
                                                       'template_head_for_acpc']), 
                        name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['acpc_aligned_head',
                                                         'acpc_brain_mask']),
                         name='outputspec')

    robust_fov = pe.Node(interface=fsl_utils.RobustFOV(),
                            name='anat_acpc_1_robustfov')
    robust_fov.inputs.brainsize = config.acpc_brainsize
    robust_fov.inputs.out_transform = 'fov_xfm.mat'

    # align head-to-head to get acpc.mat (for human)
    if acpc_target == 'whole-head':
        preproc.connect(inputnode, 'anat_leaf', robust_fov, 'in_file')
    
    # align brain-to-brain to get acpc.mat (for monkey)
    if acpc_target == 'brain':
        initial_skullstrip =  skullstrip_anatomical(method=skullstrip_tool, config=config, 
                                                    wf_name="anat_acpc_0_pre")
        preproc.connect(inputnode, 'anat_leaf', 
                        initial_skullstrip, 'inputspec.anat_data')
        if skullstrip_tool == 'mask':
            preproc.connect(inputnode, 'brain_mask',
                            initial_skullstrip, 'inputspec.brain_mask') 
        elif skullstrip_tool == 'unet':
            preproc.connect(inputnode, 'template_brain_only_for_anat',
                            initial_skullstrip, 'inputspec.template_brain_only_for_anat')   
            preproc.connect(inputnode, 'template_skull_for_anat',
                            initial_skullstrip, 'inputspec.template_skull_for_anat') 
        preproc.connect(initial_skullstrip, 'outputspec.brain', 
                        robust_fov, 'in_file')

    convert_fov_xfm = pe.Node(interface=fsl_utils.ConvertXFM(),
                                name='anat_acpc_2_fov_convertxfm')
    convert_fov_xfm.inputs.invert_xfm = True

    preproc.connect(robust_fov, 'out_transform',
                    convert_fov_xfm, 'in_file')

    align = pe.Node(interface=fsl.FLIRT(),
                    name='anat_acpc_3_flirt')
    align.inputs.interp = 'spline'
    align.inputs.searchr_x = [30, 30]
    align.inputs.searchr_y = [30, 30]
    align.inputs.searchr_z = [30, 30]

    preproc.connect(robust_fov, 'out_roi', align, 'in_file')

    # align head-to-head to get acpc.mat (for human)
    if acpc_target == 'whole-head':
        preproc.connect(inputnode, 'template_head_for_acpc', align, 'reference')
    
    # align brain-to-brain to get acpc.mat (for monkey)
    if acpc_target=='brain':
        preproc.connect(inputnode, 'template_brain_for_acpc', align, 'reference')

    concat_xfm = pe.Node(interface=fsl_utils.ConvertXFM(),
                            name='anat_acpc_4_concatxfm')
    concat_xfm.inputs.concat_xfm = True

    preproc.connect(convert_fov_xfm, 'out_file', concat_xfm, 'in_file')
    preproc.connect(align, 'out_matrix_file', concat_xfm, 'in_file2')

    aff_to_rig_imports = ['import os', 'from numpy import *']
    aff_to_rig = pe.Node(util.Function(input_names=['in_xfm', 'out_name'],
                                        output_names=['out_mat'],
                                        function=fsl_aff_to_rigid,
                                        imports=aff_to_rig_imports),
                            name='anat_acpc_5_aff2rigid')
    aff_to_rig.inputs.out_name = 'acpc.mat'

    preproc.connect(concat_xfm, 'out_file', aff_to_rig, 'in_xfm')

    apply_xfm = pe.Node(interface=fsl.ApplyWarp(),
                        name='anat_acpc_6_applywarp')
    apply_xfm.inputs.interp = 'spline'
    apply_xfm.inputs.relwarp = True

    preproc.connect(inputnode, 'anat_leaf', apply_xfm, 'in_file')
    preproc.connect(inputnode, 'template_head_for_acpc', apply_xfm, 'ref_file')
    preproc.connect(aff_to_rig, 'out_mat', apply_xfm, 'premat')
    preproc.connect(apply_xfm, 'out_file', output_node, 'acpc_aligned_head')

    if skullstrip_tool == 'mask':
        apply_xfm_mask = pe.Node(interface=fsl.ApplyWarp(),
                            name='anat_mask_acpc_7_applywarp')
        apply_xfm_mask.inputs.interp = 'nn'
        apply_xfm_mask.inputs.relwarp = True

        preproc.connect(inputnode, 'brain_mask', apply_xfm_mask, 'in_file')
        preproc.connect(inputnode, 'template_brain_for_acpc', apply_xfm_mask, 'ref_file')
        preproc.connect(aff_to_rig, 'out_mat', apply_xfm_mask, 'premat')
        preproc.connect(apply_xfm_mask, 'out_file', output_node, 'acpc_brain_mask')


    return preproc


def skullstrip_anatomical(method='afni', config=None, wf_name='skullstrip_anatomical'):

    method = method.lower()

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_data', 
                                                       'brain_mask',
                                                       'rawavg',
                                                       'template_brain_only_for_anat',
                                                       'template_skull_for_anat',
                                                       'freesurfer_subject_dir',
                                                       'wmparc']), 
                         name='inputspec')
    outputnode = pe.Node(util.IdentityInterface(fields=['skullstrip',
                                                        'brain',
                                                        'brain_mask']),
                        name='outputspec')

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
        
        inputnode_afni.inputs.set(
                mask_vol=config.skullstrip_mask_vol,
                shrink_factor=config.skullstrip_shrink_factor,
                var_shrink_fac=config.skullstrip_var_shrink_fac,
                shrink_fac_bot_lim=config.skullstrip_shrink_factor_bot_lim,
                avoid_vent=config.skullstrip_avoid_vent,
                niter=config.skullstrip_n_iterations,
                pushout=config.skullstrip_pushout,
                touchup=config.skullstrip_touchup,
                fill_hole=config.skullstrip_fill_hole,
                avoid_eyes=config.skullstrip_avoid_eyes,
                use_edge=config.skullstrip_use_edge,
                exp_frac=config.skullstrip_exp_frac,
                smooth_final=config.skullstrip_smooth_final,
                push_to_edge=config.skullstrip_push_to_edge,
                use_skull=config.skullstrip_use_skull,
                perc_int=config.skullstrip_perc_int,
                max_inter_iter=config.skullstrip_max_inter_iter,
                blur_fwhm=config.skullstrip_blur_fwhm,
                fac=config.skullstrip_fac,
                monkey=config.skullstrip_monkey,
            )

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

        preproc.connect(inputnode, 'anat_data',
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

        preproc.connect(inputnode, 'anat_data',
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
        
        inputnode_bet.inputs.set(
                frac=config.bet_frac,
                mask_boolean=config.bet_mask_boolean,
                mesh_boolean=config.bet_mesh_boolean,
                outline=config.bet_outline,
                padding=config.bet_padding,
                radius=config.bet_radius,
                reduce_bias=config.bet_reduce_bias,
                remove_eyes=config.bet_remove_eyes,
                robust=config.bet_robust,
                skull=config.bet_skull,
                surfaces=config.bet_surfaces,
                threshold=config.bet_threshold,
                vertical_gradient=config.bet_vertical_gradient,
            )

        preproc.connect(inputnode, 'anat_data',
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

        preproc.connect(inputnode, 'anat_data',
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

        preproc.connect(inputnode, 'anat_data',
                        anat_skullstrip_ants, 'inputnode.in_files')

        preproc.connect(anat_skullstrip_ants, 'copy_xform.out_file',
                        outputnode, 'skullstrip')

        preproc.connect(anat_skullstrip_ants, 'copy_xform.out_file',
                        outputnode, 'brain')

        preproc.connect(anat_skullstrip_ants, 'atropos_wf.copy_xform.out_mask',
                        outputnode, 'brain_mask')

    elif method == 'freesurfer':

        # register FS brain mask to native space
        fs_brain_mask_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                        name='fs_brain_mask_to_native')
        fs_brain_mask_to_native.inputs.reg_header = True
        preproc.connect(inputnode, 'brain_mask',
                        fs_brain_mask_to_native, 'source_file')
        preproc.connect(inputnode, 'rawavg',
                        fs_brain_mask_to_native, 'target_file')
        preproc.connect(inputnode, 'freesurfer_subject_dir',
                        fs_brain_mask_to_native, 'subjects_dir')

        # convert brain mask file from .mgz to .nii.gz
        fs_brain_mask_to_nifti = pe.Node(util.Function(input_names=['in_file'], 
                                            output_names=['out_file'],
                                            function=mri_convert),                        
                                name='fs_brainmask_to_nifti')                                
        preproc.connect(fs_brain_mask_to_native, 'transformed_file',
                        fs_brain_mask_to_nifti, 'in_file')

        # binarize the brain mask                
        binarize_fs_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                        name='binarize_fs_brainmask')
        binarize_fs_brain_mask.inputs.args = '-bin'
        preproc.connect(fs_brain_mask_to_nifti, 'out_file',
                        binarize_fs_brain_mask, 'in_file')

        # fill holes
        fill_fs_brain_mask = pe.Node(interface=afni.MaskTool(),
                            name='fill_fs_brainmask')
        fill_fs_brain_mask.inputs.fill_holes = True
        fill_fs_brain_mask.inputs.outputtype = 'NIFTI_GZ'
        preproc.connect(binarize_fs_brain_mask, 'out_file',
                    fill_fs_brain_mask, 'in_file')
        preproc.connect(fill_fs_brain_mask, 'out_file',
                    outputnode, 'brain_mask')

        # apply mask
        fs_brain = pe.Node(interface=fsl.ApplyMask(),
                        name='anat_skullstrip')
        preproc.connect(inputnode, 'anat_data',
                        fs_brain, 'in_file')
        preproc.connect(fill_fs_brain_mask, 'out_file',
                        fs_brain, 'mask_file')
        preproc.connect(fs_brain, 'out_file',
                        outputnode, 'brain')
    
    elif method == 'freesurfer-abcd':
        ### ABCD harmonization - anatomical brain mask generation ###
        # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh#L151-L159
        
        wmparc_to_nifti = pe.Node(util.Function(input_names=['in_file','reslice_like','args'],
                                            output_names=['out_file'],
                                            function=mri_convert),
                                name='wmparc_to_nifti')
        wmparc_to_nifti.inputs.args = '-rt nearest'

        preproc.connect(inputnode, 'wmparc',
                        wmparc_to_nifti, 'in_file')
        preproc.connect(inputnode, 'anat_data',
                        wmparc_to_nifti, 'reslice_like')

        binary_mask = pe.Node(interface=fsl.maths.MathsCommand(), 
                                name='binarize_wmparc')
        binary_mask.inputs.args = '-bin -dilD -dilD -dilD -ero -ero'

        preproc.connect(wmparc_to_nifti, 'out_file',
                        binary_mask, 'in_file')

        wb_command_fill_holes = pe.Node(util.Function(input_names=['in_file'],
                                            output_names=['out_file'],
                                            function=wb_command),
                                name='wb_command_fill_holes')

        preproc.connect(binary_mask, 'out_file',
                        wb_command_fill_holes, 'in_file')

        binary_mask2 = pe.Node(interface=fsl.maths.MathsCommand(),
                                name='binarize_wmparc2')
        binary_mask2.inputs.args = '-bin'

        preproc.connect(wb_command_fill_holes, 'out_file',
                        binary_mask2, 'in_file')

        brain_mask_to_t1_restore = pe.Node(interface=fsl.ApplyWarp(),
                    name='brain_mask_to_t1_restore')
        brain_mask_to_t1_restore.inputs.interp = 'nn'
        brain_mask_to_t1_restore.inputs.premat = config.identityMatrix

        preproc.connect(binary_mask2, 'out_file',
                        brain_mask_to_t1_restore, 'in_file')

        preproc.connect(inputnode, 'anat_data',
                        brain_mask_to_t1_restore, 'ref_file')

        preproc.connect(brain_mask_to_t1_restore, 'out_file',
                        outputnode, 'brain_mask')

        fs_brain = pe.Node(interface=fsl.MultiImageMaths(), name='fs_brain')
        fs_brain.inputs.op_string = "-mul %s"

        preproc.connect(brain_mask_to_t1_restore, 'out_file', 
                        fs_brain, 'in_file')

        preproc.connect(inputnode, 'anat_data', 
                        fs_brain, 'operand_files')

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

        preproc.connect(inputnode, 'anat_data',
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
        preproc.connect(inputnode, 'anat_data', unet_mask, 'cimg_in')

        """
        Revised mask with ANTs
        """
        # fslmaths <whole head> -mul <mask> brain.nii.gz
        unet_masked_brain = pe.Node(interface=fsl.MultiImageMaths(), name='unet_masked_brain')
        unet_masked_brain.inputs.op_string = "-mul %s"
        preproc.connect(inputnode, 'anat_data', unet_masked_brain, 'in_file')
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
        preproc.connect(inputnode, 'anat_data', native_head_to_template_head, 'in_file')
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
        preproc.connect(inputnode, 'anat_data', refined_brain, 'in_file')
        preproc.connect(refined_mask, 'out_file', refined_brain, 'operand_files')
        
        preproc.connect(refined_mask, 'out_file', outputnode, 'brain_mask')
        preproc.connect(refined_brain, 'out_file', outputnode, 'brain')

    return preproc


# TODO function node for fnirt-based brain extraction
def fnirt_based_brain_extraction(wf_name='fnirt_based_brain_extraction', config=None):

    ### ABCD Harmonization - FNIRT-based brain extraction ###
    # Ref: https://github.com/Washington-University/HCPpipelines/blob/master/PreFreeSurfer/scripts/BrainExtraction_FNIRTbased.sh

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_data',
                                                       'ref_mask_2mm',
                                                       'template_skull_for_anat_2mm',
                                                       'template_brain_mask_for_anat']), 
                         name='inputspec')
    outputnode = pe.Node(util.IdentityInterface(fields=['anat_brain',
                                                        'anat_brain_mask']),
                        name='outputspec')

    # Register to 2mm reference image (linear then non-linear)
    # linear registration to 2mm reference
    # flirt -interp spline -dof 12 -in "$Input" -ref "$Reference2mm" -omat "$WD"/roughlin.mat -out "$WD"/"$BaseName"_to_MNI_roughlin.nii.gz -nosearch
    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_reg')
    linear_reg.inputs.dof = 12
    linear_reg.inputs.interp = 'spline'
    linear_reg.inputs.no_search = True

    preproc.connect(inputnode, 'anat_data',
                    linear_reg, 'in_file')
    
    preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                    linear_reg, 'reference')

    # non-linear registration to 2mm reference
    # fnirt --in="$Input" --ref="$Reference2mm" --aff="$WD"/roughlin.mat --refmask="$Reference2mmMask" \
    # --fout="$WD"/str2standard.nii.gz --jout="$WD"/NonlinearRegJacobians.nii.gz \
    # --refout="$WD"/IntensityModulatedT1.nii.gz --iout="$WD"/"$BaseName"_to_MNI_nonlin.nii.gz \
    # --logout="$WD"/NonlinearReg.txt --intout="$WD"/NonlinearIntensities.nii.gz \
    # --cout="$WD"/NonlinearReg.nii.gz --config="$FNIRTConfig"
    non_linear_reg = pe.Node(interface=fsl.FNIRT(),
                         name='non_linear_reg')
    non_linear_reg.inputs.field_file = True # --fout
    non_linear_reg.inputs.jacobian_file = True # --jout
    non_linear_reg.inputs.modulatedref_file = True # --refout
    # non_linear_reg.inputs.warped_file = True # --iout
    # non_linear_reg.inputs.log_file = True # --logout
    non_linear_reg.inputs.out_intensitymap_file = True # --intout
    non_linear_reg.inputs.fieldcoeff_file = True # --cout

    preproc.connect(inputnode, 'anat_data',
                    non_linear_reg, 'in_file')

    preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                    non_linear_reg, 'ref_file')

    preproc.connect(linear_reg, 'out_matrix_file',
                    non_linear_reg, 'affine_file')

    preproc.connect(inputnode, 'ref_mask_2mm',
                    non_linear_reg, 'refmask_file')

    # # Overwrite the image output from FNIRT with a spline interpolated highres version
    # creating spline interpolated hires version
    # applywarp --rel --interp=spline --in="$Input" --ref="$Reference" -w "$WD"/str2standard.nii.gz --out="$WD"/"$BaseName"_to_MNI_nonlin.nii.gz
    apply_warp = pe.Node(interface=fsl.ApplyWarp(),
                        name='apply_warp')
    apply_warp.inputs.interp = 'spline'
    apply_warp.inputs.premat = config.identityMatrix

    preproc.connect(inputnode, 'anat_data',
                    apply_warp, 'in_file')

    preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                    apply_warp, 'ref_file')

    preproc.connect(non_linear_reg, 'field_file',
                    apply_warp, 'field_file')

    # # Invert warp and transform dilated brain mask back into native space, and use it to mask input image
    # # Input and reference spaces are the same, using 2mm reference to save time
    # computing inverse warp
    # invwarp --ref="$Reference2mm" -w "$WD"/str2standard.nii.gz -o "$WD"/standard2str.nii.gz
    inverse_warp = pe.Node(interface=fsl.InvWarp(), name='inverse_warp')
    inverse_warp.inputs.output_type = 'NIFTI_GZ'

    preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                    inverse_warp, 'reference')

    preproc.connect(non_linear_reg, 'field_file',
                    inverse_warp, 'warp')

    # applying inverse warp
    # applywarp --rel --interp=nn --in="$ReferenceMask" --ref="$Input" -w "$WD"/standard2str.nii.gz -o "$OutputBrainMask"
    apply_inv_warp = pe.Node(interface=fsl.ApplyWarp(),
                        name='apply_inv_warp')
    apply_inv_warp.inputs.interp = 'nn'
    apply_inv_warp.inputs.relwarp = True

    preproc.connect(inputnode, 'template_brain_mask_for_anat',
                    apply_inv_warp, 'in_file')

    preproc.connect(inputnode, 'anat_data',
                    apply_inv_warp, 'ref_file')

    preproc.connect(inverse_warp, 'inverse_warp',
                    apply_inv_warp, 'field_file')

    preproc.connect(apply_inv_warp, 'out_file',
                    outputnode, 'anat_brain_mask')
    
    # creating mask
    # fslmaths "$Input" -mas "$OutputBrainMask" "$OutputBrainExtractedImage"
    apply_mask = pe.Node(interface=fsl.MultiImageMaths(),
                            name='apply_mask')
    apply_mask.inputs.op_string = '-mas %s'

    preproc.connect(inputnode, 'anat_data', 
                    apply_mask, 'in_file')

    preproc.connect(apply_inv_warp, 'out_file',
                    apply_mask, 'operand_files')

    preproc.connect(apply_mask, 'out_file',
                    outputnode, 'anat_brain')
    
    return preproc


# TODO function node for FAST bias field correction
def fast_bias_field_correction(wf_name='fast_bias_field_correction', config=None):

    ### ABCD Harmonization - FAST bias field correction ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/92913242419d492aee733a45d454ea319fbaac35/PreFreeSurfer/PreFreeSurferPipeline.sh#L688-L694

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_data',
                                                       'anat_brain',
                                                       'anat_brain_mask']), 
                         name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['anat_restore',
                                                        'anat_brain_restore']),
                        name='outputspec')

    # fast -b -B -o ${T1wFolder}/T1w_fast -t 1 ${T1wFolder}/T1w_acpc_dc_brain.nii.gz
    fast_bias_field_correction = pe.Node(interface=fsl.FAST(),
                            name='fast_bias_field_correction')
    fast_bias_field_correction.inputs.img_type = 1
    fast_bias_field_correction.inputs.output_biasfield = True
    fast_bias_field_correction.inputs.output_biascorrected = True

    preproc.connect(inputnode, 'anat_brain',
                    fast_bias_field_correction, 'in_files')

    preproc.connect(fast_bias_field_correction, 'restored_image',
                    outputnode, 'anat_brain_restore')

    # FAST does not output a non-brain extracted image so create an inverse mask, 
    # apply it to T1w_acpc_dc.nii.gz, insert the T1w_fast_restore to the skull of 
    # the T1w_acpc_dc.nii.gz and use that for the T1w_acpc_dc_restore head

    # fslmaths ${T1wFolder}/T1w_acpc_brain_mask.nii.gz -mul -1 -add 1 ${T1wFolder}/T1w_acpc_inverse_brain_mask.nii.gz
    inverse_brain_mask = pe.Node(interface=fsl.ImageMaths(),
                                    name='inverse_brain_mask')
    inverse_brain_mask.inputs.op_string = '-mul -1 -add 1'

    preproc.connect(inputnode, 'anat_brain_mask',
                    inverse_brain_mask, 'in_file')

    # fslmaths ${T1wFolder}/T1w_acpc_dc.nii.gz -mul ${T1wFolder}/T1w_acpc_inverse_brain_mask.nii.gz ${T1wFolder}/T1w_acpc_dc_skull.nii.gz
    # TODO what's this step doing?
    apply_mask = pe.Node(interface=fsl.MultiImageMaths(),
                            name='apply_mask')
    apply_mask.inputs.op_string = '-mul %s'

    preproc.connect(inputnode, 'anat_data',
                    apply_mask, 'in_file')

    preproc.connect(inverse_brain_mask, 'out_file',
                    apply_mask, 'operand_files')

    # fslmaths ${T1wFolder}/T1w_fast_restore.nii.gz -add ${T1wFolder}/T1w_acpc_dc_skull.nii.gz ${T1wFolder}/${T1wImage}_acpc_dc_restore
    # TODO what's this step doing?
    anat_restore = pe.Node(interface=fsl.MultiImageMaths(),
                            name='anat_restore')
    anat_restore.inputs.op_string = '-add %s'

    preproc.connect(fast_bias_field_correction, 'restored_image',
                    anat_restore, 'in_file')

    preproc.connect(apply_mask, 'out_file',
                    anat_restore, 'operand_files')

    preproc.connect(anat_restore, 'out_file',
                    outputnode, 'anat_restore')

    return preproc



def create_anat_preproc(method='afni', already_skullstripped=False,
                        config=None, acpc_target='whole-head', wf_name='anat_preproc', sub_dir=None):
    """The main purpose of this workflow is to process T1 scans. Raw mprage
    file is deobliqued, reoriented into RPI and skullstripped. Also, a whole
    brain only mask is generated from the skull stripped image for later use
    in registration.

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
            Path to skull stripped RPI oriented mprage file with normalized
            intensities.

        outputspec.brain : string
            Path to skull stripped RPI brain image with original intensity
            values and not normalized or scaled.

    Order of commands:
    - Deobliqing the scans. ::
        3drefit -deoblique mprage.nii.gz

    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior
      Inferior-to-Superior  (RPI) orientation ::
        3dresample -orient RPI
                   -prefix mprage_RPI.nii.gz
                   -inset mprage.nii.gz

    - Skull-Stripping the image ::
        Using AFNI ::
            3dSkullStrip -input mprage_RPI.nii.gz
                         -o_ply mprage_RPI_3dT.nii.gz
        or using BET ::
            bet mprage_RPI.nii.gz

    - The skull-stripping step modifies the intensity values. To get back the
      original intensity values, we do an element wise product of RPI data
      with step function of skull-stripped data ::
        3dcalc -a mprage_RPI.nii.gz
               -b mprage_RPI_3dT.nii.gz
               -expr 'a*step(b)'
               -prefix mprage_RPI_3dc.nii.gz

    .. exec::
        import yaml
        from urllib.request import urlopen
        from CPAC.anat_preproc import create_anat_preproc
        from CPAC.utils.configuration import Configuration

        wf = create_anat_preproc(
            config=Configuration(yaml.safe_load(urlopen(
                'https://raw.githubusercontent.com/FCP-INDI/C-PAC/develop/'
                'dev/docker_data/default_pipeline.yml'
            )))
        )
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/anatpreproc.dot'
        )

    High Level Workflow Graph:

    .. image:: ../../images/generated/anatpreproc.png
       :width: 500

    Detailed Workflow Graph:

    .. image:: ../../images/generated/anatpreproc_detailed.png
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
                                                       'template_brain_only_for_acpc',
                                                       'template_skull_for_acpc',
                                                       'ref_mask_2mm',
                                                       'template_skull_for_anat_2mm',
                                                       'template_brain_mask_for_anat',
                                                       'template_cmass']), 
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['refit',
                                                        'reorient',
                                                        'skullstrip',
                                                        'brain',
                                                        'brain_mask',
                                                        'anat_skull_leaf',
                                                        'freesurfer_subject_dir',
                                                        'center_of_mass',
                                                        'wm_mask',
                                                        'gm_mask',
                                                        'csf_mask',
                                                        'curv',
                                                        'pial',
                                                        'smoothwm',
                                                        'sphere',
                                                        'sulc',
                                                        'thickness',
                                                        'volume',
                                                        'white']),
                        name='outputspec')

    anat_deoblique = pe.Node(interface=afni.Refit(),
                             name='anat_deoblique')
    anat_deoblique.inputs.deoblique = True

    preproc.connect(inputnode, 'anat', anat_deoblique, 'in_file')
    preproc.connect(anat_deoblique, 'out_file', outputnode, 'refit')

    # Anatomical reorientation
    anat_reorient = pe.Node(interface=afni.Resample(),
                            name='anat_reorient')
    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(anat_deoblique, 'out_file', anat_reorient, 'in_file')
    preproc.connect(anat_reorient, 'out_file', outputnode, 'reorient')

    anat_leaf = pe.Node(util.IdentityInterface(fields=['anat_data']),
                        name='anat_leaf')

    if config.non_local_means_filtering and config.n4_bias_field_correction and config.acpc_run_preprocessing=='before':
        denoise_before_acpc = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise_before_acpc')
        # align ABCD pipeline
        denoise_before_acpc.inputs.dimension = 3
        denoise_before_acpc.inputs.noise_model = 'Rician'
        preproc.connect(anat_reorient, 'out_file', denoise_before_acpc, 'input_image')

        n4_before_acpc = pe.Node(interface = ants.N4BiasFieldCorrection(), name='anat_n4_before_acpc')
        n4_before_acpc.inputs.dimension = 3
        preproc.connect(denoise_before_acpc, 'output_image', n4_before_acpc, 'input_image')
    
    if not config.acpc_align:
        preproc.connect(anat_reorient, 'out_file', anat_leaf, 'anat_data')

    # ACPC alignment     
    if config.acpc_align:
        acpc_align = acpc_alignment(skullstrip_tool=method, config=config, acpc_target=acpc_target, wf_name='acpc_align')

        if config.acpc_run_preprocessing=='before':
            preproc.connect(n4_before_acpc, 'output_image', acpc_align, 'inputspec.anat_leaf')
        else:
            preproc.connect(anat_reorient, 'out_file', acpc_align, 'inputspec.anat_leaf')

        preproc.connect(inputnode, 'brain_mask', acpc_align, 'inputspec.brain_mask')
        preproc.connect(inputnode, 'template_brain_only_for_acpc', acpc_align, 'inputspec.template_brain_for_acpc')
        preproc.connect(inputnode, 'template_skull_for_acpc', acpc_align, 'inputspec.template_head_for_acpc')
        preproc.connect(acpc_align, 'outputspec.acpc_aligned_head', anat_leaf, 'anat_data')
        if method == 'unet':
            preproc.connect(inputnode, 'template_brain_only_for_anat', acpc_align, 'inputspec.template_brain_only_for_anat')
            preproc.connect(inputnode, 'template_skull_for_anat', acpc_align, 'inputspec.template_skull_for_anat')
    
    # Disable non_local_means_filtering and n4_bias_field_correction when run niworkflows-ants
    if method == 'niworkflows-ants':
        config.non_local_means_filtering = False
        config.n4_bias_field_correction = False

    if config.acpc_run_preprocessing=='after':

        if config.non_local_means_filtering and config.n4_bias_field_correction:
            denoise = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise')
            preproc.connect(anat_leaf, 'anat_data', denoise, 'input_image')

            n4 = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, shrink_factor=2, copy_header=True),
                name='anat_n4')
            preproc.connect(denoise, 'output_image', n4, 'input_image')

        elif config.non_local_means_filtering and not config.n4_bias_field_correction:
            denoise = pe.Node(interface = ants.DenoiseImage(), name = 'anat_denoise')
            preproc.connect(anat_leaf, 'anat_data', denoise, 'input_image')

        elif not config.non_local_means_filtering and config.n4_bias_field_correction:
            n4 = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, shrink_factor=2, copy_header=True),
                name='anat_n4')
            preproc.connect(anat_leaf, 'anat_data', n4, 'input_image')

    anat_leaf2 = pe.Node(util.IdentityInterface(fields=['anat_data']),
                         name='anat_leaf2')

    if config.acpc_run_preprocessing=='after':
        if config.n4_bias_field_correction:
            preproc.connect(n4, 'output_image', anat_leaf2, 'anat_data')
        elif config.non_local_means_filtering and not config.n4_bias_field_correction:
            preproc.connect(denoise, 'output_image', anat_leaf2, 'anat_data')
        else:
            preproc.connect(anat_leaf, 'anat_data', anat_leaf2, 'anat_data')
    else:
        preproc.connect(anat_leaf, 'anat_data', anat_leaf2, 'anat_data')


    if already_skullstripped:
        anat_skullstrip = pe.Node(interface=util.IdentityInterface(fields=['out_file']),
                                    name='anat_skullstrip')

        preproc.connect(anat_leaf2, 'anat_data',
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

        if method == 'freesurfer':
            ### recon-all -all ###
            reconall = pe.Node(interface=freesurfer.ReconAll(),
                name='anat_freesurfer')
            
            num_strat = len(config.skullstrip_option)-1 # the number of strategy that FreeSurfer will be

            freesurfer_subject_dir = os.path.join(sub_dir, f'anat_preproc_freesurfer_{num_strat}', 'anat_freesurfer')
            
            # create the node dir
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
            if not os.path.exists(os.path.join(sub_dir, f'anat_preproc_freesurfer_{num_strat}')):
                os.mkdir(os.path.join(sub_dir, f'anat_preproc_freesurfer_{num_strat}'))
            if not os.path.exists(freesurfer_subject_dir):
                os.mkdir(freesurfer_subject_dir)

            reconall.inputs.directive = 'all'
            reconall.inputs.subjects_dir = freesurfer_subject_dir
            reconall.inputs.openmp = config.num_omp_threads

            preproc.connect(anat_leaf2, 'anat_data',
                            reconall, 'T1_files')

            preproc.connect(reconall, 'subjects_dir',
                            outputnode, 'freesurfer_subject_dir')

            ### segmentation output ###
            # register FS segmentations (aseg.mgz) to native space
            fs_aseg_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                            name='fs_aseg_to_native')
            fs_aseg_to_native.inputs.reg_header = True
            fs_aseg_to_native.inputs.interp = 'nearest'
            fs_aseg_to_native.inputs.subjects_dir = freesurfer_subject_dir

            preproc.connect(reconall, 'aseg',
                            fs_aseg_to_native, 'source_file')
            preproc.connect(reconall, 'rawavg',
                            fs_aseg_to_native, 'target_file')

            # convert registered FS segmentations from .mgz to .nii.gz
            fs_aseg_to_nifti = pe.Node(util.Function(input_names=['in_file'], 
                                                output_names=['out_file'],
                                                function=mri_convert),                        
                                                name='fs_aseg_to_nifti')
            fs_aseg_to_nifti.inputs.args = '-rt nearest'

            preproc.connect(fs_aseg_to_native, 'transformed_file',
                            fs_aseg_to_nifti, 'in_file')

            pick_tissue = pe.Node(util.Function(input_names=['multiatlas_Labels'], 
                                                output_names=['csf_mask', 'gm_mask', 'wm_mask'],
                                                function=pick_tissue_from_labels_file), 
                                                name=f'anat_preproc_freesurfer_tissue_mask')

            pick_tissue.inputs.include_ventricles = True

            preproc.connect(fs_aseg_to_nifti, 'out_file',
                            pick_tissue, 'multiatlas_Labels')
            
            preproc.connect(pick_tissue, 'wm_mask',
                            outputnode, 'wm_mask')

            preproc.connect(pick_tissue, 'gm_mask',
                            outputnode, 'gm_mask')

            preproc.connect(pick_tissue, 'csf_mask',
                            outputnode, 'csf_mask')

            preproc.connect(reconall, 'curv',
                            outputnode, 'curv')

            preproc.connect(reconall, 'pial',
                            outputnode, 'pial')

            preproc.connect(reconall, 'smoothwm',
                            outputnode, 'smoothwm')

            preproc.connect(reconall, 'sphere',
                            outputnode, 'sphere')

            preproc.connect(reconall, 'sulc',
                            outputnode, 'sulc')

            preproc.connect(reconall, 'thickness',
                            outputnode, 'thickness')

            preproc.connect(reconall, 'volume',
                            outputnode, 'volume')

            preproc.connect(reconall, 'white',
                            outputnode, 'white')


        elif method == 'freesurfer-abcd':
            
            brain_extraction = fnirt_based_brain_extraction(config=config)

            preproc.connect(anat_leaf2, 'anat_data',
                            brain_extraction, 'inputspec.anat_data')

            preproc.connect(inputnode, 'ref_mask_2mm',
                            brain_extraction, 'inputspec.ref_mask_2mm')

            preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                            brain_extraction, 'inputspec.template_skull_for_anat_2mm')

            preproc.connect(inputnode, 'template_brain_mask_for_anat',
                            brain_extraction, 'inputspec.template_brain_mask_for_anat')

            fast_correction = fast_bias_field_correction(config=config)

            preproc.connect(anat_leaf2, 'anat_data',
                            fast_correction, 'inputspec.anat_data')

            preproc.connect(brain_extraction, 'outputspec.anat_brain',
                            fast_correction, 'inputspec.anat_brain')

            preproc.connect(brain_extraction, 'outputspec.anat_brain_mask',
                            fast_correction, 'inputspec.anat_brain_mask')

            ### ABCD Harmonization ###
            # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/92913242419d492aee733a45d454ea319fbaac35/FreeSurfer/FreeSurferPipeline.sh#L140-L144

            # flirt -interp spline -in "$T1wImage" -ref "$T1wImage" -applyisoxfm 1 -out "$T1wImageFile"_1mm.nii.gz
            resample_head_1mm = pe.Node(interface=fsl.FLIRT(),
                            name='resample_anat_head_1mm')
            resample_head_1mm.inputs.apply_isoxfm = 1

            preproc.connect(anat_leaf2, 'anat_data',
                            resample_head_1mm, 'in_file')
            
            preproc.connect(anat_leaf2, 'anat_data',
                            resample_head_1mm, 'reference')

            # applywarp --rel --interp=spline -i "$T1wImage" -r "$T1wImageFile"_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wImageFile"_1mm.nii.gz
            applywarp_head_to_head_1mm = pe.Node(interface=fsl.ApplyWarp(),
                        name='applywarp_head_to_head_1mm')
            applywarp_head_to_head_1mm.inputs.interp = 'spline'
            applywarp_head_to_head_1mm.inputs.premat = config.identityMatrix

            preproc.connect(anat_leaf2, 'anat_data',
                            applywarp_head_to_head_1mm, 'in_file')

            preproc.connect(resample_head_1mm, 'out_file',
                            applywarp_head_to_head_1mm, 'ref_file')

            # applywarp --rel --interp=nn -i "$T1wImageBrain" -r "$T1wImageFile"_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wImageBrainFile"_1mm.nii.gz
            applywarp_brain_to_head_1mm = pe.Node(interface=fsl.ApplyWarp(),
                        name='applywarp_brain_to_head_1mm')
            applywarp_brain_to_head_1mm.inputs.interp = 'nn'
            applywarp_brain_to_head_1mm.inputs.premat = config.identityMatrix

            preproc.connect(fast_correction, 'outputspec.anat_brain_restore',
                            applywarp_brain_to_head_1mm, 'in_file')

            preproc.connect(resample_head_1mm, 'out_file',
                            applywarp_brain_to_head_1mm, 'ref_file')

            # fslstats $T1wImageBrain -M
            average_brain = pe.Node(interface=fsl.ImageStats(),
                        name='average_brain')
            average_brain.inputs.op_string = '-M'
            average_brain.inputs.output_type = 'NIFTI_GZ'

            preproc.connect(fast_correction, 'outputspec.anat_brain_restore',
                            average_brain, 'in_file')

            # fslmaths "$T1wImageFile"_1mm.nii.gz -div $Mean -mul 150 -abs "$T1wImageFile"_1mm.nii.gz
            normalize_head = pe.Node(util.Function(input_names=['in_file', 'number', 'out_file_suffix'],
                                                    output_names=['out_file'],
                                                    function=fslmaths_command),
                                    name='normalize_head')
            normalize_head.inputs.out_file_suffix = '_norm'

            preproc.connect(applywarp_head_to_head_1mm, 'out_file', 
                            normalize_head, 'in_file')

            preproc.connect(average_brain, 'out_stat',
                            normalize_head, 'number')

            ### recon-all -all step ###
            reconall = pe.Node(interface=freesurfer.ReconAll(),
                name='anat_freesurfer')

            # TODO split recon-all -all into 4 steps as ABCD does
            # ABCD performs autorecon1 without brain extraction
            num_strat = len(config.skullstrip_option)-1 # the number of strategy that FreeSurfer will be

            freesurfer_subject_dir = os.path.join(sub_dir, f'anat_preproc_freesurfer_abcd_{num_strat}', 'anat_freesurfer')
            
            # create the node dir
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
            if not os.path.exists(os.path.join(sub_dir, f'anat_preproc_freesurfer_abcd_{num_strat}')):
                os.mkdir(os.path.join(sub_dir, f'anat_preproc_freesurfer_abcd_{num_strat}'))
            if not os.path.exists(freesurfer_subject_dir):
                os.mkdir(freesurfer_subject_dir)

            reconall.inputs.directive = 'all'
            reconall.inputs.subjects_dir = freesurfer_subject_dir
            reconall.inputs.openmp = config.num_omp_threads

            # preproc.connect(anat_leaf2, 'anat_data',
            #                 reconall, 'T1_files')

            preproc.connect(normalize_head, 'out_file',
                            reconall, 'T1_files')

            preproc.connect(reconall, 'subjects_dir',
                            outputnode, 'freesurfer_subject_dir')
            
            '''
            ### ABCD Harmonization - FreeSurfer brain mask refinement ###
            # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/92913242419d492aee733a45d454ea319fbaac35/FreeSurfer/FreeSurferPipeline.sh#L169-L172

            # mri_convert "$T1wImageBrainFile"_1mm.nii.gz "$SubjectDIR"/"$SubjectID"/mri/brainmask.mgz --conform
            t1_brain_to_mgz = pe.Node(util.Function(input_names=['in_file'],
                                                    output_names=['out_file'],
                                                    function=mri_convert),
                                        name='t1_brain_to_mgz')

            preproc.connect(applywarp_brain_to_head_1mm, 'out_file',
                            t1_brain_to_mgz, 'in_file')

            # mri_em_register -mask "$SubjectDIR"/"$SubjectID"/mri/brainmask.mgz "$SubjectDIR"/"$SubjectID"/mri/nu.mgz $FREESURFER_HOME/average/RB_all_2008-03-26.gca "$SubjectDIR"/"$SubjectID"/mri/transforms/talairach_with_skull.lta
            # TODO write a function node

            # preproc.connect(t1_brain_to_mgz, 'out_file',
            #                 mri_em_register, 'brain_mask')

            # mri_watershed -T1 -brain_atlas $FREESURFER_HOME/average/RB_all_withskull_2008-03-26.gca "$SubjectDIR"/"$SubjectID"/mri/transforms/talairach_with_skull.lta "$SubjectDIR"/"$SubjectID"/mri/T1.mgz "$SubjectDIR"/"$SubjectID"/mri/brainmask.auto.mgz 
            '''

        anat_skullstrip = skullstrip_anatomical(method=method, config=config,
                                                wf_name="{0}_skullstrip".format(wf_name))

        preproc.connect(anat_leaf2, 'anat_data',
                        anat_skullstrip, 'inputspec.anat_data')
        preproc.connect(inputnode, 'template_brain_only_for_anat',
                        anat_skullstrip, 'inputspec.template_brain_only_for_anat') 
        preproc.connect(inputnode, 'template_skull_for_anat', 
                        anat_skullstrip, 'inputspec.template_skull_for_anat')
            
        if method == 'mask' and config.acpc_align:
            preproc.connect(acpc_align, 'outputspec.acpc_brain_mask', 
                            anat_skullstrip, 'inputspec.brain_mask') 
        elif method == 'freesurfer':
            preproc.connect(reconall, 'brainmask',
                            anat_skullstrip, 'inputspec.brain_mask')

            preproc.connect(reconall, 'rawavg',
                            anat_skullstrip, 'inputspec.rawavg')

            preproc.connect(reconall, 'subjects_dir',
                            anat_skullstrip, 'inputspec.freesurfer_subject_dir')
        elif method == 'freesurfer-abcd':
            preproc.connect(reconall, 'wmparc',
                            anat_skullstrip, 'inputspec.wmparc')

            preproc.connect(reconall, 'subjects_dir',
                            anat_skullstrip, 'inputspec.freesurfer_subject_dir')                    
        else:
            preproc.connect(inputnode, 'brain_mask',
                            anat_skullstrip, 'inputspec.brain_mask')
        
        preproc.connect(anat_skullstrip, 'outputspec.brain_mask',
                        outputnode, 'brain_mask')
        preproc.connect(anat_skullstrip, 'outputspec.brain', 
                        outputnode, 'brain')

    preproc.connect(anat_leaf2, 'anat_data', outputnode, 'anat_skull_leaf')

    return preproc

