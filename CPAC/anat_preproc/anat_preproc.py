# -*- coding: utf-8 -*-
import os
from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces.fsl import utils as fsl_utils
from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.utility as util
from CPAC.anat_preproc.ants import init_brain_extraction_wf
from nipype.interfaces import freesurfer
from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string, \
    fsl_aff_to_rigid, \
    mri_convert, \
    wb_command, \
    fslmaths_command
from CPAC.unet.function import predict_volumes

from CPAC.seg_preproc.utils import pick_tissue_from_labels_file


def acpc_alignment(config=None, acpc_target='whole-head', mask=False,
                   wf_name='acpc_align'):
    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_leaf',
                                                       'anat_brain',
                                                       'brain_mask',
                                                       'template_brain_for_acpc',
                                                       'template_head_for_acpc']),
                        name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['acpc_aligned_head',
                                                         'acpc_brain_mask',
                                                         'from-affine_to-rigid_mode-image_xfm']),
                          name='outputspec')

    robust_fov = pe.Node(interface=fsl_utils.RobustFOV(),
                         name='anat_acpc_1_robustfov')
    robust_fov.inputs.brainsize = config.anatomical_preproc['acpc_alignment'][
        'brain_size']
    robust_fov.inputs.out_transform = 'fov_xfm.mat'

    # align head-to-head to get acpc.mat (for human)
    if acpc_target == 'whole-head':
        preproc.connect(inputnode, 'anat_leaf', robust_fov, 'in_file')

    # align brain-to-brain to get acpc.mat (for monkey)
    if acpc_target == 'brain':
        preproc.connect(inputnode, 'anat_brain', robust_fov, 'in_file')

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
        preproc.connect(inputnode, 'template_head_for_acpc', align,
                        'reference')

    # align brain-to-brain to get acpc.mat (for monkey)
    if acpc_target == 'brain':
        preproc.connect(inputnode, 'template_brain_for_acpc', align,
                        'reference')

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
    preproc.connect(aff_to_rig, 'out_mat', output_node, 'from-affine_to-rigid_mode-image_xfm')

    apply_xfm = pe.Node(interface=fsl.ApplyWarp(),
                        name='anat_acpc_6_applywarp')
    apply_xfm.inputs.interp = 'spline'
    apply_xfm.inputs.relwarp = True

    preproc.connect(inputnode, 'anat_leaf', apply_xfm, 'in_file')
    preproc.connect(inputnode, 'template_head_for_acpc', apply_xfm,
                    'ref_file')
    preproc.connect(aff_to_rig, 'out_mat', apply_xfm, 'premat')
    preproc.connect(apply_xfm, 'out_file', output_node, 'acpc_aligned_head')

    if mask:
        apply_xfm_mask = pe.Node(interface=fsl.ApplyWarp(),
                                 name='anat_mask_acpc_7_applywarp')
        apply_xfm_mask.inputs.interp = 'nn'
        apply_xfm_mask.inputs.relwarp = True

        preproc.connect(inputnode, 'brain_mask', apply_xfm_mask, 'in_file')
        preproc.connect(inputnode, 'template_brain_for_acpc', apply_xfm_mask,
                        'ref_file')
        preproc.connect(aff_to_rig, 'out_mat', apply_xfm_mask, 'premat')
        preproc.connect(apply_xfm_mask, 'out_file', output_node,
                        'acpc_brain_mask')

    return preproc


def afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt):
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
                                       'NN_smooth',
                                       'smooth_final',
                                       'avoid_eyes',
                                       'use_edge',
                                       'exp_frac',
                                       'push_to_edge',
                                       'use_skull',
                                       'perc_int',
                                       'max_inter_iter',
                                       'blur_fwhm',
                                       'fac',
                                       'monkey']),
        name=f'AFNI_options_{pipe_num}')

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
                                                         'NN_smooth',
                                                         'smooth_final',
                                                         'avoid_eyes',
                                                         'use_edge',
                                                         'exp_frac',
                                                         'push_to_edge',
                                                         'use_skull',
                                                         'perc_int',
                                                         'max_inter_iter',
                                                         'blur_fwhm',
                                                         'fac',
                                                         'monkey'],
                                            output_names=['expr'],
                                            function=create_3dskullstrip_arg_string),
                              name=f'anat_skullstrip_args_{pipe_num}')

    inputnode_afni.inputs.set(
        mask_vol=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['mask_vol'],
        shrink_factor=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['shrink_factor'],
        var_shrink_fac=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['var_shrink_fac'],
        shrink_fac_bot_lim=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['shrink_factor_bot_lim'],
        avoid_vent=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['avoid_vent'],
        niter=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['n_iterations'],
        pushout=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['pushout'],
        touchup=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['touchup'],
        fill_hole=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['fill_hole'],
        NN_smooth=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['NN_smooth'],
        smooth_final=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['smooth_final'],
        avoid_eyes=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['avoid_eyes'],
        use_edge=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['use_edge'],
        exp_frac=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['exp_frac'],
        push_to_edge=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['push_to_edge'],
        use_skull=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['use_skull'],
        perc_int=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['perc_int'],
        max_inter_iter=
        cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['max_inter_iter'],
        fac=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['fac'],
        blur_fwhm=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['blur_fwhm'],
        monkey=cfg.anatomical_preproc['brain_extraction'][
            'AFNI-3dSkullStrip']['monkey'],
    )

    wf.connect([
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
            ('NN_smooth', 'NN_smooth'),
            ('smooth_final', 'smooth_final'),
            ('push_to_edge', 'push_to_edge'),
            ('use_skull', 'use_skull'),
            ('perc_int', 'perc_int'),
            ('max_inter_iter', 'max_inter_iter'),
            ('blur_fwhm', 'blur_fwhm'),
            ('fac', 'fac'),
            ('monkey', 'monkey')
        ])
    ])

    anat_skullstrip = pe.Node(interface=afni.SkullStrip(),
                              name=f'anat_skullstrip_{pipe_num}')
    anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, anat_skullstrip, 'in_file')
    wf.connect(skullstrip_args, 'expr', anat_skullstrip, 'args')

    # Generate anatomical brain mask
    anat_brain_mask = pe.Node(interface=afni.Calc(),
                              name=f'anat_brain_mask_{pipe_num}')

    anat_brain_mask.inputs.expr = 'step(a)'
    anat_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(anat_skullstrip, 'out_file',
               anat_brain_mask, 'in_file_a')

    outputs = {
        'space-T1w_desc-brain_mask': (anat_brain_mask, 'out_file')
    }

    return (wf, outputs)


def fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt):
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
        name=f'BET_options_{pipe_num}')

    anat_skullstrip = pe.Node(
        interface=fsl.BET(), name=f'anat_BET_skullstrip_{pipe_num}')
    anat_skullstrip.inputs.output_type = 'NIFTI_GZ'

    inputnode_bet.inputs.set(
        frac=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['frac'],
        mask_boolean=
        cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['mask_boolean'],
        mesh_boolean=
        cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['mesh_boolean'],
        outline=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['outline'],
        padding=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['padding'],
        radius=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['radius'],
        reduce_bias=
        cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['reduce_bias'],
        remove_eyes=
        cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['remove_eyes'],
        robust=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['robust'],
        skull=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['skull'],
        surfaces=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['surfaces'],
        threshold=cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['threshold'],
        vertical_gradient=
        cfg.anatomical_preproc['brain_extraction'][
            'FSL-BET']['vertical_gradient'],
    )

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, anat_skullstrip, 'in_file')

    wf.connect([
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

    outputs = {
        'space-T1w_desc-brain_mask': (anat_skullstrip, 'mask_file')
    }

    return (wf, outputs)


def niworkflows_ants_brain_connector(wf, cfg, strat_pool, pipe_num, opt):
    # Skull-stripping using niworkflows-ants
    anat_skullstrip_ants = init_brain_extraction_wf(tpl_target_path=
                                                    cfg.anatomical_preproc[
                                                        'brain_extraction'][
                                                        'niworkflows-ants'][
                                                        'template_path'],
                                                    tpl_mask_path=
                                                    cfg.anatomical_preproc[
                                                        'brain_extraction'][
                                                        'niworkflows-ants'][
                                                        'mask_path'],
                                                    tpl_regmask_path=
                                                    cfg.anatomical_preproc[
                                                        'brain_extraction'][
                                                        'niworkflows-ants'][
                                                        'regmask_path'],
                                                    name='anat_skullstrip_ants')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, anat_skullstrip_ants, 'inputnode.in_files')

    outputs = {
        'space-T1w_desc-brain_mask':
            (anat_skullstrip_ants, 'atropos_wf.copy_xform.out_mask')
    }

    return (wf, outputs)


def unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt):
    """
    UNet
    options (following numbers are default):
    input_slice: 3
    conv_block: 5
    kernel_root: 16
    rescale_dim: 256
    """

    unet_mask = pe.Node(util.Function(input_names=['model_path', 'cimg_in'],
                                      output_names=['out_path'],
                                      function=predict_volumes),
                        name=f'unet_mask_{pipe_num}')

    node, out = strat_pool.get_data('unet_model')
    wf.connect(node, out, unet_mask, 'model_path')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, unet_mask, 'cimg_in')

    """
    Revised mask with ANTs
    """
    # fslmaths <whole head> -mul <mask> brain.nii.gz
    unet_masked_brain = pe.Node(interface=fsl.MultiImageMaths(),
                                name=f'unet_masked_brain_{pipe_num}')
    unet_masked_brain.inputs.op_string = "-mul %s"

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, unet_masked_brain, 'in_file')
    wf.connect(unet_mask, 'out_path', unet_masked_brain, 'operand_files')

    # flirt -v -dof 6 -in brain.nii.gz -ref NMT_SS_0.5mm.nii.gz -o brain_rot2atl -omat brain_rot2atl.mat -interp sinc
    native_brain_to_template_brain = pe.Node(interface=fsl.FLIRT(),
                                             name=f'native_brain_to_template_'
                                                  f'brain_{pipe_num}')
    native_brain_to_template_brain.inputs.dof = 6
    native_brain_to_template_brain.inputs.interp = 'sinc'
    wf.connect(unet_masked_brain, 'out_file',
               native_brain_to_template_brain, 'in_file')

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, native_brain_to_template_brain, 'reference')

    # flirt -in head.nii.gz -ref NMT_0.5mm.nii.gz -o head_rot2atl -applyxfm -init brain_rot2atl.mat
    native_head_to_template_head = pe.Node(interface=fsl.FLIRT(),
                                           name=f'native_head_to_template_'
                                                f'head_{pipe_num}')
    native_head_to_template_head.inputs.apply_xfm = True

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, native_head_to_template_head, 'in_file')

    wf.connect(native_brain_to_template_brain, 'out_matrix_file',
               native_head_to_template_head, 'in_matrix_file')

    node, out = strat_pool.get_data('T1w_template')
    wf.connect(node, out, native_head_to_template_head, 'reference')

    # fslmaths NMT_SS_0.5mm.nii.gz -bin templateMask.nii.gz
    template_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                  name=f'template_brain_mask_{pipe_num}')
    template_brain_mask.inputs.args = '-bin'

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, template_brain_mask, 'in_file')

    # ANTS 3 -m  CC[head_rot2atl.nii.gz,NMT_0.5mm.nii.gz,1,5] -t SyN[0.25] -r Gauss[3,0] -o atl2T1rot -i 60x50x20 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
    ants_template_head_to_template = pe.Node(interface=ants.Registration(),
                                             name=f'template_head_to_'
                                                  f'template_{pipe_num}')
    ants_template_head_to_template.inputs.metric = ['CC']
    ants_template_head_to_template.inputs.metric_weight = [1, 5]
    ants_template_head_to_template.inputs.transforms = ['SyN']
    ants_template_head_to_template.inputs.transform_parameters = [(0.25,)]
    ants_template_head_to_template.inputs.interpolation = 'NearestNeighbor'
    ants_template_head_to_template.inputs.number_of_iterations = [
        [60, 50, 20]]
    ants_template_head_to_template.inputs.smoothing_sigmas = [[0.6, 0.2, 0.0]]
    ants_template_head_to_template.inputs.shrink_factors = [[4, 2, 1]]
    ants_template_head_to_template.inputs.convergence_threshold = [1.e-8]
    wf.connect(native_head_to_template_head, 'out_file',
               ants_template_head_to_template, 'fixed_image')

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, ants_template_head_to_template, 'moving_image')

    # antsApplyTransforms -d 3 -i templateMask.nii.gz -t atl2T1rotWarp.nii.gz atl2T1rotAffine.txt -r brain_rot2atl.nii.gz -o brain_rot2atl_mask.nii.gz
    template_head_transform_to_template = pe.Node(
        interface=ants.ApplyTransforms(),
        name=f'template_head_transform_to_template_{pipe_num}')
    template_head_transform_to_template.inputs.dimension = 3

    wf.connect(template_brain_mask, 'out_file',
               template_head_transform_to_template, 'input_image')
    wf.connect(native_brain_to_template_brain, 'out_file',
               template_head_transform_to_template, 'reference_image')
    wf.connect(ants_template_head_to_template, 'forward_transforms',
               template_head_transform_to_template, 'transforms')

    # convert_xfm -omat brain_rot2native.mat -inverse brain_rot2atl.mat 
    invt = pe.Node(interface=fsl.ConvertXFM(), name='convert_xfm')
    invt.inputs.invert_xfm = True
    wf.connect(native_brain_to_template_brain, 'out_matrix_file', invt,
               'in_file')

    # flirt -in brain_rot2atl_mask.nii.gz -ref brain.nii.gz -o brain_mask.nii.gz -applyxfm -init brain_rot2native.mat
    template_brain_to_native_brain = pe.Node(interface=fsl.FLIRT(),
                                             name=f'template_brain_to_native_'
                                                  f'brain_{pipe_num}')
    template_brain_to_native_brain.inputs.apply_xfm = True
    wf.connect(template_head_transform_to_template, 'output_image',
               template_brain_to_native_brain, 'in_file')
    wf.connect(unet_masked_brain, 'out_file', template_brain_to_native_brain,
               'reference')
    wf.connect(invt, 'out_file', template_brain_to_native_brain,
               'in_matrix_file')

    # fslmaths brain_mask.nii.gz -thr .5 -bin brain_mask_thr.nii.gz
    refined_mask = pe.Node(interface=fsl.Threshold(), name=f'refined_mask'
                                                           f'_{pipe_num}')
    refined_mask.inputs.thresh = 0.5
    refined_mask.inputs.args = '-bin'
    wf.connect(template_brain_to_native_brain, 'out_file', refined_mask,
               'in_file')

    outputs = {
        'space-T1w_desc-brain_mask': (refined_mask, 'out_file')
    }

    return (wf, outputs)


def freesurfer_brain_connector(wf, cfg, strat_pool, pipe_num, opt):
    # register FS brain mask to native space
    fs_brain_mask_to_native = pe.Node(
        interface=freesurfer.ApplyVolTransform(),
        name='fs_brain_mask_to_native')
    fs_brain_mask_to_native.inputs.reg_header = True

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, fs_brain_mask_to_native, 'source_file')

    node, out = strat_pool.get_data('raw_average')
    wf.connect(node, out, fs_brain_mask_to_native, 'target_file')

    node, out = strat_pool.get_data('freesurfer_subject_dir')
    wf.connect(node, out, fs_brain_mask_to_native, 'subjects_dir')

    # convert brain mask file from .mgz to .nii.gz
    fs_brain_mask_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                                   output_names=['out_file'],
                                                   function=mri_convert),
                                     name=f'fs_brainmask_to_nifti_{pipe_num}')
    wf.connect(fs_brain_mask_to_native, 'transformed_file',
               fs_brain_mask_to_nifti, 'in_file')

    # binarize the brain mask
    binarize_fs_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                     name=f'binarize_fs_brainmask_{pipe_num}')
    binarize_fs_brain_mask.inputs.args = '-bin'
    wf.connect(fs_brain_mask_to_nifti, 'out_file',
               binarize_fs_brain_mask, 'in_file')

    # fill holes
    fill_fs_brain_mask = pe.Node(interface=afni.MaskTool(),
                                 name=f'fill_fs_brainmask_{pipe_num}')
    fill_fs_brain_mask.inputs.fill_holes = True
    fill_fs_brain_mask.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(binarize_fs_brain_mask, 'out_file',
               fill_fs_brain_mask, 'in_file')

    outputs = {
        'space-T1w_desc-brain_mask': (fill_fs_brain_mask, 'out_file')
    }

    return (wf, outputs)


def freesurfer_abcd_brain_connector(wf, cfg, strat_pool, pipe_num, opt):

    ### ABCD harmonization - anatomical brain mask generation ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh#L151-L156

    wmparc_to_nifti = pe.Node(util.Function(input_names=['in_file',
                                                         'reslice_like',
                                                         'args'],
                                            output_names=['out_file'],
                                            function=mri_convert),
                              name=f'wmparc_to_nifti_{pipe_num}')
    wmparc_to_nifti.inputs.args = '-rt nearest'

    node, out = strat_pool.get_data('wmparc')
    wf.connect(node, out, wmparc_to_nifti, 'in_file')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, wmparc_to_nifti, 'reslice_like')

    binary_mask = pe.Node(interface=fsl.maths.MathsCommand(), 
                          name=f'binarize_wmparc_{pipe_num}')
    binary_mask.inputs.args = '-bin -dilD -dilD -dilD -ero -ero'

    wf.connect(wmparc_to_nifti, 'out_file', binary_mask, 'in_file')

    wb_command_fill_holes = pe.Node(util.Function(input_names=['in_file'],
                                                  output_names=['out_file'],
                                                  function=wb_command),
                                    name=f'wb_command_fill_holes_{pipe_num}')

    wf.connect(binary_mask, 'out_file', wb_command_fill_holes, 'in_file')

    binary_filled_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                 name=f'binarize_filled_wmparc_{pipe_num}')
    binary_filled_mask.inputs.args = '-bin'

    wf.connect(wb_command_fill_holes, 'out_file',
               binary_filled_mask, 'in_file')

    brain_mask_to_t1_restore = pe.Node(interface=fsl.ApplyWarp(),
                                       name=f'brain_mask_to_t1_restore_{pipe_num}')
    brain_mask_to_t1_restore.inputs.interp = 'nn'
    brain_mask_to_t1_restore.inputs.premat = cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    wf.connect(binary_filled_mask, 'out_file',
               brain_mask_to_t1_restore, 'in_file')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, brain_mask_to_t1_restore, 'ref_file')

    # TODO check if it's necessary
    # fs_brain = pe.Node(interface=fsl.MultiImageMaths(), name='fs_brain')
    # fs_brain.inputs.op_string = "-mul %s"

    # wf.connect(brain_mask_to_t1_restore, 'out_file', 
    #            fs_brain, 'in_file')

    # node, out = strat_pool.get_data('desc-preproc_T1w')
    # wf.connect(node, out, fs_brain, 'operand_files')

    outputs = {
        'space-T1w_desc-brain_mask': (brain_mask_to_t1_restore, 'out_file')
        # 'desc-brain_T1w': (fs_brain, 'out_file')
    }

    return (wf, outputs)


def anatomical_init(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "anatomical_init",
     "config": "None",
     "switch": "None",
     "option_key": "None",
     "option_val": "None",
     "inputs": ["T1w"],
     "outputs": ["desc-preproc_T1w",
                 "desc-reorient_T1w"]}
    '''

    anat_deoblique = pe.Node(interface=afni.Refit(),
                             name=f'anat_deoblique_{pipe_num}')
    anat_deoblique.inputs.deoblique = True

    node, out = strat_pool.get_data('T1w')
    wf.connect(node, out, anat_deoblique, 'in_file')

    anat_reorient = pe.Node(interface=afni.Resample(),
                            name=f'anat_reorient_{pipe_num}')
    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(anat_deoblique, 'out_file', anat_reorient, 'in_file')

    outputs = {'desc-preproc_T1w': (anat_reorient, 'out_file'),
               'desc-reorient_T1w': (anat_reorient, 'out_file')}

    return (wf, outputs)


def acpc_align_head(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_head",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "T1w_ACPC_template"],
     "outputs": ["desc-preproc_T1w",
                 "from-affine_to-rigid_mode-image_xfm"]}
    '''

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=False,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T1w_ACPC_template')
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'from-affine_to-rigid_mode-image_xfm': (
            acpc_align, 'outputspec.from-affine_to-rigid_mode-image_xfm')
    }

    return (wf, outputs)


def acpc_align_head_with_mask(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_head_with_mask",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 "space-T1w_desc-brain_mask"),
                "T1w_ACPC_template"],
     "outputs": ["desc-preproc_T1w",
                 "space-T1w_desc-brain_mask",
                 "from-affine_to-rigid_mode-image_xfm"]}
    '''

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=True,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T1w_ACPC_template')
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'space-T1w_desc-brain_mask': (
            acpc_align, 'outputspec.acpc_brain_mask'),
        'from-affine_to-rigid_mode-image_xfm': (
            acpc_align, 'outputspec.from-affine_to-rigid_mode-image_xfm')
    }

    return (wf, outputs)


def acpc_align_brain(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_brain",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 "T1w_brain_ACPC_template")],
     "outputs": ["desc-preproc_T1w",
                 "from-affine_to-rigid_mode-image_xfm"]}
    '''

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=False,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T1w_brain_ACPC_template')
    wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'from-affine_to-rigid_mode-image_xfm': (
            acpc_align, 'outputspec.from-affine_to-rigid_mode-image_xfm')
    }

    return (wf, outputs)


def acpc_align_brain_with_mask(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_brain_with_mask",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 "space-T1w_desc-brain_mask"),
                "T1w_brain_ACPC_template"],
     "outputs": ["desc-preproc_T1w",
                 "space-T1w_desc-brain_mask"]}
    '''

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=True,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, acpc_align, 'inputspec.brain_mask')

    node, out = strat_pool.get_data('T1w_brain_ACPC_template')
    wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'space-T1w_desc-brain_mask': (
        acpc_align, 'outputspec.acpc_brain_mask')
    }

    return (wf, outputs)


def non_local_means(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "nlm_filtering",
     "config": ["anatomical_preproc", "non_local_means_filtering"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["desc-preproc_T1w"]}
    '''

    denoise = pe.Node(interface=ants.DenoiseImage(),
                      name=f'anat_denoise_{pipe_num}')

    denoise.inputs.noise_model = cfg.anatomical_preproc['non_local_means_filtering']['noise_model']

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, denoise, 'input_image')

    outputs = {
        'desc-preproc_T1w': (denoise, 'output_image')
    }

    return (wf, outputs)


def n4_bias_correction(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "n4_bias_correction",
     "config": ["anatomical_preproc", "n4_bias_field_correction"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["desc-preproc_T1w",
                 "desc-n4_T1w"]}
    '''

    n4 = pe.Node(interface=ants.N4BiasFieldCorrection(dimension=3,
                                                      copy_header=True),
                 name=f'anat_n4_{pipe_num}')
    n4.inputs.shrink_factor = cfg.anatomical_preproc['n4_bias_field_correction']['shrink_factor']

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, n4, 'input_image')

    outputs = {
        'desc-preproc_T1w': (n4, 'output_image'),
        'desc-n4_T1w': (n4, 'output_image')
    }

    return (wf, outputs)


def brain_mask_afni(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_afni",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "3dSkullStrip",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    wf, outputs = afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


def brain_mask_acpc_afni(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_acpc_afni",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "3dSkullStrip",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-T1w_desc-acpcbrain_mask"]}
    '''

    wf, wf_outputs = afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


def brain_mask_fsl(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_fsl",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "BET",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    wf, outputs = fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


def brain_mask_acpc_fsl(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_acpc_fsl",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "BET",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-T1w_desc-acpcbrain_mask"]}
    '''

    wf, wf_outputs = fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


def brain_mask_niworkflows_ants(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_niworkflows_ants",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "niworkflows-ants",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    wf, outputs = niworkflows_ants_brain_connector(wf, cfg, strat_pool,
                                                   pipe_num, opt)

    return (wf, outputs)


def brain_mask_acpc_niworkflows_ants(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_acpc_niworkflows_ants",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "niworkflows-ants",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-T1w_desc-acpcbrain_mask"]}
    '''

    wf, wf_outputs = niworkflows_ants_brain_connector(wf, cfg, strat_pool,
                                                      pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


def brain_mask_unet(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_unet",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "UNet",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "T1w_brain_template",
                "T1w_template",
                "unet_model"],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    wf, outputs = unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


def brain_mask_acpc_unet(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_acpc_unet",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "UNet",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "T1w_brain_template",
                "T1w_template",
                "unet_model"],
     "outputs": ["space-T1w_desc-acpcbrain_mask"]}
    '''

    wf, wf_outputs = unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


def brain_mask_freesurfer(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_freesurfer",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "Freesurfer",
     "inputs": ["space-T1w_desc-brain_mask",
                "raw_average",
                "freesurfer_subject_dir"],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    wf, outputs = freesurfer_brain_connector(wf, cfg, strat_pool, pipe_num,
                                             opt)

    return (wf, outputs)


def brain_mask_acpc_freesurfer(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_acpc_freesurfer",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "Freesurfer",
     "inputs": ["space-T1w_desc-brain_mask",
                "raw_average",
                "freesurfer_subject_dir"],
     "outputs": ["space-T1w_desc-acpcbrain_mask"]}
    '''

    wf, wf_outputs = freesurfer_brain_connector(wf, cfg, strat_pool, pipe_num,
                                                opt)

    outputs = {'space-T1w_desc-acpcbrain_mask':
                   wf_outputs['space-T1w_desc-brain_mask']}

    return (wf, outputs)


def brain_mask_freesurfer_abcd(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_freesurfer_abcd",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "FreeSurfer-ABCD",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "wmparc",
                "freesurfer_subject_dir"],
     "outputs": ["space-T1w_desc-brain_mask"]}
    '''

    wf, outputs = freesurfer_abcd_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


def brain_mask_acpc_freesurfer_abcd(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_mask_acpc_freesurfer_abcd",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "FreeSurfer-ABCD",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "wmparc",
                "freesurfer_subject_dir"],
     "outputs": ["space-T1w_desc-acpcbrain_mask"]}
    '''

    wf, wf_outputs = freesurfer_abcd_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {'space-T1w_desc-acpcbrain_mask':
                wf_outputs['space-T1w_desc-brain_mask']}

    return (wf, outputs)


def brain_extraction(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_extraction",
     "config": "None",
     "switch": "None",
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 ["space-T1w_desc-brain_mask", "space-T1w_desc-acpcbrain_mask"])],
     "outputs": ["desc-brain_T1w"]}
    '''

    '''
    brain_mask_deoblique = pe.Node(interface=afni.Refit(),
                                   name='brain_mask_deoblique')
    brain_mask_deoblique.inputs.deoblique = True
    wf.connect(inputnode, 'brain_mask',
                    brain_mask_deoblique, 'in_file')

    brain_mask_reorient = pe.Node(interface=afni.Resample(),
                                  name='brain_mask_reorient')
    brain_mask_reorient.inputs.orientation = 'RPI'
    brain_mask_reorient.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(brain_mask_deoblique, 'out_file',
                    brain_mask_reorient, 'in_file')
    '''

    anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                       name=f'brain_extraction_{pipe_num}')

    anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
    anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_a')

    node, out = strat_pool.get_data(['space-T1w_desc-brain_mask',
                                     'space-T1w_desc-acpcbrain_mask'])
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_b')

    outputs = {
        'desc-brain_T1w': (anat_skullstrip_orig_vol, 'out_file')
    }

    return (wf, outputs)


def freesurfer_preproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "freesurfer_preproc",
     "config": ["surface_analysis"],
     "switch": ["run_freesurfer"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 "freesurfer_subject_dir"],
     "outputs": ["space-T1w_desc-brain_mask",
                 "freesurfer_subject_dir",
                 "label-CSF_mask",
                 "label-WM_mask",
                 "label-GM_mask",
                 "surface_curvature",
                 "pial_surface_mesh",
                 "smoothed_surface_mesh",
                 "spherical_surface_mesh",
                 "sulcal_depth_surface_maps",
                 "cortical_thickness_surface_maps",
                 "cortical_volume_surface_maps",
                 "white_matter_surface_mesh",
                 "raw_average"]}
    '''

    # TODO check if freesurfer_subject_dir exists
    # if strat_pool.check_rpool('freesurfer_subject_dir')
    reconall = pe.Node(interface=freesurfer.ReconAll(),
                       name=f'anat_freesurfer_{pipe_num}')

    freesurfer_subject_dir = os.path.join(cfg.pipeline_setup['working_directory']['path'],
        f'anat_preproc_freesurfer_{pipe_num}',
        'anat_freesurfer')

    if not os.path.exists(freesurfer_subject_dir):
        os.makedirs(freesurfer_subject_dir)

    reconall.inputs.directive = 'all'
    reconall.inputs.subjects_dir = freesurfer_subject_dir
    reconall.inputs.openmp = cfg.pipeline_setup['system_config'][
        'num_OMP_threads']

    node, out = strat_pool.get_data(["desc-preproc_T1w", "desc-reorient_T1w",
                                     "T1w"])
    wf.connect(node, out, reconall, 'T1_files')

    # register FS brain mask to native space
    fs_brain_mask_to_native = pe.Node(
        interface=freesurfer.ApplyVolTransform(),
        name='fs_brain_mask_to_native')
    fs_brain_mask_to_native.inputs.reg_header = True

    wf.connect(reconall, 'brainmask', fs_brain_mask_to_native, 'source_file')
    wf.connect(reconall, 'rawavg', fs_brain_mask_to_native, 'target_file')
    wf.connect(reconall, 'subjects_dir',
               fs_brain_mask_to_native, 'subjects_dir')

    # convert brain mask file from .mgz to .nii.gz
    fs_brain_mask_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                                   output_names=['out_file'],
                                                   function=mri_convert),
                                     name='fs_brainmask_to_nifti')
    wf.connect(fs_brain_mask_to_native, 'transformed_file',
               fs_brain_mask_to_nifti, 'in_file')

    # binarize the brain mask
    binarize_fs_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                     name='binarize_fs_brainmask')
    binarize_fs_brain_mask.inputs.args = '-bin'
    wf.connect(fs_brain_mask_to_nifti, 'out_file',
               binarize_fs_brain_mask, 'in_file')

    # fill holes
    fill_fs_brain_mask = pe.Node(interface=afni.MaskTool(),
                                 name='fill_fs_brainmask')
    fill_fs_brain_mask.inputs.fill_holes = True
    fill_fs_brain_mask.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(binarize_fs_brain_mask, 'out_file',
               fill_fs_brain_mask, 'in_file')

    # register FS segmentations (aseg.mgz) to native space
    fs_aseg_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                                name='fs_aseg_to_native')
    fs_aseg_to_native.inputs.reg_header = True
    fs_aseg_to_native.inputs.interp = 'nearest'
    fs_aseg_to_native.inputs.subjects_dir = freesurfer_subject_dir

    wf.connect(reconall, 'aseg', fs_aseg_to_native, 'source_file')
    wf.connect(reconall, 'rawavg', fs_aseg_to_native, 'target_file')

    # convert registered FS segmentations from .mgz to .nii.gz
    fs_aseg_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=mri_convert),
                               name='fs_aseg_to_nifti')
    fs_aseg_to_nifti.inputs.args = '-rt nearest'

    wf.connect(fs_aseg_to_native, 'transformed_file',
               fs_aseg_to_nifti, 'in_file')

    pick_tissue = pe.Node(util.Function(input_names=['multiatlas_Labels'],
                                        output_names=['csf_mask', 'gm_mask',
                                                      'wm_mask'],
                                        function=pick_tissue_from_labels_file),
                          name=f'anat_preproc_freesurfer_tissue_mask')
    pick_tissue.inputs.include_ventricles = True

    wf.connect(fs_aseg_to_nifti, 'out_file', pick_tissue, 'multiatlas_Labels')

    outputs = {
        'space-T1w_desc-brain_mask': (fill_fs_brain_mask, 'out_file'),
        'freesurfer_subject_dir': (reconall, 'subjects_dir'),
        'label-CSF_mask': (pick_tissue, 'csf_mask'),
        'label-WM_mask': (pick_tissue, 'wm_mask'),
        'label-GM_mask': (pick_tissue, 'gm_mask'),
        'surface_curvature': (reconall, 'curv'),
        'pial_surface_mesh': (reconall, 'pial'),
        'smoothed_surface_mesh': (reconall, 'smoothwm'),
        'spherical_surface_mesh': (reconall, 'sphere'),
        'sulcal_depth_surface_maps': (reconall, 'sulc'),
        'cortical_thickness_surface_maps': (reconall, 'thickness'),
        'cortical_volume_surface_maps': (reconall, 'volume'),
        'white_matter_surface_mesh': (reconall, 'white'),
        'raw_average': (reconall, 'rawavg')
    }

    return (wf, outputs)


def fnirt_based_brain_extraction(config=None, wf_name='fnirt_based_brain_extraction'):

    ### ABCD Harmonization - FNIRT-based brain extraction ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PreFreeSurfer/scripts/BrainExtraction_FNIRTbased.sh

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_data',
                                                       'ref_mask_res-2',
                                                       'template_skull_for_anat',
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
    # non_linear_reg.inputs.warped_file = 'T1w_acpc_to_MNI_nonlin.nii.gz' # --iout
    # non_linear_reg.inputs.log_file = 'NonlinearReg.txt' # --logout
    non_linear_reg.inputs.out_intensitymap_file = True # --intout
    non_linear_reg.inputs.fieldcoeff_file = True # --cout
    non_linear_reg.inputs.config_file = config.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT']['fnirt_config']

    preproc.connect(inputnode, 'anat_data',
                    non_linear_reg, 'in_file')

    preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                    non_linear_reg, 'ref_file')

    preproc.connect(linear_reg, 'out_matrix_file',
                    non_linear_reg, 'affine_file')

    preproc.connect(inputnode, 'ref_mask_res-2',
                    non_linear_reg, 'refmask_file')

    # Overwrite the image output from FNIRT with a spline interpolated highres version
    # creating spline interpolated hires version
    # applywarp --rel --interp=spline --in="$Input" --ref="$Reference" -w "$WD"/str2standard.nii.gz --out="$WD"/"$BaseName"_to_MNI_nonlin.nii.gz
    apply_warp = pe.Node(interface=fsl.ApplyWarp(),
                        name='apply_warp')

    apply_warp.inputs.interp = 'spline'
    apply_warp.inputs.relwarp = True

    preproc.connect(inputnode, 'anat_data',
                    apply_warp, 'in_file')

    preproc.connect(inputnode, 'template_skull_for_anat',
                    apply_warp, 'ref_file')

    preproc.connect(non_linear_reg, 'field_file',
                    apply_warp, 'field_file')

    # Invert warp and transform dilated brain mask back into native space, and use it to mask input image
    # Input and reference spaces are the same, using 2mm reference to save time
    # invwarp --ref="$Reference2mm" -w "$WD"/str2standard.nii.gz -o "$WD"/standard2str.nii.gz
    inverse_warp = pe.Node(interface=fsl.InvWarp(), name='inverse_warp')
    inverse_warp.inputs.output_type = 'NIFTI_GZ'

    preproc.connect(inputnode, 'template_skull_for_anat_2mm',
                    inverse_warp, 'reference')

    preproc.connect(non_linear_reg, 'field_file',
                    inverse_warp, 'warp')

    # Apply inverse warp
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
    
    # Apply mask to create brain
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


def fast_bias_field_correction(config=None, wf_name='fast_bias_field_correction'):

    ### ABCD Harmonization - FAST bias field correction ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PreFreeSurfer/PreFreeSurferPipeline.sh#L688-L694

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_data',
                                                       'anat_brain',
                                                       'anat_brain_mask']), 
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['anat_restore',
                                                        'anat_brain_restore',
                                                        'bias_field']),
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

    preproc.connect(fast_bias_field_correction, 'bias_field',
                    outputnode, 'bias_field')

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
    apply_mask = pe.Node(interface=fsl.MultiImageMaths(),
                            name='apply_mask')
    apply_mask.inputs.op_string = '-mul %s'

    preproc.connect(inputnode, 'anat_data',
                    apply_mask, 'in_file')

    preproc.connect(inverse_brain_mask, 'out_file',
                    apply_mask, 'operand_files')

    # fslmaths ${T1wFolder}/T1w_fast_restore.nii.gz -add ${T1wFolder}/T1w_acpc_dc_skull.nii.gz ${T1wFolder}/${T1wImage}_acpc_dc_restore
    anat_restore = pe.Node(interface=fsl.MultiImageMaths(),
                            name='get_anat_restore')
    anat_restore.inputs.op_string = '-add %s'

    preproc.connect(fast_bias_field_correction, 'restored_image',
                    anat_restore, 'in_file')

    preproc.connect(apply_mask, 'out_file',
                    anat_restore, 'operand_files')

    preproc.connect(anat_restore, 'out_file',
                    outputnode, 'anat_restore')

    return preproc


def freesurfer_abcd_preproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "freesurfer_abcd_preproc",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "FreeSurfer-ABCD",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 "T1w_template",
                 "T1w_brain_template_mask",
                 "ref_mask_res-2",
                 "T1w_template_res-2"],
     "outputs": ["desc-restore_T1w",
                 "desc-restore-brain_T1w",
                 "desc-fast_biasfield",
                 "wmparc",
                 "freesurfer_subject_dir"]}
    '''

    # fnirt-based brain extraction
    brain_extraction = fnirt_based_brain_extraction(config=cfg,
                                                    wf_name=f'fnirt_based_brain_extraction_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, brain_extraction, 'inputspec.anat_data')

    node, out = strat_pool.get_data('ref_mask_res-2')
    wf.connect(node, out, brain_extraction, 'inputspec.ref_mask_res-2')

    node, out = strat_pool.get_data('T1w_template')
    wf.connect(node, out, brain_extraction, 'inputspec.template_skull_for_anat')

    node, out = strat_pool.get_data('T1w_template_res-2')
    wf.connect(node, out, brain_extraction, 'inputspec.template_skull_for_anat_2mm')

    node, out = strat_pool.get_data('T1w_brain_template_mask')
    wf.connect(node, out, brain_extraction, 'inputspec.template_brain_mask_for_anat')

    # fast bias field correction
    fast_correction = fast_bias_field_correction(config=cfg,
                                                 wf_name=f'fast_bias_field_correction_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, fast_correction, 'inputspec.anat_data')

    wf.connect(brain_extraction, 'outputspec.anat_brain', fast_correction, 'inputspec.anat_brain')

    wf.connect(brain_extraction, 'outputspec.anat_brain_mask', fast_correction, 'inputspec.anat_brain_mask')

    ### ABCD Harmonization ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/FreeSurfer/FreeSurferPipeline.sh#L140-L144

    # flirt -interp spline -in "$T1wImage" -ref "$T1wImage" -applyisoxfm 1 -out "$T1wImageFile"_1mm.nii.gz
    resample_head_1mm = pe.Node(interface=fsl.FLIRT(),
                                name=f'resample_anat_head_1mm_{pipe_num}')
    resample_head_1mm.inputs.interp = 'spline'
    resample_head_1mm.inputs.apply_isoxfm = 1

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, resample_head_1mm, 'in_file')

    wf.connect(node, out, resample_head_1mm, 'reference')

    # applywarp --rel --interp=spline -i "$T1wImage" -r "$T1wImageFile"_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wImageFile"_1mm.nii.gz
    applywarp_head_to_head_1mm = pe.Node(interface=fsl.ApplyWarp(),
                                         name=f'applywarp_head_to_head_1mm_{pipe_num}')
    applywarp_head_to_head_1mm.inputs.relwarp = True
    applywarp_head_to_head_1mm.inputs.interp = 'spline'
    applywarp_head_to_head_1mm.inputs.premat = cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    wf.connect(node, out, applywarp_head_to_head_1mm, 'in_file')

    wf.connect(resample_head_1mm, 'out_file',
               applywarp_head_to_head_1mm, 'ref_file')

    # applywarp --rel --interp=nn -i "$T1wImageBrain" -r "$T1wImageFile"_1mm.nii.gz --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wImageBrainFile"_1mm.nii.gz
    applywarp_brain_to_head_1mm = pe.Node(interface=fsl.ApplyWarp(),
                name=f'applywarp_brain_to_head_1mm_{pipe_num}')
    applywarp_brain_to_head_1mm.inputs.relwarp = True
    applywarp_brain_to_head_1mm.inputs.interp = 'nn'
    applywarp_brain_to_head_1mm.inputs.premat = cfg.registration_workflows['anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    wf.connect(fast_correction, 'outputspec.anat_brain_restore',
                    applywarp_brain_to_head_1mm, 'in_file')

    wf.connect(resample_head_1mm, 'out_file',
                    applywarp_brain_to_head_1mm, 'ref_file')

    # fslstats $T1wImageBrain -M
    average_brain = pe.Node(interface=fsl.ImageStats(),
                name=f'average_brain_{pipe_num}')
    average_brain.inputs.op_string = '-M'
    average_brain.inputs.output_type = 'NIFTI_GZ'

    wf.connect(fast_correction, 'outputspec.anat_brain_restore',
               average_brain, 'in_file')

    # fslmaths "$T1wImageFile"_1mm.nii.gz -div $Mean -mul 150 -abs "$T1wImageFile"_1mm.nii.gz
    normalize_head = pe.Node(util.Function(input_names=['in_file', 'number', 'out_file_suffix'],
                                           output_names=['out_file'],
                                           function=fslmaths_command),
                             name=f'normalize_head_{pipe_num}')
    normalize_head.inputs.out_file_suffix = '_norm'

    wf.connect(applywarp_head_to_head_1mm, 'out_file', 
               normalize_head, 'in_file')

    wf.connect(average_brain, 'out_stat',
               normalize_head, 'number')

    ### recon-all -all step ###
    reconall = pe.Node(interface=freesurfer.ReconAll(),
                name=f'anat_freesurfer_{pipe_num}')

    sub_dir = cfg.pipeline_setup['working_directory']['path']
    freesurfer_subject_dir = os.path.join(sub_dir,
                                          f'anat_preproc_freesurfer_{pipe_num}',
                                          'anat_freesurfer')

    # create the directory for FreeSurfer node
    if not os.path.exists(freesurfer_subject_dir):
        os.makedirs(freesurfer_subject_dir)

    reconall.inputs.directive = 'all'
    reconall.inputs.subjects_dir = freesurfer_subject_dir
    reconall.inputs.openmp = cfg.pipeline_setup['system_config']['num_OMP_threads']

    wf.connect(normalize_head, 'out_file',
               reconall, 'T1_files')

    outputs = {
        'desc-restore_T1w': (fast_correction, 'outputspec.anat_restore'),
        'desc-restore-brain_T1w': (fast_correction, 'outputspec.anat_brain_restore'),
        'desc-fast_biasfield': (fast_correction, 'outputspec.bias_field'),
        'wmparc': (reconall, 'wmparc'),
        'freesurfer_subject_dir': (reconall, 'subjects_dir'),
    }

    return (wf, outputs)


def correct_restore_brain_intensity_abcd(wf, cfg, strat_pool, pipe_num, opt=None):
    # TODO check docstring? function name?
    '''
    {"name": "correct_restore_brain_intensity_abcd",
     "config": ["anatomical_preproc", "brain_extraction"],
     "switch": "None",
     "option_key": "using",
     "option_val": "FreeSurfer-ABCD",
     "inputs": [(["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                 "desc-n4_T1w",
                 "desc-restore-brain_T1w",
                 "space-T1w_desc-brain_mask",
                 "desc-fast_biasfield",
                 "from-affine_to-rigid_mode-image_xfm",
                 "from-T1w_to-template_mode-image_xfm")],
     "outputs": ["desc-restore-brain_T1w"]}
    '''

    ### ABCD Harmonization - Myelin Map ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PreFreeSurfer/PreFreeSurferPipeline.sh#L655-L656
    # fslmerge -t ${T1wFolder}/xfms/${T1wImage}_dc ${T1wFolder}/${T1wImage}_acpc ${T1wFolder}/${T1wImage}_acpc ${T1wFolder}/${T1wImage}_acpc
    merge_t1_acpc_to_list = pe.Node(util.Merge(3),
                                    name=f'merge_t1_acpc_to_list_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, merge_t1_acpc_to_list, 'in1')
    wf.connect(node, out, merge_t1_acpc_to_list, 'in2')
    wf.connect(node, out, merge_t1_acpc_to_list, 'in3')

    merge_t1_acpc = pe.Node(interface=fsl.Merge(),
                            name='merge_t1_acpc')

    merge_t1_acpc.inputs.dimension = 't'

    wf.connect(merge_t1_acpc_to_list, 'out',
        merge_t1_acpc, 'in_files')

    # fslmaths ${T1wFolder}/xfms/${T1wImage}_dc -mul 0 ${T1wFolder}/xfms/${T1wImage}_dc
    multiply_t1_acpc_by_zero = pe.Node(interface=fsl.ImageMaths(),
                                       name=f'multiply_t1_acpc_by_zero_{pipe_num}')
    
    multiply_t1_acpc_by_zero.inputs.op_string = '-mul 0'

    wf.connect(merge_t1_acpc, 'merged_file', 
        multiply_t1_acpc_by_zero, 'in_file')

    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh#L157
    # convertwarp --relout --rel --ref="$T1wFolder"/"$T1wImageBrainMask" --premat="$T1wFolder"/xfms/"$InitialT1wTransform" \
    # --warp1="$T1wFolder"/xfms/"$dcT1wTransform" --out="$T1wFolder"/xfms/"$OutputOrigT1wToT1w"
    convertwarp_orig_t1_to_t1 = pe.Node(interface=fsl.ConvertWarp(), 
                                        name=f'convertwarp_orig_t1_to_t1_{pipe_num}')

    convertwarp_orig_t1_to_t1.inputs.out_relwarp = True
    convertwarp_orig_t1_to_t1.inputs.relwarp = True

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, convertwarp_orig_t1_to_t1, 'reference')

    node, out = strat_pool.get_data('from-affine_to-rigid_mode-image_xfm')
    wf.connect(node, out, convertwarp_orig_t1_to_t1, 'premat')
    wf.connect(multiply_t1_acpc_by_zero, 'out_file',
        convertwarp_orig_t1_to_t1, 'warp1')

    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/scripts/CreateMyelinMaps.sh#L72-L73
    # applywarp --rel --interp=spline -i "$BiasField" -r "$T1wImageBrain" -w "$AtlasTransform" -o "$BiasFieldOutput"
    applywarp_biasfield = pe.Node(interface=fsl.ApplyWarp(), 
                                  name=f'applywarp_biasfield_{pipe_num}')

    applywarp_biasfield.inputs.relwarp = True
    applywarp_biasfield.inputs.interp = 'spline'

    node, out = strat_pool.get_data('desc-fast_biasfield')
    wf.connect(node, out, applywarp_biasfield, 'in_file')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, applywarp_biasfield, 'ref_file')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, applywarp_biasfield, 'field_file')

    # fslmaths "$BiasFieldOutput" -thr 0.1 "$BiasFieldOutput"
    threshold_biasfield = pe.Node(interface=fsl.ImageMaths(),
                                  name=f'threshold_biasfield_{pipe_num}')

    threshold_biasfield.inputs.op_string = '-thr 0.1'
    wf.connect(applywarp_biasfield, 'out_file', 
        threshold_biasfield, 'in_file')

    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/scripts/CreateMyelinMaps.sh#L67-L70
    # applywarp --rel --interp=spline -i "$OrginalT1wImage" -r "$T1wImageBrain" -w "$OutputOrigT1wToT1w" -o "$OutputT1wImage"
    applywarp_t1 = pe.Node(interface=fsl.ApplyWarp(), 
                           name=f'applywarp_t1_{pipe_num}')
    
    applywarp_t1.inputs.relwarp = True
    applywarp_t1.inputs.interp = 'spline'
    
    node, out = strat_pool.get_data('desc-n4_T1w')
    wf.connect(node, out, applywarp_t1, 'in_file')
    
    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, applywarp_t1, 'ref_file')
    
    wf.connect(convertwarp_orig_t1_to_t1, 'out_file',
        applywarp_t1, 'field_file')

    # fslmaths "$OutputT1wImage" -abs "$OutputT1wImage" -odt float
    abs_t1 = pe.Node(interface=fsl.ImageMaths(),
                     name=f'abs_t1_{pipe_num}')

    abs_t1.inputs.op_string = '-abs'
    wf.connect(applywarp_t1, 'out_file', abs_t1, 'in_file')

    # fslmaths "$OutputT1wImage" -div "$BiasField" "$OutputT1wImageRestore"
    div_t1_by_biasfield = pe.Node(interface=fsl.ImageMaths(),
                                  name=f'div_t1_by_biasfield_{pipe_num}')

    div_t1_by_biasfield.inputs.op_string = '-div'

    wf.connect(abs_t1, 'out_file', div_t1_by_biasfield, 'in_file')

    node, out = strat_pool.get_data('desc-fast_biasfield')
    wf.connect(node, out, div_t1_by_biasfield, 'in_file2')

    # fslmaths "$OutputT1wImageRestore" -mas "$T1wImageBrain" "$OutputT1wImageRestoreBrain"
    apply_mask = pe.Node(interface=fsl.maths.ApplyMask(),
                         name=f'get_restored_corrected_brain_{pipe_num}')

    wf.connect(div_t1_by_biasfield, 'out_file',
        apply_mask, 'in_file')

    node, out = strat_pool.get_data('space-T1w_desc-brain_mask')
    wf.connect(node, out, apply_mask, 'mask_file')

    outputs = {
        'desc-restore-brain_T1w': (apply_mask, 'out_file')
    }

    return (wf, outputs)