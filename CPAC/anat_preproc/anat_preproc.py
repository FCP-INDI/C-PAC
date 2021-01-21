# -*- coding: utf-8 -*-
import os
from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces.fsl import utils as fsl_utils
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.anat_preproc.ants import init_brain_extraction_wf
from nipype.interfaces import freesurfer
from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string, \
    fsl_aff_to_rigid, \
    mri_convert
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


def acpc_alignment(config=None, acpc_target='whole-head', mask=False,
                   wf_name='acpc_align'):
    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_leaf',
                                                       'anat_brain',
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
                              name='anat_skullstrip_args')

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
                              name='anat_skullstrip')
    anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, anat_skullstrip, 'in_file')
    wf.connect(skullstrip_args, 'expr', anat_skullstrip, 'args')

    # Generate anatomical brain mask
    anat_brain_mask = pe.Node(interface=afni.Calc(),
                              name='anat_brain_mask')

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
        name='BET_options')

    anat_skullstrip = pe.Node(
        interface=fsl.BET(), name='anat_skullstrip')
    anat_skullstrip.inputs.output_type = 'NIFTI_GZ'

    inputnode_bet.inputs.set(
        frac=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['frac'],
        mask_boolean=
        cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['mask_boolean'],
        mesh_boolean=
        cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['mesh_boolean'],
        outline=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['outline'],
        padding=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['padding'],
        radius=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['radius'],
        reduce_bias=
        cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['reduce_bias'],
        remove_eyes=
        cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['remove_eyes'],
        robust=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['robust'],
        skull=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['skull'],
        surfaces=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['surfaces'],
        threshold=cfg.anatomical_wf['brain_extraction']['extraction'][
            'FSL-BET']['threshold'],
        vertical_gradient=
        cfg.anatomical_wf['brain_extraction']['extraction'][
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
                                                        'extraction'][
                                                        'niworkflows-ants'][
                                                        'template_path'],
                                                    tpl_mask_path=
                                                    cfg.anatomical_preproc[
                                                        'brain_extraction'][
                                                        'extraction'][
                                                        'niworkflows-ants'][
                                                        'mask_path'],
                                                    tpl_regmask_path=
                                                    cfg.anatomical_preproc[
                                                        'brain_extraction'][
                                                        'extraction'][
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
                        name='unet_mask')

    node, out = strat_pool.get_data('unet_model')
    wf.connect(node, out, unet_mask, 'model_path')

    node, out = strat_pool.get_data(['desc-wf_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, unet_mask, 'cimg_in')

    """
    Revised mask with ANTs
    """
    # fslmaths <whole head> -mul <mask> brain.nii.gz
    unet_masked_brain = pe.Node(interface=fsl.MultiImageMaths(),
                                name='unet_masked_brain')
    unet_masked_brain.inputs.op_string = "-mul %s"

    node, out = strat_pool.get_data(['desc-wf_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, unet_masked_brain, 'in_file')
    wf.connect(unet_mask, 'out_path', unet_masked_brain, 'operand_files')

    # flirt -v -dof 6 -in brain.nii.gz -ref NMT_SS_0.5mm.nii.gz -o brain_rot2atl -omat brain_rot2atl.mat -interp sinc
    # TODO: antsRegistration -z 0 -d 3 -r [NMT_SS_0.5mm.nii.gz,brain.nii.gz,0] -o [transform,brain_rot2atl.nii.gz,brain_inv_rot2atl.nii.gz] -t Rigid[0.1] -m MI[NMT_SS_0.5mm.nii.gz,brain.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-08,10] -s 3.0x2.0x1.0x0.0 -f 8x4x2x1 -u 1 -t Affine[0.1] -m MI[NMT_SS_0.5mm.nii.gz,brain.nii.gz,1,32,Regular,0.25] -c [1000x500x250x100,1e-08,10] -s 3.0x2.0x1.0x0.0 -f 8x4x2x1 -u 1
    native_brain_to_template_brain = pe.Node(interface=fsl.FLIRT(),
                                             name='native_brain_to_template_brain')
    native_brain_to_template_brain.inputs.dof = 6
    native_brain_to_template_brain.inputs.interp = 'sinc'
    wf.connect(unet_masked_brain, 'out_file',
               native_brain_to_template_brain, 'in_file')

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, native_brain_to_template_brain, 'reference')

    # flirt -in head.nii.gz -ref NMT_0.5mm.nii.gz -o head_rot2atl -applyxfm -init brain_rot2atl.mat
    # TODO: antsApplyTransforms -d 3 -i head.nii.gz -r NMT_0.5mm.nii.gz -n Linear -o head_rot2atl.nii.gz -v -t transform1Rigid.mat -t transform2Affine.mat -t transform0DerivedInitialMovingTranslation.mat
    native_head_to_template_head = pe.Node(interface=fsl.FLIRT(),
                                           name='native_head_to_template_head')
    native_head_to_template_head.inputs.apply_xfm = True

    node, out = strat_pool.get_data(['desc-wf_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, native_head_to_template_head, 'in_file')

    wf.connect(native_brain_to_template_brain, 'out_matrix_file',
               native_head_to_template_head, 'in_matrix_file')

    node, out = strat_pool.get_data('T1w_template')
    wf.connect(node, out, native_head_to_template_head, 'reference')

    # fslmaths NMT_SS_0.5mm.nii.gz -bin templateMask.nii.gz
    template_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                  name='template_brain_mask')
    template_brain_mask.inputs.args = '-bin'

    node, out = strat_pool.get_data('T1w_brain_template')
    wf.connect(node, out, template_brain_mask, 'in_file')

    # ANTS 3 -m  CC[head_rot2atl.nii.gz,NMT_0.5mm.nii.gz,1,5] -t SyN[0.25] -r Gauss[3,0] -o atl2T1rot -i 60x50x20 --use-Histogram-Matching  --number-of-affine-iterations 10000x10000x10000x10000x10000 --MI-option 32x16000
    ants_template_head_to_template = pe.Node(interface=ants.Registration(),
                                             name='template_head_to_template')
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
        name='template_head_transform_to_template')
    template_head_transform_to_template.inputs.dimension = 3

    wf.connect(template_brain_mask, 'out_file',
               template_head_transform_to_template, 'input_image')
    wf.connect(native_brain_to_template_brain, 'out_file',
               template_head_transform_to_template, 'reference_image')
    wf.connect(ants_template_head_to_template, 'forward_transforms',
               template_head_transform_to_template, 'transforms')

    # TODO: replace convert_xfm and flirt with:
    # antsApplyTransforms -d 3 -i brain_rot2atl_mask.nii.gz -r brain.nii.gz -n linear -o brain_mask.nii.gz -t [transform0DerivedInitialMovingTranslation.mat,1] -t [transform2Affine.mat,1] -t [transform1Rigid.mat,1]
    # convert_xfm -omat brain_rot2native.mat -inverse brain_rot2atl.mat 
    invt = pe.Node(interface=fsl.ConvertXFM(), name='convert_xfm')
    invt.inputs.invert_xfm = True
    wf.connect(native_brain_to_template_brain, 'out_matrix_file', invt,
               'in_file')

    # flirt -in brain_rot2atl_mask.nii.gz -ref brain.nii.gz -o brain_mask.nii.gz -applyxfm -init brain_rot2native.mat
    template_brain_to_native_brain = pe.Node(interface=fsl.FLIRT(),
                                             name='template_brain_to_native_brain')
    template_brain_to_native_brain.inputs.apply_xfm = True
    wf.connect(template_head_transform_to_template, 'output_image',
               template_brain_to_native_brain, 'in_file')
    wf.connect(unet_masked_brain, 'out_file', template_brain_to_native_brain,
               'reference')
    wf.connect(invt, 'out_file', template_brain_to_native_brain,
               'in_matrix_file')

    # fslmaths brain_mask.nii.gz -thr .5 -bin brain_mask_thr.nii.gz
    refined_mask = pe.Node(interface=fsl.Threshold(), name='refined_mask')
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

    outputs = {
        'space-T1w_desc-brain_mask': (fill_fs_brain_mask, 'out_file')
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
     "outputs": ["desc-reorient_T1w"]}
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

    outputs = {'desc-reorient_T1w': (anat_reorient, 'out_file')}

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
     "outputs": ["desc-preproc_T1w"]}
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
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head')
    }

    return (wf, outputs)


def acpc_align_head_with_mask(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_head_with_mask",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "space-T1w_desc-brain_mask",
                "T1w_ACPC_template"],
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

    node, out = strat_pool.get_data('T1w_ACPC_template')
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'space-T1w_desc-brain_mask': (
        acpc_align, 'outputspec.acpc_brain_mask')
    }

    return (wf, outputs)


def acpc_align_brain(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_brain",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "T1w_brain_ACPC_template"],
     "outputs": ["desc-preproc_T1w"]}
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
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head')
    }

    return (wf, outputs)


def acpc_align_brain_with_mask(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "acpc_alignment_brain_with_mask",
     "config": ["anatomical_preproc", "acpc_alignment"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                "space-T1w_desc-brain_mask",
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
     "config": ["anatomical_preproc"],
     "switch": ["non_local_means_filtering"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["desc-preproc_T1w"]}
    '''
    # TODO: naming, and put two of these in a row and try!
    denoise = pe.Node(interface=ants.DenoiseImage(),
                      name=f'anat_denoise_{pipe_num}')

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
     "config": ["anatomical_preproc"],
     "switch": ["n4_bias_field_correction"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["desc-preproc_T1w"]}
    '''

    n4 = pe.Node(interface=ants.N4BiasFieldCorrection(dimension=3,
                                                      shrink_factor=2,
                                                      copy_header=True),
                 name=f'anat_n4_{pipe_num}')

    node, out = strat_pool.get_data(['desc-preproc_T1w', 'desc-reorient_T1w',
                                     'T1w'])
    wf.connect(node, out, n4, 'input_image')

    outputs = {
        'desc-preproc_T1w': (n4, 'output_image')
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
                "T1w_template"],
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
                "T1w_template"],
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


def brain_extraction(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "brain_extraction",
     "config": "None",
     "switch": "None",
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_T1w", "desc-reorient_T1w", "T1w"],
                ["space-T1w_desc-brain_mask", "space-T1w_desc-acpcbrain_mask"]],
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
                "T1w_brain_template",
                "T1w_template"],
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
                 "white_matter_surface_mesh"]}
    '''

    reconall = pe.Node(interface=freesurfer.ReconAll(),
                       name=f'anat_freesurfer_{pipe_num}')

    freesurfer_subject_dir = os.path.join(
        cfg.pipeline_setup['working_directory'],
        f'anat_preproc_freesurfer_{pipe_num}',
        'anat_freesurfer')

    if not os.path.exists(freesurfer_subject_dir):
        os.mkdirs(freesurfer_subject_dir)

    reconall.inputs.directive = 'all'
    reconall.inputs.subjects_dir = freesurfer_subject_dir
    reconall.inputs.openmp = cfg.pipeline_setup['system_config'][
        'num_OMP_threads']

    node, out = strat_pool.get_data(["desc-preproc_T1w", "desc-reorient_T1w",
                                     "T1w"])
    wf.connect(node, out, reconall, 'T1_files')

    ### segmentation output ###
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
