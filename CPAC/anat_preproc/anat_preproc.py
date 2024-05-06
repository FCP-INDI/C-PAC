# -*- coding: utf-8 -*-
# Copyright (C) 2012-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
# from copy import deepcopy
import os
from CPAC.pipeline.nodeblock import nodeblock
from nipype.interfaces import afni
from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces import freesurfer
import nipype.interfaces.utility as util
from nipype.interfaces.fsl import utils as fsl_utils
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.anat_preproc.ants import init_brain_extraction_wf
from CPAC.anat_preproc.utils import create_3dskullstrip_arg_string, \
    freesurfer_hemispheres, \
    fsl_aff_to_rigid, \
    mri_convert, \
    wb_command, \
    fslmaths_command, \
    VolumeRemoveIslands, \
    normalize_wmparc, \
    pad
from CPAC.utils.interfaces.fsl import Merge as fslMerge


def acpc_alignment(config=None, acpc_target='whole-head', mask=False,
                   wf_name='acpc_align'):
    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_leaf',
                                                       'anat_brain',
                                                       'brain_mask',
                                                       'template_brain_only_for_anat',
                                                       'template_brain_for_acpc',
                                                       'template_head_for_acpc']),
                        name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['acpc_aligned_head',
                                                         'acpc_brain_mask',
                                                         'from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm',
                                                         'acpc_aligned_brain',
                                                         'acpc_brain_mask']),
                          name='outputspec')
    if config.anatomical_preproc['acpc_alignment']['FOV_crop'] == 'robustfov':
        robust_fov = pe.Node(interface=fsl_utils.RobustFOV(),
                            name='anat_acpc_1_robustfov')
        robust_fov.inputs.brainsize = config.anatomical_preproc['acpc_alignment']['brain_size']
        robust_fov.inputs.out_transform = 'fov_xfm.mat'

        fov, in_file = (robust_fov, 'in_file')
        fov, fov_mtx = (robust_fov, 'out_transform')
        fov, fov_outfile = (robust_fov, 'out_roi')
    
    elif config.anatomical_preproc['acpc_alignment']['FOV_crop'] == 'flirt':
        # robustfov doesn't work on some monkey data. prefer using flirt.
        # ${FSLDIR}/bin/flirt -in "${Input}" -applyxfm -ref "${Input}" -omat "$WD"/roi2full.mat -out "$WD"/robustroi.nii.gz
        # adopted from DCAN NHP https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/PreFreeSurfer/scripts/ACPCAlignment.sh#L80-L81
        flirt_fov = pe.Node(interface=fsl.FLIRT(),
                                name='anat_acpc_1_fov')
        flirt_fov.inputs.args = '-applyxfm'

        fov, in_file = (flirt_fov, 'in_file')
        fov, ref_file = (flirt_fov, 'reference')
        fov, fov_mtx = (flirt_fov, 'out_matrix_file')
        fov, fov_outfile = (flirt_fov, 'out_file')

    # align head-to-head to get acpc.mat (for human)
    if acpc_target == 'whole-head':
        preproc.connect(inputnode, 'anat_leaf', fov, in_file)
        if config.anatomical_preproc['acpc_alignment']['FOV_crop'] == 'flirt':
            preproc.connect(inputnode, 'anat_leaf', fov, ref_file)

    # align brain-to-brain to get acpc.mat (for monkey)
    if acpc_target == 'brain':
        preproc.connect(inputnode, 'anat_brain', fov, in_file)
        if config.anatomical_preproc['acpc_alignment']['FOV_crop'] == 'flirt':
            preproc.connect(inputnode, 'anat_brain', fov, ref_file)

    convert_fov_xfm = pe.Node(interface=fsl_utils.ConvertXFM(),
                              name='anat_acpc_2_fov_convertxfm')
    convert_fov_xfm.inputs.invert_xfm = True

    preproc.connect(fov, fov_mtx,
                    convert_fov_xfm, 'in_file')

    align = pe.Node(interface=fsl.FLIRT(),
                    name='anat_acpc_3_flirt')
    align.inputs.interp = 'spline'
    align.inputs.searchr_x = [30, 30]
    align.inputs.searchr_y = [30, 30]
    align.inputs.searchr_z = [30, 30]

    preproc.connect(fov, fov_outfile, align, 'in_file')

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
    preproc.connect(aff_to_rig, 'out_mat', output_node, 'from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm')

    apply_xfm = pe.Node(interface=fsl.ApplyWarp(),
                        name='anat_acpc_6_applywarp')
    apply_xfm.inputs.interp = 'spline'
    apply_xfm.inputs.relwarp = True

    preproc.connect(inputnode, 'anat_leaf', apply_xfm, 'in_file')
    preproc.connect(inputnode, 'template_head_for_acpc', apply_xfm,
                    'ref_file')
    preproc.connect(aff_to_rig, 'out_mat', apply_xfm, 'premat')
    preproc.connect(apply_xfm, 'out_file', output_node, 'acpc_aligned_head')

    if acpc_target == 'brain':
        apply_xfm_brain = pe.Node(interface=fsl.ApplyWarp(),
                        name='anat_acpc_brain_6_applywarp')
        apply_xfm_brain.inputs.interp = 'spline'
        apply_xfm_brain.inputs.relwarp = True

        preproc.connect(inputnode, 'anat_brain', apply_xfm_brain, 'in_file')
        preproc.connect(inputnode, 'template_brain_for_acpc', apply_xfm_brain,
                        'ref_file')
        preproc.connect(aff_to_rig, 'out_mat', apply_xfm_brain, 'premat')
        preproc.connect(apply_xfm_brain, 'out_file', output_node, 'acpc_aligned_brain')

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


def T2wToT1wReg(wf_name='T2w_to_T1w_reg'):
   
    # Adapted from DCAN lab
    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/PreFreeSurfer/scripts/T2wToT1wReg.sh

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['T1w',
                                                       'T1w_brain',
                                                       'T2w',
                                                       'T2w_brain']),
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['T2w_to_T1w']),
                          name='outputspec')

    # ${FSLDIR}/bin/epi_reg --epi="$T2wImageBrain" --t1="$T1wImage" --t1brain="$WD"/"$T1wImageBrainFile" --out="$WD"/T2w2T1w
    T2w2T1w_reg = pe.Node(interface=fsl.EpiReg(),
                                  name='T2w2T1w_reg')
    T2w2T1w_reg.inputs.out_base = 'T2w2T1w'

    preproc.connect(inputnode, 'T2w_brain', T2w2T1w_reg ,'epi')
    preproc.connect(inputnode, 'T1w', T2w2T1w_reg ,'t1_head')
    preproc.connect(inputnode, 'T1w_brain', T2w2T1w_reg ,'t1_brain')

    # ${FSLDIR}/bin/applywarp --rel --interp=spline --in="$T2wImage" --ref="$T1wImage" --premat="$WD"/T2w2T1w.mat --out="$WD"/T2w2T1w
    T2w2T1w = pe.Node(interface=fsl.ApplyWarp(),
                        name='T2w2T1w_applywarp')
    T2w2T1w.inputs.interp = 'spline'
    T2w2T1w.inputs.relwarp = True

    preproc.connect(inputnode, 'T2w', T2w2T1w, 'in_file')
    preproc.connect(inputnode, 'T1w', T2w2T1w, 'ref_file')
    preproc.connect(T2w2T1w_reg, 'epi2str_mat', T2w2T1w, 'premat')

    # ${FSLDIR}/bin/fslmaths "$WD"/T2w2T1w -add 1 "$WD"/T2w2T1w -odt float
    T2w2T1w_final = pe.Node(interface=fsl.ImageMaths(),
                                  name='T2w2T1w_final')
    T2w2T1w_final.inputs.op_string = "-add 1" 

    preproc.connect(T2w2T1w, 'out_file', T2w2T1w_final, 'in_file')
    preproc.connect(T2w2T1w_final, 'out_file', outputnode, 'T2w_to_T1w')

    return preproc


def BiasFieldCorrection_sqrtT1wXT1w(config=None, wf_name='biasfield_correction_t1t2'):
   
    # Adapted from DCAN lab
    # https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/PreFreeSurfer/scripts/BiasFieldCorrection_sqrtT1wXT1w.sh

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['T1w',
                                                       'T1w_brain',
                                                       'T2w']),
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['T1w_biascorrected',
                                                       'T1w_brain_biascorrected',
                                                       'T2w_biascorrected',
                                                       'T2w_brain_biascorrected',
                                                        'biasfield']),
                          name='outputspec')

    # 1. Form sqrt(T1w*T2w), mask this and normalise by the mean
    # ${FSLDIR}/bin/fslmaths $T1wImage -mul $T2wImage -abs -sqrt $WD/T1wmulT2w.nii.gz -odt float
    T1wmulT2w = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='T1wmulT2w')
    T1wmulT2w.inputs.op_string = "-mul %s -abs -sqrt"
    
    preproc.connect(inputnode, 'T1w', T1wmulT2w, 'in_file')
    preproc.connect(inputnode, 'T2w', T1wmulT2w, 'operand_files')

    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w.nii.gz -mas $T1wImageBrain $WD/T1wmulT2w_brain.nii.gz
    T1wmulT2w_brain = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='T1wmulT2w_brain')
    T1wmulT2w_brain.inputs.op_string = "-mas %s "

    preproc.connect(T1wmulT2w, 'out_file', T1wmulT2w_brain, 'in_file')
    preproc.connect(inputnode, 'T1w_brain', T1wmulT2w_brain, 'operand_files')

    # meanbrainval=`${FSLDIR}/bin/fslstats $WD/T1wmulT2w_brain.nii.gz -M`
    meanbrainval = pe.Node(interface=fsl.ImageStats(),
                             name='image_stats',
                             iterfield=['in_file'])
    meanbrainval.inputs.op_string = '-M'

    preproc.connect(T1wmulT2w_brain, 'out_file', meanbrainval, 'in_file')

    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w_brain.nii.gz -div $meanbrainval $WD/T1wmulT2w_brain_norm.nii.gz
    T1wmulT2w_brain_norm = pe.Node(interface=fsl.ImageMaths(),
                                  name='T1wmulT2w_brain_norm')
    
    def form_meanbrainval_string(meanbrainval):
        return '-div %f' % (meanbrainval)

    preproc.connect(T1wmulT2w_brain, 'out_file', T1wmulT2w_brain_norm, 'in_file')
    preproc.connect(meanbrainval, ('out_stat', form_meanbrainval_string), 
                    T1wmulT2w_brain_norm, 'op_string')

    # 2. Smooth the normalised sqrt image, using within-mask smoothing : s(Mask*X)/s(Mask)
    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -bin -s $BiasFieldSmoothingSigma $WD/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz
    SmoothNorm = pe.Node(interface=fsl.ImageMaths(),
                                  name='SmoothNorm')
    SmoothNorm.inputs.op_string = "-bin -s %f" % (config.anatomical_preproc['t1t2_bias_field_correction']['BiasFieldSmoothingSigma'])

    preproc.connect(T1wmulT2w_brain_norm, 'out_file', SmoothNorm, 'in_file')

    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -s $BiasFieldSmoothingSigma -div $WD/SmoothNorm_s${BiasFieldSmoothingSigma}.nii.gz $WD/T1wmulT2w_brain_norm_s${BiasFieldSmoothingSigma}.nii.gz
    def T1wmulT2w_brain_norm_s_string(sigma, in_file):
        return "-s %f -div %s" %(sigma, in_file)

    T1wmulT2w_brain_norm_s_string = pe.Node(util.Function(input_names=['sigma', 'in_file'],
                                      output_names=['out_str'],
                                      function=T1wmulT2w_brain_norm_s_string),
                                      name='T1wmulT2w_brain_norm_s_string')
    T1wmulT2w_brain_norm_s_string.inputs.sigma = config.anatomical_preproc['t1t2_bias_field_correction']['BiasFieldSmoothingSigma']

    preproc.connect(SmoothNorm, 'out_file', T1wmulT2w_brain_norm_s_string, 'in_file')
    
    T1wmulT2w_brain_norm_s = pe.Node(interface=fsl.ImageMaths(),
                                  name='T1wmulT2w_brain_norm_s')
    
    preproc.connect(T1wmulT2w_brain_norm, 'out_file', T1wmulT2w_brain_norm_s, 'in_file')
    preproc.connect(T1wmulT2w_brain_norm_s_string, 'out_str', T1wmulT2w_brain_norm_s, 'op_string')

    # 3. Divide normalised sqrt image by smoothed version (to do simple bias correction)
    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -div $WD/T1wmulT2w_brain_norm_s$BiasFieldSmoothingSigma.nii.gz $WD/T1wmulT2w_brain_norm_modulate.nii.gz
    T1wmulT2w_brain_norm_modulate = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='T1wmulT2w_brain_norm_modulate')
    T1wmulT2w_brain_norm_modulate.inputs.op_string = "-div %s" 

    preproc.connect(T1wmulT2w_brain_norm, 'out_file', T1wmulT2w_brain_norm_modulate, 'in_file')
    preproc.connect(T1wmulT2w_brain_norm_s, 'out_file', T1wmulT2w_brain_norm_modulate, 'operand_files')

    # 4. Create a mask using a threshold at Mean - 0.5*Stddev, with filling of holes to remove any non-grey/white tissue.
    # STD=`${FSLDIR}/bin/fslstats $WD/T1wmulT2w_brain_norm_modulate.nii.gz -S`
    STD = pe.Node(interface=fsl.ImageStats(),
                             name='STD',
                             iterfield=['in_file'])
    STD.inputs.op_string = '-S'

    preproc.connect(T1wmulT2w_brain_norm_modulate, 'out_file', STD, 'in_file')

    # MEAN=`${FSLDIR}/bin/fslstats $WD/T1wmulT2w_brain_norm_modulate.nii.gz -M`
    MEAN = pe.Node(interface=fsl.ImageStats(),
                             name='MEAN',
                             iterfield=['in_file'])
    MEAN.inputs.op_string = '-M'

    preproc.connect(T1wmulT2w_brain_norm_modulate, 'out_file', MEAN, 'in_file')
    
    # Lower=`echo "$MEAN - ($STD * $Factor)" | bc -l`
    def form_lower_string(mean, std):
        Factor = 0.5 #Leave this at 0.5 for now it is the number of standard deviations below the mean to threshold the non-brain tissues at
        lower = str(float(mean)-(float(std)*float(Factor)))
        return '-thr %s -bin -ero -mul 255' % (lower)

    form_lower_string = pe.Node(util.Function(input_names=['mean', 'std'],
                                      output_names=['out_str'],
                                      function=form_lower_string),
                                      name='form_lower_string')

    preproc.connect(MEAN, 'out_stat', form_lower_string, 'mean')
    preproc.connect(STD, 'out_stat', form_lower_string, 'std')

    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w_brain_norm_modulate -thr $Lower -bin -ero -mul 255 $WD/T1wmulT2w_brain_norm_modulate_mask
    T1wmulT2w_brain_norm_modulate_mask = pe.Node(interface=fsl.ImageMaths(),
                                                name='T1wmulT2w_brain_norm_modulate_mask')

    preproc.connect(T1wmulT2w_brain_norm_modulate, 'out_file', T1wmulT2w_brain_norm_modulate_mask, 'in_file')
    preproc.connect(form_lower_string, 'out_str', T1wmulT2w_brain_norm_modulate_mask, 'op_string')

    # ${CARET7DIR}/wb_command -volume-remove-islands $WD/T1wmulT2w_brain_norm_modulate_mask.nii.gz $WD/T1wmulT2w_brain_norm_modulate_mask.nii.gz
    T1wmulT2w_brain_norm_modulate_mask_roi = pe.Node(interface=VolumeRemoveIslands(),
                                                    name='remove_islands')

    preproc.connect(T1wmulT2w_brain_norm_modulate_mask, 'out_file', T1wmulT2w_brain_norm_modulate_mask_roi, 'in_file')

    # 5. Extrapolate normalised sqrt image from mask region out to whole FOV
    # ${FSLDIR}/bin/fslmaths $WD/T1wmulT2w_brain_norm.nii.gz -mas $WD/T1wmulT2w_brain_norm_modulate_mask.nii.gz -dilall $WD/bias_raw.nii.gz -odt float
    bias_raw = pe.Node(interface=fsl.MultiImageMaths(),
                        name='bias_raw')
    bias_raw.inputs.op_string = "-mas %s -dilall "

    preproc.connect(T1wmulT2w_brain_norm, 'out_file', bias_raw, 'in_file')
    preproc.connect(T1wmulT2w_brain_norm_modulate_mask_roi, 'out_file', bias_raw, 'operand_files')

    # ${FSLDIR}/bin/fslmaths $WD/bias_raw.nii.gz -s $BiasFieldSmoothingSigma $OutputBiasField
    OutputBiasField = pe.Node(interface=fsl.ImageMaths(),
                                  name='OutputBiasField')
    OutputBiasField.inputs.op_string = "-s %f " % (config.anatomical_preproc['t1t2_bias_field_correction']['BiasFieldSmoothingSigma'])

    preproc.connect(bias_raw, 'out_file', OutputBiasField, 'in_file')

    # 6. Use bias field output to create corrected images
    def file_to_a_list(infile_1, infile_2):
        return list([infile_1,infile_2])
    
    file_to_a_list = pe.Node(util.Function(input_names=['infile_1', 'infile_2'],
                                      output_names=['out_list'],
                                      function=file_to_a_list),
                                      name='file_to_a_list')

    preproc.connect(OutputBiasField, 'out_file', file_to_a_list, 'infile_1')
    preproc.connect(inputnode, 'T1w_brain', file_to_a_list, 'infile_2')

    # ${FSLDIR}/bin/fslmaths $T1wImage -div $OutputBiasField -mas $T1wImageBrain $OutputT1wRestoredBrainImage -odt float
    OutputT1wRestoredBrainImage = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='OutputT1wRestoredBrainImage')
    OutputT1wRestoredBrainImage.inputs.op_string = "-div %s -mas %s " 

    preproc.connect(inputnode, 'T1w', OutputT1wRestoredBrainImage, 'in_file')
    preproc.connect(file_to_a_list,'out_list',OutputT1wRestoredBrainImage, 'operand_files')
    
    # ${FSLDIR}/bin/fslmaths $T1wImage -div $OutputBiasField $OutputT1wRestoredImage -odt float
    OutputT1wRestoredImage = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='OutputT1wRestoredImage')
    OutputT1wRestoredImage.inputs.op_string = "-div %s "

    preproc.connect(inputnode, 'T1w', OutputT1wRestoredImage, 'in_file')
    preproc.connect(OutputBiasField, 'out_file', OutputT1wRestoredImage, 'operand_files')

    # ${FSLDIR}/bin/fslmaths $T2wImage -div $OutputBiasField -mas $T1wImageBrain $OutputT2wRestoredBrainImage -odt float
    OutputT2wRestoredBrainImage = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='OutputT2wRestoredBrainImage')
    OutputT2wRestoredBrainImage.inputs.op_string = "-div %s -mas %s " 
    
    preproc.connect(inputnode, 'T2w', OutputT2wRestoredBrainImage, 'in_file')
    preproc.connect(file_to_a_list,'out_list',OutputT2wRestoredBrainImage, 'operand_files')

    # ${FSLDIR}/bin/fslmaths $T2wImage -div $OutputBiasField $OutputT2wRestoredImage -odt float
    OutputT2wRestoredImage = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='OutputT2wRestoredImage')
    OutputT2wRestoredImage.inputs.op_string = "-div %s "

    preproc.connect(inputnode, 'T2w', OutputT2wRestoredImage, 'in_file')
    preproc.connect(OutputBiasField, 'out_file', OutputT2wRestoredImage, 'operand_files')

    preproc.connect(OutputT1wRestoredImage, 'out_file', outputnode, 'T1w_biascorrected')
    preproc.connect(OutputT1wRestoredBrainImage, 'out_file', outputnode, 'T1w_brain_biascorrected')
    preproc.connect(OutputT2wRestoredImage, 'out_file', outputnode, 'T2w_biascorrected')
    preproc.connect(OutputT2wRestoredBrainImage, 'out_file', outputnode, 'T2w_brain_biascorrected')
    preproc.connect(OutputBiasField, 'out_file', outputnode, 'biasfield')

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
    
    if strat_pool.check_rpool('desc-preproc_T1w'): 
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, anat_skullstrip, 'in_file')

    elif strat_pool.check_rpool('desc-preproc_T2w'): 
        node, out = strat_pool.get_data('desc-preproc_T2w')
        wf.connect(node, out, anat_skullstrip, 'in_file')

    wf.connect(skullstrip_args, 'expr', anat_skullstrip, 'args')

    # Generate anatomical brain mask
    anat_brain_mask = pe.Node(interface=afni.Calc(),
                              name=f'anat_brain_mask_{pipe_num}')

    anat_brain_mask.inputs.expr = 'step(a)'
    anat_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(anat_skullstrip, 'out_file',
               anat_brain_mask, 'in_file_a')

    if strat_pool.check_rpool('desc-preproc_T1w'): 
        outputs = {
            'space-T1w_desc-brain_mask': (anat_brain_mask, 'out_file')
        }

    elif strat_pool.check_rpool('desc-preproc_T2w'):
        outputs = {
            'space-T2w_desc-brain_mask': (anat_brain_mask, 'out_file')
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
        mask_boolean= True,
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
    
    anat_robustfov = pe.Node(
        interface=fsl.RobustFOV(), name=f'anat_RobustFOV_{pipe_num}')

    anat_robustfov.inputs.output_type = 'NIFTI_GZ'

    anat_pad_RobustFOV_cropped = pe.Node(util.Function(input_names=['cropped_image_path', 'target_image_path'],
                                    output_names=['padded_image_path'],
                                    function=pad),
                        name=f'anat_pad_mask_{pipe_num}'
                        )

    if strat_pool.check_rpool('desc-preproc_T1w'): 
        node, out = strat_pool.get_data('desc-preproc_T1w')
        if cfg.anatomical_preproc['brain_extraction']['FSL-BET']['Robustfov']:
           wf.connect(node, out, anat_robustfov, 'in_file')
           wf.connect(node, out, anat_pad_RobustFOV_cropped, 'target_image_path')
           wf.connect(anat_robustfov, 'out_roi', anat_pad_RobustFOV_cropped, 'cropped_image_path')
           wf.connect(anat_pad_RobustFOV_cropped, 'padded_image_path', anat_skullstrip,'in_file')
        else :
           wf.connect(node, out, anat_skullstrip, 'in_file')

    elif strat_pool.check_rpool('desc-preproc_T2w'):
        node, out = strat_pool.get_data('desc-preproc_T2w')
        if cfg.anatomical_preproc['brain_extraction']['FSL-BET']['Robustfov']:
           wf.connect(node, out, anat_robustfov, 'in_file')
           wf.connect(node, out, anat_pad_RobustFOV_cropped, 'target_image_path')
           wf.connect(anat_robustfov, 'out_roi', anat_pad_RobustFOV_cropped, 'cropped_image_path')
           wf.connect(anat_pad_RobustFOV_cropped, 'padded_image_path', anat_skullstrip,'in_file')
        else :
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

    if strat_pool.check_rpool('desc-preproc_T1w'): 
        outputs = {
            'space-T1w_desc-brain_mask': (anat_skullstrip, 'mask_file')
        }

    elif strat_pool.check_rpool('desc-preproc_T2w'):
        outputs = {
            'space-T2w_desc-brain_mask': (anat_skullstrip, 'mask_file')
        }

    return (wf, outputs)


def niworkflows_ants_brain_connector(wf, cfg, strat_pool, pipe_num, opt):
    # Skull-stripping using niworkflows-ants
    anat_skullstrip_ants = init_brain_extraction_wf(
        tpl_target_path=cfg.anatomical_preproc['brain_extraction'][
                                               'niworkflows-ants'][
                                               'template_path'],
        tpl_mask_path=cfg.anatomical_preproc['brain_extraction'][
                                             'niworkflows-ants'][
                                             'mask_path'],
        tpl_regmask_path=cfg.anatomical_preproc['brain_extraction'][
                                                'niworkflows-ants'][
                                                'regmask_path'],
        name='anat_skullstrip_ants',
        atropos_use_random_seed=cfg.pipeline_setup['system_config'][
            'random_seed'] is None)

    if strat_pool.check_rpool('desc-preproc_T1w'): 
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, anat_skullstrip_ants, 'inputnode.in_files')

    elif strat_pool.check_rpool('desc-preproc_T2w'):
        node, out = strat_pool.get_data('desc-preproc_T2w')
        wf.connect(node, out, anat_skullstrip_ants, 'inputnode.in_files')
    
    if strat_pool.check_rpool('desc-preproc_T1w'): 
        outputs = {
            'space-T1w_desc-brain_mask': (anat_skullstrip_ants, 'atropos_wf.copy_xform.out_mask'),
            'desc-preproc_T1w': (anat_skullstrip_ants, 'copy_xform.out_file')
        }

    elif strat_pool.check_rpool('desc-preproc_T2w'):
        outputs = {
            'space-T2w_desc-brain_mask': (anat_skullstrip_ants, 'atropos_wf.copy_xform.out_mask'),
            'desc-preproc_T2w': (anat_skullstrip_ants, 'copy_xform.out_file')
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
    from CPAC.unet.function import predict_volumes
    unet_mask = pe.Node(util.Function(input_names=['model_path', 'cimg_in'],
                                      output_names=['out_path'],
                                      function=predict_volumes),
                        name=f'unet_mask_{pipe_num}')

    node, out = strat_pool.get_data('unet-model')
    wf.connect(node, out, unet_mask, 'model_path')

    if strat_pool.check_rpool('desc-preproc_T1w'):
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, unet_mask, 'cimg_in')

    elif strat_pool.check_rpool('desc-preproc_T2w'):
        node, out = strat_pool.get_data('desc-preproc_T2w')
        wf.connect(node, out, unet_mask, 'cimg_in')

    """
    Revised mask with ANTs
    """
    # fslmaths <whole head> -mul <mask> brain.nii.gz
    unet_masked_brain = pe.Node(interface=fsl.MultiImageMaths(),
                                name=f'unet_masked_brain_{pipe_num}')
    unet_masked_brain.inputs.op_string = "-mul %s"

    if strat_pool.check_rpool('desc-preproc_T1w'):
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, unet_masked_brain, 'in_file')
        
    elif strat_pool.check_rpool('desc-preproc_T2w'):
        node, out = strat_pool.get_data('desc-preproc_T2w')
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

    node, out = strat_pool.get_data('T1w-brain-template')
    wf.connect(node, out, native_brain_to_template_brain, 'reference')

    # flirt -in head.nii.gz -ref NMT_0.5mm.nii.gz -o head_rot2atl -applyxfm -init brain_rot2atl.mat
    native_head_to_template_head = pe.Node(interface=fsl.FLIRT(),
                                           name=f'native_head_to_template_'
                                                f'head_{pipe_num}')
    native_head_to_template_head.inputs.apply_xfm = True

    if strat_pool.check_rpool('desc-preproc_T1w'):
        node, out = strat_pool.get_data('desc-preproc_T1w')
        wf.connect(node, out, native_head_to_template_head, 'in_file')
        
    elif strat_pool.check_rpool('desc-preproc_T2w'):
        node, out = strat_pool.get_data('desc-preproc_T2w')
        wf.connect(node, out, native_head_to_template_head, 'in_file')

    wf.connect(native_brain_to_template_brain, 'out_matrix_file',
               native_head_to_template_head, 'in_matrix_file')

    node, out = strat_pool.get_data('T1w-template')
    wf.connect(node, out, native_head_to_template_head, 'reference')

    # fslmaths NMT_SS_0.5mm.nii.gz -bin templateMask.nii.gz
    template_brain_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                  name=f'template_brain_mask_{pipe_num}')
    template_brain_mask.inputs.args = '-bin'

    node, out = strat_pool.get_data('T1w-brain-template')
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

    node, out = strat_pool.get_data('T1w-brain-template')
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

    node, out = strat_pool.get_data('pipeline-fs_brainmask') 
    wf.connect(node, out, fs_brain_mask_to_native, 'source_file')

    node, out = strat_pool.get_data('pipeline-fs_raw-average')
    wf.connect(node, out, fs_brain_mask_to_native, 'target_file')

    node, out = strat_pool.get_data('freesurfer-subject-dir')
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
    '''
    ABCD harmonization - anatomical brain mask generation

    Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh#L151-L156
    '''

    wmparc_to_nifti = pe.Node(util.Function(input_names=['in_file',
                                                         'reslice_like',
                                                         'args'],
                                            output_names=['out_file'],
                                            function=mri_convert),
                              name=f'wmparc_to_nifti_{pipe_num}')
    
    # Register wmparc file if ingressing FreeSurfer data
    if strat_pool.check_rpool('pipeline-fs_xfm'):

        wmparc_to_native = pe.Node(util.Function(input_names=['source_file',
                                                            'target_file',
                                                            'xfm',
                                                            'out_file'],
                                                output_names=['transformed_file'],
                                                function=normalize_wmparc),
                                        name=f'wmparc_to_native_{pipe_num}')
        
        wmparc_to_native.inputs.out_file = 'wmparc_warped.mgz'

        node, out = strat_pool.get_data('pipeline-fs_wmparc')
        wf.connect(node, out, wmparc_to_native, 'source_file')

        node, out = strat_pool.get_data('pipeline-fs_raw-average')
        wf.connect(node, out, wmparc_to_native, 'target_file')

        node, out = strat_pool.get_data('pipeline-fs_xfm')
        wf.connect(node, out, wmparc_to_native, 'xfm')

        wf.connect(wmparc_to_native, 'transformed_file', wmparc_to_nifti, 'in_file')
    
    else:
        
        node, out = strat_pool.get_data('pipeline-fs_wmparc')
        wf.connect(node, out, wmparc_to_nifti, 'in_file')

    wmparc_to_nifti.inputs.args = '-rt nearest'

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

    outputs = {
        'space-T1w_desc-brain_mask': (brain_mask_to_t1_restore, 'out_file')
    }

    return (wf, outputs)


def freesurfer_fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt):

    node_id = f'{opt.lower()}_{pipe_num}'

    # mri_convert -it mgz ${SUBJECTS_DIR}/${subject}/mri/brainmask.mgz -ot nii brainmask.nii.gz
    convert_fs_brainmask_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                                   output_names=['out_file'],
                                                   function=mri_convert),
                                            name=f'convert_fs_brainmask_to_nifti_{node_id}')

    node, out = strat_pool.get_data('pipeline-fs_brainmask')
    wf.connect(node, out, convert_fs_brainmask_to_nifti, 'in_file')

    # mri_convert -it mgz ${SUBJECTS_DIR}/${subject}/mri/T1.mgz -ot nii T1.nii.gz
    convert_fs_T1_to_nifti = pe.Node(util.Function(input_names=['in_file'],
                                                   output_names=['out_file'],
                                                   function=mri_convert),
                                     name=f'convert_fs_T1_to_nifti_{node_id}')

    node, out = strat_pool.get_data('pipeline-fs_T1')
    wf.connect(node, out, convert_fs_T1_to_nifti, 'in_file')

    # 3dresample -orient RPI -inset brainmask.nii.gz -prefix brain_fs.nii.gz
    reorient_fs_brainmask = pe.Node(interface=afni.Resample(),
                                    name=f'reorient_fs_brainmask_{node_id}',
                                    mem_gb=0,
                                    mem_x=(0.0115, 'in_file', 't'))
    reorient_fs_brainmask.inputs.orientation = 'RPI'
    reorient_fs_brainmask.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(convert_fs_brainmask_to_nifti, 'out_file', 
        reorient_fs_brainmask, 'in_file')

    # fslmaths brain_fs.nii.gz -abs -bin brain_fs_mask.nii.gz
    binarize_fs_brain = pe.Node(interface=fsl.maths.MathsCommand(),
                                name=f'binarize_fs_brain_{node_id}')
    binarize_fs_brain.inputs.args = '-abs -bin'

    wf.connect(reorient_fs_brainmask, 'out_file',
               binarize_fs_brain, 'in_file')

    # 3dresample -orient RPI -inset T1.nii.gz -prefix head_fs.nii.gz
    reorient_fs_T1 = pe.Node(interface=afni.Resample(),
                             name=f'reorient_fs_T1_{node_id}',
                             mem_gb=0,
                             mem_x=(0.0115, 'in_file', 't'))
    reorient_fs_T1.inputs.orientation = 'RPI'
    reorient_fs_T1.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(convert_fs_T1_to_nifti, 'out_file', 
        reorient_fs_T1, 'in_file')

    # flirt -in head_fs.nii.gz -ref ${FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz \
    # -out tmp_head_fs2standard.nii.gz -omat tmp_head_fs2standard.mat -bins 256 -cost corratio \
    # -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear
    convert_head_to_template = pe.Node(interface=fsl.FLIRT(),
                                       name=f'convert_head_to_template_{node_id}')
    convert_head_to_template.inputs.cost = 'corratio'
    convert_head_to_template.inputs.interp = 'trilinear'
    convert_head_to_template.inputs.bins = 256
    convert_head_to_template.inputs.dof = 12
    convert_head_to_template.inputs.searchr_x = [-90, 90]
    convert_head_to_template.inputs.searchr_y = [-90, 90]
    convert_head_to_template.inputs.searchr_z = [-90, 90]

    wf.connect(reorient_fs_T1, 'out_file',
        convert_head_to_template, 'in_file')

    node, out = strat_pool.get_data('T1w-ACPC-template')
    wf.connect(node, out, convert_head_to_template, 'reference')

    # convert_xfm -omat tmp_standard2head_fs.mat -inverse tmp_head_fs2standard.mat
    convert_xfm = pe.Node(interface=fsl_utils.ConvertXFM(),
                              name=f'convert_xfm_{node_id}')
    convert_xfm.inputs.invert_xfm = True

    wf.connect(convert_head_to_template, 'out_matrix_file',
        convert_xfm, 'in_file')

    # bet tmp_head_fs2standard.nii.gz tmp.nii.gz -f ${bet_thr_tight} -m
    skullstrip = pe.Node(interface=fsl.BET(), 
                         name=f'anat_BET_skullstrip_{node_id}')
    skullstrip.inputs.output_type = 'NIFTI_GZ'
    skullstrip.inputs.mask=True

    if opt == 'FreeSurfer-BET-Tight':
        skullstrip.inputs.frac=0.3
    elif opt == 'FreeSurfer-BET-Loose':
        skullstrip.inputs.frac=0.1

    wf.connect(convert_head_to_template, 'out_file', 
        skullstrip, 'in_file')
    
    # fslmaths tmp_mask.nii.gz -mas ${CCSDIR}/templates/MNI152_T1_1mm_first_brain_mask.nii.gz tmp_mask.nii.gz
    apply_mask = pe.Node(interface=fsl.maths.ApplyMask(),
                         name=f'apply_mask_{node_id}')

    wf.connect(skullstrip, 'out_file',
        apply_mask, 'in_file')

    node, out = strat_pool.get_data('T1w-brain-template-mask-ccs')
    wf.connect(node, out, apply_mask, 'mask_file')

    # flirt -in tmp_mask.nii.gz -applyxfm -init tmp_standard2head_fs.mat -out brain_fsl_mask_tight.nii.gz \
    # -paddingsize 0.0 -interp nearestneighbour -ref head_fs.nii.gz
    convert_template_mask_to_native = pe.Node(interface=fsl.FLIRT(),
                                              name=f'convert_template_mask_to_native_{node_id}')
    convert_template_mask_to_native.inputs.apply_xfm = True
    convert_template_mask_to_native.inputs.padding_size = 0
    convert_template_mask_to_native.inputs.interp = 'nearestneighbour'

    wf.connect(apply_mask, 'out_file',
        convert_template_mask_to_native, 'in_file')

    wf.connect(convert_xfm, 'out_file',
        convert_template_mask_to_native, 'in_matrix_file')

    wf.connect(reorient_fs_T1, 'out_file',
        convert_template_mask_to_native, 'reference')

    # fslmaths brain_fs_mask.nii.gz -add brain_fsl_mask_tight.nii.gz -bin brain_mask_tight.nii.gz
    # BinaryMaths doesn't use -bin! 
    combine_mask = pe.Node(interface=fsl.BinaryMaths(),
                           name=f'combine_mask_{node_id}')

    if opt == 'FreeSurfer-BET-Tight':
        combine_mask.inputs.operation = 'add'
    elif opt == 'FreeSurfer-BET-Loose':
        combine_mask.inputs.operation = 'mul'

    wf.connect(binarize_fs_brain, 'out_file',
        combine_mask, 'in_file')

    wf.connect(convert_template_mask_to_native, 'out_file',
        combine_mask, 'operand_file')

    binarize_combined_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                                     name=f'binarize_combined_mask_{node_id}')
    binarize_combined_mask.inputs.args = '-bin'

    wf.connect(combine_mask, 'out_file',
               binarize_combined_mask, 'in_file')

    # CCS brain mask is in FS space, transfer it back to native T1 space
    fs_fsl_brain_mask_to_native = pe.Node(interface=freesurfer.ApplyVolTransform(),
                                      name=f'fs_fsl_brain_mask_to_native_{node_id}')
    fs_fsl_brain_mask_to_native.inputs.reg_header = True
    fs_fsl_brain_mask_to_native.inputs.interp = 'nearest'

    wf.connect(binarize_combined_mask, 'out_file', 
        fs_fsl_brain_mask_to_native, 'source_file')

    node, out = strat_pool.get_data('pipeline-fs_raw-average')
    wf.connect(node, out, fs_fsl_brain_mask_to_native, 'target_file')

    node, out = strat_pool.get_data('freesurfer-subject-dir')
    wf.connect(node, out, fs_fsl_brain_mask_to_native, 'subjects_dir')

    if opt == 'FreeSurfer-BET-Tight':
        outputs = {
            'space-T1w_desc-tight_brain_mask': (fs_fsl_brain_mask_to_native, 'transformed_file')
        }
    elif opt == 'FreeSurfer-BET-Loose':
        outputs = {
            'space-T1w_desc-loose_brain_mask': (fs_fsl_brain_mask_to_native, 'transformed_file')
        }

    return (wf, outputs)


def mask_T2(wf_name='mask_T2'):
    # create T2 mask based on T1 mask
    # reference https://github.com/DCAN-Labs/dcan-macaque-pipeline/blob/master/PreliminaryMasking/macaque_masking.py
    
    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['T1w',
                                                       'T1w_mask',
                                                       'T2w']),
                        name='inputspec')

    outputnode = pe.Node(util.IdentityInterface(fields=['T1w_brain',
                                                        'T2w_mask',
                                                       'T2w_brain']),
                          name='outputspec')

    # mask_t1w = 'fslmaths {t1w} -mas {t1w_mask_edit} {t1w_brain}'.format(**kwargs)
    mask_t1w = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='mask_t1w')
    mask_t1w.inputs.op_string = "-mas %s "

    preproc.connect(inputnode, 'T1w', mask_t1w, 'in_file')
    preproc.connect(inputnode, 'T1w_mask', mask_t1w, 'operand_files')


    # t1w2t2w_rigid = 'flirt -dof 6 -cost mutualinfo -in {t1w} -ref {t2w} ' \
    #                     '-omat {t1w2t2w}'.format(**kwargs)

    t1w2t2w_rigid = pe.Node(interface=fsl.FLIRT(),
                            name='t1w2t2w_rigid')

    t1w2t2w_rigid.inputs.dof = 6
    t1w2t2w_rigid.inputs.cost = 'mutualinfo'
    preproc.connect(inputnode, 'T1w', t1w2t2w_rigid, 'in_file')
    preproc.connect(inputnode, 'T2w', t1w2t2w_rigid, 'reference')

    # t1w2t2w_mask = 'flirt -in {t1w_mask_edit} -interp nearestneighbour -ref {' \
    #                 't2w} -o {t2w_brain_mask} -applyxfm -init {' \
    #                 't1w2t2w}'.format(**kwargs)
    t1w2t2w_mask = pe.Node(interface=fsl.FLIRT(),
                                    name='t1w2t2w_mask')
    t1w2t2w_mask.inputs.apply_xfm = True
    t1w2t2w_mask.inputs.interp = 'nearestneighbour'

    preproc.connect(inputnode, 'T1w_mask', t1w2t2w_mask, 'in_file')
    preproc.connect(inputnode, 'T2w', t1w2t2w_mask, 'reference')
    preproc.connect(t1w2t2w_rigid, 'out_matrix_file', t1w2t2w_mask, 'in_matrix_file')

    # mask_t2w = 'fslmaths {t2w} -mas {t2w_brain_mask} ' \
    #         '{t2w_brain}'.format(**kwargs)
    mask_t2w = pe.Node(interface=fsl.MultiImageMaths(),
                                  name='mask_t2w')
    mask_t2w.inputs.op_string = "-mas %s "

    preproc.connect(inputnode, 'T2w', mask_t2w, 'in_file')
    preproc.connect(t1w2t2w_mask, 'out_file', mask_t2w, 'operand_files')

    preproc.connect(mask_t1w, 'out_file', outputnode, 'T1w_brain')
    preproc.connect(mask_t2w, 'out_file', outputnode, 'T2w_brain')
    preproc.connect(t1w2t2w_mask, 'out_file', outputnode, 'T2w_mask')

    return preproc


@nodeblock(
    name="anatomical_init",
    config=["anatomical_preproc"],
    switch=["run"],
    inputs=["T1w"],
    outputs=["desc-preproc_T1w", "desc-reorient_T1w", "desc-head_T1w"],
)
def anatomical_init(wf, cfg, strat_pool, pipe_num, opt=None):

    anat_deoblique = pe.Node(interface=afni.Refit(),
                             name=f'anat_deoblique_{pipe_num}')
    anat_deoblique.inputs.deoblique = True

    node, out = strat_pool.get_data('T1w')
    wf.connect(node, out, anat_deoblique, 'in_file')

    anat_reorient = pe.Node(interface=afni.Resample(),
                            name=f'anat_reorient_{pipe_num}',
                            mem_gb=0,
                            mem_x=(0.0115, 'in_file', 't'))
    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(anat_deoblique, 'out_file', anat_reorient, 'in_file')

    outputs = {'desc-preproc_T1w': (anat_reorient, 'out_file'),
               'desc-reorient_T1w': (anat_reorient, 'out_file'),
               'desc-head_T1w': (anat_reorient, 'out_file')}

    return (wf, outputs)

@nodeblock(
    name="acpc_alignment_head",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=["desc-head_T1w", "desc-preproc_T1w", "T1w-ACPC-template"],
    outputs=[
        "desc-head_T1w",
        "desc-preproc_T1w",
        "from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm",
    ],
)
def acpc_align_head(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=False,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data(['desc-preproc_T1w','desc-head_T1w'])
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T1w-ACPC-template')
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    outputs = {
        'desc-head_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm': (
            acpc_align,
            'outputspec.from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm')
    }

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_head_with_mask",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=[
        (
            "desc-head_T1w",
            "desc-preproc_T1w",
            ["space-T1w_desc-brain_mask", "space-T1w_desc-brain_mask"],
        ),
        "T1w-ACPC-template",
        "T1w-brain-ACPC-template",
    ],
    outputs=[
        "desc-head_T1w",
        "desc-preproc_T1w",
        ["space-T1w_desc-brain_mask", "space-T1w_desc-brain_mask"],
        "from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm",
    ],
)
def acpc_align_head_with_mask(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=True,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data(['desc-head_T1w', 'desc-preproc_T1w'])
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T1w-ACPC-template')
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    if strat_pool.check_rpool("space-T1w_desc-brain_mask"):
        node, out = strat_pool.get_data("space-T1w_desc-brain_mask")
        wf.connect(node, out, acpc_align, 'inputspec.brain_mask')

        node, out = strat_pool.get_data('T1w-brain-ACPC-template')
        wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-head_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'space-T1w_desc-brain_mask': (
            acpc_align, 'outputspec.acpc_brain_mask'),
        'from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm': (
            acpc_align, 'outputspec.from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm')
    }

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_brain",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=[
        (
            "desc-preproc_T1w",
            "desc-tempbrain_T1w",
            "T1w-ACPC-template",
            "T1w-brain-ACPC-template",
        )
    ],
    outputs=[
        "desc-preproc_T1w",
        "desc-acpcbrain_T1w",
        "from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm",
    ],
)
def acpc_align_brain(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=False,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('desc-tempbrain_T1w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_brain')

    node, out = strat_pool.get_data('T1w-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    node, out = strat_pool.get_data('T1w-brain-ACPC-template')
    wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'desc-acpcbrain_T1w': (acpc_align, 'outputspec.acpc_aligned_brain'),
        'from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm': (
            acpc_align, 'outputspec.from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm')
    }

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_brain_with_mask",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=[
        ("desc-preproc_T1w", "desc-tempbrain_T1w", "space-T1w_desc-brain_mask"),
        "T1w-ACPC-template",
        "T1w-brain-ACPC-template",
    ],
    outputs=[
        "desc-preproc_T1w",
        "desc-acpcbrain_T1w",
        "space-T1w_desc-brain_mask",
        "space-T1w_desc-prebrain_mask",
    ],
)
def acpc_align_brain_with_mask(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=True,
                                wf_name=f'acpc_align_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('desc-tempbrain_T1w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_brain')

    node, out = strat_pool.get_data("space-T1w_desc-brain_mask")
    wf.connect(node, out, acpc_align, 'inputspec.brain_mask')

    node, out = strat_pool.get_data('T1w-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    node, out = strat_pool.get_data('T1w-brain-ACPC-template')
    wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-preproc_T1w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'desc-acpcbrain_T1w': (acpc_align, 'outputspec.acpc_aligned_brain'),
        'space-T1w_desc-brain_mask': (acpc_align, 'outputspec.acpc_brain_mask'),
        'space-T1w_desc-prebrain_mask': (strat_pool.get_data('space-T1_desc-brain_mask'))
    }

    return (wf, outputs)


@nodeblock(
    name="registration_T2w_to_T1w",
    config=["anatomical_preproc"],
    switch=["run_t2"],
    inputs=[
        (
            "desc-preproc_T1w",
            "desc-preproc_T2w",
            "desc-acpcbrain_T1w",
            "desc-acpcbrain_T2w",
        )
    ],
    outputs=["desc-preproc_T2w"],
)
def registration_T2w_to_T1w(wf, cfg, strat_pool, pipe_num, opt=None):

    T2_to_T1_reg = T2wToT1wReg(wf_name=f'T2w_to_T1w_Reg_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, T2_to_T1_reg, 'inputspec.T1w')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, T2_to_T1_reg, 'inputspec.T2w')

    node, out = strat_pool.get_data(['desc-acpcbrain_T1w'])
    wf.connect(node, out, T2_to_T1_reg, 'inputspec.T1w_brain')

    node, out = strat_pool.get_data(['desc-acpcbrain_T2w'])
    wf.connect(node, out, T2_to_T1_reg, 'inputspec.T2w_brain')

    outputs = {
        'desc-preproc_T2w': (T2_to_T1_reg, 'outputspec.T2w_to_T1w')
    }

    return (wf, outputs)


@nodeblock(
    name="nlm_filtering",
    switch=[
        ["anatomical_preproc", "non_local_means_filtering", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=["desc-preproc_T1w"],
    outputs=["desc-preproc_T1w"],
)
def non_local_means(wf, cfg, strat_pool, pipe_num, opt=None):

    denoise = pe.Node(interface=ants.DenoiseImage(),
                      name=f'anat_denoise_{pipe_num}')

    denoise.inputs.noise_model = cfg.anatomical_preproc['non_local_means_filtering']['noise_model']

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, denoise, 'input_image')

    outputs = {
        'desc-preproc_T1w': (denoise, 'output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="n4_bias_correction",
    switch=[
        ["anatomical_preproc", "n4_bias_field_correction", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=["desc-preproc_T1w"],
    outputs={
        "desc-preproc_T1w": {
            "Description": "T1w image that has been N4-bias-field corrected using ANTs N4BiasFieldCorrection."
        },
        "desc-n4_T1w": {
            "Description": "T1w image that has been N4-bias-field corrected using ANTs N4BiasFieldCorrection."
        },
    },
)
def n4_bias_correction(wf, cfg, strat_pool, pipe_num, opt=None):

    n4 = pe.Node(interface=ants.N4BiasFieldCorrection(dimension=3,
                                                      copy_header=True),
                 name=f'anat_n4_{pipe_num}')
    n4.inputs.shrink_factor = cfg.anatomical_preproc['n4_bias_field_correction']['shrink_factor']

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, n4, 'input_image')

    outputs = {
        'desc-preproc_T1w': (n4, 'output_image'),
        'desc-n4_T1w': (n4, 'output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="t1t2_bias_correction",
    config=["anatomical_preproc", "t1t2_bias_field_correction"],
    switch=["run"],
    inputs=[("desc-preproc_T1w", "desc-preproc_T2w", "desc-acpcbrain_T1w")],
    outputs=[
        "desc-preproc_T1w",
        "desc-brain_T1w",
        "desc-preproc_T2w",
        "desc-brain_T2w",
        "desc-biasfield_T1wT2w",
    ],
)
def t1t2_bias_correction(wf, cfg, strat_pool, pipe_num, opt=None):

    t1t2_bias_correction = BiasFieldCorrection_sqrtT1wXT1w(config=cfg, wf_name=f't1t2_bias_correction_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, t1t2_bias_correction, 'inputspec.T1w')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, t1t2_bias_correction, 'inputspec.T2w')

    node, out = strat_pool.get_data("desc-acpcbrain_T1w")
    wf.connect(node, out, t1t2_bias_correction, 'inputspec.T1w_brain')

    outputs = {
        'desc-preproc_T1w': (t1t2_bias_correction, 'outputspec.T1w_biascorrected'),
        'desc-brain_T1w': (t1t2_bias_correction, 'outputspec.T1w_brain_biascorrected'),
        'desc-preproc_T2w': (t1t2_bias_correction, 'outputspec.T2w_biascorrected'),
        'desc-brain_T2w': (t1t2_bias_correction, 'outputspec.T2w_brain_biascorrected'),
        'desc-biasfield_T1wT2w': (t1t2_bias_correction, 'outputspec.biasfield'),
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_afni",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="3dSkullStrip",
    inputs=["desc-preproc_T1w"],
    outputs=["space-T1w_desc-brain_mask"],
)
def brain_mask_afni(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_afni",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="3dSkullStrip",
    inputs=["desc-preproc_T1w"],
    outputs=["space-T1w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_afni(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_fsl",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="BET",
    inputs=["desc-preproc_T1w"],
    outputs=["space-T1w_desc-brain_mask"],
)
def brain_mask_fsl(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_fsl",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="BET",
    inputs=["desc-preproc_T1w"],
    outputs=["space-T1w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_fsl(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_niworkflows_ants",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="niworkflows-ants",
    inputs=["desc-preproc_T1w"],
    outputs=["space-T1w_desc-brain_mask", "desc-preproc_T1w"],
)
def brain_mask_niworkflows_ants(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = niworkflows_ants_brain_connector(wf, cfg, strat_pool,
                                                   pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_niworkflows_ants",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="niworkflows-ants",
    inputs=["desc-preproc_T1w"],
    outputs=["space-T1w_desc-acpcbrain_mask", "desc-preproc_T1w"],
)
def brain_mask_acpc_niworkflows_ants(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = niworkflows_ants_brain_connector(wf, cfg, strat_pool,
                                                      pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask'],
        'desc-preproc_T1w':
            wf_outputs['desc-preproc_T1w']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_unet",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="UNet",
    inputs=["desc-preproc_T1w", "T1w-brain-template", "T1w-template", "unet-model"],
    outputs=["space-T1w_desc-brain_mask"],
)
def brain_mask_unet(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_unet",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="UNet",
    inputs=["desc-preproc_T1w", "T1w-brain-template", "T1w-template", "unet-model"],
    outputs=["space-T1w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_unet(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T1w_desc-acpcbrain_mask':
            wf_outputs['space-T1w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_freesurfer",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-Brainmask",
    inputs=[
        "pipeline-fs_raw-average",
        "pipeline-fs_brainmask",
        "freesurfer-subject-dir",
    ],
    outputs=["space-T1w_desc-brain_mask"],
)
def brain_mask_freesurfer(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = freesurfer_brain_connector(wf, cfg, strat_pool, pipe_num,
                                             opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_freesurfer",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-Brainmask",
    inputs=[
        "space-T1w_desc-brain_mask",
        "pipeline-fs_raw-average",
        "freesurfer-subject-dir",
    ],
    outputs=["space-T1w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_freesurfer(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = freesurfer_brain_connector(wf, cfg, strat_pool, pipe_num,
                                                opt)

    outputs = {'space-T1w_desc-acpcbrain_mask':
                   wf_outputs['space-T1w_desc-brain_mask']}

    return (wf, outputs)


@nodeblock(
    name="brain_mask_freesurfer_abcd",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-ABCD",
    inputs=[
        "desc-preproc_T1w",
        "pipeline-fs_wmparc",
        "pipeline-fs_raw-average",
        "pipeline-fs_xfm",
        "freesurfer-subject-dir",
    ],
    outputs=["space-T1w_desc-brain_mask"],
)
def brain_mask_freesurfer_abcd(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = freesurfer_abcd_brain_connector(wf, cfg, strat_pool,
                                                  pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_freesurfer_fsl_tight",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-BET-Tight",
    inputs=[
        "pipeline-fs_brainmask",
        "pipeline-fs_T1",
        "pipeline-fs_raw-average",
        "freesurfer-subject-dir",
        "T1w-brain-template-mask-ccs",
        "T1w-ACPC-template",
    ],
    outputs=["space-T1w_desc-tight_brain_mask"],
)
def brain_mask_freesurfer_fsl_tight(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = freesurfer_fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_freesurfer_abcd",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-ABCD",
    inputs=[
        "desc-preproc_T1w",
        "pipeline-fs_wmparc",
        "pipeline-fs_raw-average",
        "pipeline-fs_xfm",
        "freesurfer-subject-dir",
    ],
    outputs=["space-T1w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_freesurfer_abcd(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = freesurfer_abcd_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {'space-T1w_desc-acpcbrain_mask':
                wf_outputs['space-T1w_desc-brain_mask']}

    return (wf, outputs)


@nodeblock(
    name="brain_mask_freesurfer_fsl_loose",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-BET-Loose",
    inputs=[
        "pipeline-fs_brainmask",
        "pipeline-fs_T1",
        "pipeline-fs_raw-average",
        "freesurfer-subject-dir",
        "T1w-brain-template-mask-ccs",
        "T1w-ACPC-template",
    ],
    outputs=["space-T1w_desc-loose_brain_mask"],
)
def brain_mask_freesurfer_fsl_loose(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = freesurfer_fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_freesurfer_fsl_tight",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-BET-Tight",
    inputs=[
        "pipeline-fs_brainmask",
        "pipeline-fs_T1",
        "T1w-brain-template-mask-ccs",
        "T1w-ACPC-template",
    ],
    outputs=["space-T1w_desc-tight_acpcbrain_mask"],
)
def brain_mask_acpc_freesurfer_fsl_tight(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = freesurfer_fsl_brain_connector(wf, cfg, strat_pool,
                                                    pipe_num, opt)

    outputs = {'space-T1w_desc-tight_acpcbrain_mask':
        wf_outputs['space-T1w_desc-tight_brain_mask']}

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_freesurfer_fsl_loose",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    option_key=["anatomical_preproc", "brain_extraction", "using"],
    option_val="FreeSurfer-BET-Loose",
    inputs=[
        "pipeline-fs_brainmask",
        "pipeline-fs_T1",
        "T1w-brain-template-mask-ccs",
        "T1w-ACPC-template",
    ],
    outputs=["space-T1w_desc-loose_acpcbrain_mask"],
)
def brain_mask_acpc_freesurfer_fsl_loose(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = freesurfer_fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {'space-T1w_desc-loose_acpcbrain_mask':
        wf_outputs['space-T1w_desc-loose_brain_mask']}

    return (wf, outputs)


@nodeblock(
    name="brain_extraction",
    switch=[
        ["anatomical_preproc", "brain_extraction", "run"],
        ["anatomical_preproc", "run"],
    ],
    inputs=[
        (
            "desc-head_T1w",
            "desc-preproc_T1w",
            ["space-T1w_desc-brain_mask", "space-T1w_desc-acpcbrain_mask"],
        )
    ],
    outputs={
        "desc-preproc_T1w": {"SkullStripped": "True"},
        "desc-brain_T1w": {"SkullStripped": "True"},
        "desc-head_T1w": {"SkullStripped": "False"},
    },
)
def brain_extraction(wf, cfg, strat_pool, pipe_num, opt=None):

    '''
    brain_mask_deoblique = pe.Node(interface=afni.Refit(),
                                   name='brain_mask_deoblique')
    brain_mask_deoblique.inputs.deoblique = True
    wf.connect(inputnode, 'brain_mask',
                    brain_mask_deoblique, 'in_file')

    brain_mask_reorient = pe.Node(interface=afni.Resample(),
                                  name='brain_mask_reorient',
                                  mem_gb=0,
                                  mem_x=(0.0115, 'in_file', 't'))
    brain_mask_reorient.inputs.orientation = 'RPI'
    brain_mask_reorient.inputs.outputtype = 'NIFTI_GZ'
    wf.connect(brain_mask_deoblique, 'out_file',
                    brain_mask_reorient, 'in_file')
    '''

    anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                       name=f'brain_extraction_{pipe_num}')

    anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
    anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

    node_T1w, out_T1w = strat_pool.get_data('desc-head_T1w')
    wf.connect(node_T1w, out_T1w, anat_skullstrip_orig_vol, 'in_file_a')

    node, out = strat_pool.get_data(['space-T1w_desc-brain_mask',
                                     'space-T1w_desc-acpcbrain_mask'])
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_b')

    outputs = {
        'desc-preproc_T1w': (anat_skullstrip_orig_vol, 'out_file'),
        'desc-brain_T1w': (anat_skullstrip_orig_vol, 'out_file'),
        'desc-head_T1w': (node_T1w, out_T1w)
    }

    return (wf, outputs)


@nodeblock(
    name="brain_extraction_temp",
    inputs=[
        (
            "desc-preproc_T1w",
            ["space-T1w_desc-brain_mask", "space-T1w_desc-acpcbrain_mask"],
        )
    ],
    outputs={
        "desc-preproc_T1w": {"SkullStripped": "True"},
        "desc-tempbrain_T1w": {"SkullStripped": "True"},
    },
)
def brain_extraction_temp(wf, cfg, strat_pool, pipe_num, opt=None):

    anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                       name=f'brain_extraction_temp_{pipe_num}')

    anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
    anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_a')

    node, out = strat_pool.get_data(['space-T1w_desc-brain_mask',
                                     'space-T1w_desc-acpcbrain_mask'])
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_b')

    outputs = {
        'desc-preproc_T1w': (anat_skullstrip_orig_vol, 'out_file'),
        'desc-tempbrain_T1w': (anat_skullstrip_orig_vol, 'out_file')
    }

    return (wf, outputs)


@nodeblock(
    name="anatomical_init_T2",
    config=["anatomical_preproc"],
    switch=["run_t2"],
    inputs=["T2w"],
    outputs=["desc-preproc_T2w", "desc-reorient_T2w", "desc-head_T2w"],
)
def anatomical_init_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    T2_deoblique = pe.Node(interface=afni.Refit(),
                             name=f'T2_deoblique_{pipe_num}')
    T2_deoblique.inputs.deoblique = True

    node, out = strat_pool.get_data('T2w')
    wf.connect(node, out, T2_deoblique, 'in_file')

    T2_reorient = pe.Node(interface=afni.Resample(),
                          name=f'T2_reorient_{pipe_num}',
                          mem_gb=0,
                          mem_x=(0.0115, 'in_file', 't'))
    T2_reorient.inputs.orientation = 'RPI'
    T2_reorient.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(T2_deoblique, 'out_file', T2_reorient, 'in_file')

    outputs = {'desc-preproc_T2w': (T2_reorient, 'out_file'),
               'desc-reorient_T2w': (T2_reorient, 'out_file'),
               'desc-head_T2w': (T2_reorient, 'out_file')}

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_head_T2",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run_t2"],
    ],
    inputs=["desc-preproc_T2w", "T2w-ACPC-template"],
    outputs=["desc-preproc_T2w"],
)
def acpc_align_head_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=False,
                                wf_name=f'acpc_align_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T2w-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    outputs = {
        'desc-preproc_T2w': (acpc_align, 'outputspec.acpc_aligned_head')
    }

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_head_with_mask_T2",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run_t2"],
    ],
    inputs=[("desc-preproc_T2w", "space-T2w_desc-brain_mask"), "T2w-ACPC-template"],
    outputs=["desc-preproc_T2w", "space-T2w_desc-brain_mask"],
)
def acpc_align_head_with_mask_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=True,
                                wf_name=f'acpc_align_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('T2w-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    outputs = {
        'desc-preproc_T2w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'space-T2w_desc-brain_mask': (
        acpc_align, 'outputspec.acpc_brain_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_brain_T2",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run_t2"],
    ],
    inputs=[
        (
            "desc-preproc_T2w",
            "desc-tempbrain_T2w",
            "T2w-ACPC-template",
            "T2w-brain-ACPC-template",
        )
    ],
    outputs=["desc-preproc_T2w", "desc-acpcbrain_T2w"],
)
def acpc_align_brain_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=False,
                                wf_name=f'acpc_align_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('desc-tempbrain_T2w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_brain')

    node, out = strat_pool.get_data('T2w-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    node, out = strat_pool.get_data('T2w-brain-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-preproc_T2w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'desc-acpcbrain_T2w': (acpc_align, 'outputspec.acpc_aligned_brain')
    }

    return (wf, outputs)


@nodeblock(
    name="acpc_alignment_T2_brain_with_mask",
    switch=[
        ["anatomical_preproc", "acpc_alignment", "run"],
        ["anatomical_preproc", "run_t2"],
    ],
    inputs=[
        ("desc-preproc_T2w", "desc-tempbrain_T2w", "space-T2w_desc-brain_mask"),
        "T2w-ACPC-template",
        "T2w-brain-ACPC-template",
    ],
    outputs=["desc-preproc_T2w", "desc-acpcbrain_T2w", "space-T2w_desc-brain_mask"],
)
def acpc_align_brain_with_mask_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    acpc_align = acpc_alignment(config=cfg,
                                acpc_target=cfg.anatomical_preproc[
                                    'acpc_alignment']['acpc_target'],
                                mask=True,
                                wf_name=f'acpc_align_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_leaf')

    node, out = strat_pool.get_data('desc-tempbrain_T2w')
    wf.connect(node, out, acpc_align, 'inputspec.anat_brain')

    node, out = strat_pool.get_data('space-T2w_desc-brain_mask')
    wf.connect(node, out, acpc_align, 'inputspec.brain_mask')

    node, out = strat_pool.get_data('T2w-ACPC-template') 
    wf.connect(node, out, acpc_align, 'inputspec.template_head_for_acpc')

    node, out = strat_pool.get_data('T2w-brain-ACPC-template')  
    wf.connect(node, out, acpc_align, 'inputspec.template_brain_for_acpc')

    outputs = {
        'desc-preproc_T2w': (acpc_align, 'outputspec.acpc_aligned_head'),
        'desc-acpcbrain_T2w': (acpc_align, 'outputspec.acpc_aligned_brain'),
        'space-T2w_desc-brain_mask': (
        acpc_align, 'outputspec.acpc_brain_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="nlm_filtering_T2",
    switch=[
        ["anatomical_preproc", "non_local_means_filtering", "run"],
        ["anatomical_preproc", "run_t2"],
    ],
    inputs=["desc-preproc_T2w"],
    outputs=["desc-preproc_T2w"],
)
def non_local_means_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    denoise = pe.Node(interface=ants.DenoiseImage(),
                      name=f'anat_denoise_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, denoise, 'input_image')

    outputs = {
        'desc-preproc_T2w': (denoise, 'output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="n4_bias_correction_T2",
    switch=[
        ["anatomical_preproc", "n4_bias_field_correction", "run"],
        ["anatomical_preproc", "run_t2"],
    ],
    inputs=["desc-preproc_T2w"],
    outputs=["desc-preproc_T2w"],
)
def n4_bias_correction_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    n4 = pe.Node(interface=ants.N4BiasFieldCorrection(dimension=3,
                                                      shrink_factor=2,
                                                      copy_header=True),
                 name=f'anat_n4_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, n4, 'input_image')

    outputs = {
        'desc-preproc_T2w': (n4, 'output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_afni_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="3dSkullStrip",
    inputs=["desc-preproc_T2w"],
    outputs=["space-T2w_desc-brain_mask"],
)
def brain_mask_afni_T2(wf, cfg, strat_pool, pipe_num, opt=None):
    wf, outputs = afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_afni_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="3dSkullStrip",
    inputs=["desc-preproc_T2w"],
    outputs=["space-T2w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_afni_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = afni_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T2w_desc-acpcbrain_mask':
            wf_outputs['space-T2w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_fsl_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="BET",
    inputs=["desc-preproc_T2w"],
    outputs=["space-T2w_desc-brain_mask"],
)
def brain_mask_fsl_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_fsl_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="BET",
    inputs=["desc-preproc_T2w"],
    outputs=["space-T2w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_fsl_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = fsl_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T2w_desc-acpcbrain_mask':
            wf_outputs['space-T2w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_niworkflows_ants_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="niworkflows-ants",
    inputs=["desc-preproc_T2w"],
    outputs=["space-T2w_desc-brain_mask"],
)
def brain_mask_niworkflows_ants_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = niworkflows_ants_brain_connector(wf, cfg, strat_pool,
                                                   pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_niworkflows_ants_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="niworkflows-ants",
    inputs=["desc-preproc_T2w"],
    outputs=["space-T2w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_niworkflows_ants_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = niworkflows_ants_brain_connector(wf, cfg, strat_pool,
                                                      pipe_num, opt)

    outputs = {
        'space-T2w_desc-acpcbrain_mask':
            wf_outputs['space-T2w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_unet_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="UNet",
    inputs=["desc-preproc_T2w", "T1w-brain-template", "T1w-template", "unet_model"],
    outputs=["space-T2w_desc-brain_mask"],
)
def brain_mask_unet_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, outputs = unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_unet_T2",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="UNet",
    inputs=["desc-preproc_T2w", "T1w-brain-template", "T1w-template", "unet_model"],
    outputs=["space-T2w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_unet_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    wf, wf_outputs = unet_brain_connector(wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'space-T2w_desc-acpcbrain_mask':
            wf_outputs['space-T2w_desc-brain_mask']
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_T2",
    config=["anatomical_preproc"],
    switch=["run_t2"],
    inputs=[
        (
            ["desc-reorient_T1w", "T1w", "desc-preproc_T1w"],
            ["desc-reorient_T2w", "T2w", "desc-preproc_T2w"],
            ["space-T1w_desc-brain_mask", "space-T1w_desc-acpcbrain_mask"],
            "space-T2w_desc-acpcbrain_mask",
        )
    ],
    outputs=["space-T2w_desc-brain_mask"],
)
def brain_mask_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    brain_mask_T2 = mask_T2(wf_name=f'brain_mask_T2_{pipe_num}')

    if not cfg.anatomical_preproc['acpc_alignment']['run']:
        node, out = strat_pool.get_data(['desc-reorient_T1w','T1w','desc-preproc_T1w'])
        wf.connect(node, out, brain_mask_T2, 'inputspec.T1w')

        node, out = strat_pool.get_data(['desc-reorient_T2w', 'T2w', 'desc-preproc_T2w'])
        wf.connect(node, out, brain_mask_T2, 'inputspec.T2w')

    else:
        node, out = strat_pool.get_data(['desc-preproc_T1w','desc-reorient_T1w','T1w'])
        wf.connect(node, out, brain_mask_T2, 'inputspec.T1w')

        node, out = strat_pool.get_data(['desc-preproc_T2w','desc-reorient_T2w', 'T2w'])
        wf.connect(node, out, brain_mask_T2, 'inputspec.T2w')

    node, out = strat_pool.get_data(["space-T1w_desc-brain_mask",
                                     "space-T1w_desc-acpcbrain_mask"])
    wf.connect(node, out, brain_mask_T2, 'inputspec.T1w_mask')
    
    outputs = {
        'space-T2w_desc-brain_mask': (brain_mask_T2, 'outputspec.T2w_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="brain_mask_acpc_T2",
    config=["anatomical_preproc"],
    switch=["run_t2"],
    inputs=[
        "desc-reorient_T1w",
        "desc-reorient_T2w",
        ["space-T1w_desc-acpcbrain_mask", "space-T1w_desc-prebrain_mask"],
    ],
    outputs=["space-T2w_desc-acpcbrain_mask"],
)
def brain_mask_acpc_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    brain_mask_T2 = mask_T2(wf_name=f'brain_mask_acpc_T2_{pipe_num}')

    node, out = strat_pool.get_data('desc-reorient_T1w')
    wf.connect(node, out, brain_mask_T2, 'inputspec.T1w')

    node, out = strat_pool.get_data('desc-reorient_T2w')
    wf.connect(node, out, brain_mask_T2, 'inputspec.T2w')

    node, out = strat_pool.get_data(["space-T1w_desc-acpcbrain_mask", "space-T1w_desc-prebrain_mask"])
    wf.connect(node, out, brain_mask_T2, 'inputspec.T1w_mask')

    outputs = {
        'space-T2w_desc-acpcbrain_mask': (brain_mask_T2, 'outputspec.T2w_mask')
    }

    return (wf, outputs)


@nodeblock(
    name="brain_extraction_T2",
    config=["anatomical_preproc"],
    switch=["run_t2"],
    inputs=[
        (
            "desc-acpcbrain_T2w",
            "desc-preproc_T2w",
            ["space-T2w_desc-brain_mask", "space-T2w_desc-acpcbrain_mask"],
        )
    ],
    outputs=["desc-brain_T2w"],
)
def brain_extraction_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    if cfg.anatomical_preproc['acpc_alignment']['run'] and cfg.anatomical_preproc['acpc_alignment']['acpc_target'] == 'brain':
        outputs = {
            'desc-brain_T2w': (strat_pool.get_data(["desc-acpcbrain_T2w"]))
        }
    else:
        anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                        name=f'brain_extraction_T2_{pipe_num}')

        anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
        anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data('desc-preproc_T2w')
        wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_a')

        node, out = strat_pool.get_data(['space-T2w_desc-brain_mask'])
        wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_b')

        outputs = {
            'desc-brain_T2w': (anat_skullstrip_orig_vol, 'out_file')
        }

    return (wf, outputs)


@nodeblock(
    name="brain_extraction_temp_T2",
    config=["anatomical_preproc"],
    switch=["run_t2"],
    inputs=[
        (
            "desc-preproc_T2w",
            ["space-T2w_desc-brain_mask", "space-T2w_desc-acpcbrain_mask"],
        )
    ],
    outputs=["desc-tempbrain_T2w"],
)
def brain_extraction_temp_T2(wf, cfg, strat_pool, pipe_num, opt=None):

    anat_skullstrip_orig_vol = pe.Node(interface=afni.Calc(),
                                       name=f'brain_extraction_temp_T2_{pipe_num}')

    anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
    anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data('desc-preproc_T2w')
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_a')

    node, out = strat_pool.get_data(['space-T2w_desc-brain_mask',
                                     'space-T2w_desc-acpcbrain_mask'])
    wf.connect(node, out, anat_skullstrip_orig_vol, 'in_file_b')

    outputs = {
        'desc-tempbrain_T2w': (anat_skullstrip_orig_vol, 'out_file')
    }

    return (wf, outputs)

@nodeblock(
    name="freesurfer_abcd_preproc",
    config=["surface_analysis", "abcd_prefreesurfer_prep"],
    switch=["run"],
    inputs=[
        "desc-preproc_T1w",
        "T1w-template",
        "T1w-brain-template-mask",
        "template-ref-mask-res-2",
        "T1w-template-res-2",
        "freesurfer-subject-dir",
    ],
    outputs=[
        "desc-restore_T1w",
        "desc-restore-brain_T1w",
        "desc-ABCDpreproc_T1w",
        "pipeline-fs_desc-fast_biasfield",
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
        "pipeline-fs_wmparc",
        "freesurfer-subject-dir",
    ],
)
def freesurfer_abcd_preproc(wf, cfg, strat_pool, pipe_num, opt=None):
    # fnirt-based brain extraction
    brain_extraction = fnirt_based_brain_extraction(config=cfg,
                                                    wf_name=f'fnirt_based_brain_extraction_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, brain_extraction, 'inputspec.anat_data')

    node, out = strat_pool.get_data('template-ref-mask-res-2')
    wf.connect(node, out, brain_extraction, 'inputspec.template-ref-mask-res-2')

    node, out = strat_pool.get_data('T1w-template')
    wf.connect(node, out, brain_extraction, 'inputspec.template_skull_for_anat')

    node, out = strat_pool.get_data('T1w-template-res-2')
    wf.connect(node, out, brain_extraction, 'inputspec.template_skull_for_anat_2mm')

    node, out = strat_pool.get_data('T1w-brain-template-mask')
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

    outputs = {
            'desc-restore_T1w': (fast_correction, 'outputspec.anat_restore'),
            'desc-restore-brain_T1w': (fast_correction,
                                       'outputspec.anat_brain_restore'),
            'pipeline-fs_desc-fast_biasfield': (fast_correction, 'outputspec.bias_field'),
            'desc-ABCDpreproc_T1w': (normalize_head, 'out_file')
            }
    return (wf, outputs)

@nodeblock(
    name="freesurfer_reconall",
    config=["surface_analysis", "freesurfer"],
    switch=["run_reconall"],
    inputs=[["desc-ABCDpreproc_T1w", "desc-preproc_T1w"]],
    outputs=[
        "freesurfer-subject-dir",
        "pipeline-fs_raw-average",
        "pipeline-fs_subcortical-seg",
        "pipeline-fs_brainmask",
        "pipeline-fs_wmparc",
        "pipeline-fs_T1",
        *freesurfer_abcd_preproc.outputs
        # we're grabbing the postproc outputs and appending them to
        # the reconall outputs
    ],
)
def freesurfer_reconall(wf, cfg, strat_pool, pipe_num, opt=None):

    reconall = pe.Node(interface=freesurfer.ReconAll(),
                       name=f'anat_freesurfer_{pipe_num}',
                       mem_gb=2.7)
    reconall.skip_timeout = True  # this Node could take > 24 hours

    freesurfer_subject_dir = os.path.join(
        cfg.pipeline_setup['working_directory']['path'],
        'cpac_'+cfg['subject_id'],
        f'anat_preproc_freesurfer_{pipe_num}',
        'anat_freesurfer')

    if not os.path.exists(freesurfer_subject_dir):
        os.makedirs(freesurfer_subject_dir)

    reconall.inputs.directive = 'all'
    reconall.inputs.subjects_dir = freesurfer_subject_dir
    reconall.inputs.openmp = cfg.pipeline_setup['system_config'][
        'num_OMP_threads']

    if cfg.surface_analysis['freesurfer']['reconall_args'] is not None:
        reconall.inputs.args = cfg.surface_analysis['freesurfer'][
            'reconall_args']

    node, out = strat_pool.get_data(["desc-ABCDpreproc_T1w","desc-preproc_T1w"])
    wf.connect(node, out, reconall, 'T1_files')

    wf, hemisphere_outputs = freesurfer_hemispheres(wf, reconall, pipe_num)

    outputs = {
        'freesurfer-subject-dir': (reconall, 'subjects_dir'),
        **hemisphere_outputs,
        'pipeline-fs_raw-average': (reconall, 'rawavg'),
        'pipeline-fs_subcortical-seg': (reconall, 'aseg'),
        'pipeline-fs_brainmask': (reconall, 'brainmask'),
        'pipeline-fs_wmparc': (reconall, 'wmparc'),
        'pipeline-fs_T1': (reconall, 'T1')
    }

    # for label, connection in outputs.items():
    #     strat_pool.set_data(label, connection[0], connection[1],
    #                         deepcopy(strat_pool.get('json')),
    #                         pipe_num, 'freesurfer_reconall', fork=False)
    # Run postproc if we ran reconall
    # wf, post_outputs = freesurfer_postproc(wf, cfg, strat_pool, pipe_num, opt)
    # Update the outputs to include the postproc outputs
    # outputs.update(post_outputs)

    return wf, outputs


def fnirt_based_brain_extraction(config=None,
                                 wf_name='fnirt_based_brain_extraction'):

    ### ABCD Harmonization - FNIRT-based brain extraction ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PreFreeSurfer/scripts/BrainExtraction_FNIRTbased.sh

    preproc = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(fields=['anat_data',
                                                       'template-ref-mask-res-2',
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

    preproc.connect(inputnode, 'template-ref-mask-res-2',
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



@nodeblock(
    name="correct_restore_brain_intensity_abcd",
    config=["anatomical_preproc", "brain_extraction"],
    option_key="using",
    option_val="FreeSurfer-ABCD",
    inputs=[
        (
            "desc-preproc_T1w",
            "desc-n4_T1w",
            "space-T1w_desc-brain_mask",
            "pipeline-fs_desc-fast_biasfield",
            "from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm",
            "from-T1w_to-template_mode-image_xfm",
            "desc-restore-brain_T1w",
        )
    ],
    outputs=["desc-restore-brain_T1w"],
)
def correct_restore_brain_intensity_abcd(wf, cfg, strat_pool, pipe_num,
                                         opt=None):

    ### ABCD Harmonization - Myelin Map ###
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PreFreeSurfer/PreFreeSurferPipeline.sh#L655-L656
    # fslmerge -t ${T1wFolder}/xfms/${T1wImage}_dc ${T1wFolder}/${T1wImage}_acpc ${T1wFolder}/${T1wImage}_acpc ${T1wFolder}/${T1wImage}_acpc
    merge_t1_acpc_to_list = pe.Node(util.Merge(3),
                                    name=f'merge_t1_acpc_to_list_{pipe_num}')

    node, out = strat_pool.get_data('desc-preproc_T1w')
    wf.connect(node, out, merge_t1_acpc_to_list, 'in1')
    wf.connect(node, out, merge_t1_acpc_to_list, 'in2')
    wf.connect(node, out, merge_t1_acpc_to_list, 'in3')

    merge_t1_acpc = pe.Node(interface=fslMerge(),
                            name=f'merge_t1_acpc_{pipe_num}')

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
    
    node, out = strat_pool.get_data("space-T1w_desc-brain_mask")
    wf.connect(node, out, convertwarp_orig_t1_to_t1, 'reference')

    node, out = strat_pool.get_data('from-T1w_to-ACPC_mode-image_desc-aff2rig_xfm')
    wf.connect(node, out, convertwarp_orig_t1_to_t1, 'premat')
    wf.connect(multiply_t1_acpc_by_zero, 'out_file',
        convertwarp_orig_t1_to_t1, 'warp1')

    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/scripts/CreateMyelinMaps.sh#L72-L73
    # applywarp --rel --interp=spline -i "$BiasField" -r "$T1wImageBrain" -w "$AtlasTransform" -o "$BiasFieldOutput"
    applywarp_biasfield = pe.Node(interface=fsl.ApplyWarp(), 
                                  name=f'applywarp_biasfield_{pipe_num}')

    applywarp_biasfield.inputs.relwarp = True
    applywarp_biasfield.inputs.interp = 'spline'

    node, out = strat_pool.get_data('pipeline-fs_desc-fast_biasfield')
    wf.connect(node, out, applywarp_biasfield, 'in_file')

    node, out = strat_pool.get_data("space-T1w_desc-brain_mask")
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
    
    node, out = strat_pool.get_data("space-T1w_desc-brain_mask")
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

    node, out = strat_pool.get_data('pipeline-fs_desc-fast_biasfield')
    wf.connect(node, out, div_t1_by_biasfield, 'in_file2')

    # fslmaths "$OutputT1wImageRestore" -mas "$T1wImageBrain" "$OutputT1wImageRestoreBrain"
    apply_mask = pe.Node(interface=fsl.maths.ApplyMask(),
                         name=f'get_restored_corrected_brain_{pipe_num}')

    wf.connect(div_t1_by_biasfield, 'out_file',
        apply_mask, 'in_file')

    node, out = strat_pool.get_data("space-T1w_desc-brain_mask")
    wf.connect(node, out, apply_mask, 'mask_file')

    outputs = {
        'desc-restore-brain_T1w': (apply_mask, 'out_file')
    }

    return (wf, outputs)

