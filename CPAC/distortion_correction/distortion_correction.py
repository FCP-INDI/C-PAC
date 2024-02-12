#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2017-2022  C-PAC Developers

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
import os
import subprocess
from CPAC.pipeline.nodeblock import nodeblock

import nibabel as nb

from CPAC.pipeline import nipype_pipeline_engine as pe
from nipype.interfaces import afni, fsl
import nipype.interfaces.utility as util
import nipype.interfaces.ants as ants
import nipype.interfaces.afni.preprocess as preprocess
import nipype.interfaces.afni.utils as afni_utils

from CPAC.pipeline.engine import wrap_block

from CPAC.utils import function
from CPAC.utils.interfaces.function import Function
from CPAC.utils.datasource import match_epi_fmaps

from CPAC.func_preproc.func_preproc import bold_mask_afni, bold_masking

from CPAC.distortion_correction.utils import run_convertwarp, \
                                             phase_encode, \
                                             run_fsl_topup, \
                                             z_pad, \
                                             choose_phase_image, \
                                             find_vnum_base


def create_afni_arg(shrink_fac):
    expr = '-shrink_fac {0} '.format(shrink_fac)
    return expr


@nodeblock(
    name="distcor_phasediff_fsl_fugue",
    config=["functional_preproc", "distortion_correction"],
    switch=["run"],
    option_key="using",
    option_val="PhaseDiff",
    inputs=[
        "phasediff",
        "phase1",
        "phase2",
        "magnitude",
        "magnitude1",
        "magnitude2",
        "deltaTE",
        "effectiveEchoSpacing",
        "ees-asym-ratio",
    ],
    outputs=["despiked-fieldmap", "fieldmap-mask"],
)
def distcor_phasediff_fsl_fugue(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    Fieldmap correction takes in an input magnitude image which is
    Skull Stripped (Tight).
    The magnitude images are obtained from each echo series. It also
    requires a phase image
    as an input, the phase image is a subtraction of the two phase
    images from each echo.

    Created on Thu Nov  9 10:44:47 2017
    @author: nrajamani

    Order of commands and inputs:

    -- SkullStrip:   3d-SkullStrip (or FSL-BET) is used to strip the
                     non-brain (tissue) regions from the fMRI
                     Parameters: -f, default: 0.5
                     in_file: fmap_mag
    -- fslmath_mag:  Magnitude image is eroded using the -ero option in
                     fslmath, in order to remove the non-zero voxels
                     Parameters: -ero
                     in_file:fmap_mag
    -- bet_anat:     Brain extraction of the anat file to provide as an
                     input for the epi-registration interface
                     Parameters: -f, default: 0.5
                     in_file: anat_file
    -- fast_anat:    Fast segmentation to provide partial volume files
                     of the anat file, which is further processed to
                     provide the white mater segmentation input for the
                     epi-registration interface.
                     The most important output required from this is
                     the second segment, (e.g.,'T1_brain_pve_2.nii.gz')
                     Parameters: -img_type = 1
                                 -bias_iters = 10 (-I)
                                 -bias_lowpass = 10 (-l)
                     in_file: brain_extracted anat_file
    -- fslmath_anat: The output of the FAST interface is then analyzed
                     to select all the voxels with more than 50%
                     partial volume into the binary mask
                     Parameters: -thr = 0.5
                     in_file : T1_brain_pve_2
    -- fslmath_wmseg:The selected voxels are now used to create a
                     binary mask which would can then be sent as the
                     white matter segmentation (wm_seg)
                     Parameters: -bin
                     in_file: T1_brain_pve_2
    -- Prepare:      Preparing the fieldmap.
                     Parameters: -deltaTE = default, 2.46 ms
                                 -Scanner = SIEMENS
                     in_files: fmap_phase
                               fmap_magnitude
                     For more details, check:
                     https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide
    -- FUGUE:        One of the steps in EPI-DistCorrection toolbox, it
                     unwarps the fieldmaps
                     Parameters: dwell_to_asymm ratio = \
                                     (0.77e-3 * 3)/(2.46e-3)
                                 dwell time = 0.0005 ms
                                 in_file = field map which is a 4D
                                           image (containing 2
                                           unwarped image)
    '''

    # Skull-strip, outputs a masked image file
    if cfg.functional_preproc['distortion_correction']['PhaseDiff'][
            'fmap_skullstrip_option'] == 'AFNI':

        skullstrip_args = pe.Node(util.Function(input_names=['shrink_fac'],
                                                output_names=['expr'],
                                                function=create_afni_arg),
                                  name=f'distcorr_skullstrip_arg_{pipe_num}')
        skullstrip_args.inputs.shrink_fac = cfg.functional_preproc[
            'distortion_correction']['PhaseDiff'][
            'fmap_skullstrip_AFNI_threshold']

        afni = pe.Node(interface=afni.SkullStrip(), name=f'distcor_phasediff_'
                                                         f'afni_skullstrip_'
                                                         f'{pipe_num}')
        afni.inputs.outputtype = 'NIFTI_GZ'
        wf.connect(skullstrip_args, 'expr', afni, 'args')

        if strat_pool.check_rpool('magnitude'):
            node, out = strat_pool.get_data('magnitude')
        elif strat_pool.check_rpool('magnitude1'):
            node, out = strat_pool.get_data('magnitude1')

        wf.connect(node, out, afni, 'in_file')

        brain_node, brain_out = (afni, 'out_file')

    elif cfg.functional_preproc['distortion_correction']['PhaseDiff'][
            'fmap_skullstrip_option'] == 'BET':

        bet = pe.Node(interface=fsl.BET(), name='distcor_phasediff_bet_'
                                                f'skullstrip_{pipe_num}')
        bet.inputs.output_type = 'NIFTI_GZ'
        bet.inputs.frac = cfg.functional_preproc['distortion_correction'][
            'PhaseDiff']['fmap_skullstrip_BET_frac']

        if strat_pool.check_rpool('magnitude'):
            node, out = strat_pool.get_data('magnitude')
        elif strat_pool.check_rpool('magnitude1'):
            node, out = strat_pool.get_data('magnitude1')
        wf.connect(node, out, bet, 'in_file')

        brain_node, brain_out = (bet, 'out_file')

    # Prepare Fieldmap

    # prepare the field map
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(), name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"

    node, out = strat_pool.get_data('deltaTE')
    wf.connect(node, out, prepare, 'delta_TE')

    if strat_pool.check_rpool('phase1') and strat_pool.check_rpool('phase2'):
        fslmaths_sub = pe.Node(interface=fsl.BinaryMaths(), 
                               name='fugue_phase_subtraction')
        fslmaths_sub.inputs.operation = 'sub'

        node, out = strat_pool.get_data('phase1')
        wf.connect(node, out, fslmaths_sub, 'in_file')
    
        node, out = strat_pool.get_data('phase2')
        wf.connect(node, out, fslmaths_sub, 'operand_file')

        node, out = (fslmaths_sub, 'out_file')

    elif strat_pool.check_rpool('phasediff'):
        node, out = strat_pool.get_data('phasediff')
    
    wf.connect(node, out, prepare, 'in_phase')
    wf.connect(brain_node, brain_out, prepare, 'in_magnitude')

    # erode the masked magnitude image
    fslmath_mag = pe.Node(interface=fsl.ErodeImage(), name='fslmath_mag')
    wf.connect(brain_node, brain_out, fslmath_mag, 'in_file')

    # calculate the absolute value of the eroded and masked magnitude
    # image
    fslmath_abs = pe.Node(interface=fsl.UnaryMaths(), name='fslmath_abs')
    fslmath_abs.inputs.operation = 'abs'
    wf.connect(fslmath_mag, 'out_file', fslmath_abs, 'in_file')

    # binarize the absolute value of the eroded and masked magnitude
    # image
    fslmath_bin = pe.Node(interface=fsl.UnaryMaths(), name='fslmath_bin')
    fslmath_bin.inputs.operation = 'bin'
    wf.connect(fslmath_abs, 'out_file', fslmath_bin, 'in_file')

    # take the absolute value of the fieldmap calculated in the prepare step
    fslmath_mask_1 = pe.Node(interface=fsl.UnaryMaths(),
                             name='fslmath_mask_1')
    fslmath_mask_1.inputs.operation = 'abs'
    wf.connect(prepare, 'out_fieldmap', fslmath_mask_1, 'in_file')

    # binarize the absolute value of the fieldmap calculated in the
    # prepare step
    fslmath_mask_2 = pe.Node(interface=fsl.UnaryMaths(),
                             name='fslmath_mask_2')
    fslmath_mask_2.inputs.operation = 'bin'
    wf.connect(fslmath_mask_1, 'out_file', fslmath_mask_2, 'in_file')

    # multiply together the binarized magnitude and fieldmap images
    fslmath_mask = pe.Node(interface=fsl.BinaryMaths(),
                           name='fslmath_mask')
    fslmath_mask.inputs.operation = 'mul'
    wf.connect(fslmath_mask_2, 'out_file', fslmath_mask, 'in_file')
    wf.connect(fslmath_bin, 'out_file', fslmath_mask, 'operand_file')

    # Note for the user. Ensure the phase image is within 0-4096 (upper
    # threshold is 90% of 4096), fsl_prepare_fieldmap will only work in the
    # case of the SIEMENS format. #Maybe we could use deltaTE also as an
    # option in the GUI.

    # fugue
    fugue1 = pe.Node(interface=fsl.FUGUE(), name='fsl_fugue')
    fugue1.inputs.save_fmap = True
    fugue1.outputs.fmap_out_file = 'fmap_rads'

    wf.connect(fslmath_mask, 'out_file', fugue1, 'mask_file')

    # FSL calls EffectiveEchoSpacing "dwell_time"
    node, out = strat_pool.get_data('effectiveEchoSpacing')
    wf.connect(node, out, fugue1, 'dwell_time')
    node, out = strat_pool.get_data('ees-asym-ratio')
    wf.connect(node, out, fugue1, 'dwell_to_asym_ratio')

    wf.connect(prepare, 'out_fieldmap', fugue1, 'fmap_in_file')

    outputs = {
        'despiked-fieldmap': (fugue1, 'fmap_out_file'),
        'fieldmap-mask': (fslmath_mask, 'out_file')
    }

    return (wf, outputs)


def same_pe_direction_prep(same_pe_epi, func_mean):
    """Skull-strip and align the EPI field map that has the same phase
    encoding direction as the BOLD scan.

    This function only exists to serve as a function node in the
    blip_distcor_wf sub-workflow because workflow decisions cannot be made
    dynamically at runtime."""

    if not same_pe_epi:
        qwarp_input = func_mean
    elif same_pe_epi:
        skullstrip_outfile = os.path.join(os.getcwd(),
                                          "{0}_mask.nii.gz".format(
                                              os.path.basename(same_pe_epi)
                                          ))
        skullstrip_cmd = ["3dAutomask", "-apply_prefix",
                          "{0}_masked.nii.gz".format(os.path.basename(
                              same_pe_epi
                          )), "-prefix", skullstrip_outfile, same_pe_epi]
        retcode = subprocess.check_output(skullstrip_cmd)

        extract_brain_outfile = os.path.join(os.getcwd(),
                                             "{0}_calc.nii.gz".format(
                                                 os.path.basename(same_pe_epi)
                                             ))
        extract_brain_cmd = ["3dcalc", "-a", same_pe_epi, "-b",
                             skullstrip_outfile, "-expr", "a*b", "-prefix",
                             extract_brain_outfile]
        retcode = subprocess.check_output(extract_brain_cmd)

        align_epi_outfile = os.path.join(os.getcwd(),
                                         "{0}_calc_flirt.nii.gz".format(
                                             os.path.basename(same_pe_epi)
                                         ))
        align_epi_cmd = ["flirt", "-in", extract_brain_outfile, "-ref",
                         func_mean, "-out", align_epi_outfile, "-cost",
                         "corratio"]
        retcode = subprocess.check_output(align_epi_cmd)

        qwarp_input = align_epi_outfile

    return qwarp_input


def calculate_blip_warp(opp_pe, same_pe):

    out_warp = os.path.join(os.getcwd(), "Qwarp_PLUS_WARP.nii.gz")

    cmd = ["3dQwarp", "-prefix", "Qwarp.nii.gz", "-plusminus", "-base",
           opp_pe, "-source", same_pe]

    retcode = subprocess.check_output(cmd)

    return out_warp


def convert_afni_to_ants(afni_warp):
    afni_warp = os.path.abspath(afni_warp)
    afni_warp_img = nb.load(afni_warp)
    afni_warp_hdr = afni_warp_img.header.copy()
    afni_warp_hdr.set_data_dtype('<f4')
    afni_warp_hdr.set_intent('vector', (), '')

    afni_warp_data = afni_warp_img.get_fdata().astype('<f4')

    ants_warp = os.path.join(os.getcwd(), os.path.basename(afni_warp))
    nb.Nifti1Image(afni_warp_data,
                   afni_warp_img.affine,
                   afni_warp_hdr).to_filename(ants_warp)

    return ants_warp


@nodeblock(
    name="distcor_blip_afni_qwarp",
    config=["functional_preproc", "distortion_correction"],
    switch=["run"],
    option_key="using",
    option_val="Blip",
    inputs=[
        ("sbref", "space-bold_desc-brain_mask"),
        "epi-1",
        "epi-1-scan-params",
        "epi-2",
        "epi-2-scan-params",
        "pe-direction",
    ],
    outputs=["sbref", "space-bold_desc-brain_mask", "ants-blip-warp"],
)
def distcor_blip_afni_qwarp(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Execute AFNI 3dQWarp to calculate the distortion "unwarp" for
    phase encoding direction EPI field map distortion correction.

        1. Skull-strip the opposite-direction phase encoding EPI.
        2. Transform the opposite-direction phase encoding EPI to the
           skull-stripped functional and pass this as the base_file to
           AFNI 3dQWarp (plus-minus).
        3. If there is a same-direction phase encoding EPI, also
           skull-strip this, and transform it to the skull-stripped
           functional. Then, pass this as the in_file to AFNI 3dQWarp
           (plus-minus).
        4. If there isn't a same-direction, pass the functional in as
           the in_file of AFNI 3dQWarp (plus-minus).
        5. Convert the 3dQWarp transforms to ANTs/ITK format.
        6. Use antsApplyTransforms, with the original functional as
           both the input and the reference, and apply the warp from
           3dQWarp. The output of this can then proceed to
           func_preproc.
    '''

    match_epi_imports = ['import json']
    match_epi_fmaps_node = \
        pe.Node(Function(input_names=['bold_pedir',
                                      'epi_fmap_one',
                                      'epi_fmap_params_one',
                                      'epi_fmap_two',
                                      'epi_fmap_params_two'],
                         output_names=['opposite_pe_epi',
                                       'same_pe_epi'],
                         function=match_epi_fmaps,
                         imports=match_epi_imports,
                         as_module=True),
                name=f'match_epi_fmaps_{pipe_num}')

    node, out = strat_pool.get_data('epi-1')
    wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_one')

    node, out = strat_pool.get_data('epi-1-scan-params')
    wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_params_one')

    if strat_pool.check_rpool('epi-2'):
        node, out = strat_pool.get_data('epi-2')
        wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_two')

        node, out = strat_pool.get_data('epi-2-scan-params')
        wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_params_two')

    node, out = strat_pool.get_data('pe-direction')
    wf.connect(node, out, match_epi_fmaps_node, 'bold_pedir')

    #interface = {'bold': (match_epi_fmaps_node, 'opposite_pe_epi'),
    #             'desc-brain_bold': 'opposite_pe_epi_brain'}
    #wf, strat_pool = wrap_block([bold_mask_afni, bold_masking],
    #                            interface, wf, cfg, strat_pool, pipe_num, opt)

    func_get_brain_mask = pe.Node(interface=preprocess.Automask(),
                                  name=f'afni_mask_opposite_pe_{pipe_num}')
    func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(match_epi_fmaps_node, 'opposite_pe_epi', func_get_brain_mask, 'in_file')

    func_edge_detect = pe.Node(interface=afni_utils.Calc(),
                               name=f'skullstrip_opposite_pe_{pipe_num}')

    func_edge_detect.inputs.expr = 'a*b'
    func_edge_detect.inputs.outputtype = 'NIFTI_GZ'
 
    wf.connect(match_epi_fmaps_node, 'opposite_pe_epi', func_edge_detect, 'in_file_a')
    wf.connect(func_get_brain_mask, 'out_file', func_edge_detect, 'in_file_b')

    opp_pe_to_func = pe.Node(interface=fsl.FLIRT(), name='opp_pe_to_func')
    opp_pe_to_func.inputs.cost = 'corratio'
 
    wf.connect(func_edge_detect, 'out_file',  opp_pe_to_func, 'in_file')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, opp_pe_to_func, 'reference')

    prep_qwarp_input_imports = ['import os', 'import subprocess']
    prep_qwarp_input = \
        pe.Node(function.Function(input_names=['same_pe_epi',
                                               'func_mean'],
                                  output_names=['qwarp_input'],
                                  function=same_pe_direction_prep,
                                  imports=prep_qwarp_input_imports),
                name='prep_qwarp_input')

    wf.connect(match_epi_fmaps_node, 'same_pe_epi',
               prep_qwarp_input, 'same_pe_epi')

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, prep_qwarp_input, 'func_mean')

    calculate_blip_warp_imports = ['import os', 'import subprocess']
    calc_blip_warp = pe.Node(function.Function(
        input_names=['opp_pe', 'same_pe'],
        output_names=['out_warp'],
        function=calculate_blip_warp,
        imports=calculate_blip_warp_imports),
                             name='calc_blip_warp')

    wf.connect(opp_pe_to_func, 'out_file', calc_blip_warp, 'opp_pe')
    wf.connect(prep_qwarp_input, 'qwarp_input', calc_blip_warp, 'same_pe')

    convert_afni_warp_imports = ['import os', 'import nibabel as nb']
    convert_afni_warp = \
        pe.Node(function.Function(input_names=['afni_warp'],
                                  output_names=['ants_warp'],
                                  function=convert_afni_to_ants,
                                  imports=convert_afni_warp_imports),
                name='convert_afni_warp')

    wf.connect(calc_blip_warp, 'out_warp', convert_afni_warp, 'afni_warp')

    # TODO: inverse source_warp (node:source_warp_inverse)
    #       wf.connect(###
    #                  output_node, 'blip_warp_inverse')

    undistort_func_mean = pe.Node(interface=ants.ApplyTransforms(),
                                  name='undistort_func_mean', mem_gb=.1)

    undistort_func_mean.inputs.out_postfix = '_antswarp'
    undistort_func_mean.interface.num_threads = 1
    undistort_func_mean.inputs.interpolation = "LanczosWindowedSinc"
    undistort_func_mean.inputs.dimension = 3
    undistort_func_mean.inputs.input_image_type = 0

    node, out = strat_pool.get_data('sbref')
    wf.connect(node, out, undistort_func_mean, 'input_image')
    wf.connect(node, out, undistort_func_mean, 'reference_image')
    wf.connect(convert_afni_warp, 'ants_warp',
               undistort_func_mean, 'transforms')

    #interface = {'desc-preproc_bold': (undistort_func_mean, 'output_image')}
    #wf, strat_pool = wrap_block([bold_mask_afni],
    #                            interface, wf, cfg, strat_pool, pipe_num, opt)

    remask = pe.Node(interface=preprocess.Automask(),
                     name=f'afni_remask_boldmask_{pipe_num}')
    remask.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(undistort_func_mean, 'output_image', remask, 'in_file')

    outputs = {
        'ants-blip-warp': (convert_afni_warp, 'ants_warp'),
        #'inv-blip-warp': None,  # TODO
        'sbref': (undistort_func_mean, 'output_image'),
        'space-bold_desc-brain_mask': (remask, 'out_file')
    }

    return (wf, outputs)


@nodeblock(
    name="distcor_blip_fsl_topup",
    config=["functional_preproc", "distortion_correction"],
    switch=["run"],
    option_key="using",
    option_val="Blip-FSL-TOPUP",
    inputs=[
        ("sbref", "space-bold_desc-brain_mask"),
        "pe-direction",
        "epi-1",
        "epi-1-pedir",
        "epi-1-TE",
        "epi-1-dwell",
        "epi-1-total-readout",
        "epi-2",
        "epi-2-pedir",
        "epi-2-TE",
        "epi-2-dwell",
        "epi-2-total-readout",
    ],
    outputs=["sbref", "space-bold_desc-brain_mask", "fsl-blip-warp"],
)
def distcor_blip_fsl_topup(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Execute FSL TOPUP to calculate the distortion "unwarp" for
    phase encoding direction EPI field map distortion correction.
    '''

    # TODO: re-integrate gradient distortion coefficient usage at a later date
    '''
    GradientDistortionCoeffs = None
    if GradientDistortionCoeffs != None:
        grad_outs = []
        gdc_phase_outs = []
        for idx, phase_img in enumerate(phase_imgs):
            wf, out_warpmask, out_applywarp = gradient_distortion_correction(
                wf, phase_img, f"phase-{idx}"
            )
            grad_outs.append(out_warpmask)
            gdc_phase_outs.append(out_applywarp)
            phaseone_gdc = gdc_phase_outs[0][1]
            phasetwo_gdc = gdc_phase_outs[1][1]

        # Make a conservative (eroded) intersection of the two masks
        Mask = pe.Node(interface=fsl.maths.ApplyMask(), name="Mask")
        Mask.inputs.args = "-bin -ero"

        wf.connect(grad_outs[0][0], grad_outs[0][1], Mask, "in_file")
        wf.connect(grad_outs[1][0], grad_outs[1][1], Mask, "mask_file")

        # Merge both sets of images
        create_list = pe.Node(interface=util.Merge(2), name="create_list")

        wf.connect(gdc_phase_outs[0][0], gdc_phase_outs[0][1], create_list, "in1")
        wf.connect(gdc_phase_outs[1][0], gdc_phase_outs[1][1], create_list, "in2")

        merge_image = pe.Node(interface=fsl.Merge(), name="merge_image")
        merge_image.inputs.dimension = "t"
        merge_image.inputs.output_type = "NIFTI_GZ"

        wf.connect(create_list, "out", merge_image, "in_files")

    else:
    '''

    create_list = pe.Node(interface=util.Merge(2), name="create_list")

    node, out = strat_pool.get_data('epi-1')
    wf.connect(node, out, create_list, 'in1')

    node, out = strat_pool.get_data('epi-2')
    wf.connect(node, out, create_list, 'in2')

    merge_image = pe.Node(interface=fsl.Merge(), name="merge_image")
    merge_image.inputs.dimension = "t"
    merge_image.inputs.output_type = "NIFTI_GZ"

    wf.connect(create_list, "out", merge_image, "in_files")

    Mask = pe.Node(interface=fsl.maths.BinaryMaths(), name="Mask")
    Mask.inputs.operand_value = 0
    Mask.inputs.operation = "mul"
    Mask.inputs.args = "-add 1"

    node, out = strat_pool.get_data('epi-1')
    wf.connect(node, out, Mask, 'in_file')

    # zpad_phases = z_pad("zpad_phases")
    # wf.connect(merge_image, "merged_file", zpad_phases, "inputspec.input_image")

    # zpad_mask = z_pad("zpad_mask")
    # wf.connect(Mask, "out_file", zpad_mask, "inputspec.input_image")

    # extrapolate existing values beyond the mask
    extrap_vals = pe.Node(interface=fsl.maths.BinaryMaths(),
                          name="extrap_vals")
    extrap_vals.inputs.operation = "add"
    extrap_vals.inputs.operand_value = 1
    extrap_vals.inputs.args = "-abs -dilM -dilM -dilM -dilM -dilM"

    # wf.connect(zpad_phases, "outputspec.output_image",
    #            extrap_vals,  "in_file")
    # wf.connect(zpad_mask, "outputspec.output_image",
    #            extrap_vals, "operand_file")

    wf.connect(merge_image, "merged_file", extrap_vals, "in_file")
    wf.connect(Mask, "out_file", extrap_vals, "operand_file")

    # phase encoding
    phase_encode_imports = [
        "import os",
        "import subprocess",
        "import numpy as np",
        "import nibabel",
        "import sys"
    ]
    phase_encoding = pe.Node(
        util.Function(
            input_names=[
                "unwarp_dir",
                "phase_one",
                "phase_two",
                "dwell_time_one",
                "dwell_time_two",
                "ro_time_one",
                "ro_time_two"
            ],
            output_names=["acq_params"],
            function=phase_encode,
            imports=phase_encode_imports,
        ),
        name="phase_encoding",
    )
    node, out = strat_pool.get_data('epi-1')
    wf.connect(node, out, phase_encoding, 'phase_one')

    node, out = strat_pool.get_data('epi-2')
    wf.connect(node, out, phase_encoding, 'phase_two')

    node, out = strat_pool.get_data('pe-direction')
    wf.connect(node, out, phase_encoding, 'unwarp_dir')
    
    if strat_pool.check_rpool('epi-1-dwell') and \
            strat_pool.check_rpool('epi-2-dwell'):

        node, out = strat_pool.get_data('epi-1-dwell')
        wf.connect(node, out, phase_encoding, 'dwell_time_one')

        node, out = strat_pool.get_data('epi-2-dwell')
        wf.connect(node, out, phase_encoding, 'dwell_time_two')

    if strat_pool.check_rpool('epi-1-total-readout') and \
            strat_pool.check_rpool('epi-2-total-readout'):

        node, out = strat_pool.get_data('epi-1-total-readout')
        wf.connect(node, out, phase_encoding, 'ro_time_one')
    
        node, out = strat_pool.get_data('epi-2-total-readout')
        wf.connect(node, out, phase_encoding, 'ro_time_two')

    topup_imports = ["import os",
                     "import subprocess"]
    run_topup = pe.Node(util.Function(input_names=["merged_file", 
                                                   "acqparams"],
                                      output_names=["out_fieldcoef", 
                                                    "out_movpar",
                                                    "corrected_outfile",
                                                    "field_out",
                                                    "out_jacs",
                                                    "log_out",
                                                    "out_xfms",
                                                    "out_warps"],
                                      function=run_fsl_topup,
                                      imports=topup_imports),
                        name="topup")

    wf.connect(extrap_vals, "out_file", run_topup, "merged_file")
    wf.connect(phase_encoding, "acq_params", run_topup, "acqparams")

    choose_phase = pe.Node(
        util.Function(
            input_names=["phase_imgs",
                         "unwarp_dir"],
            output_names=["out_phase_image",
                          "vnum"],
            function=choose_phase_image
        ), name="choose_phase",
    )

    wf.connect(create_list, 'out', choose_phase, 'phase_imgs')

    node, out = strat_pool.get_data("pe-direction")
    wf.connect(node, out, choose_phase, "unwarp_dir")

    vnum_base = pe.Node(
        util.Function(
            input_names=["vnum",
                         "motion_mat_list",
                         "jac_matrix_list",
                         "warp_field_list"],
            output_names=["out_motion_mat",
                          "out_jacobian",
                          "out_warp_field"],
            function=find_vnum_base
        ), name="Motion_Jac_Warp_matrices",
    )

    wf.connect(choose_phase, 'vnum', vnum_base, 'vnum')
    wf.connect(run_topup, 'out_xfms', vnum_base, 'motion_mat_list')
    wf.connect(run_topup, 'out_jacs', vnum_base, 'jac_matrix_list')
    wf.connect(run_topup, 'out_warps', vnum_base, 'warp_field_list')

    mean_bold = strat_pool.node_data("sbref")

    flirt = pe.Node(interface=fsl.FLIRT(), name="flirt")
    flirt.inputs.dof = 6
    flirt.inputs.interp = 'spline'
    flirt.inputs.out_matrix_file = 'SBRef2PhaseTwo_gdc.mat'

    wf.connect(mean_bold.node, mean_bold.out, flirt, 'in_file')
    wf.connect(choose_phase, 'out_phase_image', flirt, 'reference')

    # fsl_convert_xfm
    convert_xfm = pe.Node(interface=fsl.ConvertXFM(), name="convert_xfm")
    convert_xfm.inputs.concat_xfm = True
    convert_xfm.inputs.out_file = 'SBRef2WarpField.mat'

    wf.connect(flirt, 'out_matrix_file', convert_xfm, 'in_file')
    wf.connect(vnum_base, 'out_motion_mat', convert_xfm, 'in_file2')

    # fsl_convert_warp
    convert_warp = pe.Node(interface=fsl.ConvertWarp(),
                           name="convert_warp")
    convert_warp.inputs.relwarp = True
    convert_warp.inputs.out_relwarp = True
    convert_warp.inputs.out_file = 'WarpField.nii.gz'

    wf.connect(choose_phase, 'out_phase_image', convert_warp, 'reference')
    wf.connect(vnum_base, 'out_warp_field', convert_warp, 'warp1')
    wf.connect(convert_xfm, 'out_file', convert_warp, 'premat')

    VolumeNumber = 1 + 1
    vnum = str(VolumeNumber).zfill(2)
    name = "PhaseTwo_aw"

    vnum_base_two = pe.Node(
        util.Function(
            input_names=["vnum",
                         "motion_mat_list",
                         "jac_matrix_list",
                         "warp_field_list"],
            output_names=["out_motion_mat",
                          "out_jacobian",
                          "out_warp_field"],
            function=find_vnum_base
        ), name=f"Motion_Jac_Warp_matrices_{name}",
    )
    vnum_base_two.inputs.vnum = vnum

    wf.connect(run_topup, 'out_xfms', vnum_base_two, 'motion_mat_list')
    wf.connect(run_topup, 'out_jacs', vnum_base_two, 'jac_matrix_list')
    wf.connect(run_topup, 'out_warps', vnum_base_two, 'warp_field_list')

    # fsl_applywarp
    aw_two = pe.Node(interface=fsl.ApplyWarp(), name="aw_two")
    aw_two.inputs.relwarp = True
    aw_two.inputs.interp = 'spline'

    node, out = strat_pool.get_data('epi-2')
    wf.connect(node, out, aw_two, 'in_file')
    wf.connect(node, out, aw_two, 'ref_file')

    wf.connect(vnum_base_two, 'out_motion_mat', aw_two, 'premat')
    wf.connect(vnum_base_two, 'out_warp_field', aw_two, 'field_file')

    mul_phase_two = pe.Node(interface=fsl.BinaryMaths(),
                            name="mul_phase_two")
    mul_phase_two.inputs.operation = 'mul'

    wf.connect(aw_two, 'out_file', mul_phase_two, 'in_file')
    wf.connect(vnum_base_two, 'out_jacobian', mul_phase_two, 'operand_file')

    # PhaseOne (first vol) - warp and Jacobian modulate to get
    # distortion corrected output
    VolumeNumber = 0 + 1
    vnum = str(VolumeNumber).zfill(2)
    name = "PhaseOne_aw"

    vnum_base_one = pe.Node(
        util.Function(
            input_names=["vnum",
                         "motion_mat_list",
                         "jac_matrix_list",
                         "warp_field_list"],
            output_names=["out_motion_mat",
                          "out_jacobian",
                          "out_warp_field"],
            function=find_vnum_base
        ), name=f"Motion_Jac_Warp_matrices_{name}",
    )
    vnum_base_one.inputs.vnum = vnum

    wf.connect(run_topup, 'out_xfms', vnum_base_one, 'motion_mat_list')
    wf.connect(run_topup, 'out_jacs', vnum_base_one, 'jac_matrix_list')
    wf.connect(run_topup, 'out_warps', vnum_base_one, 'warp_field_list')

    # fsl_applywarp to phaseOne
    aw_one = pe.Node(interface=fsl.ApplyWarp(), name="aw_one")
    aw_one.inputs.relwarp = True
    aw_one.inputs.interp = 'spline'

    node, out = strat_pool.get_data('epi-1')
    wf.connect(node, out, aw_one, 'in_file')
    wf.connect(node, out, aw_one, 'ref_file')

    wf.connect(vnum_base_one, 'out_motion_mat', aw_one, 'premat')
    wf.connect(vnum_base_one, 'out_warp_field', aw_one, 'field_file')

    mul_phase_one = pe.Node(interface=fsl.BinaryMaths(), name="mul_phase_one")
    mul_phase_one.inputs.operation = 'mul'

    wf.connect(aw_one, 'out_file', mul_phase_one, 'in_file')
    wf.connect(vnum_base_one, 'out_jacobian', mul_phase_one, 'operand_file')

    # Scout - warp and Jacobian modulate to get distortion corrected output
    aw_jac = pe.Node(interface=fsl.ApplyWarp(), name="aw_jac")
    aw_jac.inputs.relwarp = True
    aw_jac.inputs.interp = 'spline'

    wf.connect(mean_bold.node, mean_bold.out, aw_jac, 'in_file') # SBRef.nii.gz
    wf.connect(mean_bold.node, mean_bold.out,
               aw_jac, 'ref_file') # SBRef.nii.gz
    wf.connect(convert_warp, 'out_file', aw_jac, 'field_file')

    mul_jac = pe.Node(interface=fsl.BinaryMaths(), name="mul_jac")
    mul_jac.inputs.operation = 'mul'
    mul_jac.inputs.out_file = "SBRef_dc_jac.nii.gz"

    wf.connect(aw_jac, 'out_file', mul_jac, 'in_file')
    wf.connect(vnum_base, 'out_jacobian', mul_jac, 'operand_file')

    # Calculate Equivalent Field Map
    tp_field_map = pe.Node(interface=fsl.BinaryMaths(), name="tp_field_map")
    tp_field_map.inputs.operation = 'mul'
    tp_field_map.inputs.operand_value = 6.283

    wf.connect(run_topup, 'out_fieldcoef', tp_field_map, 'in_file')

    mag_field_map = pe.Node(interface=fsl.MeanImage(),
                            name="mag_field_map")
    mag_field_map.inputs.dimension = 'T'
    mag_field_map.inputs.out_file = 'Magnitude.nii.gz'

    wf.connect(run_topup, 'corrected_outfile', mag_field_map, 'in_file')

    # fsl_bet
    bet = pe.Node(interface=fsl.BET(), name="bet")
    bet.inputs.frac = 0.35
    bet.inputs.mask = True

    wf.connect(mag_field_map, 'out_file', bet, 'in_file')

    outputs = {
        'sbref': (mul_jac, 'out_file'),
        'fsl-blip-warp': (convert_warp, 'out_file')
        #'space-bold_desc-brain_mask': (mask_sbref, 'out_file')
    }

    return (wf, outputs)
