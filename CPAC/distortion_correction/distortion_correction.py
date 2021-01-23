#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import subprocess

import nibabel as nb

import nipype.pipeline.engine as pe
from nipype.interfaces import afni, fsl
import nipype.interfaces.utility as util
import nipype.interfaces.ants as ants

from CPAC.pipeline.engine import wrap_block

from CPAC.utils import function
from CPAC.utils.interfaces.function import Function
from CPAC.utils.datasource import match_epi_fmaps

from CPAC.func_preproc.func_preproc import bold_mask_afni, bold_masking


def create_afni_arg(shrink_fac):
    expr = '-shrink_fac {0} '.format(shrink_fac)
    return expr


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

    Node Block:
    {"name": "distcor_phasediff_fsl_fugue",
     "config": ["functional_preproc", "distortion_correction"],
     "switch": ["run"],
     "option_key": "using",
     "option_val": "PhaseDiff",
     "inputs": ["diff_phase",
                "diff_mag_one",
                "deltaTE",
                "diff_phase_dwell",
                "dwell_asym_ratio"],
     "outputs": ["despiked_fieldmap",
                 "fieldmap_mask"]}
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

        node, out = strat_pool.get_data('diff_mag_one')
        wf.connect(node, out, afni, 'in_file')

        brain_node, brain_out = (afni, 'out_file')

    elif cfg.functional_wf['distortion_correction']['PhaseDiff'][
            'fmap_skullstrip_option'] == 'BET':

        bet = pe.Node(interface=fsl.BET(), name='distcor_phasediff_bet_'
                                                f'skullstrip_{pipe_num}')
        bet.inputs.output_type = 'NIFTI_GZ'
        bet.inputs.frac = cfg.functional_preproc['distortion_correction'][
            'PhaseDiff']['fmap_skullstrip_BET_frac']

        node, out = strat_pool.get_data('diff_mag_one')
        wf.connect(node, out, bet, 'in_file')

        brain_node, brain_out = (bet, 'out_file')

    # Prepare Fieldmap

    # prepare the field map
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(), name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"

    node, out = strat_pool.get_data('deltaTE')
    wf.connect(node, out, prepare, 'delta_TE')

    node, out = strat_pool.get_data('diff_phase')
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
    fugue1 = pe.Node(interface=fsl.FUGUE(), name='fugue1')
    fugue1.inputs.save_fmap = True
    fugue1.outputs.fmap_out_file = 'fmap_rads'

    wf.connect(fslmath_mask, 'out_file', fugue1, 'mask_file')

    node, out = strat_pool.get_data('diff_phase_dwell')
    wf.connect(node, out, fugue1, 'dwell_time')

    node, out = strat_pool.get_data('dwell_asym_ratio')
    wf.connect(node, out, fugue1, 'dwell_to_asym_ratio')

    wf.connect(prepare, 'out_fieldmap', fugue1, 'fmap_in_file')

    outputs = {
        'despiked_fieldmap': (fugue1, 'fmap_out_file'),
        'fieldmap_mask': (fslmath_mask, 'out_file')
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

    afni_warp_data = afni_warp_img.get_data().astype('<f4')

    ants_warp = os.path.join(os.getcwd(), os.path.basename(afni_warp))
    nb.Nifti1Image(afni_warp_data,
                   afni_warp_img.affine,
                   afni_warp_hdr).to_filename(ants_warp)

    return ants_warp


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

    Node Block:
    {"name": "distcor_blip_afni_qwarp",
     "config": ["functional_preproc", "distortion_correction"],
     "switch": ["run"],
     "option_key": "using",
     "option_val": "Blip",
     "inputs": ["epi_1",
                "epi_1_scan_params",
                "epi_2",
                "epi_2_scan_params",
                "pe_direction"],
     "outputs": []}
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

    node, out = strat_pool.get_data('epi_1')
    wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_one')

    node, out = strat_pool.get_data('epi_1_scan_params')
    wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_params_one')

    if strat_pool.check_rpool('epi_2'):
        node, out = strat_pool.get_data('epi_2')
        wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_two')

        node, out = strat_pool.get_data('epi_2_scan_params')
        wf.connect(node, out, match_epi_fmaps_node, 'epi_fmap_params_two')

    node, out = strat_pool.get_data('pe_direction')
    wf.connect(node, out, match_epi_fmaps_node, 'bold_pedir')

    interface = {'bold': (match_epi_fmaps_node, 'opposite_pe_epi'),
                 'desc-brain_bold': 'opposite_pe_epi_brain'}
    wf, strat_pool = wrap_block([bold_mask_afni, bold_masking],
                                interface, wf, cfg, strat_pool, pipe_num, opt)

    opp_pe_to_func = pe.Node(interface=fsl.FLIRT(), name='opp_pe_to_func')
    opp_pe_to_func.inputs.cost = 'corratio'

    node, out = strat_pool.get_data('opposite_pe_epi_brain')
    wf.connect(node, out, opp_pe_to_func, 'in_file')

    node, out = strat_pool.get_data('desc-mean_bold')
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

    node, out = strat_pool.get_data('desc-mean_bold')
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

    node, out = strat_pool.get_data('desc-mean_bold')
    wf.connect(node, out, undistort_func_mean, 'input_image')
    wf.connect(node, out, undistort_func_mean, 'reference_image')
    wf.connect(convert_afni_warp, 'ants_warp',
               undistort_func_mean, 'transforms')

    interface = {'bold': (undistort_func_mean, 'output_image'),
                 'space-bold_desc-brain_mask': 'opposite_pe_epi_brain'}
    wf, strat_pool = wrap_block([bold_mask_afni],
                                interface, wf, cfg, strat_pool, pipe_num, opt)

    outputs = {
        'blip_warp': (convert_afni_warp, 'ants_warp'),
        'blip_warp_inverse': None,  # TODO
        'desc-mean_bold': (undistort_func_mean, 'output_image'),
        'space-bold_desc-brain_mask':
            strat_pool.get_data('space-bold_desc-brain_mask')
    }

    return (wf, outputs)

