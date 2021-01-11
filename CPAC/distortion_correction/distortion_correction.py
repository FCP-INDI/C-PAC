#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import subprocess

import nibabel as nb

import nipype.pipeline.engine as pe
from nipype.interfaces import afni, fsl
import nipype.interfaces.utility as util
import nipype.interfaces.ants as ants

from CPAC.utils import function
from CPAC.utils.interfaces.function import Function
from CPAC.utils.datasource import match_epi_fmaps
from CPAC.func_preproc import skullstrip_functional


def createAFNIiterable(shrink_fac):
    expr = '-shrink_fac {0} '.format(shrink_fac)
    return expr


def create_EPI_DistCorr(use_BET, wf_name='epi_distcorr'):
    """
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
                                           unwarpped image)
    """

    preproc = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['anat_file',
                                                       'fmap_pha',
                                                       'fmap_mag']),
                        name='inputspec')

    inputNode_delTE = pe.Node(util.IdentityInterface(fields=['deltaTE']),
                              name='deltaTE_input')

    inputNode_dwellT = pe.Node(util.IdentityInterface(fields=['dwellT']),
                               name='dwellT_input')

    inputNode_dwell_asym_ratio = pe.Node(util.IdentityInterface(
        fields=['dwell_asym_ratio']
    ),
                                         name='dwell_asym_ratio_input')

    inputNode_bet_frac = pe.Node(util.IdentityInterface(fields=['bet_frac']),
                                 name='bet_frac_input')

    inputNode_afni_threshold = pe.Node(util.IdentityInterface(
        fields=['afni_threshold']
    ),
                                       name='afni_threshold_input')

    outputNode = pe.Node(util.IdentityInterface(fields=['fieldmap',
                                                        'fmap_despiked',
                                                        'fmapmagbrain',
                                                        'fieldmapmask']),
                         name='outputspec')

    # Skull-strip, outputs a masked image file
    if use_BET is False:
        skullstrip_args = pe.Node(util.Function(input_names=['shrink_fac'],
                                                output_names=['expr'],
                                                function=createAFNIiterable),
                                  name='distcorr_skullstrip_arg')

        preproc.connect(
            inputNode_afni_threshold, 'afni_threshold',
            skullstrip_args, 'shrink_fac')

        bet = pe.Node(interface=afni.SkullStrip(), name='bet')
        bet.inputs.outputtype = 'NIFTI_GZ'
        preproc.connect(skullstrip_args, 'expr', bet, 'args')
        preproc.connect(inputNode, 'fmap_mag', bet, 'in_file')
        preproc.connect(bet, 'out_file', outputNode, 'magnitude_image')
    else:
        bet = pe.Node(interface=fsl.BET(), name='bet')
        bet.inputs.output_type = 'NIFTI_GZ'
        preproc.connect(inputNode_bet_frac, 'bet_frac', bet, 'frac')
        preproc.connect(inputNode, 'fmap_mag', bet, 'in_file')
        preproc.connect(bet, 'out_file', outputNode, 'magnitude_image')

    # Prepare Fieldmap

    # prepare the field map
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(), name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"
    preproc.connect(inputNode_delTE, 'deltaTE', prepare, 'delta_TE')
    preproc.connect(inputNode, 'fmap_pha', prepare, 'in_phase')
    preproc.connect(bet, 'out_file', prepare, 'in_magnitude')
    preproc.connect(prepare, 'out_fieldmap', outputNode, 'fieldmap')

    # erode the masked magnitude image
    fslmath_mag = pe.Node(interface=fsl.ErodeImage(), name='fslmath_mag')
    preproc.connect(bet, 'out_file', fslmath_mag, 'in_file')
    preproc.connect(fslmath_mag, 'out_file', outputNode, 'fmapmagbrain')

    # calculate the absolute value of the eroded and masked magnitude
    # image
    fslmath_abs = pe.Node(interface=fsl.UnaryMaths(), name='fslmath_abs')
    fslmath_abs.inputs.operation = 'abs'
    preproc.connect(fslmath_mag, 'out_file', fslmath_abs, 'in_file')
    preproc.connect(fslmath_abs, 'out_file', outputNode, 'fmapmag_abs')

    # binarize the absolute value of the eroded and masked magnitude
    # image
    fslmath_bin = pe.Node(interface=fsl.UnaryMaths(), name='fslmath_bin')
    fslmath_bin.inputs.operation = 'bin'
    preproc.connect(fslmath_abs, 'out_file', fslmath_bin, 'in_file')
    preproc.connect(fslmath_bin, 'out_file', outputNode, 'fmapmag_bin')

    # take the absolute value of the fieldmap calculated in the prepare step
    fslmath_mask_1 = pe.Node(interface=fsl.UnaryMaths(),
                             name='fslmath_mask_1')
    fslmath_mask_1.inputs.operation = 'abs'
    preproc.connect(prepare, 'out_fieldmap', fslmath_mask_1, 'in_file')
    preproc.connect(fslmath_mask_1, 'out_file', outputNode, 'fieldmapmask_abs')

    # binarize the absolute value of the fieldmap calculated in the
    # prepare step
    fslmath_mask_2 = pe.Node(interface=fsl.UnaryMaths(),
                             name='fslmath_mask_2')
    fslmath_mask_2.inputs.operation = 'bin'
    preproc.connect(fslmath_mask_1, 'out_file', fslmath_mask_2, 'in_file')
    preproc.connect(fslmath_mask_2, 'out_file', outputNode, 'fieldmapmask_bin')

    # multiply together the binarized magnitude and fieldmap images
    fslmath_mask = pe.Node(interface=fsl.BinaryMaths(),
                           name='fslmath_mask')
    fslmath_mask.inputs.operation = 'mul'
    preproc.connect(fslmath_mask_2, 'out_file', fslmath_mask, 'in_file')
    preproc.connect(fslmath_bin, 'out_file', fslmath_mask, 'operand_file')
    preproc.connect(fslmath_mask, 'out_file', outputNode, 'fieldmapmask')

    # Note for the user. Ensure the phase image is within 0-4096 (upper
    # threshold is 90% of 4096), fsl_prepare_fieldmap will only work in the
    # case of the SIEMENS format. #Maybe we could use deltaTE also as an
    # option in the GUI.

    # fugue
    fugue1 = pe.Node(interface=fsl.FUGUE(), name='fugue1')
    fugue1.inputs.save_fmap = True
    fugue1.outputs.fmap_out_file = 'fmap_rads'
    preproc.connect(fslmath_mask, 'out_file', fugue1, 'mask_file')
    preproc.connect(inputNode_dwellT, 'dwellT', fugue1, 'dwell_time')
    preproc.connect(inputNode_dwell_asym_ratio, 'dwell_asym_ratio',
                    fugue1, 'dwell_to_asym_ratio')
    preproc.connect(prepare, 'out_fieldmap', fugue1, 'fmap_in_file')
    preproc.connect(fugue1, 'fmap_out_file', outputNode, 'fmap_despiked')

    return preproc


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


def blip_distcor_wf(wf_name='blip_distcor'):
    """Execute AFNI 3dQWarp to calculate the distortion "unwarp" for
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

    :param wf_name:
    :return:
    """

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['func_mean',
                                                        'opposite_pe_epi',
                                                        'same_pe_epi']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['blip_warp',
                                                         'blip_warp_inverse',
                                                         'new_func_mean',
                                                         'new_func_mask']),
                          name='outputspec')

    skullstrip_opposite_pe = skullstrip_functional(
        skullstrip_tool='afni', wf_name=f'{wf_name}_skullstrip_opp_pe')

    wf.connect(input_node, 'opposite_pe_epi',
               skullstrip_opposite_pe, 'inputspec.func')

    opp_pe_to_func = pe.Node(interface=fsl.FLIRT(), name='opp_pe_to_func')
    opp_pe_to_func.inputs.cost = 'corratio'

    wf.connect(skullstrip_opposite_pe, 'outputspec.func_brain',
               opp_pe_to_func, 'in_file')
    wf.connect(input_node, 'func_mean', opp_pe_to_func, 'reference')

    prep_qwarp_input_imports = ['import os', 'import subprocess']
    prep_qwarp_input = \
        pe.Node(function.Function(input_names=['same_pe_epi',
                                               'func_mean'],
                                  output_names=['qwarp_input'],
                                  function=same_pe_direction_prep,
                                  imports=prep_qwarp_input_imports),
                name='prep_qwarp_input')

    wf.connect(input_node, 'same_pe_epi', prep_qwarp_input, 'same_pe_epi')
    wf.connect(input_node, 'func_mean', prep_qwarp_input, 'func_mean')

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

    wf.connect(input_node, 'func_mean',
               undistort_func_mean, 'input_image')
    wf.connect(input_node, 'func_mean',
               undistort_func_mean, 'reference_image')
    wf.connect(convert_afni_warp, 'ants_warp',
               undistort_func_mean, 'transforms')

    create_new_mask = skullstrip_functional(skullstrip_tool='afni',
                                            wf_name="{0}_new_func_mask".format(
                                                wf_name))
    wf.connect(undistort_func_mean, 'output_image',
               create_new_mask, 'inputspec.func')

    wf.connect(convert_afni_warp, 'ants_warp',
               output_node, 'blip_warp')

    wf.connect(undistort_func_mean, 'output_image',
               output_node, 'new_func_mean')
    wf.connect(create_new_mask, 'outputspec.func_brain_mask',
               output_node, 'new_func_mask')

    return wf


def connect_distortion_correction(workflow, strat_list, c, diff, blip,
                                  fmap_rp_list, unique_id=None):

    # No need to fork if not running distortion_correction
    if True not in c.functional_preproc['distortion_correction']['run']:
        return (workflow, strat_list)

    # Distortion Correction
    new_strat_list = []

    # Distortion Correction - Field Map Phase-difference
    if "PhaseDiff" in c.functional_preproc[
        'distortion_correction'
    ]['using'] and diff:
        for num_strat, strat in enumerate(strat_list):
            if unique_id is None:
                workflow_name = f'diff_distcor_{num_strat}'
            else:
                workflow_name = f'diff_distcor_{unique_id}_{num_strat}'

            if 'BET' in c.functional_preproc[
                'distortion_correction'
            ]['PhaseDiff']['fmap_skullstrip_option']:
                epi_distcorr = create_EPI_DistCorr(
                    use_BET=True,
                    wf_name=workflow_name
                )
                epi_distcorr.inputs.bet_frac_input.bet_frac = \
                    c.functional_preproc[
                        'distortion_correction'
                    ]['PhaseDiff']['fmap_skullstrip_frac']
                epi_distcorr.get_node('bet_frac_input').iterables = \
                    ('bet_frac',
                     c.functional_preproc[
                         'distortion_correction'
                     ]['PhaseDiff']['fmap_skullstrip_frac'])
            else:
                epi_distcorr = create_EPI_DistCorr(
                    use_BET=False,
                    wf_name=workflow_name
                )
                epi_distcorr.inputs.afni_threshold_input.afni_threshold = \
                    c.functional_preproc[
                        'distortion_correction'
                    ]['PhaseDiff']['fmap_skullstrip_threshold']

            node, out_file = strat['anatomical_skull_leaf']
            workflow.connect(node, out_file, epi_distcorr,
                             'inputspec.anat_file')

            node, out_file = strat['diff_phase']
            workflow.connect(node, out_file, epi_distcorr,
                             'inputspec.fmap_pha')

            node, out_file = strat['diff_mag_one']
            workflow.connect(node, out_file, epi_distcorr,
                             'inputspec.fmap_mag')

            node, out_file = strat['deltaTE']
            workflow.connect(node, out_file, epi_distcorr,
                             'deltaTE_input.deltaTE')

            node, out_file = strat['diff_phase_dwell']
            workflow.connect(node, out_file, epi_distcorr,
                             'dwellT_input.dwellT')

            node, out_file = strat['dwell_asym_ratio']
            workflow.connect(node, out_file, epi_distcorr,
                             'dwell_asym_ratio_input.dwell_asym_ratio')

            if False in c.functional_preproc['distortion_correction']['run']:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(epi_distcorr.name)

            strat.update_resource_pool({
                'despiked_fieldmap': (epi_distcorr,
                                      'outputspec.fmap_despiked'),
                'fieldmap_mask': (epi_distcorr, 'outputspec.fieldmapmask'),
            })

    strat_list += new_strat_list

    # Distortion Correction - "Blip-Up / Blip-Down"
    if "Blip" in c.functional_preproc[
        'distortion_correction'
    ]['using'] and blip:
        for num_strat, strat in enumerate(strat_list):
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
                        name='match_epi_fmaps_{0}'.format(num_strat))

            if fmap_rp_list:
                epi_rp_key = fmap_rp_list[0]
                epi_param_rp_key = "{0}_scan_params".format(epi_rp_key)
                node, node_out = strat[epi_rp_key]
                workflow.connect(node, node_out,
                                 match_epi_fmaps_node, 'epi_fmap_one')
                node, node_out = strat[epi_param_rp_key]
                workflow.connect(node, node_out,
                                 match_epi_fmaps_node, 'epi_fmap_params_one')
                if len(epi_rp_key) > 1:
                    epi_rp_key = fmap_rp_list[1]
                    epi_param_rp_key = "{0}_scan_params".format(epi_rp_key)
                    node, node_out = strat[epi_rp_key]
                    workflow.connect(node, node_out,
                                     match_epi_fmaps_node, 'epi_fmap_two')
                    node, node_out = strat[epi_param_rp_key]
                    workflow.connect(node, node_out,
                                     match_epi_fmaps_node,
                                     'epi_fmap_params_two')

            node, node_out = strat['pe_direction']
            workflow.connect(node, node_out,
                             match_epi_fmaps_node, 'bold_pedir')

            if unique_id is None:
                workflow_name = f'blip_correct_{num_strat}'
            else:
                workflow_name = f'blip_correct_{unique_id}_{num_strat}'

            blip_correct = blip_distcor_wf(
                wf_name=workflow_name)

            node, out_file = strat["mean_functional"]
            workflow.connect(node, out_file,
                             blip_correct, 'inputspec.func_mean')

            workflow.connect(match_epi_fmaps_node, 'opposite_pe_epi',
                             blip_correct, 'inputspec.opposite_pe_epi')

            workflow.connect(match_epi_fmaps_node, 'same_pe_epi',
                             blip_correct, 'inputspec.same_pe_epi')

            if False in c.functional_preproc['distortion_correction']['run']:
                strat = strat.fork()
                new_strat_list.append(strat)

            strat.append_name(blip_correct.name)

            strat.update_resource_pool({
                'blip_warp': (blip_correct, 'outputspec.blip_warp'),
                'blip_warp_inverse': (blip_correct,
                                      'outputspec.blip_warp_inverse'),
                'mean_functional': (blip_correct, 'outputspec.new_func_mean'),
                'functional_brain_mask': (blip_correct,
                                          'outputspec.new_func_mask')
            }, override=True)

    strat_list += new_strat_list

    return (workflow, strat_list)
