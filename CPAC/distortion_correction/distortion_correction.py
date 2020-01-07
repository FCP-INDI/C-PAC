#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import os
import subprocess

import nibabel as nb

import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl
import nipype.interfaces.utility as util

from CPAC.utils import function
from CPAC.func_preproc import skullstrip_functional
from CPAC.registration.registration import create_wf_apply_ants_warp


def createAFNIiterable(shrink_fac):
    expr = '-shrink_fac {0} '.format(shrink_fac)
    return expr


def create_EPI_DistCorr(use_BET,wf_name = 'epi_distcorr'):
    """
    Fieldmap correction takes in an input magnitude image which is Skull Stripped (Tight).
    The magnitude images are obtained from each echo series. It also requires a phase image
    as an input, the phase image is a subtraction of the two phase images from each echo.

    Created on Thu Nov  9 10:44:47 2017
    @author: nrajamani

    Order of commands and inputs:

    -- SkullStrip:   3d-SkullStrip (or FSL-BET) is used to strip the non-brain (tissue) regions
                     from the fMRI
                     Parameters: -f, default: 0.5
                     in_file: fmap_mag
    -- fslmath_mag:  Magnitude image is eroded using the -ero option in fslmath, in order to remove
                     the non-zero voxels
                     Parameters: -ero
                     in_file:fmap_mag
    -- bet_anat   :  Brain extraction of the anat file to provide as an input for the epi-registration interface
                     Parameters: -f, default: 0.5
                     in_file: anat_file
    -- fast_anat  :  Fast segmentation to provide partial volume files of the anat file, which is further processed
                     to provide the white mater segmentation input for the epi-registration interface.
                     The most important output required from this is the second segment, (e.g.,'T1_brain_pve_2.nii.gz')
                     Parameters: -img_type = 1
                               -bias_iters = 10 (-I)
                               -bias_lowpass = 10 (-l)
                     in_file: brain_extracted anat_file
    -- fslmath_anat: The output of the FAST interface is then analyzed to select all the voxels with more than 50%
                     partial volume into the binary mask
                     Parameters: -thr = 0.5
                     in_file : T1_brain_pve_2
    -- fslmath_wmseg:The selected voxels are now used to create a binary mask which would can then be sent as the
                     white matter segmentation (wm_seg)
                     Parameters: -bin
                     in_file: T1_brain_pve_2
    -- Prepare      :Preparing the fieldmap.
                     Parameters: -deltaTE = default, 2.46 ms
                                 -Scanner = SIEMENS
                     in_files: fmap_phase
                               fmap_magnitude
                     For more details, check:https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide
    -- FUGUE        :One of the steps in EPI-DistCorrection toolbox, it unwarps the fieldmaps
                     Parameters: dwell_to_asymm ratio = (0.77e-3 * 3)/(2.46e-3)
                                 dwell time = 0.0005 ms
                                 in_file = field map which is a 4D image (containing 2 unwarpped image)
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
    
    inputNode_dwell_asym_ratio = pe.Node(util.IdentityInterface(fields=['dwell_asym_ratio']),
                                         name='dwell_asym_ratio_input')
    
    inputNode_bet_frac = pe.Node(util.IdentityInterface(fields=['bet_frac']),
                                 name='bet_frac_input')

    inputNode_afni_threshold = pe.Node(util.IdentityInterface(fields=['afni_threshold']),
        name='afni_threshold_input')
    
    
    outputNode = pe.Node(util.IdentityInterface(fields=['fieldmap',
                                                        'fmap_despiked',
                                                        'fmapmagbrain',
                                                        'fieldmapmask']),
                         name='outputspec')
    
    # Skull-strip
    if use_BET == False:
        skullstrip_args = pe.Node(util.Function(input_names=['shrink_fac'],output_names=['expr'],function=createAFNIiterable),name='distcorr_skullstrip_arg')
        preproc.connect(inputNode_afni_threshold,'afni_threshold',skullstrip_args,'shrink_fac')

        bet = pe.Node(interface=afni.SkullStrip(),name='bet')
        bet.inputs.outputtype = 'NIFTI_GZ'
        preproc.connect(skullstrip_args,'expr', bet, 'args')
        preproc.connect(inputNode, 'fmap_mag', bet, 'in_file')
        preproc.connect(bet, 'out_file', outputNode, 'magnitude_image')
    else:
        bet = pe.Node(interface=fsl.BET(),name='bet')
        bet.inputs.output_type='NIFTI_GZ'
        preproc.connect(inputNode_bet_frac, 'bet_frac', bet, 'frac')
        preproc.connect(inputNode, 'fmap_mag', bet, 'in_file')
        preproc.connect(bet, 'out_file', outputNode, 'magnitude_image')

    # Prepare Fieldmap
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(), name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"
    preproc.connect(inputNode_delTE, 'deltaTE', prepare, 'delta_TE')
    preproc.connect(inputNode, 'fmap_pha', prepare, 'in_phase')
    preproc.connect(bet, 'out_file', prepare, 'in_magnitude')
    preproc.connect(prepare, 'out_fieldmap', outputNode, 'fieldmap')

    fslmath_mag = pe.Node(interface=fsl.ErodeImage(),name='fslmath_mag')
    preproc.connect(bet,'out_file',fslmath_mag,'in_file')
    preproc.connect(fslmath_mag,'out_file',outputNode,'fmapmagbrain')

    ##generating mask-step 1#
    fslmath_abs = pe.Node(interface=fsl.UnaryMaths(),name = 'fslmath_abs')
    fslmath_abs.inputs.operation = 'abs'
    preproc.connect(fslmath_mag,'out_file',fslmath_abs,'in_file')
    preproc.connect(fslmath_abs,'out_file',outputNode,'fmapmag_abs')

    #generating mask-step 2#
    fslmath_bin = pe.Node(interface=fsl.UnaryMaths(),name='fslmath_bin')
    fslmath_bin.inputs.operation = 'bin'
    preproc.connect(fslmath_abs,'out_file',fslmath_bin,'in_file')
    preproc.connect(fslmath_bin,'out_file',outputNode,'fmapmag_bin')

    #generating mask-step 3#
    fslmath_mask_1 = pe.Node(interface=fsl.UnaryMaths(),name = 'fslmath_mask_1')
    fslmath_mask_1.inputs.operation = 'abs'
    preproc.connect(prepare,'out_fieldmap',fslmath_mask_1,'in_file')
    preproc.connect(fslmath_mask_1,'out_file',outputNode,'fieldmapmask_abs')

    #generating mask-step 4#
    fslmath_mask_2 = pe.Node(interface=fsl.UnaryMaths(),name = 'fslmath_mask_2')
    fslmath_mask_2.inputs.operation = 'bin'
    preproc.connect(fslmath_mask_1,'out_file',fslmath_mask_2,'in_file')
    preproc.connect(fslmath_mask_2,'out_file',outputNode,'fieldmapmask_bin')

    #generating mask-step 5#
    fslmath_mask = pe.Node(interface=fsl.BinaryMaths(),name='fslmath_mask')
    fslmath_mask.inputs.operation = 'mul'
    preproc.connect(fslmath_mask_2,'out_file',fslmath_mask,'in_file')
    preproc.connect(fslmath_bin,'out_file',fslmath_mask,'operand_file')
    preproc.connect(fslmath_mask,'out_file',outputNode,'fieldmapmask')

    # Note for the user. Ensure the phase image is within 0-4096 (upper
    # threshold is 90% of 4096), fsl_prepare_fieldmap will only work in the
    # case of the SIEMENS format. #Maybe we could use deltaTE also as an
    # option in the GUI.

    # fugue
    fugue1 = pe.Node(interface=fsl.FUGUE(), name='fugue1')
    fugue1.inputs.save_fmap = True
    fugue1.outputs.fmap_out_file = 'fmap_rads'
    preproc.connect(fslmath_mask,'out_file', fugue1, 'mask_file')
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
                                          "{0}_mask.nii.gz".format(os.path.basename(same_pe_epi)))
        skullstrip_cmd = ["3dAutomask", "-apply_prefix",
                          "{0}_masked.nii.gz".format(os.path.basename(same_pe_epi)),
                          "-prefix", skullstrip_outfile, same_pe_epi]
        retcode = subprocess.check_output(skullstrip_cmd)

        extract_brain_outfile = os.path.join(os.getcwd(),
                                             "{0}_calc.nii.gz".format(os.path.basename(same_pe_epi)))
        extract_brain_cmd = ["3dcalc", "-a", same_pe_epi, "-b",
                             skullstrip_outfile, "-expr", "a*b", "-prefix",
                             extract_brain_outfile]
        retcode = subprocess.check_output(extract_brain_cmd)

        align_epi_outfile = os.path.join(os.getcwd(),
                                         "{0}_calc_flirt.nii.gz".format(os.path.basename(same_pe_epi)))
        align_epi_cmd = ["flirt", "-in", extract_brain_outfile, "-ref",
                         func_mean, "-out", align_epi_outfile, "-cost",
                         "corratio"]
        retcode = subprocess.check_output(align_epi_cmd)

        qwarp_input = align_epi_outfile

    return qwarp_input


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
    """Execute AFNI 3dQWarp to calculate the distortion "unwarp" for phase
    encoding direction EPI field map distortion correction.

        1. Skull-strip the opposite-direction phase encoding EPI.
        2. Transform the opposite-direction phase encoding EPI to the
           skull-stripped functional and pass this as the base_file to
           AFNI 3dQWarp (plus-minus).
        3. If there is a same-direction phase encoding EPI, also skull-strip
           this, and transform it to the skull-stripped functional. Then, pass
           this as the in_file to AFNI 3dQWarp (plus-minus).
        4. If there isn't a same-direction, pass the functional in as the
           in_file of AFNI 3dQWarp (plus-minus).
        5. Convert the 3dQWarp transforms to ANTs/ITK format.
        6. Use antsApplyTransforms, with the original functional as both the
           input and the reference, and apply the warp from 3dQWarp. The
           output of this can then proceed to func_preproc.

    :param wf_name:
    :return:
    """

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['func_mean',
                                                        'opposite_pe_epi',
                                                        'same_pe_epi']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['blip_warp',
                                                         'new_func_mean',
                                                         'new_func_mask']),
                          name='outputspec')

    skullstrip_opposite_pe = skullstrip_functional("afni", False,
                                                   "{0}_skullstrip_opp_pe".format(wf_name))

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

    calc_blip_warp = pe.Node(afni.QwarpPlusMinus(), name='calc_blip_warp')
    calc_blip_warp.inputs.plusminus = True
    calc_blip_warp.inputs.outputtype = "NIFTI_GZ"
    calc_blip_warp.inputs.out_file = os.path.abspath("Qwarp.nii.gz")

    wf.connect(opp_pe_to_func, 'out_file', calc_blip_warp, 'base_file')
    wf.connect(prep_qwarp_input, 'qwarp_input', calc_blip_warp, 'in_file')

    convert_afni_warp_imports = ['import os', 'import nibabel as nb']
    convert_afni_warp = \
        pe.Node(function.Function(input_names=['afni_warp'],
                                  output_names=['ants_warp'],
                                  function=convert_afni_to_ants,
                                  imports=convert_afni_warp_imports),
                name='convert_afni_warp')

    wf.connect(calc_blip_warp, 'source_warp', convert_afni_warp, 'afni_warp')

    undistort_func_mean = create_wf_apply_ants_warp()
    undistort_func_mean.inputs.inputspec.interpolation = "LanczosWindowedSinc"

    wf.connect(input_node, 'func_mean',
               undistort_func_mean, 'inputspec.input_image')
    wf.connect(input_node, 'func_mean',
               undistort_func_mean, 'inputspec.reference_image')
    wf.connect(convert_afni_warp, 'ants_warp',
               undistort_func_mean, 'inputspec.transforms')

    create_new_mask = skullstrip_functional("afni", False,
                                            "{0}_new_func_mask".format(wf_name))
    wf.connect(undistort_func_mean, 'outputspec.output_image',
               create_new_mask, 'inputspec.func')

    wf.connect(convert_afni_warp, 'ants_warp', output_node, 'blip_warp')
    wf.connect(undistort_func_mean, 'outputspec.output_image',
               output_node, 'new_func_mean')
    wf.connect(create_new_mask, 'outputspec.func_brain_mask',
               output_node, 'new_func_mask')

    return wf