#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 10:44:47 2017

@author: nrajamani
"""

import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl
import nipype.interfaces.utility as util

def createAFNIiterable(shrink_fac):
    expr = '-shrink_fac {0} '.format(shrink_fac)
    return expr


def create_EPI_DistCorr(use_BET,wf_name = 'epi_distcorr'):
    """
    Fieldmap correction takes in an input magnitude image which is Skull Stripped (Tight).
    The magnitude images are obtained from each echo series. It also requires a phase image
    as an input, the phase image is a subtraction of the two phase images from each echo.

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
                                                       'func_file',
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


def blip_distcorr_wf():
    wf = pe.Workflow(name=wf_name)
    return wf