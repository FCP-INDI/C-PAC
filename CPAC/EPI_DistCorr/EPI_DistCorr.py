#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 10:44:47 2017

@author: nrajamani
"""

import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl
import nipype.interfaces.utility as util



"""
  Fieldmap correction takes in an input magnitude image which is Skull Stripped (Tight).
  The magnitude images are obtained from each echo series. It also requires a phase image
  as an input, the phase image is a subtraction of the two phase images from each echo.
  
  Order of commands and inputs:
    -- Extraction of brain volumes: using first three columes from the fMRI
                                    
    -- SkullStrip: 3d-SkullStrip (or BET) is used to strip the non-brain (tissue) regions
                   from the fMRI
                   
    -- PrepareFieldMap: Preparing the fieldmap. 
                           For more details, check:https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FUGUE/Guide
                           in_file = phase_file, magnitude_file
                           deltaTE = Delay in the time taken by the scanners to record the two echos.
                                     (default == 2.46 ms)
                           scanner = SIEMENS, default.
    -- FUGUE           : One of the steps in EPI-DistCorrection toolbox, it unwarps the fieldmaps 
                            in_file = field map which is a 4D image (containing 2 unwarpped image)
                            mask_file = epi_mask.nii
                            dwell_to_asymm ratio = (0.77e-3 * 3)/(2.46e-3)
  """                          
def create_EPI_DistCorr(wf_name = 'epi_distcorr'):
    preproc = pe.Workflow(name=wf_name)
                          
    inputNode = pe.Node(util.IdentityInterface(fields=['anat_file','func_file','fmap_pha','fmap_mag']),name = 'inputspec')
    
    inputNode_delTE = pe.Node(util.IdentityInterface(fields=['delTE']),name='input_delTE')
    
    inputNode_dwellT = pe.Node(util.IdentityInterface(fields=['dwellT']),name = 'input_dwellT')
    
    inputNode_dwell_asym_ratio = pe.Node(util.IdentityInterface(fields=['dwell_asym_ratio']),name='input_dwell_asym_ratio')
    
    inputNode_bet_frac = pe.Node(util.IdentityInterface(fields=['bet_frac']),name='input_bet_frac')
    
    inputNode_skullstrip_method = pe.Node(util.IdentityInterface(fields=['FSL-BET','AFNI-3dSkullStrip']),name = 'input_skullstrip_method')
    
    outputNode = pe.Node(util.IdentityInterface(fields=['func_file','fieldmap','epireg','fmap_despiked','partial_volume_files','partial_volume_map','threshold_image','fmapmagbrain','T1_wm_seg']),name='outputspec')

## Specify commands to be run
# Extract first three volumes from fmri
#fslroi = pe.Node(interface=fsl.ExtractROI(),name='fslroi')
    
    #  fslroi.inputs.t_min=0
    #fslroi.inputs.t_size=2
    #fslroi.inputs.roi_file = 'roi_file.nii.gz'
    #preproc.connect(inputNode,'anat_file',fslroi,'in_file')
    #preproc.connect(fslroi,'roi_file',outputNode,'roi_file')
##Extract the magnitude file for fmapmagbrain##
#fslroi_mag = pe.Node(interface=fsl.ExtractROI(),name='fslroi_mag')
# fslroi_mag.inputs.t_min=0
#    fslroi_mag.inputs.t_size=3
#    preproc.connect(inputNode, 'fmap_pha', fslroi_mag, 'in_file')
#    preproc.connect(fslroi_mag, 'roi_file',outputNode,'roi_file_mag')
# Skullstrip
    if (inputNode_skullstrip_method,skullstrip_method == 'AFNI-3dSkullStrip'):
        bet = pe.Node(interface=afni.SkullStrip(),name='bet')
        bet.inputs.outputtype = 'NIFTI_GZ'
        bet.inputs.args = '-shrink_fac SF'
        preproc.connect(inputNode,str('bet_frac'),bet,'args')
        preproc.connect(inputNode,'fmap_mag',bet,'in_file')
        preproc.connect(bet,'out_file',outputNode,'magnitude_image')
    
    
    elif(inputNode_skullstrip_method == 'FSL-BET'):
        bet = pe.Node(interface=fsl.BET(),name='bet')
        bet.inputs.output_type='NIFTI_GZ'
        preproc.connect(inputNode,'bet_frac',bet,'frac')
        preproc.connect(inputNode,'fmap_mag',bet,'in_file')
        preproc.connect(bet,'out_file',outputNode,'magnitude_image')
#generating fmap_mag_brain#
    fslmath_mag = pe.Node(interface=fsl.ErodeImage(),name='fslmath_mag')
    fslmath_mag.inputs.args = '-ero'
    preproc.connect(bet,'out_file',fslmath_mag,'in_file')
    preproc.connect(fslmath_mag,'out_file',outputNode,'fmapmagbrain')
#SkullStrip the anatomy file
    bet_anat = pe.Node(interface=fsl.BET(),name='bet_anat')
    bet_anat.inputs.output_type = 'NIFTI_GZ'
    bet_anat.inputs.frac = 0.5
    preproc.connect(inputNode,'anat_file',bet_anat,'in_file')
    preproc.connect(bet_anat,'out_file',outputNode,'stripped_anat')
    
    fast_anat = pe.Node(interface=fsl.FAST(),name='fast_anat')
    fast_anat.inputs.output_biascorrected=True
    fast_anat.inputs.img_type=1
    fast_anat.inputs.bias_iters=10
    fast_anat.inputs.bias_lowpass=10
    fast_anat.inputs.segments=True
    fast_anat.outputs_basename='T1'
    preproc.connect(bet_anat,'out_file',fast_anat,'in_files')
    preproc.connect(fast_anat,'partial_volume_map',outputNode,'partial_volume_map')
    preproc.connect(fast_anat,'partial_volume_files',outputNode,'partial_volume_files')
    
    fslmath_anat = pe.Node(interface=fsl.Threshold(),name='fsl_anat')
    fslmath_anat.inputs.thresh = 0.5
    preproc.connect(fast_anat,'T1_brain',fslmath_anat,'in_file')
    preproc.connect(fslmath_anat,'out_file',outputNode,'threshold_image')
    
    fslmath_wmseg = pe.Node(interface=fsl.UnaryMaths(),name='fslmath_wmseg')
    fslmath_wmseg.inputs.operation = 'bin'
    preproc.connect(fslmath_anat,'out_file',fslmath_wmseg,'in_file')
    preproc.connect(fslmath_wmseg,'out_file',outputNode,'T1_wm_seg')
#Note for the user. Ensure the phase image is within 0-4096 (upper threshold is 90% of 4096), fsl_prepare_fieldmap will only work 
#in the case of the SIEMENS format. #Maybe we could use deltaTE also as an option in the GUI. 
# Prepare Fieldmap
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(),name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"
    #prepare.inputs.delta_TE = 2.46
    preproc.connect(inputNode, 'delTE', prepare, 'delta_TE')
    preproc.connect(inputNode,'fmap_pha',prepare,'in_phase')
    preproc.connect(bet,'out_file',prepare,'in_magnitude')
    preproc.connect(prepare,'out_fieldmap',outputNode,'fieldmap')
#fugue 
    fugue1= pe.Node(interface = fsl.FUGUE(),name='fugue1')
    fugue1.inputs.save_fmap=True
    fugue1.inputs.despike_2dfilter=True
    fugue1.outputs.fmap_out_file='fmap_despiked'
    #fugue1.inputs.dwell_time = 0.0005
    #fugue1.inputs.dwell_to_asym_ratio= 0.9390243902
    preproc.connect(inputNode, 'dwellT', fugue1, 'dwell_time')
    preproc.connect(inputNode, 'dwell_asym_ratio', fugue1, 'dwell_to_asym_ratio')
    preproc.connect(prepare,'out_fieldmap',fugue1,'fmap_in_file')
    preproc.connect(fugue1,'fmap_out_file',outputNode,'fmap_despiked')
# Co-Register EPI and Correct field inhomogeniety distortions
# echospacing can also be in the GUI
    epireg = pe.Node(interface=fsl.epi.EpiReg(), name='epireg')
    #epireg.inputs.echospacing=0.0005
    epireg.inputs.pedir='-y'
    epireg.inputs.output_type='NIFTI_GZ'
    preproc.connect(inputNode,'dwellT',epireg,'echospacing')
    preproc.connect(fslmath_wmseg,'out_file',epireg,'wmseg')
    preproc.connect(inputNode,'func_file',epireg,'epi')
    preproc.connect(bet_anat,'out_file', epireg,'t1_brain')
    preproc.connect(inputNode, 'anat_file', epireg, 't1_head')
    preproc.connect(inputNode, 'fmap_mag',epireg, 'fmapmag')
    preproc.connect(fslmath_mag,'out_file',epireg,'fmapmagbrain')
    preproc.connect(fugue1, 'fmap_out_file', epireg, 'fmap')
    preproc.connect(epireg,'out_file',outputNode,'epireg')
    preproc.connect(inputNode, 'func_file',outputNode,'func_file')
    return preproc
