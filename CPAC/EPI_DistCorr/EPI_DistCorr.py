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
                                    
    -- SkullStrip: SkullStrip (or BET) is used to strip the non-brain (tissue) regions
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
    
  #  inputnode_delTE = pe.Node(util.IdentityInterface(fields = ['delTE']), name = 'input_delTE')
    
 #   inputnode_dwellT = pe.Node(util.IdentityInterface(fields = ['dwellTE']), name = 'input_dwellT')
    
 #   inputnode_asymR = pe.Node(util.IdentityInterface(fields = ['asymR']),name = 'input_asymR')
    
    outputNode = pe.Node(util.IdentityInterface(fields=['func_file','roi_file','fieldmap','epireg','fmap_despiked','shiftfile']),name='outputspec')

## Specify commands to be run
# Extract first three volumes from fmri
    fslroi = pe.Node(interface=fsl.ExtractROI(),name='fslroi')
    
    fslroi.inputs.t_min=0
    fslroi.inputs.t_size=3
    preproc.connect(inputNode,'anat_file',fslroi,'in_file') 
    preproc.connect(fslroi,'roi_file',outputNode,'fslroi_file')
# Skullstrip

    skullstrip = pe.Node(interface=fsl.BET(),name='skullstrip')
    skullstrip.inputs.output_type = 'NIFTI_GZ'
    skullstrip.inputs.frac = 1.0
    preproc.connect(inputNode,'fmap_mag',skullstrip,'in_file')
    preproc.connect(skullstrip,'out_file',outputNode,'magnitude_image')

#SkullStrip the anatomy file
    skullstrip_anat = pe.Node(interface=fsl.BET(),name='skullstrip_anat')
    skullstrip_anat.inputs.output_type = 'NIFTI_GZ'
    skullstrip_anat.inputs.frac = 0.5
    preproc.connect(inputNode,'anat_file',skullstrip_anat,'in_file')
    preproc.connect(skullstrip_anat,'out_file',outputNode,'stripped_anat')
#Note for the user. Ensure the phase image is within 0-4096 (upper threshold is 90% of 4096), fsl_prepare_fieldmap will only work 
#in the case of the SIEMENS format. #Maybe we could use deltaTE also as an option in the GUI. 
# Prepare Fieldmap
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(),name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"
    prepare.inputs.delta_TE = 2.46
   # preproc.connect(inputnode_delTE, 'input_delTE', prepare, 'delta_TE')
    preproc.connect(inputNode,'fmap_pha',prepare,'in_phase')
    preproc.connect(skullstrip,'out_file',prepare,'in_magnitude')
    preproc.connect(prepare,'out_fieldmap',outputNode,'fieldmap')
#fugue 
    fugue1= pe.Node(interface = fsl.FUGUE(),name='fugue1')
    fugue1.inputs.save_fmap=True
    fugue1.inputs.despike_2dfilter=True
    fugue1.outputs.fmap_out_file='fmap_despiked'
    fugue1.inputs.dwell_time = 0.00231
    fugue1.inputs.dwell_to_asym_ratio=0.93902439 
    #preproc.connect(inputnode_dwellT, 'input_dwelT', fugue1, 'dwell_time')
    #preproc.connect(inputnode_asymR, 'input_asymR', fugue1, 'dwell_to_asym_ratio')
    preproc.connect(prepare,'out_fieldmap',fugue1,'fmap_in_file')
    preproc.connect(fugue1,'fmap_out_file',outputNode,'fmap_despiked')
# Co-Register EPI and Correct field inhomogeniety distortions
# echospacing can also be in the GUI
    epireg = pe.Node(interface=fsl.epi.EpiReg(), name='epireg')
    epireg.inputs.echospacing=0.00046
    epireg.inputs.pedir='y'
    epireg.inputs.output_type='NIFTI_GZ'
    preproc.connect(inputNode,'func_file',epireg,'epi')
    preproc.connect(skullstrip_anat,'out_file', epireg,'t1_brain')
    preproc.connect(inputNode, 'anat_file', epireg, 't1_head')
    preproc.connect(inputNode, 'fmap_mag',epireg, 'fmapmag')
    preproc.connect(skullstrip,'out_file',epireg,'fmapmagbrain')
    preproc.connect(fugue1, 'fmap_out_file', epireg, 'fmap')
    preproc.connect(epireg,'out_file',outputNode,'epireg')
    preproc.connect(inputNode, 'func_file',outputNode,'func_file')
    return preproc
