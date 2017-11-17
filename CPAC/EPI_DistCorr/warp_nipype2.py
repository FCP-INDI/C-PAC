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
def warp_nipype():
    preproc = pe.Workflow(name='preprocflow')
                          
    inputNode = pe.Node(util.IdentityInterface(fields=['anat_file','func_file','fmap_pha','fmap_mag','t1_head','mask']),name = 'inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['roi_file','fieldmap','epireg','fmap_despiked','shiftfile']),name='outputspec')

## Specify commands to be run
# Extract first three volumes from fmri
    fslroi = pe.Node(interface=fsl.ExtractROI(),name='fslroi')
    
    fslroi.inputs.t_min=0
    fslroi.inputs.t_size=3
    preproc.connect(inputNode,'anat_file',fslroi,'in_file') 
    preproc.connect(fslroi,'roi_file',outputNode,'fslroi_file')
# Skullstrip

    skullstrip = pe.Node(interface=afni.preprocess.SkullStrip(),name='skullstrip')
    skullstrip.inputs.outputtype='NIFTI_GZ'

    preproc.connect(inputNode,'fmap_mag',skullstrip,'in_file')


#SkullStrip the anatomy file
    skullstrip_anat = pe.Node(interface=afni.preprocess.SkullStrip(),name='skullstrip_anat')
    skullstrip_anat.inputs.outputtype='NIFTI_GZ'

    preproc.connect(inputNode,'anat_file',skullstrip_anat,'in_file')
    preproc.connect(skullstrip_anat,'out_file',outputNode,'stripped_anat')
#Note for the user. Ensure the phase image is within 0-4096 (upper threshold is 90% of 4096), fsl_prepare_fieldmap will only work 
#in the case of the SIEMENS format. #Maybe we could use deltaTE also as an option in the GUI. 
# Prepare Fieldmap
    prepare = pe.Node(interface=fsl.epi.PrepareFieldmap(),name='prepare')
    prepare.inputs.output_type = "NIFTI_GZ"
    prepare.inputs.delta_TE = 2.46

    preproc.connect(inputNode,'fmap_pha',prepare,'in_phase')
    preproc.connect(inputNode,'fmap_mag',prepare,'in_magnitude')
    preproc.connect(prepare,'out_fieldmap',outputNode,'fieldmap')

#perform fugue to compute vsm 
    fugue = pe.Node(interface = fsl.FUGUE(),name='fugue')
    fugue.inputs.unwarp_direction='y'
    fugue.inputs.save_unmasked_fmap=True
    fugue.inputs.shift_out_file = 'opname_shiftout_file'
    fugue.inputs.dwell_to_asym_ratio=0.93902439
    preproc.connect(inputNode, 'fmap_pha', fugue, 'phasemap_in_file')
    preproc.connect(inputNode,'mask',fugue,'mask_file')
    preproc.connect(prepare,'out_fieldmap',fugue,'in_file')
    preproc.connect(fugue, 'shift_out_file',outputNode,'shiftfile')
    
#perform fugue to regularise fieldmap, you can also provide this as a GUI option:
    ##to perform gausian smoothing, despiking, and median filtering
    fugue1= pe.Node(interface = fsl.FUGUE(),name='fugue1')
    fugue1.inputs.save_fmap=True
    fugue1.inputs.despike_2dfilter=True
    fugue1.outputs.fmap_out_file='fmap_despiked'
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
    preproc.connect(inputNode, 't1_head', epireg, 't1_head')
    preproc.connect(inputNode, 'fmap_mag',epireg, 'fmapmag')
    preproc.connect(skullstrip,'out_file',epireg,'fmapmagbrain')
    preproc.connect(fugue1, 'fmap_out_file', epireg, 'fmap')
    preproc.connect(epireg,'out_file',outputNode,'epireg')
    return preproc
