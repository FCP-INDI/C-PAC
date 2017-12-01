#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 10:36:37 2017

@author: nrajamani
"""
import os
import glob
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe
from nipype.interfaces.fsl import ImageStats
#rom nose.tools import *
## This script runs the warp_nipype workflow to execute the interfaces and with the inputs already provided
## It accepts as it's arguments the input whihc you;d give the warp_nipe file, an output_dir which can be either be mentioned or if it is set to none will write it in the current working directory.
## The argument run can either be tset to true(default) or to false. If set to false, it should connect to the nipype workflow and return the workflow object instead
##What all should it return?: if the run had been set to true, it would generate the filepath of the output from the workflow. If the run was set to false, then it will return the file path of the 
##base directory, the workflow object, etc.
def run_warp_nipype(inputs,output_dir=None,run=True):
   import EPI_distCorr
   warp_workflow = pe.Workflow(name = 'preproc')
   if output_dir == None:
     output_dir = '/home/nrajamani'
        
   workflow_dir = os.path.join(output_dir,"Workflow_output")
   warp_workflow.base_dir = workflow_dir
   # taken from QAP files 
   #resource_pool = {}
   
   num_of_cores = 1
   #resource_pool({'epireg': (warp_nipype2.warp_nipype, 'outputspec.epireg')})
   t_node = EPI_DistCorr.create_EPI_DistCorr()####
   t_node.inputs.inputspec.anat_file=  '/home/nrajamani/sub-A00073677/ses-NFB3/anat/sub-A00073677_ses-NFB3_T1w.nii.gz'
   t_node.inputs.inputspec.func_file= '/home/nrajamani/FieldMap_SubjectExampleData/SubjectData/epi_run1/fMT0160-0014-00002-000002-01.nii' 
   t_node.inputs.inputspec.fmap_pha= '/home/nrajamani/sub-A00073677/ses-NFB3/fmap/sub-A00073677_ses-NFB3_magnitude1.nii.gz'
   t_node.inputs.inputspec.t1_head = '/usr/local/fsl/data/standard/MNI152lin_T1_2mm.nii.gz'
   #'home/nrajamani/FieldMap_SubjectExampleData/SubjectData/epi_run2/fMT0160-0015-00002-000002-01_BRAIN.nii.gz',
   t_node.inputs.inputspec.fmap_mag= '/home/nrajamani/sub-A00073677/ses-NFB3/fmap/sub-A00073677_ses-NFB3_phasediff.nii.gz'
   #'home/nrajamani/FieldMap_SubjectExampleData/SubjectData/epi_run2/fMT0160-0015-00003-000003-01_BRAIN.nii.gz',
   t_node.inputs.inputspec.mask ='/usr/local/fsl/data/standard/MNI152lin_T1_2mm_brain_mask.nii.gz'
   
#   for image in inputs:
#       if not(image.endswith('.nii') or image.endswith('.nii.gz')):
#           raise 'The input image is not the right format'
#   try:
#       for image in inputs:
#           size = image.get_shape()
#           assert len(size) == 3
#   except:
#       if len(size) < 3:
#           raise 'input image is not 3D'
#   intensity = ImageStats(in_file = t_node.inputs.inputspec.fmap_pha, op_string = '-p 90')
#   if intensity < 3686:
#      raise 'input phase image does not have the correct range values'         
   dataSink = pe.Node(nio.DataSink(), name='dataSink_file')
   dataSink.inputs.base_directory = workflow_dir
   #node, out_file = resource_pool["epireg"]
   warp_workflow.connect(t_node,'outputspec.roi_file',dataSink,'roi_file')
   warp_workflow.connect(t_node,'outputspec.fieldmap',dataSink,'fieldmap_file')
   warp_workflow.connect(t_node,'outputspec.epireg',dataSink,'epireg_file')
   if run == True:
       warp_workflow.run(plugin='MultiProc', plugin_args ={'n_procs': num_of_cores})
       #outpath = glob.glob(os.path.join(workflow_dir, "EPI_DistCorr","*"))[0]
       #return outpath
   else:
       return warp_workflow, warp_workflow.base_dir

run_warp_nipype(['anat_file','func_file','fmap_pha','t1_head', 'mask'],output_dir=None,run=True)  
    
    
    
    
    
    
