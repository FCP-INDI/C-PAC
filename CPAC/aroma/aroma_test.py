#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This script runs the warp_nipype workflow to execute the interfaces
and with the inputs already provided.

It accepts as its arguments the input whihc you'd give the warp_nipe
file, an output_dir which can be either be mentioned or if it is set to
None will write it in the current working directory.

The argument run can either be set to true(default) or to false. If
set to False, it should connect to the nipype workflow and return the
workflow object instead.

What all should it return?:
    If the run had been set to true, it would generate the filepath of
    the output from the workflow. If the run was set to false, then it
    will return the file path of the base directory, the workflow
    object, etc.

Created on Thu Nov  9 10:36:37 2017

@author: nrajamani
"""
import glob
import os
import pytest

import nipype.interfaces.io as nio
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nipype_pipeline_engine.plugins import MultiProcPlugin
from nipype.interfaces.fsl import ImageStats
# from nose.tools import *


@pytest.mark.skip(reason='needs refactoring')
def run_warp_nipype(inputs, output_dir=None, run=True):
    import EPI_DistCorr
    warp_workflow = pe.Workflow(name='preproc')
    if output_dir is None:
        output_dir = '/home/nrajamani'

    workflow_dir = os.path.join(output_dir,
                                "workflow_output_with_aroma_with_change")
    warp_workflow.base_dir = workflow_dir
    # taken from QAP files
    # resource_pool = {}

    num_of_cores = 1
    # resource_pool({
    #     'epireg': (warp_nipype2.warp_nipype, 'outputspec.epireg')})
    t_node = EPI_DistCorr.create_EPI_DistCorr()  # ###
    t_node.inputs.inputspec.anat_file = \
        '/Users/nanditharajamani/Downloads/ExBox19/T1.nii.gz'
    t_node.inputs.inputspec.func_file = \
        '/Users/nanditharajamani/Downloads/ExBox19/func.nii.gz'
    t_node.inputs.inputspec.fmap_pha = \
        '/Users/nanditharajamani/Downloads/ExBox19/fmap_phase.nii.gz'
    t_node.inputs.inputspec.fmap_mag = \
        '/Users/nanditharajamani/Downloads/ExBox19/fmap_mag.nii.gz'
    t_node.inputs.inputspec.bbr_schedule = \
        '/usr/local/fsl/etc/flirtsch/bbr.sch'
    t_node.inputs.inputspec.deltaTE = 2.46
    t_node.inputs.inputspec.dwellT = 0.0005
    t_node.inputs.inputspec.dwell_asym_ratio = 0.93902439
    t_node.inputs.inputspec.bet_frac = 0.5
    # 'home/nrajamani/FieldMap_SubjectExampleData/SubjectData/epi_run2/fMT0160-0015-00003-000003-01_BRAIN.nii.gz',
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
#   intensity = ImageStats(
#       in_file = t_node.inputs.inputspec.fmap_pha, op_string = '-p 90')
#   if intensity < 3686:
#      raise 'input phase image does not have the correct range values'
    dataSink = pe.Node(nio.DataSink(), name='dataSink_file')
    dataSink.inputs.base_directory = workflow_dir
    # node, out_file = resource_pool["epireg"]
    # warp_workflow.connect(t_node,'outputspec.roi_file',dataSink,'roi_file')
    warp_workflow.connect(t_node, 'outputspec.fieldmap',
                          dataSink, 'fieldmap_file')
    warp_workflow.connect(t_node, 'outputspec.fmapmagbrain',
                          dataSink, 'fmapmagbrain')
    warp_workflow.connect(t_node, 'outputspec.fieldmapmask',
                          dataSink, 'fieldmapmask')
    warp_workflow.connect(t_node, 'outputspec.fmap_despiked',
                          dataSink, 'fmap_despiked')
    warp_workflow.connect(t_node, 'outputspec.struct',
                          dataSink, 'epi2struct')
    warp_workflow.connect(t_node, 'outputspec.anat_func',
                          dataSink, 'anat_func')
    if run is True:
        plugin_args = {'n_procs': num_of_cores}
        warp_workflow.run(plugin=MultiProcPlugin(plugin_args),
                          plugin_args=plugin_args)
        # outpath = glob.glob(
        #     os.path.join(workflow_dir, "EPI_DistCorr","*"))[0]
        # return outpath
    else:
        return warp_workflow, warp_workflow.base_dir


if __name__ == '__main__':
    run_warp_nipype([
        'anat_file',
        'func_file',
        'fmap_pha',
        'fmap_mag',
        'deltaTE',
        'dwellT',
        'dwell_asym_ratio',
        'bet_frac',
        'bbr_schedule'
    ], output_dir=None, run=True)
