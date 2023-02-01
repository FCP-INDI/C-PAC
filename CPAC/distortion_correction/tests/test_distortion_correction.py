#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pytest

import nipype.interfaces.io as nio
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nipype_pipeline_engine.plugins import MultiProcPlugin
from CPAC.utils.test_resources import setup_test_wf
# from CPAC.func_preproc.func_preproc import create_func_preproc
# from CPAC.distortion_correction.distortion_correction import blip_distcor_wf


@pytest.mark.skip(reason='needs refactoring')
def run_warp_nipype(inputs, output_dir=None, run=True):
    """
    Created on Thu Nov  9 10:36:37 2017

    @author: nrajamani
    """
    # rom nose.tools import *
    #  This script runs the warp_nipype workflow to execute the
    #  interfaces and with the inputs already provided
    #  It accepts as it's arguments the input whihc you;d give the
    #  warp_nipe file, an output_dir which can be either be mentioned
    #  or if it is set to none will write it in the current working
    #  directory.
    #  The argument run can either be set to true(default) or to
    #  false. If set to false, it should connect to the nipype workflow
    #  and return the workflow object instead
    #  What all should it return?: if the run had been set to true, it
    #  would generate the filepath of the output from the workflow. If
    #  the run was set to false, then it will return the file path of
    #  the base directory, the workflow object, etc.
    import EPI_DistCorr
    warp_workflow = pe.Workflow(name='preproc')
    if output_dir is None:
        output_dir = '/Users/nanditharajamani'

    workflow_dir = os.path.join(
        output_dir, "workflow_output_with_fieldmap_with_change")
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
    #   intensity = ImageStats(in_file = t_node.inputs.inputspec.fmap_pha,
    #       op_string = '-p 90')
    #   if intensity < 3686:
    #      raise 'input phase image does not have the correct range values'
    dataSink = pe.Node(nio.DataSink(), name='dataSink_file')
    dataSink.inputs.base_directory = workflow_dir
    # node, out_file = resource_pool["epireg"]
    # warp_workflow.connect(t_node,'outputspec.roi_file',dataSink,'roi_file')
    warp_workflow.connect(
        t_node, 'outputspec.fieldmap', dataSink, 'fieldmap_file')
    warp_workflow.connect(
        t_node, 'outputspec.fmapmagbrain', dataSink, 'fmapmagbrain')
    warp_workflow.connect(
        t_node, 'outputspec.fieldmapmask', dataSink, 'fieldmapmask')
    warp_workflow.connect(
        t_node, 'outputspec.fmap_despiked', dataSink, 'fmap_despiked')
    warp_workflow.connect(
        t_node, 'outputspec.struct', dataSink, 'epi2struct')
    warp_workflow.connect(
        t_node, 'outputspec.anat_func', dataSink, 'anat_func')
    if run is True:
        plugin_args = {'n_procs': num_of_cores}
        warp_workflow.run(plugin=MultiProcPlugin(plugin_args),
                          plugin_args=plugin_args)
        # outpath = glob.glob(
        #     os.path.join(workflow_dir, "EPI_DistCorr","*"))[0]
        # return outpath
    else:
        return warp_workflow, warp_workflow.base_dir


@pytest.mark.skip(reason='needs refactoring')
def test_phasediff_distcor():
    run_warp_nipype([
        'anat_file', 'func_file', 'fmap_pha', 'fmap_mag', 'deltaTE', 'dwellT',
        'dwell_asym_ratio', 'bet_frac', 'bbr_schedule'
    ], output_dir=None, run=True)



@pytest.mark.skip(reason='needs refactoring')
def test_blip_distcor():

    # good data to use
    s3_prefix = \
        "s3://fcp-indi/data/Projects/HBN/MRI/Site-CBIC/sub-NDARAB708LM5"
    s3_paths = [
        "func/sub-NDARAB708LM5_task-rest_run-1_bold.nii.gz",
        "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz",
        "fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.nii.gz"
    ]

    wf, ds, local_paths = setup_test_wf(s3_prefix, s3_paths,
                                        "test_blip_distcor",
                                        workdirs_to_keep=[
                                            "func_preproc_automask"])

    func_preproc = create_func_preproc(skullstrip_tool='fsl',
                                       wf_name='func_preproc_automask')
    func_preproc.inputs.inputspec.func = local_paths[
        "func/sub-NDARAB708LM5_task-rest_run-1_bold.nii.gz"]

    blip = blip_distcor_wf("blip_distcor")
    blip.inputs.inputspec.opposite_pe_epi = local_paths[
        "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz"]
    blip.inputs.inputspec.same_pe_epi = local_paths[
        "fmap/sub-NDARAB708LM5_dir-AP_acq-fMRI_epi.nii.gz"]

    wf.connect(func_preproc, 'outputspec.func_mean',
               blip, 'inputspec.func_mean')

    wf.connect(blip, 'outputspec.blip_warp', ds, 'source_warp')
    wf.connect(blip, 'outputspec.new_func_mean', ds, 'new_func_mean')
    wf.connect(blip, 'outputspec.new_func_mask', ds, 'new_func_mask')

    wf.run()


@pytest.mark.skip(reason='needs refactoring')
def test_blip_distcor_only_opposite():

    # good data to use
    s3_prefix = \
        "s3://fcp-indi/data/Projects/HBN/MRI/Site-CBIC/sub-NDARAB708LM5"
    s3_paths = [
        "func/sub-NDARAB708LM5_task-rest_run-1_bold.nii.gz",
        "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz",
    ]

    wf, ds, local_paths = setup_test_wf(s3_prefix, s3_paths,
                                        "test_blip_distcor_only_opposite",
                                        workdirs_to_keep=[
                                            "func_preproc_automask"])

    func_preproc = create_func_preproc(skullstrip_tool='fsl',
                                       wf_name='func_preproc_automask')
    func_preproc.inputs.inputspec.func = local_paths[
        "func/sub-NDARAB708LM5_task-rest_run-1_bold.nii.gz"]

    blip = blip_distcor_wf("blip_distcor_only_opposite")
    blip.inputs.inputspec.opposite_pe_epi = local_paths[
        "fmap/sub-NDARAB708LM5_dir-PA_acq-fMRI_epi.nii.gz"]
    blip.inputs.inputspec.same_pe_epi = None

    wf.connect(func_preproc, 'outputspec.func_mean',
               blip, 'inputspec.func_mean')

    wf.connect(blip, 'outputspec.blip_warp', ds, 'source_warp')
    wf.connect(blip, 'outputspec.new_func_mean', ds, 'new_func_mean')
    wf.connect(blip, 'outputspec.new_func_mask', ds, 'new_func_mask')

    wf.run()
