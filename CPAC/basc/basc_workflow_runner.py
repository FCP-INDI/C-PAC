#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:05:19 2017

@author: aki.nikolaidis
"""
import os
import numpy as np
import nibabel as nb
import sys


import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, bootstrap_list, proc_mem, similarity_metric, cross_cluster=False, roi2_mask_file=None, blocklength=1, affinity_threshold=0.5, out_dir=None, run=True):
    """Run the 'template_workflow' function to execute the modular workflow
    with the provided inputs.
    :type input_resource: str
    :param input_resource: The filepath of the { input resource }. Can have
                           multiple.
    :type out_dir: str
    :param out_dir: (default: None) The output directory to write the results
                    to; if left as None, will write to the current directory.
    :type run: bool
    :param run: (default: True) Will run the workflow; if set to False, will
                connect the Nipype workflow and return the workflow object
                instead.
    :rtype: str
    :return: (if run=True) The filepath of the generated anatomical_reorient
             file.
    :rtype: Nipype workflow object
    :return: (if run=False) The connected Nipype workflow object.
    :rtype: str
    :return: (if run=False) The base directory of the workflow if it were to
             be run.
    """

    import os
    import glob
    import utils

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import pandas as pd
    
    from basc import create_basc

    workflow = pe.Workflow(name='basc_workflow_runner')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output")
    workflow.base_dir = workflow_dir

    resource_pool = {}
    
    basc = create_basc(proc_mem, name='basc')
    basc.inputs.inputspec.subject_file_list=subject_file_list
    basc.inputs.inputspec.roi_mask_file=roi_mask_file
    basc.inputs.inputspec.dataset_bootstraps=dataset_bootstraps
    basc.inputs.inputspec.timeseries_bootstraps=timeseries_bootstraps
    basc.inputs.inputspec.n_clusters=n_clusters
    basc.inputs.inputspec.output_size=output_size
    basc.inputs.inputspec.bootstrap_list=bootstrap_list
    basc.inputs.inputspec.proc_mem=proc_mem
    basc.inputs.inputspec.similarity_metric=similarity_metric
    basc.inputs.inputspec.cross_cluster=cross_cluster
    basc.inputs.inputspec.roi2_mask_file=roi2_mask_file
    basc.inputs.inputspec.blocklength=blocklength
    basc.inputs.inputspec.affinity_threshold=affinity_threshold
    
    
    # shitty Steve pseudo-code that is shitty
#    thing = basc.outputs.outputspec.individual_cluster_voxel_scores_imgs
#    print(thing)
#    print(thing)
#    print(thing)
#    for filepath in basc.outputs.outputspec.individual_cluster_voxel_scores_imgs:
#        filename = os.path.basename(filepath)
#        resource_pool[filename] = filepath

    resource_pool['group_stability_matrix'] = (basc, 'outputspec.group_stability_matrix')
    resource_pool['clusters_G'] = (basc, 'outputspec.clusters_G')
    resource_pool['ism_gsm_corr_file'] = (basc, 'outputspec.ism_gsm_corr_file')
    resource_pool['gsclusters_img'] = (basc, 'outputspec.gsclusters_img')
    resource_pool['cluster_voxel_scores_img'] = (basc, 'outputspec.cluster_voxel_scores_img')
    # below commented out to try out the shitty Steve pseudo-code above, that is shitty
    resource_pool['individual_cluster_voxel_scores_imgs'] = (basc, 'outputspec.individual_cluster_voxel_scores_imgs')
    resource_pool['cluster_voxel_scores'] = (basc, 'outputspec.cluster_voxel_scores')
    resource_pool['k_mask'] = (basc, 'outputspec.k_mask')
    resource_pool['ind_group_cluster_stability'] = (basc, 'outputspec.ind_group_cluster_stability')
    resource_pool['ind_group_cluster_stability_set'] = (basc, 'outputspec.ind_group_cluster_stability_set')
    #resource_pool['individualized_group_clusters_img_file'] = (basc, 'outputspec.individualized_group_clusters_img_file')

    ds = pe.Node(nio.DataSink(), name='datasink_workflow_name')
    ds.inputs.base_directory = workflow_dir
    
    for output in resource_pool.keys():
        node, out_file = resource_pool[output]
        workflow.connect(node, out_file, ds, output)

    plugin_args = { 'n_procs' : int(proc_mem[0]),'memory_gb': int(proc_mem[1])}#, 'raise_insufficient': True, 'maxtasksperchild': 1,'chunksize':1}


    if run == True:
        workflow.run(plugin='MultiProc', plugin_args= plugin_args) #
                     # {'n_procs': 1})
        outpath = glob.glob(os.path.join(workflow_dir, "*", "*"))
        return outpath
    else:
        return workflow, workflow.base_dir
