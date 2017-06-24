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

#import CPAC
#from CPAC.basc.utils import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average, individual_stability_matrix
#from CPAC.basc import group_stability_matrix, individual_group_clustered_maps, individual_stability_matrix, nifti_individual_stability, ndarray_to_vol, create_basc
#from CPAC.utils import safe_shape

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, n_clusters, output_size, cross_cluster=False, roi2_mask_file=None, affinity_threshold=0.5, out_dir=None, run=True):
   #run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, k_clusters, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)


    #subject_file, roi_mask_file, n_bootstraps, k_clusters, cross_cluster=False, roi2_mask_file=None, cbb_block_size=None, affinity_threshold=0.5
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

    #output = "{ output resource name }"

    workflow = pe.Workflow(name='basc_workflow_runner')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output")
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    basc=create_basc(name='basc')
    basc.inputs.inputspec.subject_file_list=subject_file_list
    basc.inputs.inputspec.roi_mask_file=roi_mask_file
    basc.inputs.inputspec.dataset_bootstraps=dataset_bootstraps
    basc.inputs.inputspec.timeseries_bootstraps=timeseries_bootstraps
    basc.inputs.inputspec.n_clusters=n_clusters
    basc.inputs.inputspec.output_size=output_size
    basc.inputs.inputspec.cross_cluster=cross_cluster
    basc.inputs.inputspec.roi2_mask_file=roi2_mask_file
    basc.inputs.inputspec.affinity_threshold=affinity_threshold
    

    resource_pool['gsm'] = (basc, 'outputspec.gsm')
    resource_pool['gsclusters'] = (basc, 'outputspec.gsclusters')
    resource_pool['gsmap'] = (basc, 'outputspec.gsmap')
    resource_pool['gsclusters_img'] = (basc, 'outputspec.gsclusters_img')
    resource_pool['gsmap_img'] = (basc, 'outputspec.gsmap_img')
    resource_pool['ismap_imgs'] = (basc, 'outputspec.ismap_imgs')


    ds = pe.Node(nio.DataSink(), name='datasink_workflow_name')
    ds.inputs.base_directory = workflow_dir

    for output in resource_pool.keys():
        node, out_file = resource_pool[output]
        workflow.connect(node, out_file, ds, output)



    if run == True:
        workflow.run(plugin='Linear', plugin_args= \
                         {'n_procs': 1})
        outpath = glob.glob(os.path.join(workflow_dir, "*", "*"))
        return outpath
    else:
        return workflow, workflow.base_dir