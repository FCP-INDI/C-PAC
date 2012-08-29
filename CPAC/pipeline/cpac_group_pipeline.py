import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys

def split_folders(path):
    folders = []
    
    while 1:
        path, folder = os.path.split(path)
        if folder != "":
            folders.append(folder)
        else:
            if path != "":
                folders.append(path)
            break
        
    folders.reverse()
    
    return folders
    
def prep_cwas_workflow(c, subject_infos):
    print 'Preparing CWAS workflow'
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print 'Subjects', s_ids

    wf = pe.Workflow(name='cwas_workflow')
    wf.base_dir = c.workingDirectory

    from CPAC.cwas import create_cwas
    import numpy as np
    regressor = np.loadtxt(c.cwasRegressorFile)

    cw = create_cwas()
    cw.inputs.inputspec.roi = c.cwasROIFile
    cw.inputs.inputspec.subjects = s_paths
    cw.inputs.inputspec.regressor = regressor
    cw.inputs.inputspec.f_samples = c.cwasFSamples
    cw.inputs.inputspec.parallel_nodes = c.cwasParallelNodes
    
    ds = pe.Node(nio.DataSink(), name='cwas_sink')
    out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'basc_results')
    ds.inputs.base_directory = out_dir
    ds.inputs.container = ''

    wf.connect(cw, 'outputspec.F_map',
               ds, 'F_map')
    wf.connect(cw, 'outputspec.p_map',
               ds, 'p_map')
    
    return wf

def prep_basc_workflow(c, subject_infos):
    print 'Preparing BASC workflow'
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print 'Subjects', s_ids
    
    wf = pe.Workflow(name='basc_workflow')
    wf.base_dir = c.workingDirectory
    
    from CPAC.basc import create_basc
    
    b = create_basc()
    b.inputs.inputspec.roi = c.bascROIFile
    b.inputs.inputspec.subjects = s_paths
    b.inputs.inputspec.k_clusters = c.bascClusters
    b.inputs.inputspec.dataset_bootstraps = c.bascDatasetBootstraps
    b.inputs.inputspec.timeseries_bootstraps = c.bascTimeseriesBootstraps
    
    aff_list = open(c.bascAffinityThresholdFile, 'r').readlines()
    aff_list = [ float(aff.rstrip('\r\n')) for aff in aff_list]
    
    b.inputs.inputspec.affinity_threshold = aff_list
    
    ds = pe.Node(nio.DataSink(), name='basc_sink')
    out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'basc_results')
    ds.inputs.base_directory = out_dir
    ds.inputs.container = ''
    
#    wf.connect(b, 'outputspec.gsm',
#               ds, 'gsm')
#    wf.connect(b, 'outputspec.gsclusters',
#               ds, 'gsclusters')
#    wf.connect(b, 'outputspec.gsmap',
#               ds, 'gsmap')
    wf.connect(b, 'outputspec.gsclusters_img',
               ds, 'gsclusters_img')
    wf.connect(b, 'outputspec.ismap_imgs',
               ds, 'ismap_imgs')
    
    return wf

def prep_group_analysis_workflow(c, subject_infos):
    print 'Preparing Group Analysis workflow'
    print 'subjects', subject_infos
    
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-c', '--config',
                        dest='config',
                        required=True,
                        help='location of config file'
                        )
    
    parser.add_argument('-s', '--subjects',
                        dest='subjects',
                        required=True,
                        help='location of subjects file'
                        )
    
    args = parser.parse_args()
    path, fname = os.path.split(os.path.realpath(args.config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])
    
    subject_paths = open(args.subjects, 'r').readlines()
    subject_paths = [s.rstrip('\r\n') for s in subject_paths]
    
    #base_path = os.path.dirname(os.path.commonprefix(subject_paths))
    base_path = c.sinkDirectory
    
    from collections import defaultdict
    analysis_map = defaultdict(list)
    
    for subject_path in subject_paths:
        #Remove the base bath offset
        rs_path = subject_path.replace(base_path, "", 1)
        folders = split_folders(rs_path)
        
        pipeline_id = folders[1]
        subject_id = folders[2]
        resource_id = folders[3]
        scan_id = folders[4]
        
        key = subject_path.replace(subject_id, '*')
        analysis_map[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))
    
    w_list = []
    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':
            w_list.append(prep_basc_workflow(c, analysis_map[(resource, glob_key)]))
#            w_list.append(prep_cwas_workflow(c, analysis_map[(resource, glob_key)]))
        elif resource == 'alff':
            w_list.append(prep_group_analysis_workflow(c, analysis_map[(resource, glob_key)]))
    
    for w in w_list:
        w.run()

