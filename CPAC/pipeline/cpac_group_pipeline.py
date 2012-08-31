import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob

from CPAC.utils.datasource import create_gpa_dataflow

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

def prep_group_analysis_workflow(c, resource, subject_infos):
    print 'Preparing Group Analysis workflow'
    print 'subjects', subject_infos
    
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    
    if c.mixed_scan_analysis == True:
        wf = pe.Workflow(name = 'group_analysis_%s'%resource)
    else:
        wf = pe.Workflow(name = 'group_analysis_%s_%s'%resource%scan_ids[0]) 
    
    wf.base_dir = c.workingDirectory
    
    #extract model files
    model_list = [line.rstrip('\r\n') for line in open(c.modelFile, 'r')]
    
    if not model_list:
        raise Exception("mode_list is empty. Please provide" \
                         "a model file with full paths of the" \
                         "folder containing models for group analysis")
    
    from collections import defaultdict
    model_map = defaultdict(list)
    
    #create a map of model as key and its sub files as values
    for model in model_list:
        if os.path.exists(model):
            files = glob.glob(os.path.join(model, '*'))
            model_map[os.path.basename(model)] = files
        else:
            raise Exception ("Path to the model %s doesn't exist"%model)
    
    print model_map
    gp_flow = create_gpa_dataflow(model_map, c.ftest)
    
    from CPAC.group_analysis import create_group_analysis
    
    gpa_wf = create_group_analysis(c.ftest)
    gpa_wf.inputs.inputspec.zmap_files = s_paths
    
    wf.connect(gp_flow, 'outputspec.mat',
               gpa_wf, 'inputspec.mat_file')
    wf.connect(gp_flow, 'outputspec.con',
               gpa_wf, 'inputspec.con_file')
    wf.connect(gp_flow, 'outputspec.grp',
                gpa_wf, 'inputspec.grp_file')
        
    if c.fTest:
        workflow.connect(gp_flow, 'outputspec.fts',
                         gpa_wf, 'inputspec.fts_file') 
    
    ds = pe.Node(nio.DataSink(), name='gpa_sink')
    out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'group_analyis_results')
    ds.inputs.base_directory = out_dir
    ds.inputs.container = resource
    
    return wf

    
if __name__ == "__main__":
    import argparse
    import re
    
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
    analysis_map_gp = defaultdict(list)
    
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
        
        # separate map for group analysis
        if c.mixed_scan_analysis == True:
            key = key.replace(scan_id, '*')
            
        analysis_map_gp[(resource_id, key)].append((pipeline_id, subject_id, scan_id, subject_path))
            
    w_list = []
    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':
            w_list.append(prep_basc_workflow(c, analysis_map[(resource, glob_key)]))
#            w_list.append(prep_cwas_workflow(c, analysis_map[(resource, glob_key)]))
    
    for resource, glob_key in analysis_map_gp.keys():
        if resource in c.derivative_list:
            w_list.append(prep_group_analysis_workflow(c, resource, analysis_map_gp[(resource, glob_key)]))
        #if re.match(r"^falff", resource):
            
    
    for w in w_list:
        w.run()

