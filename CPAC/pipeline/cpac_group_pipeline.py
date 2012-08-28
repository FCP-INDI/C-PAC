import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys




def prep_cwas_workflow(c, subject_paths):
    print 'Preparing CWAS workflow'
    print 'subjects', subject_paths

    from CPAC.cwas import create_cwas

def prep_basc_workflow(c, subject_paths):
    print 'Preparing BASC workflow'
    s_ids, scan_ids, s_paths = *zip(subject_paths)
    print 'Subjects', s_ids
    
    wf = pe.Workflow(name='basc_workflow')
    
    from CPAC.basc import create_basc
    b = basc.create_basc()
    b.inputs.inputspec.roi = c.bascROIFile
    b.inputs.inputspec.subjects = s_paths
    b.inputs.inputspec.k_clusters = c.bascClusters
    b.inputs.inputspec.dataset_bootstraps = c.bascDatasetBootstraps
    b.inputs.inputspec.timeseries_bootstraps = c.bascTimeseriesBootstraps
    b.inputs.inputspec.affinity_threshold = c.bascAffinityTheshold
    
    ds = pe.Node(nio.DataSink(), name='basc_sink')
    ds.inputs.base_directory = c.sinkDirectory
    ds.inputs.container = 'basc_results'
    
    wf.connect(b, 'outputspec.gsm',
               ds, 'gsm')
    wf.connect(b, 'outputspec.gsclusters',
               ds, 'gsclusters')
    wf.connect(b, 'outputspec.gsmap',
               ds, 'gsmap')
    wf.connect(b, 'outputspec.gsclusters_img',
               ds, 'gsclusters_img')
    wf.connect(b, 'outputspec.ismap_imgs',
               ds, 'ismap_imgs')
    
    return wf

def prep_group_analysis_workflow(c, subject_paths):
    print 'Preparing Group Analysis workflow'
    print 'subjects', subject_paths
    
    
    
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
    
    spr_pattern = re.compile(base_path  + "/(\w+)/(\w+)/(\w+)/(\w+)")
    for subject_path in subject_paths:
        m = spr_pattern.search(subject_path)
        pipeline_id = m.group(1)
        subject_id = m.group(2)
        resource_id = m.group(3)
        scan_id = m.group(4)
        
        key = subject_path.replace(subject_id, '*')
        analysis_map[(resource_id, key)].append((subject_id, scan_id, subject_path))
        
    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(analysis_map)
    
    w_list = []
    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':
            w_list.append(prep_cwas_workflow(c, analysis_map[(resource, glob_key)]))
        elif resource == 'alff':
            w_list.append(prep_group_analysis_workflow(c, analysis_map[(resource, glob_key)]))
    
    return w_list

