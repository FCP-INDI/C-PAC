import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

import re
import os
import sys

from CPAC.cwas import create_cwas
from CPAC.basc import create_basc

def prep_cwas_workflow(c, subject_paths):
    print 'Preparing CWAS workflow'
    print 'subjects', subject_paths

def prep_basc_workflow(c, subject_paths):
    print 'Preparing BASC workflow'
    print 'subjects', subject_paths
    
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
    
    spr_pattern = re.compile(base_path  + "/(\w+)/(\w+)/(\w+)")
    for subject_path in subject_paths:
        m = spr_pattern.search(subject_path)
        pipeline_id = m.group(1)
        subject_id = m.group(2)
        resource_id = m.group(3)

        key = subject_path.replace(subject_id, '*')
        analysis_map[(resource_id, key)].append((subject_id,subject_path))
        
    import pprint
    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(analysis_map)
    
    for resource, glob_key in analysis_map.keys():
        if resource == 'functional_mni':
            prep_cwas_workflow(c, analysis_map[(resource, glob_key)])
        elif resource == 'alff':
            prep_group_analysis_workflow(c, analysis_map[(resource, glob_key)])

