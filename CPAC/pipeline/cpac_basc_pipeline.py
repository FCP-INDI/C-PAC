import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
from CPAC.utils import Configuration

import re
import os
import sys
import glob

def prep_basc_workflow(c, subject_infos):
    print('Preparing BASC workflow')
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print('Subjects', s_ids)
    
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

    wf.run(plugin='MultiProc',
                         plugin_args={'n_procs': c.numCoresPerSubject})


def run(config, subject_infos):
    import re
    import subprocess
    subprocess.getoutput('source ~/.bashrc')
    import os
    import sys
    import pickle
    import yaml
    import yamlordereddictloader

    c = Configuration(yaml.safe_load(open(os.path.realpath(config), 'r')))


    prep_basc_workflow(c, pickle.load(open(subject_infos, 'r') ))

