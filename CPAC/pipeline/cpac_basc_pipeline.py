# Copyright (C) 2022-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
import re
import os
import sys
import glob
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nipype_pipeline_engine.plugins import MultiProcPlugin
from CPAC.utils import Configuration


def prep_basc_workflow(c, subject_infos):
    print('Preparing BASC workflow')
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print('Subjects', s_ids)
    
    wf = pe.Workflow(name='basc_workflow')
    wf.base_dir = c.pipeline_setup['working_directory']['path']
    
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

    plugin_args = {'n_procs': c.numCoresPerSubject}
    wf.run(plugin=MultiProcPlugin(plugin_args),
           plugin_args=plugin_args)


def run(config, subject_infos):
    import re
    import subprocess
    subprocess.getoutput('source ~/.bashrc')
    import os
    import sys
    import pickle
    import yaml

    c = Configuration(yaml.safe_load(open(os.path.realpath(config), 'r')))


    prep_basc_workflow(c, pickle.load(open(subject_infos, 'r') ))

