import os

import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe

from CPAC.utils.configuration import Configuration


def prep_cwas_workflow(c, subject_infos):
    print('Preparing CWAS workflow')
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print('Subjects', s_ids)

    wf = pe.Workflow(name='cwas_workflow')
    wf.base_dir = c.workingDirectory

    from CPAC.cwas import create_cwas
    import numpy as np
    regressor = np.loadtxt(c.cwasRegressorFile)

    cw = create_cwas()
    cw.inputs.inputspec.roi         = c.cwasROIFile
    cw.inputs.inputspec.subjects    = s_paths
    cw.inputs.inputspec.regressor   = regressor
    cw.inputs.inputspec.cols        = c.cwasRegressorCols
    cw.inputs.inputspec.f_samples   = c.cwasFSamples
    cw.inputs.inputspec.strata      = c.cwasRegressorStrata # will stay None?
    cw.inputs.inputspec.parallel_nodes = c.cwasParallelNodes

    ds = pe.Node(nio.DataSink(), name='cwas_sink')
    out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'cwas_results')
    ds.inputs.base_directory = out_dir
    ds.inputs.container = ''

    wf.connect(cw, 'outputspec.F_map',
               ds, 'F_map')
    wf.connect(cw, 'outputspec.p_map',
               ds, 'p_map')

    wf.run(plugin='MultiProc',
                         plugin_args={'n_procs': c.numCoresPerSubject})


def run(config, subject_infos):
    import subprocess
    subprocess.getoutput('source ~/.bashrc')
    import os
    import pickle
    import yaml
    import yamlordereddictloader

    c = Configuration(yaml.safe_load(open(os.path.realpath(config), 'r')))

    prep_cwas_workflow(c, pickle.load(open(subject_infos, 'r') ))
