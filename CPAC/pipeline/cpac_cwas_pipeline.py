import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob


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


    w.run(plugin='MultiProc',
                         plugin_args={'n_procs': c.numCoresPerSubject})



if __name__ == "__main__":
    import argparse
    import re
    import commands
    commands.getoutput('source ~/.bashrc')
    import os
    import sys
    import pickle

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--config',
                        dest='config',
                        required=True,
                        help='location of config file'
                        )

    parser.add_argument('-i', '--subject_infos',
                        dest='subject_infos',
                        required=True,
                        help='subject_info'
                        )




    args = parser.parse_args()
    path, fname = os.path.split(os.path.realpath(args.config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])


    prep_cwas_workflow(c, pickle.load(open(args.subject_infos, 'r') ))


