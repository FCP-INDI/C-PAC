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
import os

import nipype.interfaces.io as nio
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nipype_pipeline_engine.plugins import MultiProcPlugin
from CPAC.utils.configuration import Configuration


def prep_cwas_workflow(c, subject_infos):
    print('Preparing CWAS workflow')
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print('Subjects', s_ids)

    wf = pe.Workflow(name='cwas_workflow')
    wf.base_dir = c.pipeline_setup['working_directory']['path']

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

    plugin_args = {'n_procs': c.numCoresPerSubject}
    wf.run(plugin=MultiProcPlugin(plugin_args),
           plugin_args=plugin_args)


def run(config, subject_infos):
    import subprocess
    subprocess.getoutput('source ~/.bashrc')
    import os
    import pickle
    import yaml

    c = Configuration(yaml.safe_load(open(os.path.realpath(config), 'r')))

    prep_cwas_workflow(c, pickle.load(open(subject_infos, 'r') ))
