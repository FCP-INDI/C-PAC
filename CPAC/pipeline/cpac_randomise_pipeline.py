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
from CPAC.pipeline.cpac_group_runner import load_config_yml
from CPAC.pipeline.nipype_pipeline_engine.plugins import MultiProcPlugin
from CPAC.utils.interfaces.fsl import Merge as fslMerge
from CPAC.utils.monitoring import log_nodes_cb


def load_subject_file(group_config_path):
    group_config_obj = load_config_yml(group_config_path)
    pipeline_output_folder = group_config_obj.pipeline_dir
    
    if not group_config_obj.participant_list == None:
        s_paths = group_config_obj.participant_list
    else:
        s_paths = [x for x in os.listdir(pipeline_output_folder) if os.path.isdir(x)]

def randomise_merged_file(s_paths):
    
    merge = pe.Node(interface=fslMerge(), name='fsl_merge')
    merge.inputs.dimension = 't'
    merge.inputs.merged_file = "randomise_merged.nii.gz"
    merge.inputs.in_files = s_paths   

def randomise_merged_mask(s_paths):

    mask = pe.Node(interface=fsl.maths.MathsCommand(), name='fsl_maths')
    mask.inputs.args = '-abs -Tmin -bin'
    mask.inputs.out_file = "randomise_mask.nii.gz"
    mask.inputs.in_file = s_paths

def prep_randomise_workflow(c, subject_infos):
    print('Preparing Randomise workflow')
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    print('Subjects', s_ids)

    wf = pe.Workflow(name='randomise_workflow')
    wf.base_dir = c.pipeline_setup['working_directory']['path']

    from CPAC.randomise import create_randomise

    rw = create_randomise()

    rw.inputs.inputspec.permutations = c.randopermutations
    rw.inputs.inputspec.subjects    = s_paths
    #rw.inputs.inputspec.pipeline_ouput_folder = c.os.path.join(c.outputDirectory,
     #                                         'pipeline_{0}'.format(c.pipelineName))
    rw.inputs.inputspec.mask_boolean   = c.mask_boolean #TODO pipe from output dir, not the user input
    rw.inputs.inputspec.tfce           = c.tfce # will stay None?
    rw.inputs.inputspec.demean         = c.demean
    rw.inputs.inputspec.c_thresh       = c.c_thresh

    ds = pe.Node(nio.DataSink(), name='randomise_sink')
    out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'randomise_results')
    ds.inputs.base_directory = out_dir
    ds.inputs.container = ''
#'tstat_files' ,'t_corrected_p_files','index_file','threshold_file','localmax_txt_file','localmax_vol_file','max_file','mean_file','pval_file','size_file'
    wf.connect(rw, 'outputspec.tstat_files',
               ds, 'tstat_files')
    wf.connect(rw, 'outputspec.t_corrected_p_files',
               ds, 't_corrected_p_files')
    wf.connect(rw, 'outputspec.index_file', ds, 'index_file')
    wf.connect(rw, 'outputspec.threshold_file', ds, 'threshold_file')
    wf.connect(rw, 'outputspec.localmax_vol_file', ds, 'localmax_vol_file')
    wf.connect(rw, 'outputspec.localmax_txt_file', ds, 'localmax_txt_file')
    wf.connect(rw, 'outputspec.max_file', ds, 'max_file')
    wf.connect(rw, 'outputspec.mean_file', ds, 'mean_file')
    wf.connect(rw, 'outputspec.max_file', ds, 'max_file')
    wf.connect(rw, 'outputspec.pval_file', ds, 'pval_file')
    wf.connect(rw, 'outputspec.size_file', ds, 'size_file')

    plugin_args = {'n_procs': c.numCoresPerSubject,
                   'status_callback': log_nodes_cb}
    wf.run(plugin=MultiProcPlugin(plugin_args),
           plugin_args=plugin_args)

    return wf

  
#def run(config, subject_infos):
#    import re
#    import commands
#    commands.getoutput('source ~/.bashrc')
#    import os
#    import sys
#    import pickle
#    import yaml

#    c = Configuration(yaml.safe_load(open(os.path.realpath(config), 'r')))


#    prep_randomise_workflow(c, pickle.load(open(subject_infos, 'r') ))

