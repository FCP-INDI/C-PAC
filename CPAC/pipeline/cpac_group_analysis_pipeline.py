import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob

from CPAC.utils.datasource import create_gpa_dataflow
from CPAC.utils import Configuration
from CPAC.utils.utils import prepare_gp_links

def prep_group_analysis_workflow(c, resource, subject_infos):
    print 'Preparing Group Analysis workflow for resource', resource
    print 'subjects', subject_infos
    
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    
    if c.mixedScanAnalysis == True:
        wf = pe.Workflow(name = 'group_analysis/%s'%resource)
    else:
        wf = pe.Workflow(name = 'group_analysis/%s/%s'%(resource,scan_ids[0])) 
    
    wf.base_dir = c.workingDirectory
    
    #extract model files
    model_list = [line.rstrip('\r\n') for line in open(c.modelFile, 'r') if not (line == '\n') and not line.startswith('#')]
    
    if not model_list:
        raise Exception("mode_list is empty. Please provide" \
                         "a model file with full paths of the" \
                         "folder containing models for group analysis")
    
    from collections import defaultdict
    model_map = defaultdict(list)
    
    #create a map of model as key and its sub files as values
    import os 
    import glob
    for model in model_list:
        if os.path.exists(model):
            files = glob.glob(os.path.join(model, '*'))
            model_map[os.path.basename(model)] = files
        else:
            raise Exception ("Path to the model %s doesn't exist"%model)
    
    #print model_map
    
    input_subject_list = [line.rstrip('\r\n') for line in open(c.groupAnalysisSubjectList, 'r') if not (line == '\n') and not line.startswith('#')]
    
    ordered_paths=[]
    for sub in input_subject_list :
       for path in s_paths:
           if sub in path:
               ordered_paths.append(path)
    print "input_subject_list", input_subject_list
    print "ordered_paths", ordered_paths
    
    strgy_path = os.path.dirname(s_paths[0]).split(scan_ids[0])[1]
    for ch in ['.']:
        if ch in strgy_path:
            strgy_path = strgy_path.replace(ch, '_')
    
    gp_flow = create_gpa_dataflow(model_map, c.fTest, "gp_dataflow%s"%strgy_path)
    gp_flow.inputs.inputspec.input_sublist = input_subject_list 
    gp_flow.inputs.inputspec.output_sublist = s_ids
    
    from CPAC.group_analysis import create_group_analysis
    
    gpa_wf = create_group_analysis(c.fTest, "gp_analysis%s"%strgy_path)
    gpa_wf.inputs.inputspec.zmap_files = ordered_paths
    gpa_wf.inputs.inputspec.z_threshold = c.zThreshold
    gpa_wf.inputs.inputspec.p_threshold = c.pThreshold
    gpa_wf.inputs.inputspec.parameters = (c.FSLDIR,
                                               'MNI152')
    
    wf.connect(gp_flow, 'outputspec.mat',
               gpa_wf, 'inputspec.mat_file')
    wf.connect(gp_flow, 'outputspec.con',
               gpa_wf, 'inputspec.con_file')
    wf.connect(gp_flow, 'outputspec.grp',
                gpa_wf, 'inputspec.grp_file')
        
    if c.fTest:
        wf.connect(gp_flow, 'outputspec.fts',
                   gpa_wf, 'inputspec.fts_file') 
    
    ds = pe.Node(nio.DataSink(), name='gpa_sink')
    #out_dir = os.path.join('group_analysis_results', resource)
    out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'group_analysis_results')
    if 'sca_roi' in resource:
        out_dir = os.path.join(out_dir, \
          re.search('ROI_number_(\d)+',os.path.splitext(os.path.splitext(os.path.basename(s_paths[0]))[0])[0]).group(0))
    if 'centrality' in resource:
         names = ['degree_centrality_binarize', 'degree_centrality_weighted', \
                  'eigenvector_centrality_binarize', 'eigenvector_centrality_weighted']
         for name in names:
             if name in os.path.basename(s_paths[0]):
                 out_dir = os.path.join(out_dir, name)
                 break
    if c.mixedScanAnalysis == True:
        out_dir = re.sub(r'(\w)*scan_(\w)*(\d)*(\w)*[/]', '', out_dir)
        
        
    ds.inputs.base_directory = out_dir
    ds.inputs.container = ''
    
    ds.inputs.regexp_substitutions = [(r'(?<=rendered)(.)*_grp_model_','/_grp_model_'),
                                      (r'(?<=model_files)(.)*_grp_model_','/_grp_model_'),
                                      (r'(?<=merged)(.)*[/]','/'),
                                      (r'(?<=stats/clusterMap)(.)*_grp_model_','/_grp_model_'),
                                      (r'(?<=stats/unthreshold)(.)*_grp_model_','/_grp_model_'),
                                      (r'(?<=stats/threshold)(.)*_grp_model_','/_grp_model_'),
                                      (r'_cluster(.)*[/]',''),
                                      (r'_slicer(.)*[/]',''),
                                      (r'_overlay(.)*[/]','')]

    if 1 in c.runSymbolicLinks:


        link_node = pe.MapNode(interface=util.Function(
                            input_names=['in_file',
                                        'resource'],
                                output_names=[],
                                function=prepare_gp_links),
                                name='link_gp_', iterfield=['in_file'])
        link_node.inputs.resource = resource
        wf.connect(ds, 'out_file', link_node, 'in_file')


    ########datasink connections#########
    
    wf.connect(gp_flow, 'outputspec.mat',
               ds, 'model_files')
    wf.connect(gp_flow, 'outputspec.grp',
               ds, 'model_files.@02')
    wf.connect(gp_flow, 'outputspec.sublist',
               ds, 'model_files.@03')
    wf.connect(gpa_wf, 'outputspec.merged',
               ds, 'merged')
    wf.connect(gpa_wf, 'outputspec.zstats',
               ds, 'stats.unthreshold')
    wf.connect(gpa_wf, 'outputspec.zfstats',
               ds,'stats.unthreshold.@01')
    wf.connect(gpa_wf, 'outputspec.fstats',
               ds,'stats.unthreshold.@02')
    wf.connect(gpa_wf, 'outputspec.cluster_threshold_zf',
               ds, 'stats.threshold')
    wf.connect(gpa_wf, 'outputspec.cluster_index_zf',
               ds,'stats.clusterMap')
    wf.connect(gpa_wf, 'outputspec.cluster_localmax_txt_zf',
               ds, 'stats.clusterMap.@01')
    wf.connect(gpa_wf, 'outputspec.overlay_threshold_zf',
               ds, 'rendered')
    wf.connect(gpa_wf, 'outputspec.rendered_image_zf',
               ds, 'rendered.@01')   
    wf.connect(gpa_wf, 'outputspec.cluster_threshold',
               ds,  'stats.threshold.@01')
    wf.connect(gpa_wf, 'outputspec.cluster_index',
               ds, 'stats.clusterMap.@02')
    wf.connect(gpa_wf, 'outputspec.cluster_localmax_txt',
               ds, 'stats.clusterMap.@03')
    wf.connect(gpa_wf, 'outputspec.overlay_threshold',
               ds, 'rendered.@02')
    wf.connect(gpa_wf, 'outputspec.rendered_image',
               ds, 'rendered.@03')
    
    ######################################
    
    wf.run(plugin='MultiProc',
                         plugin_args={'n_procs': c.numCoresPerSubject})

    



def run(config, subject_infos, resource):
    import re
    import commands
    commands.getoutput('source ~/.bashrc')
    import os
    import sys
    import pickle
    import yaml
    
    c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))

    prep_group_analysis_workflow(c, pickle.load(open(resource, 'r') ), pickle.load(open(subject_infos, 'r')))

