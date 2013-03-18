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
    
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))    
    
    if type(c.modelFile) is list:
        model_sub_list = c.modelFile
    elif os.path.exists(c.modelFile):    
        model_sub_list = [line.rstrip('\r\n').split() \
                          for line in open(c.modelFile, 'r') \
                          if not (line == '\n') and not line.startswith('#')]
    else:
        raise ValueError("modelFile %s has invalid entry" %(c.modelFile)) 

    for model_sub in model_sub_list:
        
        model, subject_list = model_sub
        
        print "running for model %s and resource %s..."%(os.path.basename(model), resource)
        
        if not os.path.exists(model):
            raise Exception("path to model %s doesn't exit"%model)
        
        if not os.path.exists(subject_list):
            raise Exception("path to input subject list %s is invalid"%subject_list)
        
        if c.mixedScanAnalysis == True:
            wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s'%(resource, os.path.basename(model)))
        else:
            wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s/%s'%(resource, os.path.basename(model), scan_ids[0])) 

        wf.base_dir = c.workingDirectory
    
        input_subject_list = [line.rstrip('\r\n') for line in open(subject_list, 'r') \
                              if not (line == '\n') and not line.startswith('#')]
    
        ordered_paths=[]
        for sub in input_subject_list :
           for path in s_paths:
               if sub in path:
                   ordered_paths.append(path)
        
        print "input_subject_list ->", input_subject_list
        #print "ordered_paths ->", ordered_paths
    
        strgy_path = os.path.dirname(s_paths[0]).split(scan_ids[0])[1]
        for ch in ['.']:
            if ch in strgy_path:
                strgy_path = strgy_path.replace(ch, '_')
        
        gp_flow = create_gpa_dataflow("gp_dataflow%s"%strgy_path)
        gp_flow.inputs.inputspec.input_sublist = input_subject_list 
        gp_flow.inputs.inputspec.output_sublist = s_ids
        gp_flow.inputs.inputspec.grp_model = model
        gp_flow.inputs.inputspec.ftest = c.fTest
        
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
        out_dir = os.path.dirname(s_paths[0]).replace(s_ids[0], 'group_analysis_results/_grp_model_%s'%(os.path.basename(model)))
        
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
        if 'tempreg_maps_z_files' in resource:
            out_dir = os.path.join(out_dir, \
                re.search('\w*[#]*\d+', os.path.splitext(os.path.splitext(os.path.basename(s_paths[0]))[0])[0]).group(0))
        
        if c.mixedScanAnalysis == True:
            out_dir = re.sub(r'(\w)*scan_(\w)*(\d)*(\w)*[/]', '', out_dir)
            
            
        ds.inputs.base_directory = out_dir
        ds.inputs.container = ''
        
        ds.inputs.regexp_substitutions = [(r'(?<=rendered)(.)*[/]','/'),
                                          (r'(?<=model_files)(.)*[/]','/'),
                                          (r'(?<=merged)(.)*[/]','/'),
                                          (r'(?<=stats/clusterMap)(.)*[/]','/'),
                                          (r'(?<=stats/unthreshold)(.)*[/]','/'),
                                          (r'(?<=stats/threshold)(.)*[/]','/'),
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
    
        print "**Workflow finished for model %s and resource %s"%(os.path.basename(model), resource)

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

