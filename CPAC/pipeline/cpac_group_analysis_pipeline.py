import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

import re
import os
import sys
import glob

from CPAC.utils.datasource import create_grp_analysis_dataflow
from CPAC.utils import Configuration
from CPAC.utils.utils import prepare_gp_links
from CPAC.group_analysis import create_group_analysis


def prep_group_analysis_workflow(c, resource, subject_infos):
    
    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))
    #print "p_id -%s, s_ids -%s, scan_ids -%s, s_paths -%s" %(p_id, s_ids, scan_ids, s_paths) 

    '''
    #diag = open(os.path.join('/home/data/Projects/CPAC_Regression_Test/2013-08-19-20_v0-3-1/fsl-model/2013-09-03', 'group_analysis_diagnostic.txt'), 'wt')

    #for tup in subject_infos:
    #    print >>diag, list(tup)

    #print >>diag, ""

    #for tup in zip(*subject_infos):
    #    print >>diag, list(tup)

    #print >>diag, ""


    print >>diag, "Working variables passed from cpac_group_runner: "
    print >>diag, ""
    print >>diag, "Pipeline ID (p_id): ", p_id
    print >>diag, "Subject IDs (s_ids): ", s_ids
    print >>diag, "Scan IDs (scan_ids): ", scan_ids
    print >>diag, "(s_paths): ", s_paths
    print >>diag, ""
    '''

    def get_phenotypic_file(phenotypic_file, m_dict, m_list, mod_path, sub_id):
        
        #print "phenotypic_file, m_dict", phenotypic_file, m_dict
        import csv
        reader = csv.reader(open(phenotypic_file, 'rU'))
        columns = {}
        order = {}
        count = 0
        headers = reader.next()
                
        for h in headers:
            columns[h] =[]
            order[h] = count
            count+=1
            
        for r in reader:
            for h, v in zip(headers, r):
                if v:
                    columns[h].append(str(v))

        if m_dict:
            for measure in m_list:
                if measure in headers:
                    #check if 'MeanFD  is present'
                    if len(columns[measure]) < 1:
                        for sub in columns[sub_id]:
                            if m_dict.get(sub):
                                if m_dict.get(sub).get(measure):
                                    columns[measure].append(m_dict[sub][measure])
                                else:
                                    raise Exception("Couldn't find %s value for subject %s"%(measure,sub))
                            else:
                                raise Exception("Couldn't find subject %s in the parameter file"%sub)
        
        b = zip(*([k] + columns[k] for k in sorted(columns, key=order.get)))
        
        
        try:
            os.makedirs(mod_path)
        except:
            print "%s already exist"%(mod_path)
            
        new_phenotypic_file = os.path.join(mod_path, os.path.basename(phenotypic_file))
                
        a = csv.writer(open(new_phenotypic_file, 'w'))
        
        for col in b:
            a.writerow(list(col))
          
        return new_phenotypic_file

    threshold_val = None
    measure_dict = None
    measure_list = ['MeanFD', 'MeanFD_Jenkinson', 'MeanDVARS']
    model_sub_list = []
    

    if c.runScrubbing == 1:

        #get scrubbing threshold
    
        if re.search('(?<=/_threshold_)\d+.\d+',s_paths[0]):

            threshold_val = re.search('(?<=/_threshold_)\d+.\d+',s_paths[0]).group(0)

        elif len(c.scrubbingThreshold) == 1:

            threshold_val = c.scrubbingThreshold[0]

        else:
            print ("Found Multiple threshold value ")
    

        print "scrubbing threshold_val -->", threshold_val

    else:

        print "No scrubbing enabled."
        print "\n"


    #pick the right parameter file from the pipeline folder
    #create a dictionary of subject and measures in measure_list
    if c.runScrubbing == 1:
  
        try:
            parameter_file = os.path.join(c.outputDirectory, p_id[0], '%s_threshold_%s_all_params.csv'%(scan_ids[0].strip('_'),threshold_val))

            if os.path.exists(parameter_file):
                import csv
                measure_dict = {}
                f = csv.DictReader(open(parameter_file,'r'))

                for line in f:
                    measure_map = {}
                    for m in measure_list:
                        if line.get(m):
                            measure_map[m] = line[m]

                    measure_dict[line['Subject']] = measure_map
            else:
                print "No file name %s found"%parameter_file
                
        except Exception:
            print "Exception while extracting parameters from movement file - %s"%(parameter_file)


    #print >>diag, "Begins to iterate over each config file listed here: ", c.modelConfigs
    #print >>diag, ""
    
    for config in c.modelConfigs:
        
        import yaml
        
        try:
            conf = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
        except:
            raise Exception("Error in reading %s configuration file" % config)

        #print >>diag, "Starting iteration for config: ", config
        #print >>diag, ""


        group_sublist = open(conf.subjectListFile, 'r')

        sublist_items = group_sublist.readlines()

        subject_list = [line.rstrip('\n') for line in sublist_items \
                              if not (line == '\n') and not line.startswith('#')]


        #print >>diag, "Subject list run-through #1: ", subject_list
        #print >>diag, ""

    
        #subject_list = [line.rstrip('\r\n') for line in open(conf.subjectListFile, 'r') \
        #                      if not (line == '\n') and not line.startswith('#')]

        # list of subject paths which DO exist
        exist_paths = []
        
        # check for missing subject for the derivative


        #print >>diag, "> Iterates over subject_list - for each subject in the subject list, it iterates over the paths in s_paths."
        #print >>diag, "> For each path, it checks if the current subject exists in this path, and then appends this subject to 'exist_paths' list."
        #print >>diag, ""

        for sub in subject_list :

            for path in s_paths:

                if sub in path:
                    exist_paths.append(sub)


        #print >>diag, "Current status of exist_paths list: "
        #print >>diag, exist_paths
        #print >>diag, ""



        # check to see if any derivatives of subjects are missing
        if len(list(set(subject_list) - set(exist_paths))) >0:
            print "-------------------------------------------"
            print "List of outputs missing for subjects:"
            print list(set(subject_list) - set(exist_paths))
            print "\n"
            print "..for derivatives:"
            print resource
            print "\n"
            print "..at paths:"
            print os.path.dirname(s_paths[0]).replace(s_ids[0], '*')
            print "-------------------------------------------"

            print '\n'

            #import warnings
            #warnings.warn(msg)
        

        mod_path = os.path.join(os.path.dirname(s_paths[0]).replace(s_ids[0], 'group_analysis_results/_grp_model_%s'%(conf.modelName)),
                                'model_files')


        #print >>diag, "> Created mod_path variable: ", mod_path
        #print >>diag, ""

                
        print "basename: ", os.path.basename(conf.subjectListFile)

        try:

            os.makedirs(mod_path)
            print "Creating directory:"
            print mod_path
            print "\n"

        except:

            print "Attempted to create directory, but path already exists:"
            print mod_path
            print '\n'
        

        new_sub_file = os.path.join(mod_path, os.path.basename(conf.subjectListFile))

        try:

            f = open(new_sub_file, 'w')
         
            for sub in exist_paths:
                print >>f, sub
        
            f.close()

            #print >>diag, "> Created new subject list file: ", new_sub_file
            #print >>diag, ""

            #print >>diag, "> ..which is filled with the subjects from exist_paths"
            #print >>diag, ""

        except:

            print "Error: Could not open subject list file: ", new_sub_file
            print ""
            raise Exception


        #print >>diag, "> Updates the FSL model config's subject list file parameter from: ", conf.subjectListFile

        conf.update('subjectListFile',new_sub_file)

        #print >>diag, "> ..to new subject list file: ", conf.subjectListFile
        #print >>diag, ""
        
        sub_id = conf.subjectColumn
        
        '''
        print >>diag, "> If measure_dict is not empty, it updates the phenotypic file with these parameters: "
        print >>diag, ""

        print >>diag, "measure_dict: ", measure_dict
        print >>diag, "measure_list: ", measure_list
        print >>diag, "mod_path: ", mod_path
        print >>diag, "sub_id: ", sub_id
        print >>diag, ""
        '''

        if measure_dict != None:
            conf.update('phenotypicFile',get_phenotypic_file(conf.phenotypicFile, measure_dict, measure_list, mod_path, sub_id))
            
            
        print "Model config dictionary ->"
        print conf.__dict__
        print '\n'



        # Run 'create_fsl_model' script to extract phenotypic data from
        # the phenotypic file for each of the subjects in the subject list

        try:

            from CPAC.utils import create_fsl_model
            create_fsl_model.run(conf, True)

            #print >>diag, "> Runs create_fsl_model."
            #print >>diag, ""

        except Exception, e:

            print "FSL Group Analysis model not successfully created - error in create_fsl_model script"
            #print "Error ->", e
            raise


            
        model_sub_list.append((conf.outputModelFilesDirectory, conf.subjectListFile))

        print "model_sub_list ->", model_sub_list

        '''
        print >>diag, "> Appends FSL config's outputModelFilesDirectory: ", conf.outputModelFilesDirectory
        print >>diag, "> and FSL config's subjectListFile: ", conf.subjectListFile
        print >>diag, "> ..to model_sub_list, which is now: ", model_sub_list
        print >>diag, ""


    print >>diag, "> Ending modelConfigs iteration."
    print >>diag, ""
    '''
    
    if len(model_sub_list) == 0:
        raise Exception("no model found")





    #start group analysis

    #print >>diag, "> Starting iteration over 'model_sub's in model_sub_list: ", model_sub_list
    #print >>diag, ""

    for model_sub in model_sub_list:

        #print >>diag, "Current model_sub: ", model_sub
        #print >>diag, ""
        
        model, subject_list = model_sub
        
        print "running for model %s and resource %s..."%(os.path.basename(model), resource)
        

        if not os.path.exists(model):
            raise Exception("path to model %s doesn't exit"%model)
        
        if not os.path.exists(subject_list):
            raise Exception("path to input subject list %s is invalid"%subject_list)
        
        #if c.mixedScanAnalysis == True:
        #    wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s'%(resource, os.path.basename(model)))
        #else:

        wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s/%s'%(resource, os.path.basename(model), scan_ids[0])) 

        wf.base_dir = c.workingDirectory
        wf.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}
        log_dir = os.path.join(c.outputDirectory, 'logs', 'group_analysis', resource, 'model_%s'%(os.path.basename(model)))
        try:
            os.makedirs(log_dir)
        except:
            print "log_dir already exist"
        



        # enable logging
    
        from nipype import config
        from nipype import logging
        
        config.update_config({'logging': {'log_directory': log_dir,
                              'log_to_file': True}})
        
        # Temporarily disable until solved
        #logging.update_logging(config)

        iflogger = logging.getLogger('interface')
    
    
        group_sublist = open(subject_list, 'r')

        #print >>diag, "> Opened subject list: ", subject_list
        #print >>diag, ""

        sublist_items = group_sublist.readlines()

        input_subject_list = [line.rstrip('\n') for line in sublist_items \
                              if not (line == '\n') and not line.startswith('#')]

        '''
        print >>diag, "Subject list run-through #2: ", input_subject_list
        print >>diag, ""    


        print >>diag, "> Now iterating over each subject in this subject list."
        print >>diag, "> For each subject, for each path in s_paths, it checks if the current"
        print >>diag, "> subject is in path, and then appends it into ordered_paths."
        print >>diag, ""
        '''

        ordered_paths = []
        pathcount = 0
        subcount = 0
        for sub in input_subject_list :
            subcount += 1
            for path in s_paths:

                if sub in path:
                    pathcount += 1
                    ordered_paths.append(path)

        '''
        print >>diag, "Ordered paths: ", ordered_paths
        print >>diag, ""

        print >>diag, "> The list ordered_paths is then fed into the group analysis"
        print >>diag, "> workflow input 'zmap_files'."
        print >>diag, ""
        '''

        print "Ordered paths length (number of subjects): ", len(ordered_paths)
        print ""
      
        iflogger.info("input_subject_list -> %s"%input_subject_list)
        print ""
    
        strgy_path = os.path.dirname(s_paths[0]).split(scan_ids[0])[1]

        for ch in ['.']:
            if ch in strgy_path:
                strgy_path = strgy_path.replace(ch, '_')

        print "strgy_path: ", strgy_path
        print ""
        


        # gp_flow
        # Extracts the model files (.con, .grp, .mat, .fts) from the model
        # directory and sends them to the create_group_analysis workflow gpa_wf

        gp_flow = create_grp_analysis_dataflow("gp_dataflow%s"%strgy_path)
        gp_flow.inputs.inputspec.grp_model = model
        gp_flow.inputs.inputspec.ftest = c.fTest
        
        
        
        # gpa_wf
        # Creates the actual group analysis workflow

        gpa_wf = create_group_analysis(c.fTest, "gp_analysis%s"%strgy_path)

        gpa_wf.inputs.inputspec.zmap_files = ordered_paths
        gpa_wf.inputs.inputspec.z_threshold = c.zThreshold
        gpa_wf.inputs.inputspec.p_threshold = c.pThreshold
        gpa_wf.inputs.inputspec.parameters = (c.FSLDIR,
                                                   'MNI152')
    
        print "group model: ", model
        print "f test: ", c.fTest
        print "z threshold: ", c.zThreshold
        print "p threshold: ", c.pThreshold
        print "parameters: ", (c.FSLDIR, 'MNI152')

    
        wf.connect(gp_flow, 'outputspec.mat',
                   gpa_wf, 'inputspec.mat_file')
        wf.connect(gp_flow, 'outputspec.con',
                   gpa_wf, 'inputspec.con_file')
        wf.connect(gp_flow, 'outputspec.grp',
                    gpa_wf, 'inputspec.grp_file')

            
        if c.fTest:
            wf.connect(gp_flow, 'outputspec.fts',
                       gpa_wf, 'inputspec.fts_file')
        


        # ds
        # Creates the datasink node for group analysis
        
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
        
#         if c.mixedScanAnalysis == True:
#             out_dir = re.sub(r'(\w)*scan_(\w)*(\d)*(\w)*[/]', '', out_dir)
              
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
        if c.fTest:
            wf.connect(gp_flow, 'outputspec.fts',
                       ds, 'model_files.@0') 
        
        wf.connect(gp_flow, 'outputspec.mat',
                   ds, 'model_files.@1' )
        wf.connect(gp_flow, 'outputspec.con',
                   ds, 'model_files.@2')
        wf.connect(gp_flow, 'outputspec.grp',
                   ds, 'model_files.@3')
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



        # Run the actual group analysis workflow with the amount
        # of processors determined in the pipeline_config file

        try:

            wf.run()#(plugin='MultiProc',
                     #            plugin_args={'n_procs': c.numCoresPerSubject})

        except:

            print "Error: Group analysis workflow run command did not complete successfully."
            print "subcount: ", subcount
            print "pathcount: ", pathcount
            print "sublist: ", sublist_items
            print "input subject list: ", input_subject_list
            print "conf: ", conf.subjectListFile
            
            raise Exception
    
        print "**Workflow finished for model %s and resource %s"%(os.path.basename(model), resource)


    #diag.close()



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


