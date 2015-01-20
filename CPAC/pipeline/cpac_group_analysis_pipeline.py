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
    
    #
    # this function runs once per output file during group analysis
    #

    # p_id = a list of pipeline IDs, i.e. the name of the output folder for
    #        the strat
    
    # s_ids = a list of all the subject IDs

    # scan_ids = a list of scan IDs

    # s_paths = a list of all of the filepaths of this particular output
    #           file that prep_group_analysis_workflow is being called for

    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))


    # set this to False for now
    fTest = False

    def get_phenotypic_file(phenotypic_file, m_dict, m_list, mod_path, sub_id):
        
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

                print '\n\nMeasure: ', measure, '\n\n'

                if measure in headers:
                    #check if 'MeanFD  is present'
                    if len(columns[measure]) < 1:

                        print '\n\ncolumns[sub_id]: ', columns[sub_id], '\n\n'

                        for sub in columns[sub_id]:

                            if m_dict.get(sub):
                                if m_dict.get(sub).get(measure):
                                    columns[measure].append(m_dict[sub][measure])
                                else:
                                    raise Exception("Couldn't find %s value for subject %s"%(measure,sub))
                            else:
                                raise Exception("Couldn't find subject %s in the parameter file"%sub)


        print '\n\ncolumns[measure]: ', columns, '\n\n'
        
        b = zip(*([k] + columns[k] for k in sorted(columns, key=order.get)))
        
        
        try:
            os.makedirs(mod_path)
        except:
            print "%s already exists"%(mod_path)
            
        new_phenotypic_file = os.path.join(mod_path, os.path.basename(phenotypic_file))
                
        a = csv.writer(open(new_phenotypic_file, 'w'))
        
        for col in b:
            a.writerow(list(col))
          
        return new_phenotypic_file

    # END get_phenotypic_file function



    threshold_val = None
    measure_dict = None
    measure_list = ['MeanFD', 'MeanFD_Jenkinson', 'MeanDVARS']
    model_sub_list = []
    

    if 1 in c.runScrubbing:

        #get scrubbing threshold
    
        if re.search('(?<=/_threshold_)\d+.\d+',s_paths[0]):

            threshold_val = re.search('(?<=/_threshold_)\d+.\d+',s_paths[0]).group(0)

        elif len(c.scrubbingThreshold) == 1:

            threshold_val = c.scrubbingThreshold[0]

        else:
            print "Found Multiple threshold value "


        print "scrubbing threshold_val -->", threshold_val

    else:

        print "No scrubbing enabled."

        if len(c.scrubbingThreshold) == 1:
            threshold_val = c.scrubbingThreshold[0]




    import yaml    

    for config in c.modelConfigs:

        print c.modelConfigs
        print config
        
        try:
            conf = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
        except:
            raise Exception("Error in reading %s configuration file" % config)

        
        group_sublist = open(conf.subject_list, 'r')

        sublist_items = group_sublist.readlines()

        subject_list = [line.rstrip('\n') for line in sublist_items \
                              if not (line == '\n') and not line.startswith('#')]

        # list of subject paths which DO exist
        exist_paths = []




        ''' begin iteration through group subject list for processing '''

        for sub in subject_list:

            # let's check to make sure the subject list is formatted for
            # repeated measures properly if repeated measures is enabled and
            # vice versa
            if (c.repeatedMeasures == True) and (',' not in sub):
                print '\n\n'
                print '[!] CPAC says: The group analysis subject list is ' \
                        'not inthe appropriate format for repeated ' \
                        'measures.\n'
                print 'Please use the appropriate format as described in ' \
                        'the CPAC User Guide or turn off Repeated Measures ' \
                        'in the CPAC pipeline configuration editor, found ' \
                        'in the \'Group Analysis Settings\' tab of the ' \
                        'pipeline configuration editor.\n'
                print 'NOTE: CPAC generates a properly-formatted group ' \
                        'analysis subject list meant for running repeated ' \
                        'measures when you create your original subject ' \
                        'list. Look for \'subject_list_group_analysis_' \
                        'repeated_measures.txt\' in the directory where ' \
                        'you created your subject list.\n\n'
                raise Exception

            elif (c.repeatedMeasures == False) and (',' in sub):
                print '\n\n'
                print '[!] CPAC says: It looks like your group analysis ' \
                        'subject list is formatted for running repeated ' \
                        'measures, but \'Run Repeated Measures\' is not ' \
                        'enabled in the pipeline configuration, found in ' \
                        'the \'Group Analysis Settings\' tab of the ' \
                        'pipeline configuration editor.\n'
                print 'Double-check your pipeline configuration?\n\n'
                raise Exception



            ''' process subject ids for repeated measures, if it is on '''
            # if repeated measures is being run and the subject list
            # is a list of subject IDs and scan IDs concatenated
            if (c.repeatedMeasures == True):

                # sub.count(',') equals 1 when there is either multiple scans
                # or multiple sessions but not both, for repeated measures

                # sub.count(',') equals 2 when there are multiple sessions
                # AND scans, for repeated measures

                if sub.count(',') == 1:
                    sub_id = sub.split(',',1)[0]
                    other_id = sub.split(',',1)[1]

                elif sub.count(',') == 2:
                    sub_id = sub.split(',',2)[0]
                    scan_id = sub.split(',',2)[1]
                    session_id = sub.split(',',2)[2]



            ''' drop subjects from the group subject list '''
            # check the path files in path_files_here folder in the subject's
            # output folder - and drop any subjects from the group analysis
            # subject list which do not exist in the paths to the output files

            for path in s_paths:

                if (c.repeatedMeasures == True):

                    if sub.count(',') == 1:
                        if (sub_id in path) and (other_id in path):
                            exist_paths.append(sub)

                    elif sub.count(',') == 2:
                        if (sub_id in path) and (scan_id in path) and \
                                (session_id in path):
                            exist_paths.append(sub)

                else:

                    if sub in path:
                        exist_paths.append(sub)
 




        # check to see if any derivatives of subjects are missing
        if len(list(set(subject_list) - set(exist_paths))) >0:
            print "List of outputs missing for subjects:"
            print list(set(subject_list) - set(exist_paths))
            print "..for derivatives:"
            print resource
            print "..at paths:"
            print os.path.dirname(s_paths[0]).replace(s_ids[0], '*')

        

        # create the path string for the group analysis output
        out_dir = os.path.dirname(s_paths[0]).split(p_id[0] + '/')
        out_dir = os.path.join(conf.output_dir, out_dir[1])
        out_dir = out_dir.replace(s_ids[0], 'group_analysis_results_%s/_grp_model_%s'%(p_id[0],conf.model_name))

        mod_path = os.path.join(out_dir, 'model_files')


        if not os.path.isdir(mod_path):
            os.makedirs(mod_path)

        


        ''' write the new subject list '''
        new_sub_file = os.path.join(mod_path, os.path.basename(conf.subject_list))

        try:

            f = open(new_sub_file, 'w')
         
            for sub in exist_paths:
                print >>f, sub
        
            f.close()

        except:

            print "Error: Could not open subject list file: ", new_sub_file
            raise Exception


        conf.update('subject_list',new_sub_file)

        sub_id = conf.subject_id_label
        


        if measure_dict != None:
            conf.update('pheno_file',get_phenotypic_file(conf.pheno_file, measure_dict, measure_list, mod_path, sub_id))
        
        print 'conf updated pheno: ', conf.pheno_file, '\n\n'

            
        print "Model config dictionary ->"
        print conf.__dict__



        # Run 'create_fsl_model' script to extract phenotypic data from
        # the phenotypic file for each of the subjects in the subject list



        ''' get the motion statistics parameter file, if present '''
        # get the parameter file so it can be passed to create_fsl_model.py
        # so MeanFD or other measures can be included in the design matrix
        parameter_file = os.path.join(c.outputDirectory, p_id[0], '%s_threshold_%s_all_params.csv'%(scan_ids[0].strip('_'),threshold_val))

        if 1 in c.runGenerateMotionStatistics:

            if not os.path.exists(parameter_file):
                print '\n\n[!] CPAC says: Could not open the parameter file. ' \
                      'If Generate Motion Statistics is enabled, this can ' \
                      'usually be found in the output directory of your ' \
                      'individual-level analysis runs.\n'
                print 'Path not found: ', parameter_file, '\n\n'
                raise Exception

        elif (1 not in c.runGenerateMotionStatistics) and (os.path.exists(parameter_file)):

            if not os.path.exists(parameter_file):
                print '\n\n[!] CPAC says: Could not open the parameter file. ' \
                      'If Generate Motion Statistics is enabled, this can ' \
                      'usually be found in the output directory of your ' \
                      'individual-level analysis runs.\n'
                print 'Path not found: ', parameter_file, '\n\n'
                raise Exception

        else:

            def no_measures_error(measure):
                print '\n\n[!] CPAC says: The measure %s was included in ' \
                      'your group analysis design matrix formula, but ' \
                      'Generate Motion Statistics was not run during ' \
                      'individual-level analysis.\n' % measure
                print 'Please run Generate Motion Statistics if you wish ' \
                      'to include this measure in your model.\n'
                print 'If you HAVE completed a run with this option ' \
                      'enabled, then you are seeing this error because ' \
                      'the motion parameter file normally created by this ' \
                      'option is missing.\n\n'
                raise Exception

            for measure in measure_list:
                if (measure in conf.design_formula):
                    no_measures_error(measure)

            parameter_file = None



        ''' run create_fsl_model.py to generate the group analysis models '''
        # path to the pipeline folder to be passed to create_fsl_model.py
        # so that certain files like output_means.csv can be accessed
        pipeline_path = os.path.join(c.outputDirectory, p_id[0])

        # the current output that cpac_group_analysis_pipeline.py and
        # create_fsl_model.py is currently being run for
        current_output = s_paths[0].replace(pipeline_path, '').split('/')[2]


        try:

            from CPAC.utils import create_fsl_model

            create_fsl_model.run(conf, fTest, parameter_file, pipeline_path, current_output, True)

            #print >>diag, "> Runs create_fsl_model."
            #print >>diag, ""

        except Exception, e:

            print "FSL Group Analysis model not successfully created - error in create_fsl_model script"
            #print "Error ->", e
            raise


            
        model_sub_list.append((conf.output_dir, conf.subject_list))


    
    if len(model_sub_list) == 0:
        raise Exception("no model found")





    ''' start group analysis '''

    print '\n\nPreparing the group analysis workflow..\n\n'

    for model_sub in model_sub_list:

        #print >>diag, "Current model_sub: ", model_sub
        #print >>diag, ""
        
        model, subject_list = model_sub
   

        if not os.path.exists(model):
            raise Exception("path to model %s doesn't exist"%model)
        
        if not os.path.exists(subject_list):
            raise Exception("path to input subject list %s is invalid" % subject_list)
        
        #if c.mixedScanAnalysis == True:
        #    wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s'%(resource, os.path.basename(model)))
        #else:
        
        
        # s_paths is a list of paths to each subject's derivative (of the current
        # derivative gpa is being run on) - s_paths_dirList is a list of each directory
        # in this path separated into list elements
        s_paths_dirList = s_paths[0].split('/')
        
        currentDerivativeFile = s_paths_dirList[-1]
        
        currentDerivative = currentDerivativeFile.split('.')[0]
        
        currentDerivative = currentDerivative.replace('#', '_')
        
        
        strgy_path = os.path.dirname(s_paths[0]).split(scan_ids[0])[1]

        for ch in ['.']:
            if ch in strgy_path:
                strgy_path = strgy_path.replace(ch, '_')
                
        # create nipype-workflow-name-friendly strgy_path
        # (remove special characters)
        strgy_path_name = strgy_path.replace('/', '__')
        
        

        wf = pe.Workflow(name = currentDerivative) 

        workDir = c.workingDirectory + '/group_analysis__%s__grp_model_%s__%s' % (resource, conf.model_name, scan_ids[0])
        workDir = workDir + '/' + strgy_path_name

        wf.base_dir = workDir
        wf.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}
        log_dir = os.path.join(conf.output_dir, 'logs', 'group_analysis', resource, 'model_%s' % (conf.model_name))
        

        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        else:
            print "log_dir already exist"
        



        # enable logging
    
        from nipype import config
        from nipype import logging
        
        config.update_config({'logging': {'log_directory': log_dir,
                              'log_to_file': True}})
        
        # Temporarily disable until solved
        #logging.update_logging(config)

        iflogger = logging.getLogger('interface')




        ''' create the list of paths to all output files to go to model '''
        # create the 'ordered_paths' list, which is a list of all of the
        # output paths of the output files being included in the current
        # group-level analysis model
        #     'ordered_paths' is later connected to the 'zmap_files' input
        #     of the group analysis workflow - the files listed in this list
        #     are merged into the merged 4D file that goes into group analysis
      
        group_sublist = open(subject_list, 'r')
        sublist_items = group_sublist.readlines()

        input_subject_list = [line.rstrip('\n') for line in sublist_items \
                              if not (line == '\n') and not line.startswith('#')]

        ordered_paths = []
        pathcount = 0
        subcount = 0
        for sub in input_subject_list:

            subcount += 1

            if (c.repeatedMeasures == True):

                # sub.count(',') equals 1 when there is either multiple scans
                # or multiple sessions but not both, for repeated measures

                # sub.count(',') equals 2 when there are multiple sessions
                # AND scans, for repeated measures

                if sub.count(',') == 1:
                    sub_id = sub.split(',',1)[0]
                    other_id = sub.split(',',1)[1]

                elif sub.count(',') == 2:
                    sub_id = sub.split(',',2)[0]
                    scan_id = sub.split(',',2)[1]
                    session_id = sub.split(',',2)[2]


            for path in s_paths:

                if (c.repeatedMeasures == True):

                    # if repeated measures is enabled, make sure all of the
                    # relevant indicators are in the path before adding it
                    # to 'ordered_paths', i.e. the session and/or scan IDs

                    if sub.count(',') == 1:
                        if (sub_id in path) and (other_id in path):
                            pathcount += 1
                            ordered_paths.append(path)

                    elif sub.count(',') == 2:
                        if (sub_id in path) and (scan_id in path) and \
                                (session_id in path):
                            pathcount += 1
                            ordered_paths.append(path)

                else:
                    if sub in path:
                        pathcount += 1
                        ordered_paths.append(path)




        print 'S_paths length: ', len(s_paths)

        print "Ordered paths length (number of subjects): ", len(ordered_paths)
      
        print "input_subject_list -> %s" % input_subject_list

        print "strgy_path: ", strgy_path


        if len(ordered_paths) == 0:
            print '\n\n\n[!] CPAC says: None of the subjects listed in the ' \
                  'group analysis subject list were found to have outputs ' \
                  'produced by individual-level analysis.\n\nEnsure that ' \
                  'the subjects listed in your group analysis subject list ' \
                  'are the same as the ones included in the individual-' \
                  'level analysis you are running group-level analysis for.' \
                  '\n\n\n'
            raise Exception



        # gp_flow
        # Extracts the model files (.con, .grp, .mat, .fts) from the model
        # directory and sends them to the create_group_analysis workflow gpa_wf

        gp_flow = create_grp_analysis_dataflow("gp_dataflow_%s" % currentDerivative)
        gp_flow.inputs.inputspec.grp_model = model
        gp_flow.inputs.inputspec.fTest = fTest
  


        # gpa_wf
        # Creates the actual group analysis workflow

        gpa_wf = create_group_analysis(fTest, "gp_analysis_%s" % currentDerivative)

        gpa_wf.inputs.inputspec.zmap_files = ordered_paths
        gpa_wf.inputs.inputspec.z_threshold = c.zThreshold
        gpa_wf.inputs.inputspec.p_threshold = c.pThreshold
        gpa_wf.inputs.inputspec.parameters = (c.FSLDIR, 'MNI152')
    
        print "group model: ", model
        print "f test: ", fTest
        print "z threshold: ", c.zThreshold
        print "p threshold: ", c.pThreshold
        print "parameters: ", (c.FSLDIR, 'MNI152')

    
        wf.connect(gp_flow, 'outputspec.mat',
                   gpa_wf, 'inputspec.mat_file')
        wf.connect(gp_flow, 'outputspec.con',
                   gpa_wf, 'inputspec.con_file')
        wf.connect(gp_flow, 'outputspec.grp',
                    gpa_wf, 'inputspec.grp_file')

            
        if fTest:
            wf.connect(gp_flow, 'outputspec.fts',
                       gpa_wf, 'inputspec.fts_file')
        


        # ds
        # Creates the datasink node for group analysis
        
        ds = pe.Node(nio.DataSink(), name='gpa_sink')
     
        if 'sca_roi' in resource:
            out_dir = os.path.join(out_dir, \
              re.search('ROI_number_(\d)+',os.path.splitext(os.path.splitext(os.path.basename(s_paths[0]))[0])[0]).group(0))
            
        if 'centrality' in resource:
            names = ['degree_centrality_binarize', 'degree_centrality_weighted', \
                     'eigenvector_centrality_binarize', 'eigenvector_centrality_weighted', \
                     'lfcd_binarize', 'lfcd_weighted']
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
    
        '''
        if 1 in c.runSymbolicLinks:
    
    
            link_node = pe.MapNode(interface=util.Function(
                                input_names=['in_file',
                                            'resource'],
                                    output_names=[],
                                    function=prepare_gp_links),
                                    name='link_gp_', iterfield=['in_file'])
            link_node.inputs.resource = resource
            wf.connect(ds, 'out_file', link_node, 'in_file')
        '''
    


        ########datasink connections#########
        if fTest:
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

        # Run the actual group analysis workflow
        wf.run()

        '''
        except:

            print "Error: Group analysis workflow run command did not complete successfully."
            print "subcount: ", subcount
            print "pathcount: ", pathcount
            print "sublist: ", sublist_items
            print "input subject list: "
            print "conf: ", conf.subjectListFile
            
            raise Exception
        '''
    
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



