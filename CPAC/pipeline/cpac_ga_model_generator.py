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



def prep_group_analysis_workflow(c, group_config_file, resource, subject_infos, threshold_val):
    
    #
    # this function runs once per output file during group analysis
    #

    import yaml
    import commands

    # p_id = a list of pipeline IDs, i.e. the name of the output folder for
    #        the strat
    
    # s_ids = a list of all the subject IDs

    # scan_ids = a list of scan IDs

    # s_paths = a list of all of the filepaths of this particular output
    #           file that prep_group_analysis_workflow is being called for

    p_id, s_ids, scan_ids, s_paths = (list(tup) for tup in zip(*subject_infos))

    try:
        group_conf = Configuration(yaml.load(open(os.path.realpath(group_config_file), 'r')))
    except Exception as e:
        err_string = "\n\n[!] CPAC says: Could not read group model " \
                     "configuration YML file. Ensure you have read access " \
                     "for the file and that it is formatted properly.\n\n" \
                     "Configuration file: %s\n\nError details: %s" \
                     % (group_config_file, e)
        raise Exception(err_string)

     
    group_sublist_file = open(group_conf.subject_list, 'r')

    group_sublist_items = group_sublist_file.readlines()

    group_sublist = [line.rstrip('\n') for line in group_sublist_items \
                          if not (line == '\n') and not line.startswith('#')]

    # list of subjects for which paths which DO exist
    exist_paths = []

    # paths to the actual derivatives for those subjects
    derivative_paths = []


    z_threshold = float(group_conf.z_threshold[0])

    p_threshold = float(group_conf.p_threshold[0])


    custom_confile = group_conf.custom_contrasts

    if ((custom_confile == None) or (custom_confile == '') or \
            ("None" in custom_confile)):

        if (len(group_conf.f_tests) == 0) or (group_conf.f_tests == None):
            fTest = False
        else:
            fTest = True

    else:

        if not os.path.exists(custom_confile):
            errmsg = "\n[!] CPAC says: You've specified a custom contrasts " \
                     ".CSV file for your group model, but this file cannot " \
                     "be found. Please double-check the filepath you have " \
                     "entered.\n\nFilepath: %s\n\n" % custom_confile
            raise Exception(errmsg)

        evs = open(custom_confile, 'r').readline()
        evs = evs.rstrip('\r\n').split(',')
        count_ftests = 0

        fTest = False

        for ev in evs:
            if "f_test" in ev:
                count_ftests += 1

        if count_ftests > 0:
            fTest = True



    ''' begin iteration through group subject list for processing '''

    print "Sorting through subject list to check for missing outputs " \
          "for %s..\n" % resource

    for ga_sub in group_sublist:
        # Strip out carriage-return character if it is there
        
        if ga_sub.endswith('\r'):
            ga_sub = ga_sub.rstrip('\r')

        # ga_sub = subject ID taken off the group analysis subject list

        # let's check to make sure the subject list is formatted for
        # repeated measures properly if repeated measures is enabled
        # and vice versa
        if (group_conf.repeated_measures == True) and (',' not in ga_sub):
            print '\n\n'
            print '[!] CPAC says: The group analysis subject list ' \
                  'is not in the appropriate format for repeated ' \
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

        elif (group_conf.repeated_measures == False) and (',' in ga_sub):
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
        if (group_conf.repeated_measures == True):

            # sub.count(',') equals 1 when there is either multiple scans
            # or multiple sessions but not both, for repeated measures

            # sub.count(',') equals 2 when there are multiple sessions
            # AND scans, for repeated measures

            if ga_sub.count(',') == 1:
                sub_id = ga_sub.split(',',1)[0]
                other_id = ga_sub.split(',',1)[1]

            elif ga_sub.count(',') == 2:
                sub_id = ga_sub.split(',',2)[0]
                scan_id = ga_sub.split(',',2)[1]
                session_id = ga_sub.split(',',2)[2]



        ''' drop subjects from the group subject list '''
        # check the path files in path_files_here folder in the
        # subject's output folder - and drop any subjects from the
        # group analysis subject list which do not exist in the paths
        # to the output files

        '''
        REVISIT THIS LATER to establish a potentially better way to
        pull output paths (instead of path_files_here)
        '''

        for path in s_paths:

            if (group_conf.repeated_measures == True):

                if ga_sub.count(',') == 1:
                    if (sub_id in path) and (other_id in path):
                        exist_paths.append(ga_sub)
                        derivative_paths.append(path)

                elif ga_sub.count(',') == 2:
                    if (sub_id in path) and (scan_id in path) and \
                            (session_id in path):
                        exist_paths.append(ga_sub)
                        derivative_paths.append(path)

            else:
                if ga_sub in path:
                    exist_paths.append(ga_sub)
                    derivative_paths.append(path)


        # END subject-dropping!

        if len(derivative_paths) == 0:
            print '\n\n\n[!] CPAC says: None of the subjects listed in the ' \
                  'group analysis subject list were found to have outputs ' \
                  'produced by individual-level analysis.\n\nEnsure that ' \
                  'the subjects listed in your group analysis subject list ' \
                  'are the same as the ones included in the individual-' \
                  'level analysis you are running group-level analysis for.' \
                  '\n\n\n'
            raise Exception

    ''' END subject list iteration '''
 

    # check to see if any derivatives of subjects are missing
    if len(list(set(group_sublist) - set(exist_paths))) >0:
        print "List of outputs missing for subjects:"
        print list(set(group_sublist) - set(exist_paths))
        print "..for derivatives:"
        print resource
        print "..at paths:"
        print os.path.dirname(s_paths[0]).replace(s_ids[0], '*')

        

    # create the path string for the group analysis output
    out_dir = os.path.dirname(s_paths[0]).split(p_id[0] + '/')
    out_dir = os.path.join(group_conf.output_dir, out_dir[1])
    out_dir = out_dir.replace(s_ids[0], 'group_analysis_results_%s/_grp_model_%s'%(p_id[0],group_conf.model_name))

    model_out_dir = os.path.join(group_conf.output_dir, 'group_analysis_results_%s/_grp_model_%s'%(p_id[0],group_conf.model_name))

    mod_path = os.path.join(out_dir, 'model_files')


    if not os.path.isdir(mod_path):
        os.makedirs(mod_path)

        
    ''' write the new subject list '''
    new_sub_file = os.path.join(mod_path, os.path.basename(group_conf.subject_list))

    try:

        f = open(new_sub_file, 'w')
         
        for sub in exist_paths:
            print >>f, sub
        
        f.close()

    except:

        print "Error: Could not open subject list file: ", new_sub_file
        raise Exception


    group_conf.update('subject_list',new_sub_file)

    sub_id_label = group_conf.subject_id_label


    # Run 'create_fsl_model' script to extract phenotypic data from
    # the phenotypic file for each of the subjects in the subject list

    ''' get the motion statistics parameter file, if present '''
    # get the parameter file so it can be passed to create_fsl_model.py
    # so MeanFD or other measures can be included in the design matrix

    measure_list = ['MeanFD', 'MeanFD_Jenkinson', 'MeanDVARS']

    for measure in measure_list:
    
        if (measure in group_conf.design_formula):    

            parameter_file = os.path.join(c.outputDirectory, p_id[0], '%s%s_all_params.csv'%(scan_ids[0].strip('_'),threshold_val))

            if 1 in c.runGenerateMotionStatistics:

                if not os.path.exists(parameter_file):
                    print '\n\n[!] CPAC says: Could not find or open the motion ' \
                          'parameter file. This is necessary if you have included ' \
                          'any of the MeanFD measures in your group model.\n\n' \
                          'If Generate Motion Statistics is enabled, this file can ' \
                          'usually be found in the output directory of your ' \
                          'individual-level analysis runs. If it is not there, ' \
                          'double-check to see if individual-level analysis had ' \
                          'completed successfully.\n'
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
                    if (measure in group_conf.design_formula):
                        no_measures_error(measure)

                parameter_file = None
                
            break
            
    else:
    
        parameter_file = None



    # path to the pipeline folder to be passed to create_fsl_model.py
    # so that certain files like output_means.csv can be accessed
    pipeline_path = os.path.join(c.outputDirectory, p_id[0])

    # the current output that cpac_group_analysis_pipeline.py and
    # create_fsl_model.py is currently being run for
    current_output = resource #s_paths[0].replace(pipeline_path, '').split('/')[2]

    # generate working directory for this output's group analysis run
    workDir = '%s/group_analysis/%s/%s_%s' % (c.workingDirectory, group_conf.model_name, resource, scan_ids[0])

    # s_paths is a list of paths to each subject's derivative (of the current
    # derivative gpa is being run on) - s_paths_dirList is a list of each directory
    # in this path separated into list elements
             
    # this makes strgy_path basically the directory path of the folders after
    # the scan ID folder level         
    strgy_path = os.path.dirname(s_paths[0]).split(scan_ids[0])[1]

    # get rid of periods in the path
    for ch in ['.']:
        if ch in strgy_path:
            strgy_path = strgy_path.replace(ch, "")
                
    # create nipype-workflow-name-friendly strgy_path
    # (remove special characters)
    strgy_path_name = strgy_path.replace('/', "_")

    workDir = workDir + '/' + strgy_path_name



    ''' merge the remaining subjects for this current output '''
    # then, take the group mask, and iterate over the list of subjects
    # remaining to extract the mean of each subject using the group
    # mask

    merge_input = " "

    merge_output_dir = workDir + "/merged_files"

    if not os.path.exists(merge_output_dir):
        os.makedirs(merge_output_dir)

    merge_output = merge_output_dir + "/" + current_output + "_merged.nii.gz"
    merge_mask_output = merge_output_dir + "/" + current_output + "_merged_mask.nii.gz"

    # create a string per derivative filled with every subject's path to the
    # derivative output file
    for derivative_path in derivative_paths:
        merge_input = merge_input + " " + derivative_path
        
    merge_string = "fslmerge -t %s %s" % (merge_output, merge_input)

    # MERGE the remaining outputs
    try:
        commands.getoutput(merge_string)
    except Exception as e:
        print "[!] CPAC says: FSL Merge failed for output: %s" % current_output
        print "Error details: %s\n\n" % e
        raise

    merge_mask_string = "fslmaths %s -abs -Tmin -bin %s" % (merge_output, merge_mask_output)

    # CREATE A MASK of the merged file
    try:
        commands.getoutput(merge_mask_string)
    except Exception as e:
        print "[!] CPAC says: FSL Mask failed for output: %s" % current_output
        print "Error details: %s\n\n" % e
        raise


    derivative_means_dict = {}
    roi_means_dict = {}

    
    # CALCULATE THE MEANS of each remaining output using the group mask
    for derivative_path in derivative_paths:

        try:
            if "Group Mask" in group_conf.mean_mask:
                maskave_output = commands.getoutput("3dmaskave -mask %s %s" % (merge_mask_output, derivative_path))
            elif "Individual Mask" in group_conf.mean_mask:
                maskave_output = commands.getoutput("3dmaskave -mask %s %s" % (derivative_path, derivative_path))
        except Exception as e:
            print "[!] CPAC says: AFNI 3dmaskave failed for output: %s\n" \
                  "(Measure Mean calculation)" % current_output
            print "Error details: %s\n\n" % e
            raise

        # get the subject ID of the current derivative path reliably
        derivative_path_subID = derivative_path.replace(pipeline_path,"").strip("/").split("/")[0]

        # this crazy-looking command simply extracts the mean from the
        # verbose AFNI 3dmaskave output string
        derivative_means_dict[derivative_path_subID] = maskave_output.split("\n")[-1].split(" ")[0]

        # derivative_means_dict is now something like this:
        # { 'sub001': 0.3124, 'sub002': 0.2981, .. }

 
        # if custom ROI means are included in the model, do the same for those
        if "Custom_ROI_Mean" in group_conf.design_formula:

            try:
            
                if "centrality" in derivative_path:
                
                    # resample custom roi mask to 3mm, then use that
                    resampled_roi_mask = merge_output_dir + "/" + current_output + "_resampled_roi_mask.nii.gz"
                    
                    commands.getoutput("flirt -in %s -ref %s -o %s -applyxfm -init %s -interp nearestneighbour" % (group_conf.custom_roi_mask, derivative_path, resampled_roi_mask, c.identityMatrix))
                    
                    ROIstats_output = commands.getoutput("3dROIstats -mask %s %s" % (resampled_roi_mask, derivative_path))       
                    
                else:    
                        
                    ROIstats_output = commands.getoutput("3dROIstats -mask %s %s" % (group_conf.custom_roi_mask, derivative_path))
                    
            except Exception as e:
                print "[!] CPAC says: AFNI 3dROIstats failed for output: %s" \
                      "\n(Custom ROI Mean calculation)" % current_output
                print "Error details: %s\n\n" % e
                raise

            ROIstats_list = ROIstats_output.split("\t")

            # calculate the number of ROIs - 3dROIstats output can be split
            # into a list, and the actual ROI means begin at a certain point
            num_rois = (len(ROIstats_list)-3)/2

            roi_means = []

            # create a list of the ROI means - each derivative of each subject
            # will have N number of ROIs depending on how many ROIs were
            # specified in the custom ROI mask
            for num in range(num_rois+3,len(ROIstats_list)):

                roi_means.append(ROIstats_list[num])


            roi_means_dict[derivative_path_subID] = roi_means

        else:

            roi_means_dict = None



    if len(derivative_means_dict.keys()) == 0:
        err_string = "[!] CPAC says: Something went wrong with the " \
                     "calculation of the output means via the group mask.\n\n"
        raise Exception(err_string)
                     


    ''' run create_fsl_model.py to generate the group analysis models '''
    
    from CPAC.utils import create_fsl_model
    create_fsl_model.run(group_conf, fTest, parameter_file, derivative_means_dict, pipeline_path, current_output, model_out_dir, roi_means_dict, True)



    ''' begin GA workflow setup '''

    if not os.path.exists(new_sub_file):
        raise Exception("path to input subject list %s is invalid" % new_sub_file)
        
    #if c.mixedScanAnalysis == True:
    #    wf = pe.Workflow(name = 'group_analysis/%s/grp_model_%s'%(resource, os.path.basename(model)))
    #else:

    wf = pe.Workflow(name = resource)

    wf.base_dir = workDir
    wf.config['execution'] = {'hash_method': 'timestamp', 'crashdump_dir': os.path.abspath(c.crashLogDirectory)}
    log_dir = os.path.join(group_conf.output_dir, 'logs', 'group_analysis', resource, 'model_%s' % (group_conf.model_name))
        

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    else:
        pass


    # gp_flow
    # Extracts the model files (.con, .grp, .mat, .fts) from the model
    # directory and sends them to the create_group_analysis workflow gpa_wf

    gp_flow = create_grp_analysis_dataflow("gp_dataflow_%s" % resource)
    gp_flow.inputs.inputspec.grp_model = os.path.join(model_out_dir, "model_files", current_output)
    gp_flow.inputs.inputspec.model_name = group_conf.model_name
    gp_flow.inputs.inputspec.ftest = fTest
  

    # gpa_wf
    # Creates the actual group analysis workflow

    gpa_wf = create_group_analysis(fTest, "gp_analysis_%s" % resource)

    gpa_wf.inputs.inputspec.merged_file = merge_output
    gpa_wf.inputs.inputspec.merge_mask = merge_mask_output

    gpa_wf.inputs.inputspec.z_threshold = z_threshold
    gpa_wf.inputs.inputspec.p_threshold = p_threshold
    gpa_wf.inputs.inputspec.parameters = (c.FSLDIR, 'MNI152')
    
   
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
            re.search('sca_roi_(\d)+',os.path.splitext(os.path.splitext(os.path.basename(s_paths[0]))[0])[0]).group(0))
            
            
    if 'dr_tempreg_maps_zstat_files_to_standard_smooth' in resource:
        out_dir = os.path.join(out_dir, \
            re.search('temp_reg_map_z_(\d)+',os.path.splitext(os.path.splitext(os.path.basename(s_paths[0]))[0])[0]).group(0))
            
            
    if 'centrality' in resource:
        names = ['degree_centrality_binarize', 'degree_centrality_weighted', \
                 'eigenvector_centrality_binarize', 'eigenvector_centrality_weighted', \
                 'lfcd_binarize', 'lfcd_weighted']

        for name in names:
            if name in os.path.basename(s_paths[0]):
                out_dir = os.path.join(out_dir, name)
                break

    if 'tempreg_maps' in resource:
        out_dir = os.path.join(out_dir, \
            re.search('\w*[#]*\d+', os.path.splitext(os.path.splitext(os.path.basename(s_paths[0]))[0])[0]).group(0))
        
#     if c.mixedScanAnalysis == True:
#         out_dir = re.sub(r'(\w)*scan_(\w)*(\d)*(\w)*[/]', '', out_dir)
              
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
    
    print "**Workflow finished for model %s and resource %s"%(os.path.basename(group_conf.output_dir), resource)
        
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



