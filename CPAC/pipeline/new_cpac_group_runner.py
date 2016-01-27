

def load_config_yml(config_file):

	# loads a configuration YAML file
    #
    # input
    #   config_file: full filepath to YAML (.yml) file
    #
    # output
    #   config: Configuration object

    import os
    import yaml
    from CPAC.utils import Configuration

    try:

        config_path = os.path.realpath(config_file)

        with open(config_path,"r") as f:
            config_dict = yaml.load(f)

        config = Configuration(config_dict)

    except Exception as e:
        err = "\n\n[!] CPAC says: Could not load or read the configuration " \
        	  "YAML file:\n%s\nDetails: %s\n\n" % (config_file, e)
        raise Exception(err)


    return config



def load_group_participant_list(ga_part_list_filepath):

    # loads the group-level analysis participant list file (.txt)
    #
    # input
    #   ga_part_list_filepath: full filepath to the group-level analysis
    #                          participant list text file
    #
    # output
    #   ga_part_list: list of participant IDs to be included in group analysis

    # get the group participant list
    if not ga_part_list_filepath.endswith(".txt"):
        err = "\n\n[!] CPAC says: The group-level analysis participant " \
              "list should be a text file (.txt).\nPath provided: %s\n\n" \
              % ga_part_list_filepath
        raise Exception(err)

    try:
        with open(ga_part_list_filepath,"r") as f:
            ga_part_list = f.readlines()
    except Exception as e:
    	err = "\n\n[!] CPAC says: Could not load or read the group-level " \
    	      "analysis participant list:\n%s\nDetails: %s\n\n" \
    	      % (ga_part_list_filepath, e)
    	raise Exception(err)

    # get rid of those \n's that love to show up everywhere
    ga_part_list = [i.rstrip("\n") for i in ga_part_list]


    return ga_part_list



def collect_group_config_info(pipeline_config):

    ga_configs = {}
    ga_partlists = {}

    for group_config_file in pipeline_config.modelConfigs:

        ga_config = load_config_yml(group_config_file)

        if len(ga_config.derivative_list) == 0:
            err = "\n\n[!] CPAC says: You do not have any derivatives " \
                  "selected to run for group-level analysis. Return to " \
                  "your group-analysis configuration file and select at " \
                  "least one.\nGroup analysis configuration file: %s\n\n" \
                   % group_config_file
            raise Exception(err)

        if ga_config.model_name in ga_configs.keys():
            err = "\n\n[!] CPAC says: You have more than one group-level " \
                  "analysis model configuration with the same name!\n\n" \
                  "Duplicate model name: %s\n\n" % ga_config.model_name
            raise Exception(err)

        ga_configs[ga_config.model_name] = (ga_config, group_config_file)


    # now the group participant lists
    #     NOTE: it is assumed the participant IDs listed in these text files
    #           are in the format of {participant}_{site}_{session} !!!
    for ga_config_tuple in ga_configs.values():

        ga_config = ga_config_tuple[0]

        ga_partlist = load_group_participant_list(ga_config.subject_list)

        ga_partlists[ga_config.model_name] = ga_partlist


    return ga_configs, ga_partlists



def assign_output_name(fullpath, deriv_folder_name):

    '''
    if ("_mask_" in fullpath) and (("sca_roi" in fullpath) or \
        ("sca_tempreg" in fullpath) or ("centrality" in fullpath)):
            
        for dirname in split_fullpath:
            if "_mask_" in dirname:
                maskname = dirname
                    
        filename = split_fullpath[-1]
            
        if ".nii.gz" in filename:
            filename = filename.replace(".nii.gz","")
        elif ".nii" in filename:
            filename = filename.replace(".nii","")
            
        output_name = deriv_folder_name + "_%s_%s" % (maskname, filename)

            
    elif ("_spatial_map_" in fullpath) and ("dr_tempreg" in fullpath):
            
        for dirname in split_fullpath:
            if "_spatial_map_" in dirname:
                mapname = dirname
                    
        filename = split_fullpath[-1]
            
        if ".nii.gz" in filename:
            filename = filename.replace(".nii.gz","")
        elif ".nii" in filename:
            filename = filename.replace(".nii","")
            
        output_name = deriv_folder_name + "_%s_%s" % (mapname, filename)
                       
            
    else:
        
        output_name = deriv_folder_name
    '''

    # let's make sure the output file paths are stored in bins separately
    # for each combination of series ID, nuisance regression strategy,
    # bandpass filter cutoffs, etc.

    output_name = fullpath.split(deriv_folder_name)[1]
    output_name = deriv_folder_name + output_name

    if ".nii.gz" in output_name:
        output_name = output_name.replace(".nii.gz","")
    elif ".nii" in output_name:
        output_name = output_name.replace(".nii","") 


    return output_name



def collect_output_paths(output_path_file, ga_configs, ga_partlists):

    import os

    output_paths = {}
    matched_parts = {}
    matched_ders = []

    # collect relevant output paths from individual-level analysis   
    for root, folders, files in os.walk(output_path_file):
    
        split_output_dir_path = output_path_file.split("/")
    
        for filename in files:
        
            if ".nii" in filename:
    
                fullpath = os.path.join(root, filename)
            
                split_fullpath = fullpath.split("/")
                
                # unique_part_ID format: {participant}_{site}_{session} ID
                unique_part_ID = \
                    split_fullpath[len(split_output_dir_path)]

                deriv_folder_name = \
                    split_fullpath[len(split_output_dir_path)+1]

                series_id = \
                    split_fullpath[len(split_output_dir_path)+2]

                # for repeated measures (potentially)
                full_ID = unique_part_ID + "," + series_id
                  
                for group_model in ga_configs.keys():

                    ga_config = ga_configs[group_model][0]
                    ga_config_file = ga_configs[group_model][1]

                    # include output path if the participant is in the group
                    # model's participant list, and if the derivative is one
                    # of the group model's selected derivatives
                    if ((unique_part_ID in ga_partlists[group_model]) or \
                    	(full_ID in ga_partlists[group_model])) and \
                        (deriv_folder_name in ga_config.derivative_list):

                        # let's keep track of which derivatives have at least
                        # one participant completed for it
                        if deriv_folder_name not in matched_ders:
                        	matched_ders.append(deriv_folder_name)

                        # create "output_name" for analysis_map_gp keying
                        output_name = assign_output_name(fullpath, \
                        	                             deriv_folder_name)
        
                        key = (group_model, ga_config_file)

                        if key not in output_paths.keys():
                        	output_paths[key] = {}

                        if series_id not in output_paths[key].keys():
                        	output_paths[key][series_id] = {}

                        if output_name not in output_paths[key][series_id].keys():
                        	output_paths[key][series_id][output_name] = []

                        output_paths[key][series_id][output_name].append(fullpath)

                        # keep track of which participants from the
                        # participant list have been found in the output paths
                        if key not in matched_parts.keys():
                            matched_parts[key] = {}

                        if series_id not in matched_parts[key].keys():
                        	matched_parts[key][series_id] = {}

                        if output_name not in matched_parts[key][series_id].keys():
                        	matched_parts[key][series_id][output_name] = []

                        if unique_part_ID not in \
                            matched_parts[key][series_id][output_name]:

                            matched_parts[key][series_id][output_name].append(unique_part_ID)


    # see if there were no output files found for the derivative
    if len(ga_config.derivative_list) != len(matched_ders):

        empty_ders = set(ga_config.derivative_list) - set(matched_ders)

        err = "[!] CPAC says: No individual-level analysis outputs " \
              "were found for the following selected derivatives within " \
              "the pipeline output directory path you provided.\n\n" \
              "Pipeline Output Directory provided: %s\nDerivatives with no " \
              "completed participants:\n%s\n\nEither make sure your " \
              "selections are correct, or that individual-level analysis " \
              "completed successfully for the derivative in " \
              "question.\n\n" % (output_path_file, str(empty_ders))
        raise Exception(err)


    return output_paths, matched_parts



def report_missing_participants(matched_parts, ga_partlists):

    missing_parts = {}
    split_ga_partlists = {}

    for group_model in ga_partlists.keys():
        for part in ga_partlists[group_model]:
   	        if "," in part:
   	            unique_id = part.split(",")[0]
   	            series_id = part.split(",")[1]
   	            if group_model not in split_ga_partlists.keys():
   	                split_ga_partlists[group_model] = {}
   	            if series_id not in split_ga_partlists[group_model].keys():
   	                split_ga_partlists[group_model][series_id] = []
                split_ga_partlists[group_model][series_id].append(unique_id)
            else:
                # if not repeated measures, OR if repeated measures but only
                # for multiple sessions
                pass


    # compare what's been matched to what's in the participant list
    for key_tuple in matched_parts.keys():

        group_model = key_tuple[0]
        matched_list = matched_parts[key_tuple]

        for series_id in matched_parts[key_tuple].keys():

            for output_name in matched_parts[key_tuple][series_id].keys():

                missing_list = []

                if len(split_ga_partlists) > 0:
                	# if repeated measures with repeated series
                    for part in split_ga_partlists[group_model][series_id]:
                        if part not in matched_parts[key_tuple][series_id][output_name]:
                	        missing_list.append(part)
                else:
                    for part in ga_partlists[group_model]:
                        if part not in matched_parts[key_tuple][series_id][output_name]:
                	        missing_list.append(part)

                if len(missing_list) > 0:

                    if key_tuple not in missing_parts.keys():
                        missing_parts[key_tuple] = {}

                    if series_id not in missing_parts[key_tuple].keys():
                    	missing_parts[key_tuple][series_id] = {}

                    missing_parts[key_tuple][series_id][output_name] = missing_list


    # report missing participants
    if len(missing_parts) > 0:

        miss_msg = "\n\n[!] CPAC warns: Output files missing for " \
                   "participants included in group-level analysis:\n"

        for key_tuple in missing_parts.keys():
        	for series_id in missing_parts[key_tuple].keys():
                for output_name in missing_parts[key_tuple][series_id].keys():
            
                miss_msg = miss_msg + "\n%s - %s\n" % (output_name, series_id)
            
                miss_msg = miss_msg + str(missing_parts[key_tuple][series_id][output_name]) + "\n"

        print miss_msg



def run(config_file, output_path_file):
    
    # Runs group analysis

    import os 

    # Load the MAIN PIPELINE config file into 'c'
    c = load_config_yml(config_file)

    if len(c.modelConfigs) == 0:
        print '[!] CPAC says: You do not have any models selected ' \
              'to run for group-level analysis. Return to your pipeline ' \
              'configuration file and create or select at least one.\n\n'
        raise Exception

    # sammin sammin mmmm samin in gray v

    # load the group model configs and the group participant lists into
    # dictionaries, with the model names as keys
    ga_configs, ga_partlists = collect_group_config_info(c)


    # gather output file paths from individual level analysis
    # return only the ones appropriate for each model and participant list
    output_paths, matched_parts = \
        collect_output_paths(output_path_file, ga_configs, ga_partlists)


    # warn the user of any missing participants
    report_missing_participants(matched_parts, ga_partlists)


    # let's get the show on the road
    procss = []
    
    for group_model_config in output_paths.keys():

        # group_model_config is a tuple:
        #    (group model name, group analysis config YAML file)
        
        group_model_name = group_model_config[0]
        group_config_file = group_model_config[1]

        key = (group_model_name, group_config_file)

        for key in output_paths.keys():

        	for series_id in output_paths[key].keys():

        		for output_name in output_paths[key][series_id].keys():

        			resource = output_name.split("/")[0]

                    # just make sure the .split above worked - ideally we
                    # should never see this error
                    # if this error trips, something has changed in the output
                    # directory structure and this script was not adapted
                    # appropriately
                    ga_config = ga_configs[group_model_name][0]
        			if resource not in ga_config.derivative_list:
        				err = "\n\n[!] CPAC says: The individual-level " \
        				      "output file paths have not been parsed " \
        				      "correctly.\n\n%s\n\n" % str(locals())
        				raise Exception(err)

                    ''' how are we handling sending in combined series for repeated measures ????? '''

                    #get all the motion parameters across subjects
                    from CPAC.utils import extract_parameters
                    scrub_threshold = \
                        extract_parameters.run(c.outputDirectory, \
                        	                   c.runScrubbing)

                    if not c.runOnGrid:

                        from CPAC.pipeline.cpac_ga_model_generator import \
                            prep_group_analysis_workflow

                        procss.append(Process(target = \
                        	prep_group_analysis_workflow, args = \
                        	    (c, group_config_file, resource, \
            		             output_paths[group_model_config][SERIES???][output_name], \
            		             output_path_file, scrub_threshold)))
            
                    else:
        
                        print "\n\n[!] CPAC says: Group-level analysis has " \
                              "not yet been implemented to handle runs on a "\
                              "cluster or grid.\n\nPlease turn off 'Run " \
                              "CPAC On A Cluster/Grid' in order to continue "\
                              "with group-level analysis. This will submit " \
                              "the job to only one node, however.\n\nWe " \
                              "will update users on when this feature will " \
                              "be available through release note " \
                              "announcements.\n\n"

    
          
    pid = open(os.path.join(c.outputDirectory, 'pid_group.txt'), 'w')
                        
    jobQueue = []
    if len(procss) <= c.numGPAModelsAtOnce:
        """
        Stream all the subjects as sublist is
        less than or equal to the number of 
        subjects that need to run
        """
        for p in procss:
            p.start()
            print >>pid,p.pid
                
    else:
        """
        Stream the subject workflows for preprocessing.
        At Any time in the pipeline c.numSubjectsAtOnce
        will run, unless the number remaining is less than
        the value of the parameter stated above
        """
        idx = 0
        while(idx < len(procss)):
                
            if len(jobQueue) == 0 and idx == 0:
                
                idc = idx
                    
                for p in procss[idc: idc + c.numGPAModelsAtOnce]:
                
                    p.start()
                    print >>pid,p.pid
                    jobQueue.append(p)
                    idx += 1
                
            else:
                
                for job in jobQueue:
                
                    if not job.is_alive():
                        print 'found dead job ', job
                        loc = jobQueue.index(job)
                        del jobQueue[loc]
                        procss[idx].start()
                
                        jobQueue.append(procss[idx])
                        idx += 1
                
    pid.close()
    
    
