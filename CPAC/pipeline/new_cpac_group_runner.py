

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
    full_deriv_list = []

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

        # create a list of all derivatives being run for group analysis
        # across ALL group models - this is used to prune the search width
        # in the output path collection
        for deriv in ga_config.derivative_list:
            if deriv not in full_deriv_list:
                full_deriv_list.append(deriv)


    # now the group participant lists
    #     NOTE: it is assumed the participant IDs listed in these text files
    #           are in the format of {participant}_{site}_{session} !!!
    for ga_config_tuple in ga_configs.values():

        ga_config = ga_config_tuple[0]

        ga_partlist = load_group_participant_list(ga_config.subject_list)

        ga_partlists[ga_config.model_name] = ga_partlist


    return ga_configs, ga_partlists, full_deriv_list



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



def gather_nifti_globs(pipeline_output_folder, resource_list):

    # the number of directory levels under each participant's output folder
    # can vary depending on what preprocessing strategies were chosen, and
    # there may be several output filepaths with varying numbers of directory
    # levels

    # this parses them quickly while also catching each preprocessing strategy

    import os
    import glob
    from __builtin__ import all as b_all

    ext = ".nii"
    nifti_globs = []

    if len(resource_list) == 0:
        err = "\n\n[!] No derivatives selected!\n\n"
        raise Exception(err)

    # remove any extra /'s
    pipeline_output_folder = pipeline_output_folder.rstrip("/")

    # grab MeanFD_Jenkinson just in case
    resource_list.append("power_params")

    print "\n\nGathering the output file paths from %s..." \
          % pipeline_output_folder

    for resource_name in resource_list:

        glob_string = os.path.join(pipeline_output_folder, "*", \
                                       resource_name, "*", "*")

        # get all glob strings that result in a list of paths where every path
        # ends with a NIFTI file
        
        prog_string = ".."

        while len(glob.glob(glob_string)) != 0:

            if b_all(ext in x for x in glob.glob(glob_string)) == True:
                nifti_globs.append(glob_string)
        
            glob_string = os.path.join(glob_string, "*")
            prog_string = prog_string + "."
            print prog_string

    if len(nifti_globs) == 0:
        err = "\n\n[!] No output filepaths found in the pipeline output " \
              "directory provided for the derivatives selected!\n\nPipeline "\
              "output directory provided: %s\nDerivatives selected:\s\n\n" \
              % (pipeline_output_folder, resource_list)
        raise Exception(err)

    return nifti_globs



def create_output_dict_list(nifti_globs, pipeline_output_folder, \
    session_list=None, get_raw_score=False):

    import glob

    # parse each result of each "valid" glob string
    output_dict_list = {}
    output_df_dict = {}

    for nifti_glob_string in nifti_globs:

        nifti_paths = glob.glob(nifti_glob_string)   

        for filepath in nifti_paths:
        
            second_half_filepath = filepath.split(pipeline_output_folder)[1]
            filename = filepath.split("/")[-1]
            
            resource_id = second_half_filepath.split("/")[2]
            series_id_string = second_half_filepath.split("/")[3]
            strat_info = second_half_filepath.split(series_id_string)[1]
            
            unique_resource_id = (resource_id,strat_info)
                        
            if unique_resource_id not in output_dict_list.keys():
                output_dict_list[unique_resource_id] = []
            
            unique_id = second_half_filepath.split("/")[1]

            series_id = series_id_string.replace("_scan_","")
            series_id = series_id.replace("_rest","")
            
            new_row_dict = {}
            
            new_row_dict["Participant"] = unique_id
            new_row_dict["Series"] = series_id
            
            if session_list:
                for session in session_list:
                    if session in second_half_filepath:
                        new_row_dict["Session"] = session
                        break
                        
            new_row_dict["Filepath"] = filepath
                        
            if get_raw_score:
                # grab raw score for measure mean just in case
                if "vmhc" in resource_id:
                    raw_score_path = filepath.replace(resource_id,"vmhc_raw_score")
                    raw_score_path = raw_score_path.replace(raw_score_path.split("/")[-1],"")
                    raw_score_path = glob.glob(os.path.join(raw_score_path,"*"))[0]
                else:                   
                    raw_score_path = filepath.replace("_zstd","")
                    raw_score_path = raw_score_path.replace("_fisher","")
                    raw_score_path = raw_score_path.replace("_zstat","")
                    
                    if "sca_roi_files_to_standard" in resource_id:
                        sub_folder = raw_score_path.split("/")[-2] + "/"
                        if "z_score" in sub_folder:
                            raw_score_path = raw_score_path.replace(sub_folder,"")
                    elif "sca_tempreg_maps_zstat" in resource_id:
                        sca_filename = raw_score_path.split("/")[-1]
                        globpath = raw_score_path.replace(sca_filename, "*")
                        globpath = os.path.join(globpath, sca_filename)
                        raw_score_path = glob.glob(globpath)[0]     
                    elif "dr_tempreg_maps" in resource_id:
                        raw_score_path = raw_score_path.replace("map_z_","map_")
                        raw_filename = raw_score_path.split("/")[-1]
                        raw_score_path = raw_score_path.replace(raw_filename,"")
                        raw_score_path = glob.glob(os.path.join(raw_score_path,"*",raw_filename))[0]       
                    else:
                        # in case filenames are different between z-standardized and raw
                        raw_score_path = raw_score_path.replace(raw_score_path.split("/")[-1],"")
                        raw_score_path = glob.glob(os.path.join(raw_score_path,"*"))[0]
                
                if not os.path.exists(raw_score_path):
                    err = "\n\n[!] The filepath for the raw score of " \
                          "%s can not be found.\nFilepath: %s\n\nThis " \
                          "is needed for the Measure Mean calculation." \
                          "\n\n" % (resource_id, raw_score_path)
                    raise Exception(err)
                    
                new_row_dict["Raw_Filepath"] = raw_score_path
                       
            # unique_resource_id is tuple (resource_id,strat_info)
            output_dict_list[unique_resource_id].append(new_row_dict)

    return output_dict_list



def create_output_df_dict(output_dict_list, inclusion_list=None, \
    session_list=None, series_list=None):

    import pandas as pd

    for unique_resource_id in output_dict_list.keys():
    
        new_df = pd.DataFrame(output_dict_list[unique_resource_id])
        
        # drop whatever is not in the inclusion lists
        if inclusion_list:
            new_df = new_df[new_df.Participant.isin(inclusion_list)]
        
        if series_list:
            new_df = new_df[new_df.Series.isin(series_list)]
            
        if session_list:
            new_df = new_df[new_df.Session.isin(session_list)]
    
        # unique_resource_id is tuple (resource_id,strat_info)
        if unique_resource_id not in output_df_dict.keys():
            output_df_dict[unique_resource_id] = new_df
            
            
    return output_df_dict



def load_pheno_csv_into_df(pheno_file):

    import os
    import pandas as pd

    if not os.path.isfile(pheno_file):
        err = "\n\n[!] CPAC says: The group-level analysis phenotype file "\
              "provided does not exist!\nPath provided: %s\n\n" \
              % pheno_file
        raise Exception(err)

    with open(os.path.abspath(pheno_file),"r") as f:
        pheno_dataframe = pd.read_csv(f)


    return pheno_dataframe



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

    # Load the MAIN PIPELINE config file into 'c' as a CONFIGURATION OBJECT
    c = load_config_yml(config_file)

    if len(c.modelConfigs) == 0:
        print '[!] CPAC says: You do not have any models selected ' \
              'to run for group-level analysis. Return to your pipeline ' \
              'configuration file and create or select at least one.\n\n'
        raise Exception

    # sammin sammin mmmm samin in gray v

    # load the group model configs and the group participant lists into
    # dictionaries, with the model names as keys

    # BREAK THIS OUT INTO A MAP!!!!!!!!
    ga_configs, ga_partlists, full_deriv_list = collect_group_config_info(c)


    # gather output file paths from individual level analysis
    # return only the ones appropriate for each model and participant list
    output_paths, matched_parts = \
        collect_output_paths(output_path_file, ga_configs, ga_partlists, \
                                 full_deriv_list)


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

                    # resource = general derivative name, like "sca_roi_" etc.
                    # output_name = derivative name and details (mask, ROI,
                    #               type of centrality, etc.)
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
            		             output_paths[group_model_config][series_id][output_name], \
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
    
