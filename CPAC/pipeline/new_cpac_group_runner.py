




def load_config_yml(config_file):

	# loads a configuration YAML file
    #
    # input
    #   config_file: full filepath to YAML (.yml) file
    #
    # output
    #   config: Configuration object

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
              "list should be a text file (.text).\nPath provided: %s\n\n" \
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



def assign_output_name(fullpath, split_fullpath, deriv_folder_name):

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


    return output_name



def run(config_file, output_path_file):
    
    # Runs group analysis

    import os

    ind_outputs = []
    ga_configs = []
    ga_partlists = []
    output_paths = {}
    matched_parts = {}
    missing_parts = {}

    analysis_map_gp = {}

    # Load the MAIN PIPELINE config file into 'c'
    c = load_config_yml(config_file)

    if len(c.modelConfigs) == 0:
        print '[!] CPAC says: You do not have any models selected ' \
              'to run for group-level analysis. Return to your pipeline ' \
              'configuration file and create or select at least one.\n\n'
        raise Exception


    # load the group model configs and the group participant lists into
    # dictionaries, with the model names as keys
    for group_config_file in c.modelConfigs:

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

        # gather which outputs the user has selected to run group analysis for
        #     this cuts down on parsing time for output paths
        for output_type in ga_config.derivative_list:
            if output_type not in ind_outputs:
                ind_outputs.append(output_type)

        ga_configs[ga_config.model_name] = (ga_config, group_config_file)


    # now the group participant lists
    #     NOTE: it is assumed the participant IDs listed in these text files
    #           are in the format of {participant}_{site}_{session} !!!
    for ga_config_tuple in ga_configs.values():

        ga_config = ga_config_tuple[0]

        ga_partlist = load_group_participant_list(ga_config.subject_list)

        ga_partlists[ga_config.model_name] = ga_partlist


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
                  
                for group_model in ga_configs.keys():

                    ga_config = ga_configs[group_model][0]
                    ga_config_file = ga_configs[group_model][1]

                    # include output path if the participant is in the group
                    # model's participant list, and if the derivative is one
                    # of the group model's selected derivatives
                    if (unique_part_ID in ga_partlists[group_model]) and \
                        (deriv_folder_name in ga_config.derivative_list):

                        # create "output_name" for analysis_map_gp keying
                        output_name = assign_output_name(fullpath, \
                        	              split_fullpath, deriv_folder_name)
        
                        key = (group_model, ga_config_file)

                        if key not in output_paths.keys():
                        	output_paths[key] = {}

                        if output_name not in output_paths[key].keys():
                        	output_paths[key][output_name] = []

                        output_paths[key][output_name].append(fullpath)

                        # keep track of which participants from the
                        # participant list have been found in the output paths
                        if key not in matched_parts.keys():
                            matched_parts[key] = {}

                        if output_name not in matched_parts[key].keys():
                        	matched_parts[key][output_name] = []

                        # WARNING: losing information here - if there are
                        #          multiple series/scans, and one of them did
                        #          not complete during individual level, then
                        #          it will not be reported as missing below -
                        #          however, the missing participant/scan will
                        #          be obvious in the updated participant lists
                        #          in the model_files outputs
                        if unique_part_ID not in \
                            matched_parts[key][output_name]:

                            matched_parts[key][output_name].append( \
                            	unique_part_ID)


    for key_tuple in output_paths:

        for output_name in output_paths[key_tuple]:

            if len(output_paths[key_tuple]) == 0:
                err = "[!] CPAC says: No individual-level analysis outputs were "\
                      "found given the selected derivatives and pipeline output "\
                      "directory path you provided.\n\nPipeline Output " \
                      "Directory provided: %s\nSelected derivative: %s\n\n" \
                      "Either make sure your selections are correct, or that " \
                      "individual-level analysis completed successfully for the "\
                      "derivative in question.\n\n" \
                  % (output_path_file, output_name)
            raise Exception(err)


    # compare what's been matched to what's in the participant list
    for key_tuple in matched_parts.keys():

        group_model = key_tuple[0]
        matched_list = matched_parts[key_tuple]

        if key_tuple not in missing_parts.keys():
            missing_parts[key_tuple] = {}

        for output_name in matched_parts[key_tuple].keys():

            #if output_name not in missing_parts[key_tuple].keys():
            # 	missing_parts[key_tuple][output_name] = []

            missing_list = ga_partlists[group_model] - matched_list

        	missing_parts[key_tuple][output_name] = missing_list


    # report missing participants
    if len(missing_parts) > 0:

        miss_msg = "\n\n[!] CPAC warns: Output files missing for " \
                   "participants included in group-level analysis:\n"

        for key_tuple in missing_parts.keys():
        	for output_name in missing_parts[key_tuple].keys():
            
            miss_msg = miss_msg + "\n" + output_name + "\n"
            
            miss_msg = miss_msg + missing_parts[key_tuple][output_name] + "\n"

        print miss_msg


    ''' output_paths can now function as analysis_map_gp '''
    ''' subject IDs etc. can be gleaned on the other side '''
    ''' all it sends over is: group model name, group config file path, list of output paths relevant '''
    '''     already reduced from available subjects and whats in the subject list '''

        
