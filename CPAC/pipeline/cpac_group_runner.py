

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



def load_text_file(filepath, label="file"):

    # loads a text file and returns the lines in a list
    #
    # input
    #   filepath: full filepath to the text file
    #
    # output
    #   lines_list: list of lines from text file

    if not filepath.endswith(".txt"):
        err = "\n\n[!] CPAC says: The %s should be a text file (.txt).\n" \
              "Path provided: %s\n\n" % (label, filepath)
        raise Exception(err)

    try:
        with open(filepath,"r") as f:
            lines_list = f.readlines()
    except Exception as e:
    	err = "\n\n[!] CPAC says: Could not load or read the %s:\n%s\n" \
              "Details: %s\n\n" % (label, filepath, e)
    	raise Exception(err)

    # get rid of those \n's that love to show up everywhere
    lines_list = [i.rstrip("\n") for i in lines_list]

    return lines_list



def load_pheno_csv_into_df(pheno_file):

    import os
    import pandas as pd

    if not os.path.isfile(pheno_file):
        err = "\n\n[!] CPAC says: The group-level analysis phenotype file "\
              "provided does not exist!\nPath provided: %s\n\n" \
              % pheno_file
        raise Exception(err)

    if not pheno_file.endswith(".csv"):
        err = "\n\n[!] CPAC says: The group-level analysis phenotype " \
              "file should be a CSV file (.csv).\nPath provided: %s\n\n" \
              % pheno_file
        raise Exception(err)

    with open(os.path.abspath(pheno_file),"r") as f:
        pheno_dataframe = pd.read_csv(f)


    return pheno_dataframe



def gather_nifti_globs(pipeline_output_folder, resource_list):

    # the number of directory levels under each participant's output folder
    # can vary depending on what preprocessing strategies were chosen, and
    # there may be several output filepaths with varying numbers of directory
    # levels

    # this parses them quickly while also catching each preprocessing strategy

    import os
    import glob
    from __builtin__ import any as b_any

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

            if b_any(ext in x for x in glob.glob(glob_string)) == True:
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



def grab_raw_score_filepath(filepath, resource_id):

    # this lives in the output path collector

    import os
    import glob

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
            try:
                raw_score_path = glob.glob(os.path.join(raw_score_path,"*"))[0]
            except:
                raw_score_path = os.path.join(raw_score_path,"*")
                
    if (raw_score_path is None) or (not os.path.exists(raw_score_path)):
        err = "\n\n[!] The filepath for the raw score of " \
              "%s can not be found.\nFilepath: %s\n\nThis " \
              "is needed for the Measure Mean calculation." \
              "\n\n" % (resource_id, raw_score_path)
        raise Exception(err)

    return raw_score_path



def find_power_params_file(filepath, resource_id, series_id):

    import os

    try:
        power_path = filepath.replace(resource_id, "power_params", 1)
        series_id_string = "_scan_%s" % series_id
        power_first_half = power_path.split(series_id_string)[0]
        power_first_half = os.path.join(power_first_half, series_id_string)
        participant_id = power_first_half.split("/")[-3]
    except Exception as e:
        err = "\n\n[!] Something went wrong with finding the power " \
              "parameters file for at least one of the participants.\n\n" \
              "Error details: %s\n\n" % e
        raise Exception(err)
    
    power_params_file = None
    for root, dirs, files in os.walk(power_first_half):
        for filename in files:
            filepath = os.path.join(root, filename)
            if "pow_params.txt" in filepath:
                power_params_file = filepath

    if not power_params_file:
        err = "\n\n[!] Could not find the power parameters file for the " \
              "following participant and series..\nParticipant: %s\n" \
              "Series: %s\n\nIt should be available here: %s\n\n" \
              % (participant_id, series_id, power_first_half)
        raise Exception(err)

    return power_params_file



def extract_power_params(power_params_lines, power_params_filepath):

    # check formatting
    if len(power_params_lines) != 2:
        err = "\n\n[!] There is something wrong with the formatting of the " \
              "power parameters file.\nFilepath: %s\n\n" \
              % power_params_filepath
        raise Exception(err)

    names_list = power_params_lines[0].split(",")
    values_list = power_params_lines[1].split(",")

    # let's make extra sure
    if (values_list[0] not in power_params_filepath) or \
        (values_list[1] not in power_params_filepath):
        err = "\n\n[!] There is a mismatch between the contents of the " \
              "power parameters file and where it is located!\n" \
              "Filepath: %s\n\n" % power_params_filepath
        raise Exception(err)

    if (names_list[2] != "MeanFD_Power") or \
        (names_list[3] != "MeanFD_Jenkinson") or \
            (names_list[-1] != "MeanDVARS"):
        err = "\n\n[!] There is a mismatch between the power parameters " \
              "format and what is expected!!\nFilepath: %s\n\n" \
              % power_params_filepath
        raise Exception(err)

    meanfd_power = values_list[2]
    meanfd_jenk = values_list[3]
    meandvars = values_list[-1]

    return meanfd_power, meanfd_jenk, meandvars
 


def create_output_dict_list(nifti_globs, pipeline_output_folder, \
                                get_motion=False, get_raw_score=False):

    import os
    import glob

    ext = ".nii"

    # parse each result of each "valid" glob string
    output_dict_list = {}
    output_df_dict = {}

    for nifti_glob_string in nifti_globs:

        nifti_paths = glob.glob(nifti_glob_string + ext + "*")   

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
                                   
            new_row_dict["Filepath"] = filepath
                        
            if get_motion:
                # if we're including motion measures
                power_params_file = find_power_params_file(filepath, \
                    resource_id, series_id)
                power_params_lines = load_text_file(power_params_file, \
                    "power parameters file")
                meanfd_p, meanfd_j, meandvars = \
                    extract_power_params(power_params_lines, \
                                         power_params_file)
                new_row_dict["MeanFD_Power"] = meanfd_p
                new_row_dict["MeanFD_Jenkinson"] = meanfd_j
                new_row_dict["MeanDVARS"] = meandvars

            if get_raw_score:
                # grab raw score for measure mean just in case
                raw_score_path = grab_raw_score_filepath(filepath, \
                                                         resource_id)                    
                new_row_dict["Raw_Filepath"] = raw_score_path
                       
            # unique_resource_id is tuple (resource_id,strat_info)
            output_dict_list[unique_resource_id].append(new_row_dict)

    return output_dict_list



def create_output_df_dict(output_dict_list, inclusion_list=None):

    import pandas as pd

    output_df_dict = {}

    # unique_resource_id is tuple (resource_id,strat_info)
    for unique_resource_id in output_dict_list.keys():
    
        new_df = pd.DataFrame(output_dict_list[unique_resource_id])
        
        # drop whatever is not in the inclusion lists
        if inclusion_list:
            new_df = new_df[new_df.Participant.isin(inclusion_list)]
                   
        # unique_resource_id is tuple (resource_id,strat_info)
        if unique_resource_id not in output_df_dict.keys():
            output_df_dict[unique_resource_id] = new_df
            
    return output_df_dict



def gather_outputs(pipeline_folder, resource_list, inclusion_list, \
                       get_motion, get_raw_score):

    # probably won't have a session list due to subject ID format!

    nifti_globs = gather_nifti_globs(pipeline_folder, resource_list)
    output_dict_list = create_output_dict_list(nifti_globs, pipeline_folder, \
                           get_motion, get_raw_score)
    output_df_dict = create_output_df_dict(output_dict_list, inclusion_list)

    return output_df_dict



def pheno_sessions_to_repeated_measures(pheno_df, sessions_list):

    # take in the selected sessions, and match them to the participant
    # unique IDs appropriately

    import pandas as pd

    # there are no new rows, since the phenotype file will have all of the
    # subject_site_session combo unique IDs on each row!!!
    sessions_col = []

    # participant IDs new columns
    participant_id_cols = {}
    i = 0

    for participant_unique_id in list(pheno_df["Participant"]):
        part_col = [0] * len(pheno_df["Participant"])
        for session in sessions_list:
            if session in participant_unique_id:
                # generate/update sessions categorical column
                sessions_col.append(session)
                part_id = participant_unique_id.replace(session, "")
                # generate/update participant ID column (1's or 0's)
                if part_id not in participant_id_cols.keys():
                    part_col[i] = 1
                    participant_id_cols["participant_%s" % part_id] = part_col
                else:
                    participant_id_cols["participant_%s" % part_id][i] = 1
        i += 1

    pheno_df["Session"] = sessions_col

    for new_col in participant_id_cols.keys():
        pheno_df[new_col] = participant_id_cols[new_col]

    return pheno_df



def pheno_series_to_repeated_measures(pheno_df, series_list, \
    repeated_sessions=False):

    # take in the selected series/scans, and create all of the permutations
    # of unique participant IDs (participant_site_session) and series/scans
    # and populate the pheno
    #   this is so the user does not have to have a specially-formatted
    #   version of the phenotype CSV for repeated measures; they can just
    #   enter the regular one

    import pandas as pd
   
    new_rows = []
    for series in series_list:
        sub_pheno_df = pheno_df.copy()
        sub_pheno_df["Series"] = series
        new_rows.append(sub_pheno_df)
    pheno_df = pd.concat(new_rows)

    if not repeated_sessions:

        # participant IDs new columns
        participant_id_cols = {}
        i = 0

        for participant_unique_id in pheno_df["Participant"]:
            part_col = [0] * len(pheno_df["Participant"])
            if participant_unique_id not in participant_id_cols.keys():
                part_col[i] = 1
                participant_id_cols["participant_%s" % participant_unique_id] = part_col
            else:
                participant_id_cols["participant_%s" % participant_unique_id][i] = 1
            i += 1

        for new_col in participant_id_cols.keys():
            pheno_df[new_col] = participant_id_cols[new_col]
        
    return pheno_df



def prep_analysis_df_dict(config_file, pipeline_output_folder):
    
    # Preps group analysis run

    # config_file: filepath to the main CPAC pipeline configuration YAML
    #              (not the group analysis config YAML)
    # pipeline_output_folder: filepath to the CPAC pipeline individual-level
    #                   analysis output directory
    #                   example:
    #                     /home/cpac_run_1/output/pipeline_040_ANTS

    import os
    import pandas as pd

    # Load the MAIN PIPELINE config file into 'c' as a CONFIGURATION OBJECT
    c = load_config_yml(config_file)

    if len(c.modelConfigs) == 0:
        print '[!] CPAC says: You do not have any models selected ' \
              'to run for group-level analysis. Return to your pipeline ' \
              'configuration file and create or select at least one.\n\n'
        raise Exception

    # load the group model configs
    group_models = []

    for group_config_file in c.modelConfigs:
        group_models.append((group_config_file, \
                             load_config_yml(group_config_file)))


    # get the lowest common denominator of group model config choices
    #   - create full participant list
    #   - create full output measure list
    #   - see if any of the models will require the raw scores
    full_inclusion_list = []
    full_output_measure_list = []
    get_motion = False
    get_raw_score = False
    
    for group_model_tuple in group_models:

        group_model = group_model_tuple[1]

        inclusion = load_text_file(group_model.participant_list, \
            "group-level analysis participant list")
        full_inclusion_list = full_inclusion_list + inclusion

        full_output_measure_list = full_output_measure_list + \
                                       group_model.derivative_list

        # if any of the models will require motion parameters
        if ("MeanFD" in group_model.design_formula) or \
            ("MeanDVARS" in group_model.design_formula):
            get_motion = True

        # make sure "None" gets processed properly here...
        if (group_model.custom_roi_mask == "None") or \
            (group_model.custom_roi_mask == "none"):
            custom_roi_mask = None
        else:
            custom_roi_mask = group_model.custom_roi_mask

        if ("Measure_Mean" in group_model.design_formula) or \
            (custom_roi_mask != None):
            get_raw_score = True

    full_inclusion_list = list(set(full_inclusion_list))
    full_output_measure_list = list(set(full_output_measure_list))

    # sammin sammin mmmm samin in gray v

    # create encompassing output dataframe dictionary
    #     note, it is still limited to the lowest common denominator of all
    #     group model choices- it does not pull in the entire output directory
    # - there will be a dataframe for each combination of output measure
    #   type and preprocessing strategy
    # - each dataframe will contain output filepaths and their associated
    #   information, and each dataframe will include ALL SERIES/SCANS
    # - the dataframes will be pruned for each model LATER
    output_df_dict = gather_outputs(pipeline_output_folder, \
                                        full_output_measure_list, \
                                        full_inclusion_list, \
                                        get_motion, \
                                        get_raw_score)


    # alright, group model processing time
    #   going to merge the phenotype DFs with the output file DF
    analysis_dict = {}

    group_model_names = []

    for group_model_tuple in group_models:

        group_config_file = group_model_tuple[0]
        group_model = group_model_tuple[1]

        model_name = group_model.model_name

        if model_name in group_model_names:
            err = "\n\n[!] You have two group analysis models with the same "\
                  "name!\n\nDuplicate name: %s\n\n" % model_name
            raise Exception(err)
        else:
            group_model_names.append(model_name)

        # load original phenotype CSV into a dataframe
        pheno_df = load_pheno_csv_into_df(group_model.pheno_file)

        # enforce the sub ID label to "Participant"
        pheno_df.rename(columns={group_model.participant_id_label:"Participant"}, \
                        inplace=True)   
        pheno_df["Participant"] = pheno_df["Participant"].astype(str)

        # unique_resource = (output_measure_type, preprocessing strategy)
        # output_df_dict[unique_resource] = dataframe
        for unique_resource in output_df_dict.keys():

            resource_id = unique_resource[0]

            if resource_id not in group_model.derivative_list:
                continue

            strat_info = unique_resource[1]
        
            # output_df has the information for ALL of the output files for
            # this unique_resource_id- all series, and if applicable, motion
            # params numbers, and paths to raw outputs (for measure mean or
            # custom ROI means)
            #   then cut it down and merge with the phenotype DF as needed
            #   depending on the analysis
            output_df = output_df_dict[unique_resource]

            # prune the output_df for this specific group model and output +
            # preprocessing strategy
            inclusion_list = load_text_file(group_model.participant_list, \
                "group-level analysis participant list")
            output_df = \
                output_df[output_df["Participant"].isin(inclusion_list)]

            new_pheno_df = pheno_df.copy()
            
            if group_model.repeated_measures == True:

                if group_model.repeated_sessions == True:
                    new_pheno_df = pheno_sessions_to_repeated_measures( \
                                       new_pheno_df, \
                                       group_model.sessions_list)

                # create new rows for all of the series, if applicable
                #   ex. if 10 subjects and two sessions, 10 rows -> 20 rows
                if group_model.repeated_series == True:
                    new_pheno_df = pheno_series_to_repeated_measures( \
                                       new_pheno_df, \
                                       group_model.series_list, \
                                       group_model.repeated_sessions)


                # drop the pheno rows - if there are participants missing in
                # the output files (ex. if ReHo did not complete for 2 of the
                # participants, etc.), then drop these rows from the phenotype
                #   we are dropping all instances of this participant, all
                #   sessions and all series, because in repeated measures/
                #   within-subject, if one goes, they all have to go    
                new_pheno_df = \
                    new_pheno_df[pheno_df["Participant"].isin(output_df["Participant"])]

                if len(new_pheno_df) == 0:
                    err = "\n\n[!] There is a mis-match between the "\
                          "participant IDs in the output directory/particip" \
                          "ant list and the phenotype file.\n\n"
                    raise Exception(err)

                join_columns = ["Participant"]

                # if Series is one of the categorically-encoded covariates,
                # make sure we only are including the series the user has
                # selected to include in the repeated measures analysis
                if "Series" in new_pheno_df:
                    # check in case the pheno has series IDs that doesn't
                    # exist in the output directory, first
                    new_pheno_df = \
                        new_pheno_df[new_pheno_df["Series"].isin(output_df["Series"])]
                    # okay, now check against the user-specified series list
                    new_pheno_df = \
                        new_pheno_df[new_pheno_df["Series"].isin(group_model.series_list)]
                    join_columns.append("Series")
                    # pull together the pheno DF and the output files DF!
                    new_pheno_df = pd.merge(new_pheno_df, output_df, how="inner",\
                        on=join_columns)

                    analysis_dict[(model_name, group_config_file, resource_id, strat_info, "repeated_measures_multiple_series")] = \
                        new_pheno_df

                else:
                    # split up the series here
                    # iterate over the Series/Scans
                    for series_df_tuple in output_df.groupby("Series"):
                        series = series_df_tuple[0]
                        # series_df = output_df but with only one of the Series
                        series_df = series_df_tuple[1]
                        # trim down the pheno DF to match the output DF and merge
                        newer_pheno_df = new_pheno_df[pheno_df["Participant"].isin(series_df["Participant"])]
                        newer_pheno_df = pd.merge(new_pheno_df, series_df, how="inner", on=["Participant"])
                        # unique_resource =
                        #              (output_measure_type, preprocessing strategy)
                        analysis_dict[(model_name, group_config_file, resource_id, strat_info, "repeated_measures_%s" % series)] = newer_pheno_df

            else:
                # no repeated measures

                # split up the output files list DataFrame by series, then
                # merge with the pheno DataFrame and send it off for analysis
            
                # essentially, make sure each series combination goes into its
                # own model (and dataframe) for this unique_resource_id

                # iterate over the Series/Scans
                for series_df_tuple in output_df.groupby("Series"):
                    series = series_df_tuple[0]
                    # series_df = output_df but with only one of the Series
                    series_df = series_df_tuple[1]
                    # trim down the pheno DF to match the output DF and merge
                    newer_pheno_df = new_pheno_df[pheno_df["Participant"].isin(series_df["Participant"])]
                    newer_pheno_df = pd.merge(new_pheno_df, series_df, how="inner", on=["Participant"])
                    # send it in
                    analysis_dict[(model_name, group_config_file, resource_id, strat_info, series)] = newer_pheno_df

    return analysis_dict



def run(config_file, pipeline_output_folder):

    import os
    from multiprocessing import Process

    # create the analysis DF dictionary
    analysis_dict = prep_analysis_df_dict(config_file, pipeline_output_folder)

    # get MAIN pipeline config loaded
    c = load_config_yml(config_file)

    # let's get the show on the road   
    procss = []
    
    for unique_resource_id in analysis_dict.keys():

        # unique_resource_id is a 5-long tuple:
        #    ( model name, group model config file, output measure name,
        #          preprocessing strategy string,
        #          series_id or "repeated_measures" )
        
        model_name = unique_resource_id[0]
        group_config_file = unique_resource_id[1]
        resource_id = unique_resource_id[2]
        preproc_strat = unique_resource_id[3]
        series_or_repeated = unique_resource_id[4]

        model_df = analysis_dict[unique_resource_id]

        if not c.runOnGrid:

            from CPAC.pipeline.cpac_ga_model_generator import \
                prep_group_analysis_workflow

            procss.append(Process(target=prep_group_analysis_workflow, \
                                  args = (model_df, config_file, model_name, \
                                          group_config_file, resource_id, \
                                          preproc_strat, series_or_repeated)))
            
        else:
            print "\n\n[!] CPAC says: Group-level analysis has not yet been "\
                  "implemented to handle runs on a cluster or grid.\n\n" \
                  "Please turn off 'Run CPAC On A Cluster/Grid' in order to "\
                  "continue with group-level analysis. This will submit " \
                  "the job to only one node, however.\n\nWe will update " \
                  "users on when this feature will be available through " \
                  "release note announcements.\n\n"  
          
    # start kicking it off
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
