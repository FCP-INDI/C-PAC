

def load_config_yml(config_file, individual=False):

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

    if individual:
        config.logDirectory = os.path.abspath(config.logDirectory)
        config.workingDirectory = os.path.abspath(config.workingDirectory)
        config.outputDirectory = os.path.abspath(config.outputDirectory)
        config.crashLogDirectory = os.path.abspath(config.crashLogDirectory)

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


def read_pheno_csv_into_df(pheno_csv, id_label=None):
    """Read the phenotypic file CSV or TSV into a Pandas DataFrame."""

    import pandas as pd

    with open(pheno_csv, "r") as f:
        if id_label:
            if '.tsv' in pheno_csv or '.TSV' in pheno_csv:
                pheno_df = pd.read_table(f, dtype={id_label: object})
            else:
                pheno_df = pd.read_csv(f, dtype={id_label: object})
        else:
            if '.tsv' in pheno_csv or '.TSV' in pheno_csv:
                pheno_df = pd.read_table(f)
            else:
                pheno_df = pd.read_csv(f)

    return pheno_df


def gather_nifti_globs(pipeline_output_folder, resource_list):

    # the number of directory levels under each participant's output folder
    # can vary depending on what preprocessing strategies were chosen, and
    # there may be several output filepaths with varying numbers of directory
    # levels

    # this parses them quickly while also catching each preprocessing strategy

    import os
    import glob
    import pandas as pd
    import pkg_resources as p
    from __builtin__ import any as b_any

    ext = ".nii"
    nifti_globs = []

    keys_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')
    try:
        keys = pd.read_csv(keys_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(keys_csv, e)
        raise Exception(err)

    derivative_list = list(
        keys[keys['Derivative'] == 'yes'][keys['Space'] == 'template'][
            keys['Values'] == 'z-score']['Resource'])
    derivative_list = derivative_list + list(
        keys[keys['Derivative'] == 'yes'][keys['Space'] == 'template'][
            keys['Values'] == 'z-stat']['Resource'])

    if len(resource_list) == 0:
        err = "\n\n[!] No derivatives selected!\n\n"
        raise Exception(err)

    # remove any extra /'s
    pipeline_output_folder = pipeline_output_folder.rstrip("/")

    print "\n\nGathering the output file paths from %s..." \
          % pipeline_output_folder

    # this is just to keep the fsl feat config file derivative_list entries
    # nice and lean
    dirs_to_grab = []
    for derivative_name in derivative_list:
        for resource_name in resource_list:
            if resource_name in derivative_name:
                dirs_to_grab.append(derivative_name)

    # grab MeanFD_Jenkinson just in case
    dirs_to_grab.append("power_params")

    for resource_name in dirs_to_grab:
        glob_string = os.path.join(pipeline_output_folder, "*",
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
              "output directory provided: %s\nDerivatives selected:%s\n\n" \
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
    if (values_list[0].replace(" ", "") not in power_params_filepath) or \
        (values_list[1].replace(" ", "") not in power_params_filepath):
        err = "\n\n[!] There is a mismatch between the contents of the " \
              "power parameters file and where it is located!\n" \
              "Filepath: %s\n\n" % power_params_filepath
        raise Exception(err)

    if (names_list[2].replace(" ", "") != "MeanFD_Power") or \
        (names_list[3].replace(" ", "") != "MeanFD_Jenkinson") or \
            (names_list[-1].replace(" ", "") != "MeanDVARS"):
        err = "\n\n[!] There is a mismatch between the power parameters " \
              "format and what is expected!!\nFilepath: %s\n\n" \
              % power_params_filepath
        raise Exception(err)

    meanfd_power = values_list[2]
    meanfd_jenk = values_list[3]
    meandvars = values_list[-1]

    return meanfd_power, meanfd_jenk, meandvars
 

def create_output_dict_list(nifti_globs, pipeline_output_folder,
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
            
            new_row_dict["participant_session_id"] = unique_id
            new_row_dict["participant_id"] = unique_id.split('_')[0]
            new_row_dict["session_id"] = unique_id.split('_')[1]
            new_row_dict["Series"] = series_id
                                   
            new_row_dict["Filepath"] = filepath
                        
            if get_motion:
                # if we're including motion measures
                power_params_file = find_power_params_file(filepath,
                    resource_id, series_id)
                power_params_lines = load_text_file(power_params_file,
                    "power parameters file")
                meanfd_p, meanfd_j, meandvars = \
                    extract_power_params(power_params_lines,
                                         power_params_file)
                new_row_dict["MeanFD_Power"] = meanfd_p
                new_row_dict["MeanFD_Jenkinson"] = meanfd_j
                new_row_dict["MeanDVARS"] = meandvars

            if get_raw_score:
                # grab raw score for measure mean just in case
                raw_score_path = grab_raw_score_filepath(filepath,
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
            new_df = new_df[new_df.participant_session_id.isin(inclusion_list)]

        if new_df.empty:
            raise Exception("the group analysis participant list you used "
                            "resulted in no outputs")
                   
        # unique_resource_id is tuple (resource_id,strat_info)
        if unique_resource_id not in output_df_dict.keys():
            output_df_dict[unique_resource_id] = new_df
            
    return output_df_dict


def gather_outputs(pipeline_folder, resource_list, inclusion_list,
                       get_motion, get_raw_score):

    # probably won't have a session list due to subject ID format!

    nifti_globs = gather_nifti_globs(pipeline_folder, resource_list)
    output_dict_list = create_output_dict_list(nifti_globs, pipeline_folder,
                           get_motion, get_raw_score)
    output_df_dict = create_output_df_dict(output_dict_list, inclusion_list)

    return output_df_dict


def pheno_sessions_to_repeated_measures(pheno_df, sessions_list):
    """Take in the selected session names, and match them to the unique
    participant-session IDs appropriately for an FSL FEAT repeated measures
    analysis.

    More info:
      https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/
          UserGuide#Paired_Two-Group_Difference_.28Two-Sample_Paired_T-Test.29

    Sample input:
        pheno_df
          sub01_ses01
          sub01_ses02
          sub02_ses01
          sub02_ses02
        sessions_list
          [ses01, ses02]

    Expected output:
        pheno_df       Session  participant_sub01  participant_sub02
          sub01_ses01    ses01                  1                  0
          sub02_ses01    ses01                  0                  1
          sub01_ses02    ses02                  1                  0
          sub02_ses02    ses02                  0                  1
    """

    # first, check to see if this design matrix setup has already been done
    # in the pheno CSV file
    #     NOTE: this is mainly for PRESET GROUP ANALYSIS MODELS!!!
    num_partic_cols = 0
    for col_names in pheno_df.columns:
        if "participant" in col_names:
            num_partic_cols += 1
    if num_partic_cols > 1 and ("session" in pheno_df.columns or "session_column_one" in pheno_df.columns):
        for part_ses_id in pheno_df["participant_session_id"]:
            if "participant_{0}".format(part_ses_id.split("_")[0]) in pheno_df.columns:
                continue
            break
        else:
            # if it's already set up properly, then just send the pheno_df
            # back and bypass all the machinery below
            return pheno_df

    # there are no new rows, since the phenotype file will have all of the
    # subject_site_session combo unique IDs on each row!!!
    sessions_col = []
    part_ids_col = []

    # participant IDs new columns
    participant_id_cols = {}
    i = 0
    for participant_unique_id in list(pheno_df["participant_session_id"]):
        part_col = [0] * len(pheno_df["participant_session_id"])
        for session in sessions_list:
            if session in participant_unique_id.split("_")[1]:
                # generate/update sessions categorical column
                part_id = participant_unique_id.split("_")[0]
                part_ids_col.append(part_id)
                sessions_col.append(session)
                header_title = "participant_%s" % part_id
                # generate/update participant ID column (1's or 0's)
                if header_title not in participant_id_cols.keys():
                    part_col[i] = 1
                    participant_id_cols[header_title] = part_col
                else:
                    participant_id_cols[header_title][i] = 1
        i += 1

    pheno_df["Session"] = sessions_col
    pheno_df["participant"] = part_ids_col

    # add new participant ID columns
    for new_col in participant_id_cols.keys():
        pheno_df[new_col] = participant_id_cols[new_col]

    return pheno_df


def pheno_series_to_repeated_measures(pheno_df, series_list,
                                      repeated_sessions=False):

    # take in the selected series/scans, and create all of the permutations
    # of unique participant IDs (participant_site_session) and series/scans
    # and populate the pheno
    #   this is so the user does not have to have a specially-formatted
    #   version of the phenotype CSV for repeated measures; they can just
    #   enter the regular one

    import pandas as pd

    # first, check to see if this design matrix setup has already been done
    # in the pheno CSV file
    num_partic_cols = 0
    for col_names in pheno_df.columns:
        if "participant" in col_names:
            num_partic_cols += 1
    if num_partic_cols > 1 and "scan" in pheno_df.columns:
        for part_ses_id in pheno_df["participant_session_id"]:
            if "participant_{0}".format(part_ses_id.split("_")[0]) in pheno_df.columns:
                continue
            break
        else:
            # if it's already set up properly, then just send the pheno_df
            # back and bypass all the machinery below
            return pheno_df

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

        for participant_unique_id in pheno_df["participant_session_id"]:

            part_col = [0] * len(pheno_df["participant_session_id"])
            header_title = "participant_%s" % participant_unique_id

            if header_title not in participant_id_cols.keys():
                part_col[i] = 1
                participant_id_cols[header_title] = part_col
            else:
                participant_id_cols[header_title][i] = 1

            i += 1

        for new_col in participant_id_cols.keys():
            pheno_df[new_col] = participant_id_cols[new_col]
        
    return pheno_df


def balance_repeated_measures(pheno_df, sessions_list, series_list=None):

    # this is for repeated measures only.
    # if the user selects a participant list like this:
    #    sub001_session_1
    #    sub001_session_2
    #    sub002_session_1
    #    sub002_session_2
    #    sub003_session_1
    # then have this drop "sub003_session_1", because repeated measures
    # requires a uniform balance of repeats

    from collections import Counter

    part_ID_count = Counter(pheno_df["participant_ID"])

    if series_list:
        sessions_x_series = len(sessions_list) * len(series_list)
    else:
        sessions_x_series = len(sessions_list)

    dropped_parts = []

    for part_ID in part_ID_count.keys():
        if part_ID_count[part_ID] != sessions_x_series:
            pheno_df = pheno_df[pheno_df.participant != part_ID]
            del pheno_df["participant_%s" % part_ID]
            dropped_parts.append(part_ID)

    return pheno_df, dropped_parts


def prep_analysis_df_dict(config_file, pipeline_output_folder):
    
    # Preps group analysis run

    # config_file: filepath to the main CPAC pipeline configuration YAML
    #              (not the group analysis config YAML)
    # pipeline_output_folder: filepath to the CPAC pipeline individual-level
    #                   analysis output directory
    #                   example:
    #                     /home/cpac_run_1/output/pipeline_040_ANTS

    import pandas as pd
    import pkg_resources as p

    keys_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')
    try:
        keys = pd.read_csv(keys_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(keys_csv, e)
        raise Exception(err)

    derivatives = list(keys[keys['Derivative'] == 'yes'][keys['Space'] == 'template'][keys['Values'] == 'z-score']['Resource'])
    derivatives = derivatives + list(keys[keys['Derivative'] == 'yes'][keys['Space'] == 'template'][keys['Values'] == 'z-stat']['Resource'])

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
        group_models.append((group_config_file,
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

        inclusion = load_text_file(group_model.participant_list,
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
    output_df_dict = gather_outputs(pipeline_output_folder,
                                    full_output_measure_list,
                                    full_inclusion_list,
                                    get_motion,
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

        if len(group_model.derivative_list) == 0:
            err = "\n\n[!] There are no derivatives listed in the " \
                  "derivative_list field of your group analysis " \
                  "configuration file.\n\nConfiguration file: " \
                  "{0}\n".format(group_config_file)
            raise Exception(err)

        # removing this due to the recent change
        '''
        for deriv_name in group_model.derivative_list:
            if deriv_name not in derivatives:
                err = "\n\n[!] One of the derivative names you provided " \
                      "({0}) in the derivative_list field in your group " \
                      "analysis configuration file is not a valid CPAC " \
                      "output name.\n\nConfiguration file: {1}" \
                      "\n".format(deriv_name, group_config_file)
                raise Exception(err)
        '''

        # load original phenotype CSV into a dataframe
        pheno_df = read_pheno_csv_into_df(group_model.pheno_file)

        # enforce the sub ID label to "Participant"
        pheno_df.rename(columns={group_model.participant_id_label:"participant_id"},
                        inplace=True)   
        pheno_df["participant_id"] = pheno_df["participant_id"].astype(str)

        # unique_resource = (output_measure_type, preprocessing strategy)
        # output_df_dict[unique_resource] = dataframe
        for unique_resource in output_df_dict.keys():

            resource_id = unique_resource[0]

            # do this backwards, because the group_model.derivative_list is a
            # list of substrings that would be in a derivative name
            # for example:
            #     group_model.derivative_list = ['centrality']
            #     this would include both 'centrality_zstd' and
            #     'centrality_smooth_zstd', both of which could be the current
            #     value of 'resource_id'
            # also, 'derivatives' is a list of group-analysis eligible
            # derivatives (standard space, z-score standardized)
            for derivative in group_model.derivative_list:
                if derivative in resource_id and resource_id in derivatives:
                    break
            else:
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
            inclusion_list = load_text_file(group_model.participant_list,
                                            "group-level analysis participant list")
            output_df = \
                output_df[output_df["participant_session_id"].isin(inclusion_list)]

            new_pheno_df = pheno_df.copy()

            # check for inconsistency with leading zeroes
            # (sometimes, the sub_ids from individual will be something like
            #  '0002601' and the phenotype will have '2601')
            sublist_subs = output_df['participant_id']
            pheno_subs = list(new_pheno_df['participant_id'])
            for sub in sublist_subs:
                if sub in pheno_subs:
                    # okay, there's at least one match
                    break
            else:
                new_sublist_subs = [str(x).lstrip('0') for x in sublist_subs]
                for sub in new_sublist_subs:
                    if sub in pheno_subs:
                        # that's better
                        output_df['participant_id'] = new_sublist_subs
                        break
                else:
                    raise Exception('the participant IDs in your group '
                                    'analysis participant list and the '
                                    'participant IDs in your phenotype file '
                                    'do not match')

            repeated_measures = False
            repeated_sessions = False
            repeated_series = False

            if len(group_model.sessions_list) > 0:
                repeated_sessions = True

            if len(group_model.series_list) > 0:
                repeated_series = True

            if repeated_sessions or repeated_series:
                repeated_measures = True

            if repeated_measures:
                if repeated_sessions:
                    # IF USING FSL PRESETS: new_pheno_df will get passed
                    #                       through unchanged
                    new_pheno_df = \
                        pheno_sessions_to_repeated_measures(new_pheno_df,
                                                            group_model.sessions_list)

                # create new rows for all of the series, if applicable
                #   ex. if 10 subjects and two sessions, 10 rows -> 20 rows
                if repeated_series:
                    # IF USING FSL PRESETS: new_pheno_df will get passed
                    #                       through unchanged
                    new_pheno_df = \
                        pheno_series_to_repeated_measures(new_pheno_df,
                                                          group_model.series_list,
                                                          repeated_sessions)

                # drop the pheno rows - if there are participants missing in
                # the output files (ex. if ReHo did not complete for 2 of the
                # participants, etc.), then drop these rows from the phenotype
                #   we are dropping all instances of this participant, all
                #   sessions and all series, because in repeated measures/
                #   within-subject, if one goes, they all have to go    
                new_pheno_df = \
                    new_pheno_df[pheno_df["participant_id"].isin(output_df["participant_id"])]

                if len(new_pheno_df) == 0:
                    err = "\n\n[!] There is a mis-match between the "\
                          "participant IDs in the output directory/particip" \
                          "ant list and the phenotype file.\n\n"
                    raise Exception(err)

                join_columns = ["participant_id"]

                if "scan" in new_pheno_df:
                    # TODO: maybe come up with something more unique than
                    # TODO: "session" or "scan" for covariate names to signal
                    # TODO: when presets are being used?
                    # if we're using one of the FSL presets!
                    # IMPT: we need to match the rows with the actual scans

                    # ALSO IMPT: we're going to rely on the series_list from
                    #            the group model config to match, so always
                    #            make sure the order remains the same
                    #            example: the 1,1,1,-1,-1,-1 condition vector
                    #                     in the preset should be the first
                    #                     scan in the list for 1,1,1 and the
                    #                     second for -1,-1,-1
                    scan_label_col = []
                    for val in new_pheno_df["scan"]:
                        if len(group_model.series_list) == 2:
                            if val == 1:
                                scan_label_col.append(
                                    group_model.series_list[0])
                            elif val == -1:
                                scan_label_col.append(
                                    group_model.series_list[1])
                    new_pheno_df["Series"] = scan_label_col

                    # now make sure the 1,1,1,-1,-1,-1,...etc. matches
                    # properly with the actual scans by merging
                    join_columns.append("Series")
                    new_pheno_df = pd.merge(new_pheno_df, output_df,
                                            how="inner", on=join_columns)
                    run_label = "repeated_measures_multiple_series"

                    analysis_dict[(model_name, group_config_file, resource_id, strat_info, run_label)] = \
                        new_pheno_df

                elif "Series" in new_pheno_df:
                    # if Series is one of the categorically-encoded covariates
                    # make sure we only are including the series the user has
                    # selected to include in the repeated measures analysis

                    # check in case the pheno has series IDs that doesn't
                    # exist in the output directory, first
                    new_pheno_df = \
                        new_pheno_df[new_pheno_df["Series"].isin(output_df["Series"])]
                    # okay, now check against the user-specified series list
                    new_pheno_df = \
                        new_pheno_df[new_pheno_df["Series"].isin(group_model.series_list)]
                    join_columns.append("Series")
                    # pull together the pheno DF and the output files DF!
                    new_pheno_df = pd.merge(new_pheno_df, output_df,
                                            how="inner", on=join_columns)

                    if repeated_sessions:
                        # this can be removed/modified once sessions are no
                        # longer integrated in the full unique participant IDs
                        new_pheno_df, dropped_parts = \
                            balance_repeated_measures(new_pheno_df,
                                                      group_model.sessions_list,
                                                      group_model.series_list)

                        run_label = "repeated_measures_multiple_sessions_and_series"
                    else:
                        run_label = "repeated_measures_multiple_series"

                    analysis_dict[(model_name, group_config_file, resource_id, strat_info, run_label)] = \
                        new_pheno_df

                else:
                    # this runs if there are repeated sessions but not
                    # repeated series
                    #   split up the series here
                    #   iterate over the Series/Scans
                    for series_df_tuple in output_df.groupby("Series"):

                        series = series_df_tuple[0]

                        # series_df is output_df but with only one of the
                        # Series
                        series_df = series_df_tuple[1]

                        # TODO: is this a mistake?
                        # trim down the pheno DF to match the output DF and
                        # merge
                        newer_pheno_df = new_pheno_df[pheno_df["participant_id"].isin(series_df["participant_id"])]
                        newer_pheno_df = pd.merge(new_pheno_df, series_df, how="inner", on=["participant_id"])

                        # this can be removed/modified once sessions are no
                        # longer integrated in the full unique participant IDs
                        if "Session" in newer_pheno_df.columns:
                            # TODO: re-visit why there is a "participant_ID"
                            # TODO: column? will this still work without
                            # TODO: presets?
                            newer_pheno_df, dropped_parts = \
                                balance_repeated_measures(newer_pheno_df,
                                                          group_model.sessions_list,
                                                          None)

                        # unique_resource =
                        #        (output_measure_type, preprocessing strategy)
                        analysis_dict[(model_name, group_config_file,
                                       resource_id, strat_info,
                                       "repeated_measures_%s" % series)] = newer_pheno_df

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
                    newer_pheno_df = \
                        new_pheno_df[pheno_df["participant_id"].isin(series_df["participant_id"])]
                    newer_pheno_df = pd.merge(new_pheno_df, series_df,
                                              how="inner",
                                              on=["participant_id"])

                    # send it in
                    analysis_dict[(model_name, group_config_file, resource_id,
                                   strat_info, series)] = newer_pheno_df

    return analysis_dict


def run_cwas_group(output_dir, working_dir, roi_file,
                   regressor_file, participant_column, columns,
                   permutations, parallel_nodes, inclusion=None):

    import os
    import numpy as np
    from multiprocessing import pool
    from CPAC.cwas.pipeline import create_cwas

    output_dir = os.path.abspath(output_dir)
    working_dir = os.path.join(working_dir, 'cpac_group_analysis', 'MDMR',
                               os.path.basename(output_dir))

    inclusion_list = None
    if inclusion:
        inclusion_list = load_text_file(inclusion, "MDMR participant "
                                                   "inclusion list")

    output_df_dct = gather_outputs(output_dir,
                                   ["functional_to_standard"],
                                   inclusion_list, False, False)

    for preproc_strat in output_df_dct.keys():
        # go over each preprocessing strategy

        df_dct = {}
        strat_df = output_df_dct[preproc_strat]


        if len(set(strat_df["Series"])) > 1:
            # more than one scan/series ID
            for strat_scan in list(set(strat_df["Series"])):
                # make a list of sub-dataframes, each one with only file paths
                # from one scan ID each
                df_dct[strat_scan] = strat_df[strat_df["Series"] == strat_scan]
        else:
            df_dct[list(set(strat_df["Series"]))[0]] = strat_df

        for df_scan in df_dct.keys():
            func_paths = {
                p.split("_")[0]: f
                for p, f in
                zip(
                    df_dct[df_scan].Participant,
                    df_dct[df_scan].Filepath
                )
            }

            cwas_wf = create_cwas(name="MDMR_{0}".format(df_scan))
            cwas_wf.inputs.inputspec.subjects = func_paths
            cwas_wf.inputs.inputspec.roi = roi_file
            cwas_wf.inputs.inputspec.regressor = regressor_file
            cwas_wf.inputs.inputspec.participant_column = participant_column
            cwas_wf.inputs.inputspec.columns = columns
            cwas_wf.inputs.inputspec.permutations = permutations
            cwas_wf.inputs.inputspec.parallel_nodes = parallel_nodes
            cwas_wf.run()


def run_cwas(pipeline_config):

    import os
    import yaml

    pipeline_config = os.path.abspath(pipeline_config)

    with open(pipeline_config, "r") as f:
        pipeconfig_dct = yaml.load(f)

    output_dir = pipeconfig_dct["outputDirectory"]
    working_dir = pipeconfig_dct["workingDirectory"]

    roi_file = pipeconfig_dct["mdmr_roi_file"]
    regressor_file = pipeconfig_dct["mdmr_regressor_file"]
    participant_column = pipeconfig_dct["mdmr_regressor_participant_column"]
    columns = pipeconfig_dct["mdmr_regressor_columns"]
    permutations = pipeconfig_dct["mdmr_permutations"]
    parallel_nodes = pipeconfig_dct["mdmr_parallel_nodes"]
    inclusion = pipeconfig_dct["mdmr_inclusion"]

    if not inclusion or "None" in inclusion or "none" in inclusion:
        inclusion = None

    run_cwas_group(output_dir, working_dir, roi_file,
                   regressor_file, participant_column, columns,
                   permutations, parallel_nodes,
                   inclusion=inclusion)


def find_other_res_template(template_path, new_resolution):
    """Find the same template/standard file in another resolution, if it
    exists.

    template_path: file path to the template NIfTI file
    new_resolution: (int) the resolution of the template file you need

    NOTE: Makes an assumption regarding the filename format of the files.
    """

    # TODO: this is assuming there is a mm resolution in the file path - not
    # TODO: robust to varying templates - look into alternatives

    ref_file = None

    if "mm" in template_path:
        template_parts = template_path.rsplit('mm', 1)

        if len(template_parts) < 2:
            # TODO: better message
            raise Exception('no resolution in the file path!')

        template_parts[0] = str(new_resolution).join(template_parts[0].rsplit(template_parts[0][-1], 1))
        ref_file = "{0}{1}".format(template_parts[0], template_parts[1])

    elif "${resolution_for_func_preproc}" in template_path:
        ref_file = template_path.replace("${resolution_for_func_preproc}",
                                         "{0}mm".format(new_resolution))

    if ref_file:
        print("\n{0}mm version of the template found:\n{1}"
              "\n\n".format(new_resolution, ref_file))

    return ref_file


def check_cpac_output_image(image_path, reference_path, out_dir=None,
                            roi_file=False):

    import os
    import nibabel as nb

    if not out_dir:
        out_dir = os.getcwd()

    # we want to preserve the original directory structure of the input image,
    # but place that sub-tree into the BASC working directory (in this case,
    # 'out_dir')
    try:
        orig_dir = "pipeline_{0}".format(image_path.split('pipeline_')[1])
    except IndexError:
        if roi_file:
            orig_dir = os.path.join("ROI_files", os.path.basename(image_path))
        else:
            raise IndexError(image_path)
    out_path = os.path.join(out_dir, 'resampled_input_images', orig_dir)

    # if this was already done
    if os.path.isfile(out_path):
        image_path = out_path

    if not os.path.isdir(out_path.replace(os.path.basename(out_path), "")):
        try:
            os.makedirs(out_path.replace(os.path.basename(out_path), ""))
        except:
            # TODO: better message
            raise Exception("couldn't make the dirs!")

    resample = False

    image_nb = nb.load(image_path)
    ref_nb = nb.load(reference_path)

    # check: do we even need to resample?
    if int(image_nb.header.get_zooms()[0]) != int(ref_nb.header.get_zooms()[0]):
        print("Input image resolution is {0}mm\nTemplate image resolution "
              "is {1}mm\n".format(image_nb.header.get_zooms()[0],
                                  ref_nb.header.get_zooms()[0]))
        resample = True
    if image_nb.shape != ref_nb.shape:
        print("Input image shape is {0}\nTemplate image shape is "
              "{1}\n".format(image_nb.shape, ref_nb.shape))
        resample = True

    if resample:
        print("Resampling input image:\n{0}\n\n..to this reference:\n{1}"
              "\n\n..and writing this file here:\n{2}"
              "\n".format(image_path, reference_path, out_path))
        cmd = ['flirt', '-in', image_path, '-ref', reference_path, '-out',
               out_path]
        return cmd
    else:
        return resample


def resample_cpac_output_image(cmd_args):

    import subprocess

    print("Running:\n{0}\n\n".format(" ".join(cmd_args)))
    retcode = subprocess.check_output(cmd_args)

    return cmd_args[-1]


def run_basc_group(pipeline_dir, roi_file, roi_file_two, ref_file,
                   num_ts_bootstraps, num_ds_bootstraps, num_clusters,
                   affinity_thresh, cross_cluster, proc, memory, out_dir,
                   output_size=800, inclusion=None, scan_inclusion=None):

    import os
    import numpy as np
    from multiprocessing import pool
    from CPAC.basc.basc_workflow_runner import run_basc_workflow

    pipeline_dir = os.path.abspath(pipeline_dir)

    # TODO: this must change once PyBASC is modified (if it is) to have a
    # TODO: separate working and output directory
    out_dir = os.path.join(out_dir, 'cpac_group_analysis', 'PyBASC',
                           os.path.basename(pipeline_dir))
    working_dir = out_dir

    inclusion_list = None
    if inclusion:
        inclusion_list = load_text_file(inclusion, "BASC participant "
                                                   "inclusion list")

    if scan_inclusion:
        scan_inclusion = scan_inclusion.split(',')

    # create encompassing output dataframe dictionary
    #     note, it is still limited to the lowest common denominator of all
    #     group model choices- it does not pull in the entire output directory
    # - there will be a dataframe for each combination of output measure
    #   type and preprocessing strategy
    # - each dataframe will contain output filepaths and their associated
    #   information, and each dataframe will include ALL SERIES/SCANS
    output_df_dct = gather_outputs(pipeline_dir,
                                   ["functional_to_standard",
                                    "functional_mni"],
                                   inclusion_list, False, False)

    for preproc_strat in output_df_dct.keys():
        # go over each preprocessing strategy

        df_dct = {}
        strat_df = output_df_dct[preproc_strat]

        nuisance_string = \
            preproc_strat[1].replace(os.path.basename(preproc_strat[1]), '')

        if len(set(strat_df["Series"])) > 1:
            # more than one scan/series ID
            for strat_scan in list(set(strat_df["Series"])):
                # make a list of sub-dataframes, each one with only file paths
                # from one scan ID each
                df_dct[strat_scan] = strat_df[strat_df["Series"] == strat_scan]
        else:
            df_dct[list(set(strat_df["Series"]))[0]] = strat_df

        for df_scan in df_dct.keys():

            # do only the selected scans
            if scan_inclusion:
                if df_scan not in scan_inclusion:
                    continue

            func_paths = list(df_dct[df_scan]["Filepath"])

            # affinity threshold is an iterable, and must match the number of
            # functional file paths for the MapNodes
            affinity_thresh = [affinity_thresh] * len(func_paths)

            # resampling if necessary
            #     each run should take the file, resample it and write it
            #     into the BASC sub-dir of the working directory
            #         should end up with a new "func_paths" list with all of
            #         these file paths in it
            ref_file_iterable = [ref_file] * len(func_paths)
            working_dir_iterable = [working_dir] * len(func_paths)
            func_cmd_args_list = map(check_cpac_output_image, func_paths,
                                     ref_file_iterable, working_dir_iterable)
            roi_cmd_args = check_cpac_output_image(roi_file, ref_file,
                                                   out_dir=working_dir,
                                                   roi_file=True)
            roi_two_cmd_args = check_cpac_output_image(roi_file_two, ref_file,
                                                       out_dir=working_dir,
                                                       roi_file=True)

            # resample them now
            if func_cmd_args_list[0]:
                p = pool.Pool(int(proc))
                func_paths = p.map(resample_cpac_output_image,
                                   func_cmd_args_list)

            # and the ROI file, too
            if roi_cmd_args:
                roi_file = resample_cpac_output_image(roi_cmd_args)
            if roi_two_cmd_args:
                roi_file_two = resample_cpac_output_image(roi_two_cmd_args)

            # add scan label and nuisance regression strategy label to the
            # output directory path
            out_dir = os.path.join(out_dir, df_scan,
                                   nuisance_string.lstrip('/'))
            # TODO: change once PyBASC delineates output/working
            working_dir = out_dir

            print('Starting the PyBASC workflow...\n')

            run_basc_workflow(func_paths,
                              roi_file,
                              num_ds_bootstraps,
                              num_ts_bootstraps,
                              num_clusters,
                              output_size=output_size,
                              bootstrap_list=list(np.ones(num_ds_bootstraps,
                                                          dtype=int)*num_ds_bootstraps),
                              proc_mem=[proc, memory],
                              similarity_metric="correlation",
                              cross_cluster=cross_cluster,
                              roi2_mask_file=roi_file_two,
                              blocklength=1,
                              affinity_threshold=affinity_thresh,
                              out_dir=out_dir,
                              run=True)


def run_basc(pipeline_config):

    import os
    import yaml

    pipeline_config = os.path.abspath(pipeline_config)

    with open(pipeline_config, "r") as f:
        pipeconfig_dct = yaml.load(f)

    output_dir = pipeconfig_dct["outputDirectory"]
    func_template = pipeconfig_dct["template_brain_only_for_func"]
    basc_roi = pipeconfig_dct["basc_roi_file"]
    basc_roi_two = pipeconfig_dct["basc_roi_file_two"]
    num_ts_bootstraps = pipeconfig_dct["basc_timeseries_bootstraps"]
    num_ds_bootstraps = pipeconfig_dct["basc_dataset_bootstraps"]
    num_clusters = pipeconfig_dct["basc_clusters"]
    affinity_thresh = pipeconfig_dct["basc_affinity_threshold"]
    output_size = pipeconfig_dct["basc_output_size"]
    cross_cluster = pipeconfig_dct["basc_cross_clustering"]
    basc_resolution = pipeconfig_dct["basc_resolution"]
    basc_proc = pipeconfig_dct["basc_proc"]
    basc_memory = pipeconfig_dct["basc_memory"]
    basc_inclusion = pipeconfig_dct["basc_inclusion"]
    basc_pipeline = pipeconfig_dct["basc_pipeline"]
    basc_scan_inclusion = pipeconfig_dct["basc_scan_inclusion"]

    if "None" in basc_inclusion or "none" in basc_inclusion:
        basc_inclusion = None

    if "None" in basc_pipeline or "none" in basc_pipeline:
        basc_pipeline = None
    else:
        # turn this into a list, even if there's only one pipeline folder
        # given
        basc_pipeline = basc_pipeline.split(",")

    # we have the functional template only for potential resampling - to make
    # sure everything is the same resolution and shape (as what the user has
    # selected)
    if "mm" in basc_resolution:
        basc_resolution = basc_resolution.replace('mm', '')

    # get the functional template, but in the specified resolution for BASC
    ref_file = find_other_res_template(func_template, basc_resolution)

    # did that actually work?
    if not os.path.isfile(ref_file):
        # TODO: better message
        raise Exception('not a thing')

    pipeline_dirs = []
    if not basc_pipeline:
        for dirname in os.listdir(output_dir):
            if "pipeline_" in dirname:
                pipeline_dirs.append(os.path.join(output_dir, dirname))
    else:
        for pipeline_name in basc_pipeline:
            pipeline_dirs.append(os.path.join(output_dir, pipeline_name))

    for pipeline in pipeline_dirs:
        run_basc_group(pipeline, basc_roi, basc_roi_two, ref_file,
                       num_ts_bootstraps, num_ds_bootstraps, num_clusters,
                       affinity_thresh, cross_cluster, basc_proc, basc_memory,
                       output_dir, output_size, inclusion=basc_inclusion,
                       scan_inclusion=basc_scan_inclusion)


def run_basc_quickrun(pipeline_dir, roi_file, roi_file_two=None,
                      ref_file=None, output_size=800, output_dir=None,
                      basc_proc=2, basc_memory=4, scan=None):
    """Start a quick-run of PyBASC using default values for most
    parameters."""

    num_ts_bootstraps = 100
    num_ds_bootstraps = 100
    num_clusters = 10
    affinity_thresh = 0.0

    if not ref_file:
        import os
        try:
            fsldir = os.environ['FSLDIR']
            ref_file = os.path.join(fsldir, 'data/standard',
                                    'MNI152_T1_3mm_brain.nii.gz')
        except KeyError:
            pass

    if not output_dir:
        import os
        output_dir = os.getcwd()

    run_basc_group(pipeline_dir, roi_file, roi_file_two, ref_file,
                   num_ts_bootstraps, num_ds_bootstraps, num_clusters,
                   affinity_thresh, True, basc_proc, basc_memory, output_dir,
                   output_size, scan_inclusion=scan)


def manage_processes(procss, output_dir, num_parallel=1):

    import os

    # start kicking it off
    pid = open(os.path.join(output_dir, 'pid_group.txt'), 'w')

    jobQueue = []
    if len(procss) <= num_parallel:
        """
        Stream all the subjects as sublist is
        less than or equal to the number of
        subjects that need to run
        """
        for p in procss:
            p.start()
            print >> pid, p.pid

    else:
        """
        Stream the subject workflows for preprocessing.
        At Any time in the pipeline c.numSubjectsAtOnce
        will run, unless the number remaining is less than
        the value of the parameter stated above
        """
        idx = 0
        while idx < len(procss):
            if len(jobQueue) == 0 and idx == 0:
                idc = idx
                for p in procss[idc: idc + num_parallel]:
                    p.start()
                    print >> pid, p.pid
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


def run_feat(config_file, pipeline_output_folder=None):

    from multiprocessing import Process

    # let's get the show on the road
    procss = []

    # get MAIN pipeline config loaded
    c = load_config_yml(config_file)

    if not pipeline_output_folder:
        import os
        pipeline_output_folder = os.path.join(c.outputDirectory, 
                                              'pipeline_{0}'.format(c.pipelineName))

    # create the analysis DF dictionary
    analysis_dict = prep_analysis_df_dict(config_file,
                                          pipeline_output_folder)

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

            procss.append(Process(target=prep_group_analysis_workflow,
                                  args=(model_df, config_file, model_name,
                                        group_config_file, resource_id,
                                        preproc_strat,
                                        series_or_repeated)))
        else:
            print "\n\n[!] CPAC says: Group-level analysis has not yet " \
                  "been implemented to handle runs on a cluster or " \
                  "grid.\n\nPlease turn off 'Run CPAC On A Cluster/" \
                  "Grid' in order to continue with group-level " \
                  "analysis. This will submit the job to only one " \
                  "node, however.\n\nWe will update users on when this " \
                  "feature will be available through release note " \
                  "announcements.\n\n"

    manage_processes(procss, c.outputDirectory, c.numGPAModelsAtOnce)


def run(config_file, pipeline_output_folder):

    # this runs all group analyses, and this function only really exists for
    # the "Run Group-Level Analysis" command on the GUI

    # get MAIN pipeline config loaded
    c = load_config_yml(config_file)

    # Run MDMR, if selected
    if 1 in c.runMDMR:
        run_cwas(config_file)

    # Run PyBASC, if selected
    if 1 in c.run_basc:
        run_basc(config_file)

    # Run FSL FEAT group analysis, if selected
    if 1 in c.run_fsl_feat:
        run_feat(config_file, pipeline_output_folder)
