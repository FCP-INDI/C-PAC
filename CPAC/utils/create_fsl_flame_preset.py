
# TODO: create a function that can help easily map raw pheno files that do not
# TODO: have the participant_session id that CPAC uses


def read_group_list_text_file(group_list_text_file):
    """Read in the group-level analysis participant-session list text file."""

    with open(group_list_text_file, "r") as f:
        group_list = f.readlines()

    # each item here includes both participant and session, and this also will
    # become the main ID column in the written design matrix CSV
    group_list = [str(x).rstrip("\n") for x in group_list if x != ""]

    return group_list


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


def write_group_list_text_file(group_list, out_file=None):
    """Write out the group-level analysis participant list as a text file."""

    import os

    # prevent duplicates - depending on how the design matrix is set up, we
    # might have multiples of the sub_ses_ID's, like if we're doing repeated
    # measures with series/scans
    new_group_list = []
    for sub_ses_id in group_list:
        if sub_ses_id not in new_group_list:
            new_group_list.append(sub_ses_id)

    if not out_file:
        out_file = os.path.join(os.getcwd(), "group_analysis_participant_"
                                             "list.txt")
    else:
        out_file = os.path.abspath(out_file)
        dir_path = out_file.split(os.path.basename(out_file))[0]
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    with open(out_file, "wt") as f:
        for part_id in new_group_list:
            f.write("{0}\n".format(part_id))

    if os.path.exists(out_file):
        print("Group-level analysis participant list written:" \
              "\n{0}\n".format(out_file))

    return out_file


def write_dataframe_to_csv(matrix_df, out_file=None):
    """Write out a matrix Pandas DataFrame into a CSV file."""

    import os

    if not out_file:
        out_file = os.path.join(os.getcwd(), "matrix.csv")
    else:
        out_file = os.path.abspath(out_file)
        dir_path = out_file.split(os.path.basename(out_file))[0]
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    try:
        matrix_df = matrix_df.drop(labels='participant_session_id', axis=1)
    except ValueError:
        pass
    except KeyError:
        pass

    matrix_df.to_csv(out_file, index=False)

    if os.path.exists(out_file):
        print("CSV file written:\n{0}\n".format(out_file))


def write_config_dct_to_yaml(config_dct, out_file=None):
    """Write out a configuration dictionary into a YAML file."""

    import os
    import CPAC

    if not out_file:
        out_file = os.path.join(os.getcwd(), "group_config.yml")
    else:
        out_file = os.path.abspath(out_file)
        dir_path = out_file.split(os.path.basename(out_file))[0]
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    if not out_file.endswith(".yml"):
        out_file = "{0}.yml".format(out_file)

    field_order = ['pipeline_dir', 'participant_list', 'output_dir', 'work_dir',
                   'log_dir', 'FSLDIR', 'run_fsl_feat', 'num_models_at_once',
                   'model_name', 'preset', 'pheno_file', 'ev_selections',
                   'participant_id_label', 'design_formula', 'mean_mask',
                   'custom_roi_mask', 'derivative_list', 'coding_scheme',
                   'group_sep', 'grouping_var', 'z_threshold', 'p_threshold',
                   'sessions_list', 'series_list', 'contrasts', 'f_tests',
                   'custom_contrasts', 'run_randomise', 'randomise_permutation',
                   'randomise_thresh', 'randomise_demean', 'randomise_tfce']

    with open(out_file, "wt") as f:
        f.write("# CPAC Group-Level Analysis Configuration File\n"
                "# Version {0}\n".format(CPAC.__version__))
        f.write("#\n# http://fcp-indi.github.io for more info.\n#\n"
                "# Tip: This file can be edited manually with "
                "a text editor for quick modifications.\n\n\n")
        f.write("# General Group-Level Analysis Settings\n"
                "##############################################################"
                "################\n\n")
        for key in field_order:
            val = config_dct[key]
            f.write("{0}: {1}\n\n".format(key, val))

            if key == 'FSLDIR':
                f.write("\n# FSL-FEAT\n########################################"
                        "######################################\n\n")
            if key == 'custom_contrasts':
                f.write("\n# FSL-Randomise\n###################################"
                        "###########################################\n\n")

    if os.path.exists(out_file):
        print("Group-level analysis configuration YAML file written:\n" \
              "{0}\n".format(out_file))


def create_design_matrix_df(group_list, pheno_df=None,
                            ev_selections=None, pheno_sub_label=None,
                            pheno_ses_label=None, pheno_site_label=None,
                            ses_id=False):
    """Create the design matrix intended for group-level analysis via the FSL
    FLAME tool.

    This does NOT create the final .mat file that FSL FLAME takes in. This is
    an intermediary design matrix CSV meant for the user to review.

    If there is a phenotype CSV provided, this function will align the
    participant-session ID labels in the CPAC individual-level analysis output
    directory with the values listed in the phenotype file.
    """

    import pandas as pd

    keep_cols = ['participant_id']

    if ses_id:
        # if the group_list is participant_session_id instead of participant_id
        map_df = pd.DataFrame({'participant_session_id': group_list})
        keep_cols += ['participant_session_id']
        part_ids = []
        sess_ids = []
        for part_ses in group_list:
            part = part_ses.split('_')[0]
            sess = part_ses.split('_')[1]
            part_ids.append(part)
            sess_ids.append(sess)
        map_df['participant_id'] = part_ids
        map_df['session'] = sess_ids
        map_df = map_df.sort_values(by=['session', 'participant_id'])
    else:
        map_df = pd.DataFrame({'participant_id': group_list})

    if pheno_df is None:
        # no phenotypic matrix provided; simpler design models
        design_df = map_df[keep_cols]

    else:
        # if a phenotype CSV file is provided with the data
        pheno_df = pheno_df.drop_duplicates()

        # replace spaces and dashes with underscores, to prevent confusion with
        # the Patsy design formula
        rename_pheno_cols = {}
        for col_name in pheno_df.columns:
            if ' ' in col_name or '-' in col_name:
                rename_pheno_cols.update({col_name: col_name.replace(' ', '_').replace('-', '_')})
        pheno_df = pheno_df.rename(columns=rename_pheno_cols)

        # align the pheno's participant ID column with the group sublist text
        # file
        if not pheno_sub_label:
            # TODO: exception message
            raise Exception("there's a pheno file, but no pheno sub label")
        else:
            # rename the pheno sub label thingy
            pheno_df = pheno_df.rename(
                columns={pheno_sub_label: 'participant_id'})
            if ev_selections:
                ev_selections.insert(0, 'participant_id')
            sort_by = ['participant_id']

            if pheno_ses_label:
                # if sessions are important in the model, do this also
                pheno_df = pheno_df.rename(
                    columns={pheno_ses_label: 'session_id'})
                if ev_selections:
                    ev_selections.append(1, 'session_id')
                # again, sort by session ID first in case of repeated
                # measures, where the sessions have to be all together first
                sort_by.insert(0, 'session_id')

            if pheno_site_label:
                # and if sites are important as well, same here
                pheno_df = pheno_df.rename(
                    columns={pheno_site_label: 'site_id'})
                if ev_selections:
                    ev_selections.append(2, 'site_id')

            if ev_selections:
                # get specific covariates!
                pheno_df = pheno_df[ev_selections]

            # check for inconsistency with leading zeroes
            # (sometimes, the sub_ids from individual will be something like
            #  '0002601' and the phenotype will have '2601')
            sublist_subs = map_df['participant_id']
            pheno_subs = list(pheno_df['participant_id'])

            for index, row in pheno_df.iterrows():
                pheno_sub_id = str(row['participant_id'])
                for sub_id in sublist_subs:
                    if str(sub_id).lstrip('0') == pheno_sub_id:
                        pheno_df.at[index, 'participant_id'] = sub_id

            for sub in sublist_subs:
                if sub in pheno_subs:
                    # okay, there's at least one match
                    break
            else:
                new_sublist_subs = [str(x).lstrip('0') for x in sublist_subs]
                for sub in new_sublist_subs:
                    if sub in pheno_subs:
                        # that's better
                        map_df['participant_id'] = new_sublist_subs
                        break
                else:
                    raise Exception('the participant IDs in your group '
                                    'analysis participant list and the '
                                    'participant IDs in your phenotype file '
                                    'do not match')

            # merge
            if pheno_ses_label:
                design_df = pheno_df.merge(map_df, on=['participant_id'])
            else:
                design_df = pheno_df.merge(map_df[['participant_id']],
                                           on='participant_id')

            design_df = design_df.sort_values(sort_by)


    return design_df


def create_contrasts_template_df(design_df, contrasts_dct_list=None):
    """Create the template Pandas DataFrame for the contrasts matrix CSV.

    The headers in the contrasts matrix needs to match the headers of the
    design matrix."""

    import pandas as pd

    contrast_cols = list(design_df.columns)
    contrast_cols.remove('participant_id')

    if contrasts_dct_list:
        # if we are initializing the contrasts matrix with pre-set contrast
        # vectors - just check for accuracy here
        for contrast_dct in contrasts_dct_list:
            # contrast_dct is a dictionary with each column name mapped to its
            # contrast vector value, like this:
            #     {contrast: "Group Mean", "Group Mean": 1, "age": 0}
            if (len(contrast_dct) - 1) != len(contrast_cols):
                # it's -1 because of the "contrast" column in contrast_dct
                # TODO: message
                raise Exception("number of columns in the contrast vector "
                                "does not match the number of covariate "
                                "columns in the design matrix")

    else:
        # if default, start it up with a blank "template" contrast vector
        contrast_one = {"Contrasts": "contrast_1"}
        contrast_two = {"Contrasts": "contrast_2"}

        for col in contrast_cols:
            contrast_one.update({col: 0})
            contrast_two.update({col: 0})

        contrasts_dct_list = [contrast_one, contrast_two]

    contrast_cols.insert(0, "Contrasts")

    # now, make the actual dataframe
    contrasts_df = pd.DataFrame(contrasts_dct_list)

    # order the columns properly
    contrasts_df = contrasts_df[contrast_cols]

    return contrasts_df


def preset_single_group_avg(group_list, pheno_df=None, covariate=None,
                            pheno_sub_label=None, output_dir=None,
                            model_name="one_sample_T-test"):
    """Set up the design matrix CSV for running a single group average
    (one-sample T-test)."""

    import os

    if not output_dir:
        output_dir = os.getcwd()

    id_cols = ["participant_id", "session_id", "site_id"]

    # change spaces and dashes to underscores to prevent confusion with the
    # Patsy design formula
    if covariate:
        covariate = covariate.lstrip(' ').rstrip(' ')
        covariate = covariate.replace(' ', '_').replace('-', '_')

    ev_selections = None
    if pheno_df is not None:
        if covariate and pheno_sub_label:
            # if we're adding an additional covariate
            ev_selections = [covariate]

    design_df = create_design_matrix_df(group_list, pheno_df,
                                        ev_selections=ev_selections,
                                        pheno_sub_label=pheno_sub_label)

    design_df["Group_Mean"] = 1

    group_mean_contrast = {"Contrasts": "Group Mean"}

    # make these loops in case we expand this to handle more than one
    # covariate past the Group Mean
    for col in design_df.columns:
        if col not in id_cols:
            if col == "Group_Mean":
                group_mean_contrast.update({col: 1})
            else:
                group_mean_contrast.update({col: 0})

    contrasts = [group_mean_contrast]

    if covariate:
        covariate_contrast = {"Contrasts": covariate}

        for col in design_df.columns:
            if col not in id_cols:
                if col == covariate:
                    covariate_contrast.update({col: 1})
                else:
                    covariate_contrast.update({col: 0})

        contrasts.append(covariate_contrast)

    contrasts_df = create_contrasts_template_df(design_df, contrasts)

    # create design and contrasts matrix file paths
    design_mat_path = os.path.join(output_dir, model_name,
                                   "design_matrix_{0}.csv".format(model_name))

    contrasts_mat_path = os.path.join(output_dir, model_name,
                                      "contrasts_matrix_{0}.csv"
                                      "".format(model_name))

    # start group config yaml dictionary
    design_formula = "Group_Mean"
    if covariate:
        design_formula = "{0} + {1}".format(design_formula, covariate)

    group_config = {"pheno_file": design_mat_path,
                    "ev_selections": {"demean": [str(covariate)],
                                      "categorical": ["Group_Mean"]},
                    "design_formula": design_formula,
                    "group_sep": "Off",
                    "grouping_var": None,
                    "sessions_list": [],
                    "series_list": [],
                    "custom_contrasts": contrasts_mat_path,
                    "model_name": model_name,
                    "output_dir": os.path.join(output_dir, model_name),
                    "work_dir": os.path.join(output_dir, model_name),
                    "log_dir": os.path.join(output_dir, model_name)}

    return design_df, contrasts_df, group_config


def preset_unpaired_two_group(group_list, pheno_df, groups, pheno_sub_label,
                              output_dir=None,
                              model_name="two_sample_unpaired_T-test"):
    """Set up the design matrix and contrasts matrix for running an unpaired
    two-group difference (two-sample unpaired T-test).

    group_list: a list of strings- sub_ses unique IDs
    pheno_df: a Pandas DataFrame object of the phenotypic file CSV/matrix
    groups: a list of either one or two strings- design matrix EV/covariate
            labels to take from the phenotype DF and include in the model
    pheno_sub_label: a string of the label name of the column in the phenotype
                     file that holds the participant/session ID for each row
    output_dir: (optional) string of the output directory path
    model_name: (optional) name/label of the model to run

    Sets up the model described here:
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/UserGuide
            #Unpaired_Two-Group_Difference_.28Two-Sample_Unpaired_T-Test.29

    Only one "group" will be provided usually if the two groups in the
    phenotypic information you wish to compare are encoded in one covariate
    column, as categorical information. Thus, providing this one name will
    pull it from the phenotype file, and this function will break it out into
    two columns using dummy-coding.
    """

    import os

    if not output_dir:
        output_dir = os.getcwd()

    id_cols = ["participant_id", "session_id", "site_id"]

    # change spaces and dashes to underscores to prevent confusion with the
    # Patsy design formula
    old_groups = groups
    groups = []
    for group in old_groups:
        group = group.lstrip(' ').rstrip(' ')
        group = group.replace(' ', '_').replace('-', '_')
        groups.append(group)

    # if the two groups are encoded in one categorical EV/column, then we will
    # have to dummy code them out
    #     if this is the case, then "groups" will be a list with only one
    #     element in it- the one EV/column that is to be split up into two
    ev_selections = []
    for group in groups:
        ev_selections.append(group)

    design_df = create_design_matrix_df(group_list, pheno_df,
                                        ev_selections=ev_selections,
                                        pheno_sub_label=pheno_sub_label)

    if len(groups) == 1:
        # we're going to split the one categorical EV into two
        new_groups = []

        # get full range of values in one-column categorical EV
        group_set = list(set(design_df[groups[0]]))

        # run this again!
        # change spaces and dashes to underscores to prevent confusion with the
        # Patsy design formula
        new_group_set = []
        for group in group_set:
            group = group.lstrip(' ').rstrip(' ')
            group = group.replace(' ', '_').replace('-', '_')
            new_group_set.append(group)

        # this preset is for an unpaired two-group difference- should only be
        # two groups encoded in this EV!
        if len(group_set) > 2:
            # TODO: message
            raise Exception("more than two groups provided, but this is a"
                            "model for a two-group difference\n\ngroups "
                            "found in column:\n{0}".format(str(group_set)))
        elif len(group_set) == 0:
            raise Exception("no groups were found - something went wrong "
                            "with reading the phenotype information")
        elif len(group_set) == 1:
            raise Exception("only one group found in the column provided, "
                            "but this is a model for a two-group difference"
                            "\n\ngroups found in column:\n"
                            "{0}".format(str(group_set)))

        # create the two new dummy-coded columns
        # column 1
        # new column name
        new_name = "{0}_{1}".format(groups[0], new_group_set[0])
        # create new column encoded in 0's
        design_df[new_name] = 0
        # map the relevant values into 1's
        design_df[new_name] = design_df[groups[0]].map({group_set[0]: 1,
                                                        group_set[1]: 0})
        # update groups list
        new_groups.append(new_name)

        # column 2
        # new column name
        new_name = "{0}_{1}".format(groups[0], new_group_set[1])
        # create new column encoded in 0's
        design_df[new_name] = 0
        # map the relevant values into 1's
        design_df[new_name] = design_df[groups[0]].map({group_set[1]: 1,
                                                        group_set[0]: 0})
        # update groups list
        new_groups.append(new_name)

        # drop original EV/column
        del design_df[groups[0]]

        # update groups list
        groups = new_groups

    # start the contrasts
    contrast_one = {"Contrasts": "{0} - {1}".format(groups[0], groups[1])}
    contrast_two = {"Contrasts": "{0} - {1}".format(groups[1], groups[0])}

    # make these loops in case we expand this to handle additional covariates
    # past the "prescribed" ones in the model/preset
    for col in design_df.columns:
        if col not in id_cols:
            if col == groups[0]:
                contrast_one.update({col: 1})
                contrast_two.update({col: -1})
            elif col == groups[1]:
                contrast_one.update({col: -1})
                contrast_two.update({col: 1})
            else:
                contrast_one.update({col: 0})
                contrast_two.update({col: 0})

    contrasts = [contrast_one, contrast_two]

    contrasts_df = create_contrasts_template_df(design_df, contrasts)

    # create design and contrasts matrix file paths
    design_mat_path = os.path.join(output_dir, model_name,
                                   "design_matrix_{0}.csv".format(model_name))

    contrasts_mat_path = os.path.join(output_dir, model_name,
                                      "contrasts_matrix_{0}.csv"
                                      "".format(model_name))

    # start group config yaml dictionary
    design_formula = "{0} + {1}".format(groups[0], groups[1])

    group_config = {"pheno_file": design_mat_path,
                    "ev_selections": {"demean": [],
                                      "categorical": str(groups)},
                    "design_formula": design_formula,
                    "group_sep": "On",
                    "grouping_var": str(groups),
                    "sessions_list": [],
                    "series_list": [],
                    "custom_contrasts": contrasts_mat_path,
                    "model_name": model_name,
                    "output_dir": os.path.join(output_dir, model_name),
                    "work_dir": os.path.join(output_dir, model_name),
                    "log_dir": os.path.join(output_dir, model_name)}

    return design_df, contrasts_df, group_config


def preset_paired_two_group(group_list, conditions, condition_type="session",
                            output_dir=None,
                            model_name="two_sample_unpaired_T-test"):
    """Set up the design matrix and contrasts matrix for running an paired
    two-group difference (two-sample paired T-test).

    group_list: a list of strings- sub_ses unique IDs
    conditions: a two-item list of strings- session or series/scan names of
                the two sessions or two scans (per participant) you wish to
                compare
    condition_type: a string, either "session" or "scan", depending on what
                    is in "conditions"
    output_dir: (optional) string of the output directory path
    model_name: (optional) name/label of the model to run

    Sets up the model described here:
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/UserGuide
            #Paired_Two-Group_Difference_.28Two-Sample_Paired_T-Test.29
    """

    import os

    if not output_dir:
        output_dir = os.getcwd()

    if len(conditions) != 2:
        # TODO: msg
        raise Exception

    sess_conditions = ["session", "Session", "sessions", "Sessions"]
    scan_conditions = ["scan", "scans", "series", "Series/Scans", "Series"]

    sesflag = False
    if condition_type in sess_conditions:
        sesflag = True

    design_df = create_design_matrix_df(group_list, ses_id=sesflag)

    # make the "condition" EV (the 1's and -1's delineating the two
    # conditions, with the "conditions" being the two sessions or two scans)
    condition_ev = []

    if condition_type in sess_conditions:
        # note: the participant_id column in design_df should be in order, so
        #       the condition_ev should come out in order:
        #           1,1,1,1,2,2,2,2 (if the sessions are 1 and 2)
        #       later on, this has to be converted to 1,1,1,1,-1,-1,-1,-1 !
        for sub_ses_id in design_df["participant_session_id"]:
            if sub_ses_id.split("_")[-1] == conditions[0]:
                condition_ev.append(conditions[0])
            elif sub_ses_id.split("_")[-1] == conditions[1]:
                condition_ev.append(conditions[1])

        group_config = {"sessions_list": conditions, "series_list": []}

    elif condition_type in scan_conditions:
        # TODO: re-visit later, when session/scan difference in how to run
        # TODO: group-level analysis repeated measures is streamlined and
        # TODO: simplified
        # the information needed in this part is not encoded in the group
        # sublist! user inputs the two scan names, and we have a list of
        # sub_ses (which needs to be doubled), with each scan paired to each
        # half of this list (will need to ensure these scans exist for each
        # selected derivative in the output directory later on)

        # the condition_ev should come out in order:
        #     scan-1,scan-1,scan-1,scan-2,scan-2,scan-2
        #         (if the scans are scan-1 and scan-2, etc.)
        # later on, this has to be converted to 1,1,1,-1,-1,-1 !
        for sub_ses_id in design_df["participant_id"]:
            condition_ev.append(conditions[0])
        for sub_ses_id in design_df["participant_id"]:
            condition_ev.append(conditions[1])

        # NOTE: there is only one iteration of the sub_ses list in
        #       design_df["participant_id"] at this point! so use append to
        #       double that column:
        design_df = design_df.append(design_df)

        group_config = {"sessions_list": [], "series_list": conditions}

    else:
        # TODO: msg
        raise Exception

    # let's check to make sure it came out right
    #   first half
    past_val = None
    for val in condition_ev[0:(len(condition_ev) / 2) - 1]:
        if past_val:
            if val != past_val:
                raise Exception('Non-equal amount of participants for each '
                                '{0}.\n'.format(condition_type))
        past_val = val
    #   second half
    past_val = None
    for val in condition_ev[(len(condition_ev) / 2):]:
        if past_val:
            if val != past_val:
                raise Exception('Non-equal amount of participants for each '
                                '{0}.\n'.format(condition_type))
        past_val = val

    design_df[condition_type] = condition_ev

    # initalize the contrast dct's
    contrast_one = {}
    contrast_two = {}

    design_formula = "{0}".format(condition_type)

    # create the participant identity columns
    for sub_ses_id in design_df["participant_id"]:
        new_part_col = []
        sub_id = sub_ses_id.split("_")[0]
        new_part_label = "participant_{0}".format(sub_id)
        for moving_sub_ses_id in design_df["participant_id"]:
            moving_sub_id = moving_sub_ses_id.split("_")[0]
            if moving_sub_id == sub_id:
                new_part_col.append(1)
            else:
                new_part_col.append(0)
        design_df[new_part_label] = new_part_col
        contrast_one.update({new_part_label: 0})
        contrast_two.update({new_part_label: 0})
        if new_part_label not in design_formula:
            design_formula = "{0} + {1}".format(design_formula,
                                                new_part_label)

    # finish the contrasts
    #   should be something like
    #                    ses,sub,sub,sub, etc.
    #     ses-1 - ses-2:   1,  0,  0,  0, 0...
    #     ses-2 - ses-1:  -1,  0,  0,  0, etc.
    contrast_one.update({
        "Contrasts": "{0}_{1} - {2}_{3}".format(condition_type,
                                                conditions[0],
                                                condition_type,
                                                conditions[1])})
    contrast_two.update({
        "Contrasts": "{0}_{1} - {2}_{3}".format(condition_type,
                                                conditions[1],
                                                condition_type,
                                                conditions[0])})

    contrast_one.update({condition_type: 1})
    contrast_two.update({condition_type: -1})

    try:
        design_df = design_df.drop(labels=['participant_session_id'],
                                   axis='columns')
    except KeyError:
        pass

    contrasts = [contrast_one, contrast_two]
    contrasts_df = create_contrasts_template_df(design_df, contrasts)

    # create design and contrasts matrix file paths
    design_mat_path = os.path.join(output_dir, model_name,
                                   "design_matrix_{0}.csv".format(model_name))

    contrasts_mat_path = os.path.join(output_dir, model_name,
                                      "contrasts_matrix_{0}.csv"
                                      "".format(model_name))

    # start group config yaml dictionary
    group_config.update({"pheno_file": design_mat_path,
                         "ev_selections": {"demean": [],
                                           "categorical": []},
                         "design_formula": design_formula,
                         "group_sep": "Off",
                         "grouping_var": None,
                         "custom_contrasts": contrasts_mat_path,
                         "model_name": model_name,
                         "output_dir": os.path.join(output_dir, model_name),
                         "work_dir": os.path.join(output_dir, model_name),
                         "log_dir": os.path.join(output_dir, model_name)})

    return design_df, contrasts_df, group_config


def preset_tripled_two_group(group_list, conditions, condition_type="Sessions",
                             output_dir=None,
                             model_name="tripled_T-test"):
    """Set up the design matrix and contrasts matrix for running a tripled
    two-group difference ('tripled' T-test).

    group_list: a list of strings- sub_ses unique IDs
    conditions: a three-item list of strings- session or series/scan names of
                the three sessions or three scans (per participant) you wish
                to compare
    condition_type: a string, either "session" or "scan", depending on what
                    is in "conditions"
    output_dir: (optional) string of the output directory path
    model_name: (optional) name/label of the model to run

    Sets up the model described here:
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/UserGuide
            #Tripled_Two-Group_Difference_.28.22Tripled.22_T-Test.29
    """

    import os

    if not output_dir:
        output_dir = os.getcwd()

    if len(conditions) != 3:
        # TODO: msg
        raise Exception('Three conditions are required for the tripled '
                        't-test.\n')

    sess_conditions = ["session", "Session", "sessions", "Sessions"]
    scan_conditions = ["scan", "scans", "series", "Series/Scans", "Series"]

    sesflag = False
    if condition_type in sess_conditions:
        sesflag = True

    design_df = create_design_matrix_df(group_list, ses_id=sesflag)

    # make the "condition" EVs (the 1's, -1's, and 0's delineating the three
    # conditions, with the "conditions" being the three sessions or three
    # scans)
    condition_ev_one = []
    condition_ev_two = []

    if condition_type in sess_conditions:
        # note: the participant_id column in design_df should be in order, so
        #       the condition_ev's should come out in order:
        #           1,1,1,-1,-1,-1, 0, 0, 0  (this is checked further down)
        #           1,1,1, 0, 0, 0,-1,-1,-1
        for sub_ses_id in design_df["participant_session_id"]:
            if sub_ses_id.split("_")[-1] == conditions[0]:
                condition_ev_one.append(1)
                condition_ev_two.append(1)
            elif sub_ses_id.split("_")[-1] == conditions[1]:
                condition_ev_one.append(-1)
                condition_ev_two.append(0)
            elif sub_ses_id.split("_")[-1] == conditions[2]:
                condition_ev_one.append(0)
                condition_ev_two.append(-1)

        group_config = {"sessions_list": conditions, "series_list": []}

    elif condition_type in scan_conditions:
        # TODO: re-visit later, when session/scan difference in how to run
        # TODO: group-level analysis repeated measures is streamlined and
        # TODO: simplified
        # the information needed in this part is not encoded in the group
        # sublist! user inputs the two scan names, and we have a list of
        # sub_ses (which needs to be doubled), with each scan paired to each
        # half of this list (will need to ensure these scans exist for each
        # selected derivative in the output directory later on)

        for sub_ses_id in design_df["participant_id"]:
            condition_ev_one.append(1)
            condition_ev_two.append(1)
        for sub_ses_id in design_df["participant_id"]:
            condition_ev_one.append(-1)
            condition_ev_two.append(0)
        for sub_ses_id in design_df["participant_id"]:
            condition_ev_one.append(0)
            condition_ev_two.append(-1)

        # NOTE: there is only one iteration of the sub_ses list in
        #       design_df["participant_id"] at this point! so use append
        #       (twice) triple that column:
        design_df_double = design_df.append(design_df)
        design_df = design_df_double.append(design_df)

        group_config = {"sessions_list": [], "series_list": conditions}

    else:
        # TODO: msg
        raise Exception

    # let's check to make sure it came out right
    #   first third
    for val in condition_ev_one[0:(len(condition_ev_one) / 3) - 1]:
        if val != 1:
            # TODO: msg
            raise Exception
    #   second third
    for val in condition_ev_one[(len(condition_ev_one) / 3):(len(condition_ev_one)/3)*2]:
        if val != -1:
            # TODO: msg
            raise Exception
    #   third... third
    for val in condition_ev_one[((len(condition_ev_one)/3)*2 + 1):]:
        if val != 0:
            # TODO: msg
            raise Exception
    #   first third
    for val in condition_ev_two[0:(len(condition_ev_two) / 3) - 1]:
        if val != 1:
            # TODO: msg
            raise Exception
    #   second third
    for val in condition_ev_two[(len(condition_ev_two) / 3):(len(condition_ev_two)/3)*2]:
        if val != 0:
            # TODO: msg
            raise Exception
    #   third... third
    for val in condition_ev_two[((len(condition_ev_two)/3)*2 + 1):]:
        if val != -1:
            # TODO: msg
            raise Exception

    # label the two covariate columns which encode the three conditions
    column_one = "{0}_column_one".format(condition_type)
    column_two = "{0}_column_two".format(condition_type)

    design_df[column_one] = condition_ev_one
    design_df[column_two] = condition_ev_two

    # initalize the contrast dct's
    contrast_one = {}
    contrast_two = {}
    contrast_three = {}

    design_formula = "{0} + {1}".format(column_one, column_two)

    # create the participant identity columns
    for sub_id in design_df["participant_id"]:
        new_part_col = []
        new_part_label = "participant_{0}".format(sub_id)
        for moving_sub_ses_id in design_df["participant_id"]:
            moving_sub_id = moving_sub_ses_id.split("_")[0]
            if moving_sub_id == sub_id:
                new_part_col.append(1)
            else:
                new_part_col.append(0)
        design_df[new_part_label] = new_part_col
        contrast_one.update({new_part_label: 0})
        contrast_two.update({new_part_label: 0})
        contrast_three.update({new_part_label: 0})
        if new_part_label not in design_formula:
            design_formula = "{0} + {1}".format(design_formula,
                                                new_part_label)

    # finish the contrasts
    #   should be something like
    #                    ses,ses,sub,sub,sub, etc.
    #     ses-1 - ses-2:   2,  1,  0,  0,  0...
    #     ses-1 - ses-3:   1,  2,  0,  0,  0...
    #     ses-2 - ses-3:  -1,  1,  0,  0,  0, etc.
    contrast_one.update({
        "Contrasts": "{0}_{1} - {2}_{3}".format(condition_type,
                                                conditions[0],
                                                condition_type,
                                                conditions[1])})
    contrast_two.update({
        "Contrasts": "{0}_{1} - {2}_{3}".format(condition_type,
                                                conditions[0],
                                                condition_type,
                                                conditions[2])})

    contrast_three.update({
        "Contrasts": "{0}_{1} - {2}_{3}".format(condition_type,
                                                conditions[1],
                                                condition_type,
                                                conditions[2])})

    contrast_one.update({column_one: 2, column_two: 1})
    contrast_two.update({column_one: 1, column_two: 2})
    contrast_three.update({column_one: -1, column_two: 1})

    try:
        design_df = design_df.drop(labels=['participant_session_id'],
                                   axis='columns')
    except KeyError:
        pass

    contrasts = [contrast_one, contrast_two, contrast_three]
    contrasts_df = create_contrasts_template_df(design_df, contrasts)

    # create design and contrasts matrix file paths
    design_mat_path = os.path.join(output_dir, model_name,
                                   "design_matrix_{0}.csv".format(model_name))

    contrasts_mat_path = os.path.join(output_dir, model_name,
                                      "contrasts_matrix_{0}.csv"
                                      "".format(model_name))

    # start group config yaml dictionary
    group_config.update({"pheno_file": design_mat_path,
                         "ev_selections": {"demean": [],
                                           "categorical": []},
                         "design_formula": design_formula,
                         "group_sep": "Off",
                         "grouping_var": None,
                         "custom_contrasts": contrasts_mat_path,
                         "model_name": model_name,
                         "output_dir": os.path.join(output_dir, model_name),
                         "work_dir": os.path.join(output_dir, model_name),
                         "log_dir": os.path.join(output_dir, model_name)})

    return design_df, contrasts_df, group_config


def run(pipeline_dir, derivative_list, z_thresh, p_thresh, preset=None,
        group_list_text_file=None, pheno_file=None, pheno_sub_label=None,
        output_dir=None, model_name=None, covariate=None, condition_type=None,
        run=False):

    # FSL FEAT presets: run regular group analysis with no changes to its
    # original flow- use the generated pheno as the pheno, use the
    # contrasts DF as a custom contrasts matrix, and auto-generate the
    # group analysis config YAML as well

    # NOTE: the input parameters above may come in as a dictionary instead
    #       or something

    import os
    import pandas as pd
    import pkg_resources as p

    # make life easy
    keys_csv = p.resource_filename('CPAC', 'resources/cpac_outputs.csv')
    try:
        keys = pd.read_csv(keys_csv)
    except Exception as e:
        err = "\n[!] Could not access or read the cpac_outputs.csv " \
              "resource file:\n{0}\n\nError details {1}\n".format(keys_csv, e)
        raise Exception(err)

    if derivative_list == 'all':
        derivative_list = ['alff', 'falff', 'reho', 'sca_roi', 'sca_tempreg',
                           'vmhc', 'centrality', 'dr_tempreg']

    if pheno_file and not pheno_sub_label:
        # TODO: message
        raise Exception("pheno file provided, but no pheno sub label")

    if pheno_sub_label and not pheno_file:
        # TODO: message
        raise Exception("pheno sub label provided, but no pheno file")

    try:
        if "None" in group_list_text_file or "none" in group_list_text_file:
            group_list_text_file = None
    except TypeError:
        pass

    if not group_list_text_file:
        from CPAC.pipeline.cpac_group_runner import grab_pipeline_dir_subs
        if (preset == "paired_two" or preset == "tripled_two") and "Sessions" in condition_type:
            group_list = grab_pipeline_dir_subs(pipeline_dir, True)
        else:
            group_list = grab_pipeline_dir_subs(pipeline_dir)
        group_list_text_file = os.path.join(output_dir, model_name,
                                            "group_participant_list_"
                                            "{0}.txt".format(model_name))

    elif isinstance(group_list_text_file, list):
        group_list = group_list_text_file

        # write out a group analysis sublist text file so that it can be
        # linked in the group analysis config yaml
        group_list_text_file = os.path.join(output_dir, model_name,
                                            "group_participant_list_"
                                            "{0}.txt".format(model_name))
    elif os.path.isfile(group_list_text_file):
        group_list = read_group_list_text_file(group_list_text_file)

        # write out a group analysis sublist text file so that it can be
        # linked in the group analysis config yaml
        group_list_text_file = os.path.join(output_dir, model_name,
                                            "group_participant_list_"
                                            "{0}.txt".format(model_name))

    if len(group_list) == 0:
        msg = "\n\n[!] C-PAC says: No participants found in the pipeline " \
              "directory you provided. Make sure the directory is the " \
              "individual-level pipeline directory that contains the sub-" \
              "directories labeled with the participant_session IDs.\n\n" \
              "Pipeline directory provided: {0}\n\n".format(pipeline_dir)
        raise Exception(msg)

    if not preset:
        # TODO: this
        pass

    group_config = {"pipeline_dir": pipeline_dir,
                    "FSLDIR": "FSLDIR",
                    "run_fsl_feat": [1],
                    "num_models_at_once": 1,
                    "preset": preset,
                    "participant_list": group_list_text_file,
                    "participant_id_label": "participant_id",
                    "mean_mask": ["Group Mask"],
                    "custom_roi_mask": None,
                    "derivative_list": derivative_list,
                    "coding_scheme": ["Treatment"],
                    "z_threshold": [float(z_thresh)],
                    "p_threshold": [float(p_thresh)],
                    "contrasts": [],
                    "f_tests": [],
                    "run_randomise": [0],
                    "randomise_permutation": 500,
                    "randomise_thresh": 5,
                    "randomise_demean": True,
                    "randomise_tfce": True}

    if preset == "single_grp":
        design_df, contrasts_df, group_config_update = \
            preset_single_group_avg(group_list, pheno_df=None, covariate=None,
                                    pheno_sub_label=None,
                                    output_dir=output_dir,
                                    model_name=model_name)

        group_config.update(group_config_update)

    elif preset == "single_grp_cov":

        if not pheno_file:
            # TODO: message
            raise Exception("pheno file not provided")

        if not covariate:
            # TODO: message
            raise Exception("covariate not provided")

        pheno_df = read_pheno_csv_into_df(pheno_file, pheno_sub_label)

        design_df, contrasts_df, group_config_update = \
            preset_single_group_avg(group_list, pheno_df, covariate=covariate,
                                    pheno_sub_label=pheno_sub_label,
                                    output_dir=output_dir,
                                    model_name=model_name)

        group_config.update(group_config_update)

    elif preset == "unpaired_two":

        if not pheno_file:
            # TODO: message
            raise Exception("pheno file not provided")

        if not covariate:
            # TODO: message
            raise Exception("the two groups were not provided")

        # we're assuming covariate will be coming in as a string of either one
        # covariate name, or a string with two covariates separated by a comma
        #     either way, it needs to be in list form in this case, not string
        covariate = covariate.split(",")

        pheno_df = read_pheno_csv_into_df(pheno_file, pheno_sub_label)

        # in this case, "covariate" gets sent in as a list of two covariates
        design_df, contrasts_df, group_config_update = \
            preset_unpaired_two_group(group_list, pheno_df,
                                      groups=covariate,
                                      pheno_sub_label=pheno_sub_label,
                                      output_dir=output_dir,
                                      model_name=model_name)

        group_config.update(group_config_update)

    elif preset == "paired_two":
        # run a two-sample paired T-test

        # we need it as repeated measures- either session or scan
        # and the list of subs
        # also: the two session or scan names (in a list together), and
        # whether they are sessions or scans

        if not covariate:
            # TODO: message
            raise Exception("the two conditions were not provided")

        if not condition_type:
            # TODO: message
            raise Exception("you didn't specify whether the two groups are "
                            "sessions or series/scans")

        # we're assuming covariate (which in this case, is the two sessions,
        # or two scans) will be coming in as a string of either one covariate
        # name, or a string with two covariates separated by a comma
        #     either way, it needs to be in list form in this case, not string
        try:
            covariate = covariate.split(",")
        except AttributeError:
            # it's already a list- keep it that way
            pass

        design_df, contrasts_df, group_config_update = \
            preset_paired_two_group(group_list,
                                    conditions=covariate,
                                    condition_type=condition_type,
                                    output_dir=output_dir,
                                    model_name=model_name)

        group_config.update(group_config_update)

    elif preset == "tripled_two":
        # run a "tripled" T-test

        # we need it as repeated measures- either session or scan
        # and the list of subs
        # also: the two session or scan names (in a list together), and
        # whether they are sessions or scans

        if not covariate:
            # TODO: message
            raise Exception("the three conditions were not provided")

        if not condition_type:
            # TODO: message
            raise Exception("you didn't specify whether the three groups are "
                            "sessions or series/scans")

        # we're assuming covariate (which in this case, is the three sessions,
        # or three scans) will be coming in as a string of either one
        # covariate name, or a string with three covariates separated by a
        # comma
        #     either way, it needs to be in list form in this case, not string
        try:
            covariate = covariate.split(",")
        except AttributeError:
            # it's already a list- keep it that way
            pass

        design_df, contrasts_df, group_config_update = \
            preset_tripled_two_group(group_list,
                                     conditions=covariate,
                                     condition_type=condition_type,
                                     output_dir=output_dir,
                                     model_name=model_name)

        group_config.update(group_config_update)

    else:
        # TODO: not a real preset!
        raise Exception("not one of the valid presets")

    # write participant list text file
    write_group_list_text_file(design_df["participant_id"],
                               group_list_text_file)

    # write design matrix CSV
    write_dataframe_to_csv(design_df, group_config["pheno_file"])

    # write custom contrasts matrix CSV
    write_dataframe_to_csv(contrasts_df, group_config["custom_contrasts"])

    # write group-level analysis config YAML
    out_config = os.path.join(output_dir, model_name,
                              "group_config_{0}.yml".format(model_name))
    write_config_dct_to_yaml(group_config, out_config)

    if run:
        # TODO: we need to separate the individual-level pipeline config from
        # TODO: the group-level one, it's too restrictive
        pass
