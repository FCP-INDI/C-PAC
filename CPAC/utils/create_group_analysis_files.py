
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


def read_pheno_csv_into_df(pheno_csv):
    """Read the phenotypic file CSV into a Pandas DataFrame."""

    import pandas as pd

    with open(pheno_csv, "r") as f:
        pheno_df = pd.read_csv(f)

    return pheno_df


def write_group_list_text_file(group_list, out_file=None):
    """Write out the group-level analysis participant list as a text file."""

    import os

    if not out_file:
        out_file = os.path.join(os.getcwd(), "group_analysis_participant_"
                                             "list.txt")
    else:
        out_file = os.path.abspath(out_file)
        dir_path = out_file.split(os.path.basename(out_file))[0]
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    with open(out_file, "wt") as f:
        for part_id in group_list:
            f.write("{0}\n".format(part_id))

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

    matrix_df.to_csv(out_file, index=False)


def write_config_dct_to_yaml(config_dct, out_file=None):
    """Write out a configuration dictionary into a YAML file."""

    import os
    import CPAC

    if not out_file:
        out_file = os.path.join(os.getcwd(), "gpa_fsl_config.yml")
    else:
        out_file = os.path.abspath(out_file)
        dir_path = out_file.split(os.path.basename(out_file))[0]
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    if not out_file.endswith(".yml"):
        out_file = "{0}.yml".format(out_file)

    field_order = ['participant_list', 'pheno_file', 'ev_selections',
                   'participant_id_label', 'design_formula', 'mean_mask',
                   'custom_roi_mask', 'derivative_list', 'coding_scheme',
                   'group_sep', 'grouping_var', 'z_threshold', 'p_threshold',
                   'sessions_list', 'series_list', 'contrasts', 'f_tests',
                   'custom_contrasts', 'model_name', 'output_dir']

    with open(out_file, "wt") as f:
        f.write("# CPAC Group-Level Analysis Configuration File\n"
                "# Version {0}\n".format(CPAC.__version__))
        f.write("#\n# http://fcp-indi.github.io for more info.\n#\n"
                "# Tip: This file can be edited manually with "
                "a text editor for quick modifications.\n\n")
        for key in field_order:
            val = config_dct[key]
            f.write("{0}: {1}\n\n".format(key, val))


def create_design_matrix_df(group_list, pheno_df=None,
                            ev_selections=None, pheno_sub_label=None,
                            pheno_ses_label=None, pheno_site_label=None):
    """Create the design matrix intended for group-level analysis via the FSL
    FLAME tool.

    This does NOT create the final .mat file that FSL FLAME takes in. This is
    an intermediary design matrix CSV meant for the user to review.

    If there is a phenotype CSV provided, this function will align the
    participant-session ID labels in the CPAC individual-level analysis output
    directory with the values listed in the phenotype file.
    """

    import pandas as pd

    # map the participant-session IDs to just participant IDs
    group_list_map = {}
    for part_ses in group_list:
        sub_id = part_ses.split("_")[0]
        ses_id = part_ses.split("_")[1]
        group_list_map[part_ses] = [sub_id, ses_id, part_ses]

    # create a dataframe mapping the 'sub01_ses-1' CPAC-style unique IDs to
    # subject and session columns, like this:
    #     sub01_ses-1    sub01    ses-1
    #     sub02_ses-1    sub02    ses-1
    map_df = pd.DataFrame.from_dict(group_list_map, orient='index')

    # also, rename the columns to be easier
    map_df = map_df.rename(
        columns={0: 'participant', 1: 'session', 2: 'participant_id'})

    # sort by sub_id and then ses_id
    map_df = map_df.sort_values(by=['participant_id'])

    # drop unique_id column (does it ever need to really be included?)
    # was just keeping it in up until here for mental book-keeping if anything
    map_df = map_df[['participant_id', 'participant', 'session']]

    if pheno_df is None:
        # no phenotypic matrix provided; simpler design models
        design_df = map_df[['participant_id']]

    else:
        # if a phenotype CSV file is provided with the data

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

            if pheno_ses_label:
                # if sessions are important in the model, do this also
                pheno_df = pheno_df.rename(
                    columns={pheno_ses_label: 'session'})
                if ev_selections:
                    ev_selections.append(1, 'session')

            if pheno_site_label:
                # and if sites are important as well, same here
                pheno_df = pheno_df.rename(
                    columns={pheno_site_label: 'site'})
                if ev_selections:
                    ev_selections.append(2, 'site')

            if ev_selections:
                # get specific covariates!
                pheno_df = pheno_df[ev_selections]

            # merge
            if pheno_ses_label:
                design_df = pheno_df.merge(map_df, on=['participant_id'])
            else:
                design_df = pheno_df.merge(map_df[['participant_id']],
                                           on='participant_id')

    return design_df


def create_contrasts_template_df(design_df, contrasts_dct_list=None):
    """Create the template Pandas DataFrame for the contrasts matrix CSV.

    The headers in the contrasts matrix needs to match the headers of the
    design matrix."""

    import pandas as pd

    contrast_cols = list(design_df.columns)
    contrast_cols.remove('participant_id')

    # TODO:
    # if session, if site, remove

    if contrasts_dct_list:
        # if we are initializing the contrasts matrix with pre-set contrast
        # vectors - just check for accuracy here
        for contrast_dct in contrasts_dct_list:
            # contrast_dct is a dictionary with each column name mapped to its
            # contrast vector value, like this:
            #     {contrast: "Group Mean", "Group Mean": 1, "age": 0}
            if (len(contrast_dct.keys()) - 1) != len(contrast_cols):
                # it's -1 because of the "contrast" column in contrast_dct
                # TODO: message
                raise Exception("number of columns in the contrast vector "
                                "does not match the number of covariate "
                                "columns in the design matrix")

    else:
        # if default, start it up with a blank "template" contrast vector
        contrast_one = {"contrasts": "contrast_1"}
        contrast_two = {"contrasts": "contrast_2"}

        for col in contrast_cols:
            contrast_one.update({col: 0})
            contrast_two.update({col: 0})

        contrasts_dct_list = [contrast_one, contrast_two]

    contrast_cols.insert(0, "contrasts")

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

    id_cols = ["participant_id", "participant", "session", "site"]

    ev_selections = None
    if pheno_df is not None:
        if covariate and pheno_sub_label:
            # if we're adding an additional covariate
            ev_selections = [covariate]

    design_df = create_design_matrix_df(group_list, pheno_df,
                                        ev_selections=ev_selections,
                                        pheno_sub_label=pheno_sub_label)

    design_df["Group_Mean"] = 1

    group_mean_contrast = {"contrasts": "Group Mean"}

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
        covariate_contrast = {"contrasts": covariate}

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
                    "ev_selections": {"demean": [covariate],
                                      "categorical": ["Group_Mean"]},
                    "design_formula": design_formula,
                    "group_sep": "Off",
                    "grouping_var": None,
                    "sessions_list": [],
                    "series_list": [],
                    "custom_contrasts": contrasts_mat_path,
                    "model_name": model_name,
                    "output_dir": output_dir}

    return design_df, contrasts_df, group_config


def preset_unpaired_two_group(group_list, pheno_df, groups, pheno_sub_label,
                              output_dir=None,
                              model_name="two_sample_unpaired_T-test"):
    """Set up the design matrix and contrasts matrix for running an unpaired
    two-group difference (two-sample unpaired T-test)."""

    import os

    if not output_dir:
        output_dir = os.getcwd()

    id_cols = ["participant_id", "participant", "session", "site"]

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

        # this preset is for an unpaired two-group difference- should only be
        # two groups encoded in this EV!
        if len(group_set) > 2:
            # TODO: message
            raise Exception("more than two groups provided, but this is a"
                            "model for a two-group difference")

        # create the two new dummy-coded columns
        # column 1
        # new column name
        new_name = "{0}-{1}".format(groups[0], group_set[0])
        # create new column encoded in 0's
        design_df[new_name] = 0
        # map the relevant values into 1's
        design_df[new_name] = design_df[groups[0]].map({group_set[0]: 1,
                                                        group_set[1]: 0})
        # update groups list
        new_groups.append(new_name)

        # column 2
        # new column name
        new_name = "{0}-{1}".format(groups[0], group_set[1])
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
    contrast_one = {"contrasts": "{0} - {1}".format(groups[0], groups[1])}
    contrast_two = {"contrasts": "{0} - {1}".format(groups[1], groups[0])}

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
                                      "categorical": groups},
                    "design_formula": design_formula,
                    "group_sep": "On",
                    "grouping_var": groups,
                    "sessions_list": [],
                    "series_list": [],
                    "custom_contrasts": contrasts_mat_path,
                    "model_name": model_name,
                    "output_dir": os.path.join(output_dir, model_name)}

    return design_df, contrasts_df, group_config


def preset_paired_two_group(group_list, conditions, condition_type="session",
                            output_dir=None,
                            model_name="two_sample_unpaired_T-test"):
    """Set up the design matrix and contrasts matrix for running an paired
    two-group difference (two-sample paired T-test)."""

    # TODO: NEXT!!!!!!!!!!!!!!!!!
    # if the conditions are delineated by sessions, then the matrix will
    # have to have both copies of the sub_ses_id in the design matrix
    #     conversely, if they are delineated by scans, then we have to have a
    #     design matrix with only one copy of the subs, then let the internal
    #     (infernal?) machinery do the doubling

    import os

    if not output_dir:
        output_dir = os.getcwd()

    # TODO: handle conditions (sessons? scans?)
    # TODO: make sure the 1, -1 vector doesn't get clobbered by Patsy
    '''
    we need, to give in a list of sub_ses. and a list of the two sessions.
        and this needs to spit out a design df that has the sub_ses doubled
        in appropriate order, with the 1 and -1. and the other columns (avoid
        letting the cpac thing process it if you can avoid it).

    but what if it's the scans instead? now it gets ugly.
        here's a list of just sub_ses (with one ses). and a list of the two
        scans, right? now what?
            you have to send it in and let cpac handle it, unfortunately.
                TODO: would it be worth it to quickly enable a custom pheno
                      input? I don't know if cpac is going to do the subs
                      columns properly, especially into the custom contrasts..
    '''

    if len(conditions) != 2:
        # TODO: msg
        raise Exception

    design_df = create_design_matrix_df(group_list)

    if condition_type == "session":
        # make the "condition" EV (the 1's and -1's delineating the two
        # conditions)
        condition_ev = []
        for sub_ses_id in design_df["participant_id"]:
            if sub_ses_id.split("_")[-1] == conditions[0]:
                condition_ev.append(1)
            elif sub_ses_id.split("_")[-1] == conditions[1]:
                condition_ev.append(-1)

        # let's check to make sure it came out right
        # first half
        for val in condition_ev[0:(len(condition_ev)/2)-1]:
            if val != 1:
                # TODO: msg
                raise Exception
        # second half
        for val in condition_ev[(len(condition_ev)/2):]:
            if val != -1:
                # TODO: msg
                raise Exception

        design_df["condition"] = condition_ev

        # start the contrasts
        contrast_one = {"contrasts": "{0} - {1}".format(groups[0], groups[1])}
        contrast_two = {"contrasts": "{0} - {1}".format(groups[1], groups[0])}

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

        contrasts_df = create_contrasts_template_df(design_df, conditions)

    elif condition_type == "scan":

    else:
        # TODO: msg
        raise Exception

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
                                      "categorical": groups},
                    "design_formula": design_formula,
                    "group_sep": "Off",
                    "grouping_var": None,
                    "sessions_list": [],
                    "series_list": [],
                    "custom_contrasts": contrasts_mat_path,
                    "model_name": model_name,
                    "output_dir": os.path.join(output_dir, model_name)}

    return design_df, contrasts_df, group_config


def run(group_list_text_file, derivative_list, z_thresh, p_thresh,
        preset=None, pheno_file=None, pheno_sub_label=None, output_dir=None,
        model_name=None, covariate=None):

    # TODO: set this up to run regular group analysis with no changes to its
    # TODO: original flow- use the generated pheno as the pheno, use the
    # TODO: contrasts DF as a custom contrasts matrix, and auto-generate the
    # TODO: group analysis config YAML as well - but factor in how to have the
    # TODO: user easily/seamlessly decide on derivatives, output dir, and
    # TODO: model name

    # TODO: up next- create some kind of CLI to run this easily from the
    # TODO: command line- test it for single grp AVG, then fix the .grp thing
    # TODO: THEN continue expanding this

    # NOTE: the input parameters above may come in as a dictionary instead
    #       or something

    import os

    if pheno_file and not pheno_sub_label:
        # TODO: message
        raise Exception("pheno file provided, but no pheno sub label")

    if pheno_sub_label and not pheno_file:
        # TODO: message
        raise Exception("pheno sub label provided, but no pheno file")

    if isinstance(group_list_text_file, list):
        group_list = group_list_text_file

        # write out a group analysis sublist text file so that it can be
        # linked in the group analysis config yaml
        out_list = os.path.join(output_dir, model_name,
                                "gpa_participant_list_"
                                "{0}.txt".format(model_name))
        group_list_text_file = write_group_list_text_file(group_list,
                                                          out_list)
    else:
        group_list = read_group_list_text_file(group_list_text_file)

    group_config = {"participant_list": group_list_text_file,
                    "participant_id_label": "participant_id",
                    "mean_mask": ["Group Mask"],
                    "custom_roi_mask": None,
                    "derivative_list": derivative_list,
                    "coding_scheme": ["Treatment"],
                    "z_threshold": [z_thresh],
                    "p_threshold": [p_thresh],
                    "contrasts": [],
                    "f_tests": []}

    if not preset:
        # TODO: this
        pass

    elif preset == "single_grp":
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

        pheno_df = read_pheno_csv_into_df(pheno_file)

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

        pheno_df = read_pheno_csv_into_df(pheno_file)

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
        # and the list of subs. that's it.

        design_df, contrasts_df, group_config_update = \
            preset_paired_two_group(group_list,
                                    conditions=covariate,
                                    output_dir=output_dir,
                                    model_name=model_name)

        pass


    else:
        # TODO: not a real preset!
        raise Exception("not one of the valid presets")

    # write design matrix CSV
    write_dataframe_to_csv(design_df, group_config["pheno_file"])

    # write custom contrasts matrix CSV
    write_dataframe_to_csv(contrasts_df, group_config["custom_contrasts"])

    # write group-level analysis config YAML
    out_config = os.path.join(output_dir, model_name,
                              "gpa_fsl_config_{0}.yml".format(model_name))
    write_config_dct_to_yaml(group_config, out_config)
