
def load_pheno_file(pheno_file):

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


def load_group_participant_list(group_participant_list_file):

    import os
    import pandas as pd


    if not os.path.isfile(group_participant_list_file):
        err = "\n\n[!] CPAC says: The group-level analysis subject list "\
              "provided does not exist!\nPath provided: %s\n\n" \
              % group_subject_list
        raise Exception(err)


    with open(group_participant_list_file,"r") as f:

        group_subs_dataframe = pd.read_csv(f)


    if "participant" not in ga_sublist.columns:
        err = "\n\n[!] CPAC says: Your group-level analysis subject "\
                "list CSV is missing a 'participant' column.\n\n"
        raise Exception(err)


    return group_subs_dataframe


def process_pheno_file(pheno_file_dataframe, group_subs_dataframe, \
                       participant_id_label):

    # drops participants from the phenotype file if they are not in the group
    # analysis participant list
    #     also handles sessions and series appropriately for repeated measures

    # input
    #   pheno_file_dataframe: Pandas dataframe of the phenotype file
    #   group_subs_dataframe: Pandas dataframe of the group analysis
    #                         participant list
    #   participant_id_label: string of the name of the participant column in
    #                         the phenotype file
    # output
    #   pheno_file_rows: a list of dictionaries, with each dictionary being
    #                    one of the rows from the phenotype file, with the
    #                    format of {header: value, header: value, ..}

    import os
    import pandas as pd


    if not isinstance(pheno_file_dataframe, pd.DataFrame):
        err = "\n\n[!] CPAC says: The phenotype information input should " \
              "be a Python Pandas dataframe object.\n\n"
        raise Exception(err)

    if not isinstance(group_subs_dataframe, pd.DataFrame):
        err = "\n\n[!] CPAC says: The group analysis participant list input "\
              "should be a Python Pandas dataframe object.\n\n"
        raise Exception(err)

    if not isinstance(participant_id_label, str):
        err = "\n\n[!] CPAC says: The participant ID label input should be " \
              "a string.\n\n"


    pheno = pheno_file_dataframe
    pheno_file_rows = []

    df_rows = []

    # convert from dataframe to list and make those strings just in case
    subjects = list(group_subs_dataframe.participant)
    subjects = [str(i) for i in subjects]

    sessions = None
    series = None

    if "session" in group_subs_dataframe.columns:
        sessions = list(group_subs_dataframe.session)
        sessions = [str(i) for i in sessions]

    if "series" in group_subs_dataframe.columns:
        series = list(group_subs_dataframe.series)
        series = [str(i) for i in series]


    # use an integer for iteration because we're not sure if there will be
    # sessions and/or series
    for i in range(0,len(subjects)):

        full_id = []

        subject = subjects[i]
        full_id.append(subject)

        if sessions and series:

            session = sessions[i]
            scan = series[i]
            full_id.append(session)
            full_id.append(scan)

            try:
                row = pheno[(pheno[participant_id_label] == subject) & \
                          (pheno.session == session) & \
                          (pheno.series == scan)]
            except:
                row = pheno[(pheno[participant_id_label] == int(subject)) & \
                          (pheno.session == session) & \
                          (pheno.series == scan)]

        elif sessions:

            session = sessions[i]
            full_id.append(session)

            try:
                row = pheno[(pheno[participant_id_label] == subject) & \
                          (pheno.session == session)]
            except:
                row = pheno[(pheno[participant_id_label] == int(subject)) & \
                          (pheno.session == session)]


        elif series:

            scan = series[i]
            full_id.append(scan)

            try:
                row = pheno[(pheno[participant_id_label] == subject) & \
                          (pheno.series == scan)]
            except:
                row = pheno[(pheno[participant_id_label] == int(subject)) & \
                          (pheno.series == scan)]


        else:

            full_id.append(subject)

            try:
                row = pheno[(pheno[participant_id_label] == subject)]
            except:
                row = pheno[(pheno[participant_id_label] == int(subject))]


        if len(row) > 1:

            err = "\n\n[!] CPAC says: Multiple phenotype entries were " \
                  "found for these criteria:\n\n%s\n\nPlease ensure " \
                  "your group analysis participant list and phenotype " \
                  "file are configured correctly.\n\n" % str(full_id)
            raise Exception(err)

        elif len(row) == 1:

            df_rows.append(row)


    new_pheno_df = pd.concat(df_rows)

    pheno_file_rows = new_pheno_df.to_dict("records")


    return pheno_file_rows


def create_pheno_dict(pheno_file_rows, ev_selections, participant_id_label):

    # creates the phenotype data dictionary in a format Patsy requires,
    # and also demeans the continuous EVs marked for demeaning

    # input
    #   pheno_file_rows: a list of dictionaries, with each dictionary being
    #                    one of the rows from the phenotype file, with the
    #                    format of {header: value, header: value, ..}
    #   ev_selections: a dictionary with keys for "categorical" and "demean",
    #                  with the entries being lists of phenotype EV names
    #   participant_id_label: string of the name of the participant column in
    #                         the phenotype file
    # output
    #   pheno_data_dict: a dictionary with each key being a phenotype column,
    #                    and each entry being a list of values, IN ORDER
    #                        data is also in Patsy-acceptable format

    import os
    import csv
    import numpy as np

    pheno_data_dict = {}

    for line in pheno_file_rows:

        for val in line.values():

            # if there are any blank values in the pheno row, skip this
            # row. if not, continue on with the "else" clause
            if val == "":
                break

        else:

            for key in line.keys():

                # if there are blank entries because of an empty row in
                # the CSV (such as ",,,,,"), move on to the next entry
                #if len(line[key]) == 0:
                #    continue

                if key not in pheno_data_dict.keys():
                    pheno_data_dict[key] = []

                # create a list within one of the dictionary values for
                # that EV if it is categorical; formats this list into a
                # form Patsy can understand regarding categoricals:
                #     example: { ADHD: ['adhd1', 'adhd1', 'adhd0'] }
                #                instead of just [1, 1, 0], etc.
                if 'categorical' in ev_selections.keys():
                    if key in ev_selections['categorical']:
                        pheno_data_dict[key].append(key + str(line[key]))

                    elif (key == subject_id_label) or (key == "session") or \
                        (key == "series"):
                        pheno_data_dict[key].append(line[key])

                    else:
                        pheno_data_dict[key].append(float(line[key]))

                elif (key == subject_id_label) or (key == "session") or \
                    (key == "series"):
                    pheno_data_dict[key].append(line[key])

                else:
                    pheno_data_dict[key].append(float(line[key]))

    # this needs to run after each list in each key has been fully
    # populated above
    for key in pheno_data_dict.keys():

        # demean the EVs marked for demeaning
        if 'demean' in ev_selections.keys():
            if key in ev_selections['demean']:

                new_demeaned_evs = []

                mean_evs = 0.0

                # populate a dictionary, a key for each demeanable EV, with
                # the value being the sum of all the values (which need to be
                # converted to float first)
                for val in pheno_data_dict[key]:
                    mean_evs += float(val)

                # calculate the mean of the current EV in this loop
                mean_evs = mean_evs / len(pheno_data_dict[key])

                # remove the EV's mean from each value of this EV
                # (demean it!)
                for val in pheno_data_dict[key]:
                    new_demeaned_evs.append(float(val) - mean_evs)

                # replace
                pheno_data_dict[key] = new_demeaned_evs

        # converts non-categorical EV lists into NumPy arrays
        # so that Patsy may read them in properly
        if 'categorical' in ev_selections.keys():
            if key not in ev_selections['categorical']:

                pheno_data_dict[key] = np.array(pheno_data_dict[key])


    return pheno_data_dict


def get_measure_dict(param_file):

    # load the CPAC-generated power parameters file and parse it

    # input
    #   param_file: a full path to the CPAC-generated power parameters CSV
    # output
    #   measure_dict: a dictionary of dictionaries in the following format
    #                     {"MeanFD_Power": {"participant_01": 15.43,
    #                                       "participant_02": 13.22},
    #                      "MeanFD_Jenkinson": {"participant_01": 18.55,
    #                                           "participant_02": 16.27},
    #                      ...}

    import os
    import pandas as pd

    if not os.path.isfile(param_file):
        err = "\n\n[!] CPAC says: You've included a motion parameter in " \
              "your group-level analysis model design formula, but " \
              "there is no motion parameters file available.\n\n"
        raise Exception(err)

    with open(param_file,"r") as f:
        motion_params = pd.read_csv(f, index_col=False)


    measures = ['MeanFD_Power', 'MeanFD_Jenkinson', 'MeanDVARS']

    measure_dict = {}


    for m in measures:

        measure_map = {}

        if m in motion_params.columns:

            part_ids = list(motion_params["Subject"])
            part_ids = [str(i) for i in part_ids]

            scan_ids = list(motion_params["Scan"])
            scan_ids = [str(i) for i in scan_ids]

            measure_vals = list(motion_params[m])
            measure_vals = [float(i) for i in measure_vals]

            for part_id, scan_id, measure_val in \
                zip(part_ids, scan_ids, measure_vals):

                measure_map[(part_id,scan_id)] = measure_val

        measure_dict[m] = measure_map


    return measure_dict


def get_custom_roi_info(roi_means_dict):

    # check
    if roi_means_dict == None:
        err_string = "\n\n[!] CPAC says: The custom ROI means were not " \
                        "calculated properly during the group analysis " \
                        "model generation.\n\n"
        raise Exception(err_string)


    roi_num = len(roi_means_dict.values()[0])

    # this will be a dictionary matching ROI regressor header labels with
    # the actual ROI dictionaries
    roi_dict_dict = {}

    # split the roi_means_dict from { subID: [mean1,mean2,mean3,..], ..}
    # to three dictionaries of { subID: mean1, .. }, { subID: mean2, .. },
    # and so on
    for num in range(0,roi_num):

        label = "Custom_ROI_Mean_%d" % int(num+1)
        temp_roi_dict = {}

        for key in roi_means_dict.keys():

            temp_roi_dict[key] = roi_means_dict[key][num-1]

        roi_dict_dict[label] = temp_roi_dict


    return roi_dict_dict


def model_group_var_separately(grouping_var, formula, pheno_data_dict, \
                                   ev_selections, coding_scheme):

    if grouping_var == None or grouping_var not in formula:
        print('\n\n[!] CPAC says: Model group variances separately is ' \
              'enabled, but the grouping variable set is either set to ' \
              'None, or was not included in the model as one of the ' \
              'EVs.\n')
        print('Design formula: ', formula)
        print('Grouping variable: ', grouping_var, '\n\n')
        raise Exception


    # do this a little early for the grouping variable so that it doesn't
    # get in the way of doing this for the other EVs once they have the
    # grouping variable in their names
    if 'categorical' in ev_selections.keys():
         for EV_name in ev_selections['categorical']:

             if EV_name == grouping_var:

                 if coding_scheme == 'Treatment':
                     formula = formula.replace(EV_name, 'C(' + EV_name + ')')
                 elif coding_scheme == 'Sum':
                     formula = formula.replace(EV_name, 'C(' + EV_name + \
                                                   ', Sum)')


    groupvar_levels = []
    grouping_var_id_dict = {}
    idx = 0

    for cat_ev_value in pheno_data_dict[grouping_var]:

        # here, each "cat_ev_value" will be one of the Patsy-format values
        # of the categorical EV that the user has selected as the grouping
        # variable, i.e. "sex1, sex1, sex0, sex1", etc..

        # cat_ev_level is the level digit or label without the EV name
        # ex. sex1 becomes 1
        cat_ev_level = str(cat_ev_value).replace(str(grouping_var), "")

        if cat_ev_level not in groupvar_levels:
            groupvar_levels.append(cat_ev_level)

        # groupvar_levels only keeps track of how many levels there are in
        # the grouping variable

        # populate this dict for creating the .grp file:
        try:
            grouping_var_id_dict[cat_ev_level].append(idx)
        except:
            grouping_var_id_dict[cat_ev_level] = [idx]

        idx += 1

    split_EVs = {}

    for key in pheno_data_dict.keys():

        # here, "key" is the name of each EV from the phenotype file, as
        # they are labeled in the phenotype file (not Patsy format)

        if (key in formula) and (key != grouping_var):

            # for the formula edit
            new_key_string = ""

            for level in groupvar_levels:

                # for the new split EV label
                groupvar_with_level = str(grouping_var) + str(level)
                new_key = key + "__" + groupvar_with_level

                # for the formula edit
                if new_key_string == "":
                    new_key_string = new_key
                else:
                    new_key_string = new_key_string + " + " + new_key

                split_EVs[new_key] = []

                # for the formula as well
                if key in ev_selections["categorical"]:
                    ev_selections["categorical"].append(new_key)

                for val, groupvar_val in zip(pheno_data_dict[key], \
                                             pheno_data_dict[grouping_var]):

                    if groupvar_with_level == groupvar_val:

                        split_EVs[new_key].append(val)

                    else:

                        split_EVs[new_key].append(0)

            del pheno_data_dict[key]
            if key in ev_selections["categorical"]:
                ev_selections["categorical"].remove(key)

            # formula edit
            formula = formula.replace(key, new_key_string)

    # put split EVs into pheno data dict
    pheno_data_dict.update(split_EVs)

    # parse through ev_selections, find the categorical names within the
    # design formula and insert C(<name>, Sum) into the design formula
    #     this is required for Patsy to process the categorical EVs
    #     properly when generating the design matrix (this goes into the
    #     .mat file)

    if 'categorical' in ev_selections.keys():
        for EV_name in ev_selections['categorical']:

            if EV_name != grouping_var:

                if coding_scheme == 'Treatment':
                    formula = formula.replace(EV_name, 'C(' + EV_name + ')')
                elif coding_scheme == 'Sum':
                    formula = formula.replace(EV_name, 'C(' + EV_name + \
                                                  ', Sum)')

    # remove intercept when modeling group variances separately
    formula = formula + " - 1"

    return pheno_data_dict, formula, grouping_var_id_dict


def check_multicollinearity(matrix):

    import numpy as np

    print("\nChecking for multicollinearity in the model..")

    U, s, V = np.linalg.svd(matrix)

    max_singular = np.max(s)
    min_singular = np.min(s)

    print("Max singular: ", max_singular)
    print("Min singular: ", min_singular)
    print("Rank: ", np.linalg.matrix_rank(matrix), "\n")

    if min_singular == 0:

        print('[!] CPAC warns: Detected multicollinearity in the ' \
                  'computed group-level analysis model. Please double-' \
                  'check your model design.\n\n')

    else:

        condition_number = float(max_singular)/float(min_singular)
        print("Condition number: %f\n\n" % condition_number)
        if condition_number > 30:

            print('[!] CPAC warns: Detected multicollinearity in the ' \
                      'computed group-level analysis model. Please double-' \
                      'check your model design.\n\n')


def write_mat_file(design_matrix, output_dir, model_name, \
                       depatsified_EV_names, current_output=None):

    import os
    import numpy as np

    dimx = None
    dimy = None

    if len(design_matrix.shape) == 1:
        dimy = 1
        dimx = design_matrix.shape[0]
    else:
        dimx, dimy = design_matrix.shape


    ppstring = '/PPheights'

    for i in range(0, dimy):

        ppstring += '\t' + '%1.5e' %(1.0)

    ppstring += '\n'


    filename = model_name + ".mat"

    out_file = os.path.join(output_dir, filename)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(out_file, 'wt') as f:

        print('/NumWaves\t%d' %dimy, file=f)
        print('/NumPoints\t%d' %dimx, file=f)
        print(ppstring, file=f)

        # print labels for the columns - mainly for double-checking your model
        col_string = '\n'

        for col in depatsified_EV_names:
            col_string = col_string + col + '\t'

        print(col_string, '\n', file=f)

        print('/Matrix', file=f)

        np.savetxt(f, design_matrix, fmt='%1.5e', delimiter='\t')


def create_grp_file(design_matrix, grouping_var_id_dict, output_dir, \
                        model_name, current_output=None):

    import os
    import numpy as np

    dimx = None
    dimy = None

    if len(design_matrix.shape) == 1:
        dimy = 1
        dimx = design_matrix.shape[0]
    else:
        dimx, dimy = design_matrix.shape

    design_matrix_ones = np.ones(dimx)

    if not (grouping_var_id_dict == None):
        i = 1
        for key in sorted(grouping_var_id_dict.keys()):

            for index in grouping_var_id_dict[key]:
                design_matrix_ones[index] = i

            i += 1

    filename = model_name + ".grp"

    out_file = os.path.join(output_dir, filename)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(out_file, "wt") as f:

        print('/NumWaves\t1', file=f)
        print('/NumPoints\t%d\n' %dimx, file=f)
        print('/Matrix', file=f)
        np.savetxt(f, design_matrix_ones, fmt='%d', delimiter='\t')


def create_design_matrix(pheno_file, ev_selections, formula, \
                             subject_id_label, sub_list=None, \
                             coding_scheme="Treatment", grouping_var=None, \
                             new_regressor_dict=None, roi_means_dict=None, \
                             output_dir=None, model_name="design", \
                             current_output=None):

    # this should allow the user to easily create a FLAMEO-formatted .mat file
    # and .grp file from the command line or from within CPAC

    # input
    #   pheno_file: full path to a CSV file with the phenotypic data
    #
    #   ev_selections: a Python dictionary of two lists denoting which EVs are
    #                  categorical, and which should be demeaned
    #                      format - {"categorical": ["dx_group", "sex"],
    #                                "demean": ["age"]}
    #
    #   formula: a string with the Patsy-format design matrix formula
    #            more info here:
    #                http://patsy.readthedocs.org/en/latest/formulas.html
    #
    #   subject_id_label: a string denoting the header label of the subject ID
    #                     column in the phenotype CSV file
    #
    #   sub_list: (optional) full path to a CSV file containing the
    #             participant IDs, and optionally session and/or series IDs
    #             that you want included from the phenotypic file
    #             NOTE: if not provided, all rows from phenotypic are included
    #                   in the model
    #
    #   coding_scheme: (optional) which encoding scheme Patsy should use when
    #                  creating the design matrix - "Treatment" (default) or
    #                  "Sum"
    #
    #   grouping_var: (optional) the grouping variable to use if modeling
    #                 group variances separately
    #
    #   new_regressor_dict: (optional) a Python dictionary containing other
    #                       dictionaries of subject IDs matched to the values
    #                       of each new regressor
    #                           format - {"MeanFD": {"sub001": 0.493,
    #                                                "sub002": 0.211,},
    #                                     "Measure_Mean": {"sub001": 0.193,
    #                                                      "sub002": 0.392},
    #                                     ..}
    #
    #   roi_means_dict: (optional) a Python dictionary of lists containing the
    #                   mean values of user-defined ROIs of the derivative for
    #                   each subject
    #                       format - {"sub001": [3.23, 2.11],
    #                                 "sub002": [1.79, 3.03]}
    #                       (with the length of the lists being the number of
    #                        ROIs specified)
    #
    #   output_dir: (optional) where to write the .mat file
    #
    #   model_name: (optional) name of the group analysis model
    #
    #   current_output: (optional) name of the derivative in the analysis
    #
    # output
    #   dmatrix: a Patsy object of the design matrix
    #   depatsified_EV_names: a list of the column names of the design matrix


    import os
    import patsy
    import numpy as np

    # if running this script alone outside of CPAC
    if output_dir == None:
        output_dir = os.getcwd()


    # let's process the phenotype file data and drop rows (participants) if
    # they are not listed in the participant list
    pheno_file_df = load_pheno_file(pheno_file)

    participant_list_df = load_group_participant_list(sub_list)

    pheno_file_rows = process_pheno_file(pheno_file_df, participant_list_df, \
                                         subject_id_label)


    # get number of subjects that have the derivative for this current model
    # (basically, the amount of time points, which must be greater than the
    # number of EVs)
    num_subjects = len(participant_list_df)


    # for repeated measures
    if "session" in participant_list_df.columns:
        ev_selections["categorical"].append("session")

    if "series" in participant_list_df.columns:
        ev_selections["categorical"].append("series")


    # start adding additionally created EVs
    if new_regressor_dict:

        for measure in new_regressor_dict.keys():

            if (measure in formula):

                measure_dict = new_regressor_dict[measure]

                for pheno_row_dict in pheno_file_rows:

                    participant_id = pheno_row_dict[subject_id_label]

                    if ("session" in pheno_row_dict.keys()) and \
                        ("series" in pheno_row_dict.keys()):
                        session_id = pheno_row_dict["session"]
                        series_id = pheno_row_dict["series"]
                        participant_tuple = \
                            (participant_id, session_id, series_id)

                    elif "session" in pheno_row_dict.keys():
                        session_id = pheno_row_dict["session"]
                        participant_tuple = (participant_id, session_id)

                    elif "series" in pheno_row_dict.keys():
                        series_id = pheno_row_dict["series"]
                        participant_tuple = (participant_id, series_id)

                    else:
                        participant_tuple = (participant_id)

                    pheno_row_dict[measure] = measure_dict[participant_tuple]


                ev_selections["demean"].append(measure)


    if "Custom_ROI_Mean" in formula:

        # include the means of the specified ROIs as regressors

        if roi_means_dict == None:
            err = "\n\n[!] You included 'Custom_ROI_Mean' in your model " \
                  "design, but there are no mean of ROI values provided." \
                  "\n\n"
            raise Exception(err)

        # roi_dict_dict is a dictionary of dictionaries, with each dictionary
        # holding all of the means for one ROI, with each entry being a mean
        # for a participant (the keys are the participant IDs)
        #     ex. {participant_01: 35.15, participant_02: 50.00}
        #             with the float values being all of the means of one of
        #             the ROIs specified
        # there will be a dictionary for each ROI specified
        roi_dict_dict = get_custom_roi_info(roi_means_dict)

        add_formula_string = ""

        for roi_column in roi_dict_dict.keys():

            roi_dict = roi_dict_dict[roi_column]

            for pheno_row_dict in pheno_file_rows:

                participant_id = pheno_row_dict[subject_id_label]

                if ("session" in pheno_row_dict.keys()) and \
                    ("series" in pheno_row_dict.keys()):
                    session_id = pheno_row_dict["session"]
                    series_id = pheno_row_dict["series"]
                    participant_tuple = \
                        (participant_id, session_id, series_id)

                elif "session" in pheno_row_dict.keys():
                    session_id = pheno_row_dict["session"]
                    participant_tuple = (participant_id, session_id)

                elif "series" in pheno_row_dict.keys():
                    series_id = pheno_row_dict["series"]
                    participant_tuple = (participant_id, series_id)

                else:
                    participant_tuple = (participant_id)

                pheno_row_dict[roi_column] = roi_dict[participant_tuple]


            ev_selections["demean"].append(roi_column)

            # create a string of all the new custom ROI regressor column names
            # to be inserted into the design formula, so that Patsy will
            # accept the phenotypic data dictionary that now has these columns

            if add_formula_string == "":
                add_formula_string = add_formula_string + roi_column
            else:
                add_formula_string = add_formula_string + " + " + roi_column


        # a regressor column of ROI means for each custom-specified ROI has
        # now been added to the model with appropriate column labels

        formula = formula.replace("Custom_ROI_Mean",add_formula_string)



    # return the data from the phenotype file processed properly for Patsy
    # and load it into 'pheno_data_dict'
    #     format: dictionary, each key is the name of an EV, and its value is
    #             a LIST of values in order of the subjects
    #                 - categorical EVs are already renamed from '0,1,..' to
    #                   'EV0,EV1,..' with EV being the EV name
    #                 - EVs to be demeaned are already demeaned
    #                 - numerical EVs (non-categorical) are in a list which
    #                   have been converted into a NumPy array
    pheno_data_dict = create_pheno_dict(pheno_file_rows, ev_selections, \
                                        subject_id_label)



    # handle modeling group variances separately (if enabled), then edit the
    # formula to be in Patsy language

    if grouping_var != None:

        pheno_data_dict, formula, grouping_var_id_dict = \
            model_group_var_separately(grouping_var, \
                                       formula, pheno_data_dict, \
                                       ev_selections, coding_scheme)

    else:

        grouping_var_id_dict = None

        if 'categorical' in ev_selections.keys():
            for EV_name in ev_selections['categorical']:

                if coding_scheme == 'Treatment':
                    formula = formula.replace(EV_name, 'C(' + EV_name + ')')
                elif coding_scheme == 'Sum':
                    formula = formula.replace(EV_name, 'C(' + EV_name + \
                                                  ', Sum)')



    # create the Patsy design matrix!

    try:

        dmatrix = patsy.dmatrix(formula, pheno_data_dict, NA_action='raise')

    except:
        print('\n\n[!] CPAC says: Design matrix creation wasn\'t ' \
                'successful - do the terms in your formula correctly ' \
                'correspond to the EVs listed in your phenotype file?\n')
        print('Phenotype file provided: ')
        print(pheno_file, '\n')
        print("Phenotypic data columns (regressors): ", list(pheno_data_dict.keys()))
        print("Formula: %s\n\n" % formula)
        raise Exception



    # check the model for multicollinearity - Patsy takes care of this, but
    # just in case

    check_multicollinearity(np.array(dmatrix))


    # prepare for final stages

    design_matrix = np.array(dmatrix, dtype=np.float16)

    column_names = dmatrix.design_info.column_names


    # check to make sure there are more time points than EVs!
    if len(column_names) >= num_subjects:
        err = "\n\n[!] CPAC says: There are more EVs than there are " \
              "subjects currently included in the model for %s. There must " \
              "be more subjects than EVs in the design.\n\nNumber of " \
              "subjects: %d\nNumber of EVs: %d\n\nNote: An 'Intercept' " \
              "column gets added to the design as an EV, so there will be " \
              "one more EV than you may have specified in your design. In " \
              "addition, if you specified to model group variances " \
              "separately, an Intercept column will not be included, but " \
              "the amount of EVs can nearly double once they are split " \
              "along the grouping variable.\n\n" \
              "If the number of subjects is lower than the number of " \
              "subjects in your group analysis subject list, this may be " \
              "because not every subject in the subject list has an output " \
              "for %s in the individual-level analysis output directory.\n\n"\
              % (current_output, num_subjects, len(column_names), \
              current_output)
        raise Exception(err)


    # remove the header formatting Patsy creates for categorical variables
    # because we are going to use depatsified_EV_names in the "Available EVs
    # for Contrasts" list on the next page, and also to test user-made custom
    # contrast files

    depatsified_EV_names = []

    for column in column_names:

        # if using Sum encoding, a column name may look like this:
        #     C(adhd, Sum)[S.adhd0]

        # this loop leaves it with only "adhd0" in this case, for the
        # contrasts list for the next GUI page

        column_string = column

        string_for_removal = ''

        for char in column_string:

            string_for_removal = string_for_removal + char

            if char == '.':
                column_string = column_string.replace(string_for_removal, '')
                string_for_removal = ''

        column_string = column_string.replace(']', '')

        depatsified_EV_names.append(column_string)


    # write the .mat file finally
    write_mat_file(design_matrix, output_dir, model_name, \
                       depatsified_EV_names, current_output)


    # write the .grp file also
    create_grp_file(design_matrix, grouping_var_id_dict, output_dir, \
                        model_name, current_output)


    # return the PATSY OBJECT of dmatrix, not the Numpy array "design_matrix"
    return dmatrix, depatsified_EV_names



def positive(dmat, a, coding, group_sep, grouping_var):

    import numpy as np

    # this is also where the "Intercept" column gets introduced into
    # the contrasts columns, for when the user uses the model builder's
    # contrast builder
    evs = dmat.design_info.column_name_indexes
    con = np.zeros(dmat.shape[1])

    if group_sep == True:

        if "__" in a and grouping_var in a:
            ev_desc = a.split("__")

            for ev in evs:
                count = 0
                for desc in ev_desc:
                    if desc in ev:
                        count += 1
                if count == len(ev_desc):
                    con[evs[ev]] = 1
                    break

            else:
                # it is a dropped term so make all other terms in that
                # category at -1
                term = a.split('[')[0]

                for ev in evs:
                    if ev.startswith(term):
                        con[evs[ev]]= -1

        elif len(a.split(grouping_var)) > 2:

            # this is if the current parsed contrast is the actual
            # grouping variable, as the Patsified name will have the
            # variable's name string in it twice

            for ev in evs:
                if a.split(".")[1] in ev:
                    con[evs[ev]] = 1
                    break
            else:
                # it is a dropped term so make all other terms in that
                # category at -1
                term = a.split('[')[0]

                for ev in evs:
                    if ev.startswith(term):
                        con[evs[ev]]= -1

    # else not modeling group variances separately
    else:

        if a in evs:
            con[evs[a]] = 1
        else:
            # it is a dropped term so make all other terms in that category
            # at -1
            term = a.split('[')[0]

            for ev in evs:
                if ev.startswith(term):
                    con[evs[ev]]= -1

        if coding == "Treatment":
            # make Intercept 0
            con[0] = 0
        elif coding == "Sum":
            # make Intercept 1
            con[1] = 1

    return con


def greater_than(dmat, a, b, coding, group_sep, grouping_var):
    c1 = positive(dmat, a, coding, group_sep, grouping_var)
    c2 = positive(dmat, b, coding, group_sep, grouping_var)
    return c1-c2


def negative(dmat, a, coding, group_sep, grouping_var):
    con = 0-positive(dmat, a, coding, group_sep, grouping_var)
    return con


def create_dummy_string(length):
    ppstring = ""
    for i in range(0, length):
        ppstring += '\t' + '%1.5e' %(1.0)
    ppstring += '\n'
    return ppstring


def create_con_file(con_dict, col_names, file_name, current_output, out_dir):

    import os

    print("col names: ")
    print(col_names)

    with open(os.path.join(out_dir, file_name) + ".con",'w+') as f:

        # write header
        num = 1

        for key in con_dict:
            f.write("/ContrastName%s\t%s\n" %(num,key))
            num += 1

        f.write("/NumWaves\t%d\n" %len(con_dict[key]))
        f.write("/NumContrasts\t%d\n" %len(con_dict))
        f.write("/PPheights%s" %create_dummy_string(len(con_dict[key])))
        f.write("/RequiredEffect%s" %create_dummy_string(len(con_dict[key])))
        f.write("\n\n")

        # print labels for the columns - mainly for double-checking your
        # model
        col_string = '\n'
        for col in col_names:
            col_string = col_string + col + '\t'
        print(col_string, '\n', file=f)

        # write data
        f.write("/Matrix\n")

        for key in con_dict:
            for v in con_dict[key]:
                f.write("%1.5e\t" %v)
            f.write("\n")


def create_fts_file(ftest_list, con_dict, model_name, current_output,
                    out_dir):

    import os
    import numpy as np

    try:
        print("\nFound f-tests in your model, writing f-tests file " \
              "(.fts)..\n")

        with open(os.path.join(out_dir, model_name + '.fts'), 'w') as f:

            print('/NumWaves\t', len(con_dict), file=f)
            print('/NumContrasts\t', len(ftest_list), file=f)

            # process each f-test
            ftst = []

            for ftest_string in ftest_list:

                ftest_vector = []

                cons_in_ftest = ftest_string.split(",")

                for con in con_dict.keys():
                    if con in cons_in_ftest:
                        ftest_vector.append(1)
                    else:
                        ftest_vector.append(0)

                ftst.append(ftest_vector)

            fts_n = np.array(ftst)

            # print labels for the columns - mainly for double-checking your
            # model
            col_string = '\n'
            for con in con_dict.keys():
                col_string = col_string + con + '\t'
            print(col_string, '\n', file=f)

            print('/Matrix', file=f)

            for i in range(fts_n.shape[0]):
                print(' '.join(fts_n[i].astype('str')), file=f)

    except Exception as e:

        filepath = os.path.join(out_dir, "model_files", current_output, \
                                    model_name + '.fts')

        errmsg = "\n\n[!] CPAC says: Could not create .fts file for " \
                    "FLAMEO or write it to disk.\nAttempted filepath: %s\n" \
                    "Error details: %s\n\n" % (filepath, e)

        raise Exception(errmsg)


def create_con_ftst_file(con_file, model_name, current_output, output_dir,
                         column_names, coding_scheme, group_sep):
    """
    Create the contrasts and fts file
    """

    import os
    import numpy as np

    with open(con_file, "r") as f:
        evs = f.readline()

    evs = evs.rstrip('\r\n').split(',')
    count_ftests = 0

    # TODO: this needs to be re-visited, but I think this was originally added
    # TODO: to counteract the fact that if someone was making a custom
    # TODO: contrasts matrix CSV, they wouldn't know that the design matrix
    # TODO: would have the Intercept added to it? but what if it wasn't?
    # TODO:     comment out for now... but test
    # remove "Contrasts" label and replace it with "Intercept"
    #evs[0] = "Intercept"

    fTest = False
    print("evs: ")
    print(evs)
    for ev in evs:
        if "f_test" in ev:
            count_ftests += 1

    if count_ftests > 0:
        fTest = True

    try:
        data = np.genfromtxt(con_file, names=True, delimiter=',', dtype=None)

    except:
        print("Error: Could not successfully read in contrast file: ",con_file)
        raise Exception

    lst = data.tolist()

    ftst = []
    fts_columns = []
    contrasts = []
    contrast_names = []

    length = None
    length = len(list(lst[0]))

    # lst = list of tuples, "tp"
    # tp = tuple in the format (contrast_name, 0, 0, 0, 0, ...)
    #      with the zeroes being the vector of contrasts for that contrast

    for tp in lst:
        contrast_names.append(tp[0])

        # create a list of integers that is the vector for the contrast
        # ex. [0, 1, 1, 0, ..]
        con_vector = list(tp)[1:(length-count_ftests)]

        fts_vector = list(tp)[(length-count_ftests):length]
        fts_columns.append(fts_vector)

        # TODO: see note about Intercept above
        # add Intercept column
        # if not group_sep:
        #     if coding_scheme == "Treatment":
        #         con_vector.insert(0, 0)
        #     elif coding_scheme == "Sum":
        #         con_vector.insert(0, 1)

        contrasts.append(con_vector)

    # contrast_names = list of the names of the contrasts (not regressors)
    # contrasts = list of lists with the contrast vectors
    num_EVs_in_con_file = len(contrasts[0])

    contrasts = np.array(contrasts, dtype=np.float16)

    fts_columns = np.array(fts_columns)

    # if there are f-tests, create the array for them
    if fTest:
        if len(contrast_names) < 2:
            errmsg = "\n\n[!] CPAC says: Not enough contrasts for running " \
                  "f-tests.\nTip: Do you have only one contrast in your " \
                  "contrasts file? f-tests require more than one contrast.\n"\
                  "Either remove the f-tests or include more contrasts.\n\n"

            raise Exception(errmsg)

        fts_n = fts_columns.T

    if len(column_names) != (num_EVs_in_con_file):
        err_string = "\n\n[!] CPAC says: The number of EVs in your model " \
                     "design matrix (found in the %s.mat file) does not " \
                     "match the number of EVs (columns) in your custom " \
                     "contrasts matrix CSV file.\n\nCustom contrasts matrix "\
                     "file: %s\n\nNumber of EVs in design matrix: %d\n" \
                     "Number of EVs in contrasts file: %d\n\nThe column " \
                     "labels in the design matrix should match those in " \
                     "your contrasts .CSV file.\nColumn labels in design " \
                     "matrix:\n%s" % (model_name, con_file, \
                     len(column_names), num_EVs_in_con_file, \
                             str(column_names))

        raise Exception(err_string)

    for design_mat_col, con_csv_col in zip(column_names, evs[1:]):
        if con_csv_col not in design_mat_col:
            errmsg = "\n\n[!] CPAC says: The names of the EVs in your " \
                     "custom contrasts .csv file do not match the names or " \
                     "order of the EVs in the design matrix. Please make " \
                     "sure these are consistent.\nDesign matrix EV columns: "\
                     "%s\nYour contrasts matrix columns: %s\n\n" \
                     % (column_names, evs[1:])

            raise Exception(errmsg)

    out_dir = os.path.join(output_dir, model_name + '.con')

    with open(out_dir,"wt") as f:
        idx = 1
        pp_str = '/PPheights'
        re_str = '/RequiredEffect'
        for name in contrast_names:

            print('/ContrastName%d' %idx, '\t', name, file=f)
            pp_str += '\t%1.5e' %(1)
            re_str += '\t%1.5e' %(1)
            idx += 1

        print('/NumWaves\t', (contrasts.shape)[1], file=f)
        print('/NumContrasts\t', (contrasts.shape)[0], file=f)
        print(pp_str, file=f)
        print(re_str + '\n', file=f)

        # print labels for the columns - mainly for double-checking your model
        col_string = '\n'
        for ev in evs:
            col_string = col_string + ev + '\t'
        print(col_string, '\n', file=f)

        print('/Matrix', file=f)
        np.savetxt(f, contrasts, fmt='%1.5e', delimiter='\t')

    if fTest:
        print("\nFound f-tests in your model, writing f-tests file (.fts)..\n")

        ftest_out_dir = os.path.join(output_dir, model_name + '.fts')

        with open(ftest_out_dir,"wt") as f:
            print('/NumWaves\t', (contrasts.shape)[0], file=f)
            print('/NumContrasts\t', count_ftests, file=f)

            # print labels for the columns - mainly for double-checking your
            # model
            col_string = '\n'
            for con in contrast_names:
                col_string = col_string + con + '\t'
            print(col_string, '\n', file=f)

            print('/Matrix', file=f)

            for i in range(fts_n.shape[0]):
                print(' '.join(fts_n[i].astype('str')), file=f)


def process_contrast(parsed_contrast, operator, ev_selections, group_sep, \
                         grouping_var, coding_scheme):

    # take the contrast strings and process them appropriately
    #     extract the two separate contrasts (if there are two), and then
    #     identify which are categorical - adapting the string if so

    parsed_EVs_in_contrast = []

    EVs_in_contrast = parsed_contrast.split(operator)

    if '' in EVs_in_contrast:
        EVs_in_contrast.remove('')


    for EV in EVs_in_contrast:

        skip = 0

        # they need to be put back into Patsy formatted header titles
        # because the dmatrix gets passed into the function that writes
        # out the contrast matrix
        if 'categorical' in ev_selections.keys():
            for cat_EV in ev_selections['categorical']:

                # second half of this if clause is in case group variances
                # are being modeled separately, and we don't want the EV
                # that is the grouping variable (which is now present in
                # other EV names) to confound this operation
                if group_sep == True:
                    gpvar = grouping_var
                else:
                    gpvar = "..."

                if (cat_EV in EV) and not (gpvar in EV and \
                    "__" in EV):

                    # handle interactions
                    if ":" in EV:
                        temp_split_EV = EV.split(":")
                        for interaction_EV in temp_split_EV:
                            if cat_EV in interaction_EV:
                                current_EV = interaction_EV
                    else:
                        current_EV = EV

                    if coding_scheme == 'Treatment':
                        cat_EV_contrast = EV.replace(EV, 'C(' + cat_EV + \
                                                         ')[T.' + current_EV+\
                                                             ']')

                    elif coding_scheme == 'Sum':
                        cat_EV_contrast = EV.replace(EV, 'C(' + cat_EV + \
                                                     ', Sum)[S.' + \
                                                     current_EV + ']')

                    parsed_EVs_in_contrast.append(cat_EV_contrast)
                    skip = 1

        if skip == 0:

            parsed_EVs_in_contrast.append(EV)

        # handle interactions
        if ":" in EV and len(parsed_EVs_in_contrast) == 2:

            parsed_EVs_in_contrast = [parsed_EVs_in_contrast[0] + ":" + \
                                         parsed_EVs_in_contrast[1]]

        if ":" in EV and len(parsed_EVs_in_contrast) == 3:

            parsed_EVs_in_contrast = [parsed_EVs_in_contrast[0], \
                                      parsed_EVs_in_contrast[1] + ":" + \
                                      parsed_EVs_in_contrast[2]]

    return parsed_EVs_in_contrast


def run(group_config, current_output, param_file=None, \
            derivative_means_dict=None, roi_means_dict=None, \
                model_out_dir=None, CPAC_run=True):

    import os
    import csv
    import numpy as np
    # open the GROUP ANALYSIS FSL .YML CONFIG FILE, not the main pipeline
    # config .yml file!
    if CPAC_run:
        c = group_config
    else:
        try:
            c = Configuration(yaml.safe_load(open(os.path.realpath(group_config), \
                                  'r')))
        except:
            raise Exception("Error in reading %s configuration file" \
                                % group_config)

    # pull in the gpa settings!
    ph = c.pheno_file

    sublist = c.subject_list

    ev_selections = c.ev_selections

    subject_id_label = c.subject_id_label

    formula = c.design_formula

    coding_scheme = c.coding_scheme[0]

    group_sep = c.group_sep

    grouping_var = c.grouping_var

    contrasts = c.contrasts

    f_tests = c.f_tests

    custom_contrasts = c.custom_contrasts

    model_name = c.model_name

    output_dir = c.output_dir

    # make sure the group analysis output directory exists
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    except:
        print('\n\n[!] CPAC says: Could not successfully create the group ' \
              'analysis output directory:\n', output_dir, '\n\nMake sure ' \
              'you have write access in this file structure.\n\n\n')
        raise Exception

    measure_dict = {}

    # extract motion measures from CPAC-generated power params file
    if param_file != None:
        measure_dict = get_measure_dict(param_file)

    # combine the motion measures dictionary with the measure_mean
    # dictionary (if it exists)
    if derivative_means_dict:
        measure_dict["Measure_Mean"] = derivative_means_dict

    # create the .mat and .grp files for FLAME
    design_matrix, regressor_names = create_design_matrix(ph, sublist, \
                                         ev_selections, formula, \
                                         subject_id_label, coding_scheme, \
                                         grouping_var, measure_dict, \
                                         roi_means_dict, model_out_dir, \
                                         model_name, current_output)

    # Create contrasts_dict dictionary for the .con file generation later
    contrasts_list = contrasts
    contrasts_dict = {}

    for contrast in contrasts_list:

        # each 'contrast' is a string the user input of the desired contrast

        # remove all spaces
        parsed_contrast = contrast.replace(' ', '')

        EVs_in_contrast = []
        parsed_EVs_in_contrast = []

        if '>' in parsed_contrast:

            parsed_EVs_in_contrast = \
                process_contrast(parsed_contrast, '>', ev_selections, \
                                 group_sep, grouping_var, coding_scheme)

            contrasts_dict[parsed_contrast] = \
                greater_than(design_matrix, parsed_EVs_in_contrast[0], \
                             parsed_EVs_in_contrast[1], coding_scheme, \
                             group_sep, grouping_var)


        elif '<' in parsed_contrast:

            parsed_EVs_in_contrast = \
                process_contrast(parsed_contrast, '<', ev_selections, \
                                 group_sep, grouping_var, coding_scheme)

            contrasts_dict[parsed_contrast] = \
                greater_than(design_matrix, parsed_EVs_in_contrast[1], \
                             parsed_EVs_in_contrast[0], coding_scheme, \
                             group_sep, grouping_var)


        else:

            contrast_string = parsed_contrast.replace('+',',+,')
            contrast_string = contrast_string.replace('-',',-,')
            contrast_items = contrast_string.split(',')

            if '' in contrast_items:
                contrast_items.remove('')

            if '+' in contrast_items and len(contrast_items) == 2:

                parsed_EVs_in_contrast = \
                    process_contrast(parsed_contrast, '+', ev_selections, \
                                     group_sep, grouping_var, coding_scheme)

                contrasts_dict[parsed_contrast] = \
                    positive(design_matrix, parsed_EVs_in_contrast[0], \
                             coding_scheme, group_sep, grouping_var)

            elif '-' in contrast_items and len(contrast_items) == 2:

                parsed_EVs_in_contrast = \
                    process_contrast(parsed_contrast, '-', ev_selections, \
                                     group_sep, grouping_var, coding_scheme)

                contrasts_dict[parsed_contrast] = \
                    negative(design_matrix, parsed_EVs_in_contrast[0], \
                             coding_scheme, group_sep, grouping_var)

            if len(contrast_items) > 2:

                idx = 0
                for item in contrast_items:

                    # they need to be put back into Patsy formatted header
                    # titles because the dmatrix gets passed into the function
                    # that writes out the contrast matrix
                    if 'categorical' in ev_selections.keys():
                        for cat_EV in ev_selections['categorical']:

                            if cat_EV in item:

                                if coding_scheme == 'Treatment':
                                    item = item.replace(item, \
                                          'C(' + cat_EV + ')[T.' + item + ']')

                                elif coding_scheme == 'Sum':
                                    item = item.replace(item, \
                                     'C(' + cat_EV + ', Sum)[S.' + item + ']')


                    if idx == 0:

                        if item != '+' and item != '-':

                            contrast_vector = positive(dmatrix, item)

                            if parsed_contrast not in contrasts_dict.keys():
                                contrasts_dict[parsed_contrast] = contrast_vector
                            else:
                                contrasts_dict[parsed_contrast] += contrast_vector

                    elif idx != 0:

                        if item != '+' and item != '-':

                            if contrast_items[idx-1] == '+':

                                contrast_vector = positive(dmatrix, item, \
                                                    coding_scheme, group_sep,\
                                                    grouping_var)

                                if parsed_contrast not in contrasts_dict.keys():
                                    contrasts_dict[parsed_contrast] = contrast_vector
                                else:
                                    contrasts_dict[parsed_contrast] += contrast_vector


                            if contrast_items[idx-1] == '-':

                                contrast_vector = negative(dmatrix, item, \
                                                    coding_scheme, group_sep,\
                                                    grouping_var)

                                if parsed_contrast not in contrasts_dict.keys():
                                    contrasts_dict[parsed_contrast] = contrast_vector
                                else:
                                    contrasts_dict[parsed_contrast] += contrast_vector

                    idx += 1

    # finally write out the .con file and .fts file (if f-tests)
    if (custom_contrasts == None) or (custom_contrasts == '') or \
        ("None" in custom_contrasts):

        print("Writing contrasts file (.con) based on contrasts provided " \
              "using the group analysis model builder's contrasts editor..")

        create_con_file(contrasts_dict, regressor_names, model_name, \
                            current_output, model_out_dir)

        if f_tests:
            create_fts_file(f_tests, contrasts_dict, model_name, \
                                current_output, model_out_dir)

    else:

        print("\nWriting contrasts file (.con) based on contrasts provided " \
              "with a custom contrasts matrix CSV file..\n")

        create_con_ftst_file(custom_contrasts, model_name, current_output, \
                                 model_out_dir, regressor_names, \
                                 coding_scheme, group_sep)
