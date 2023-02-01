
def create_dummy_string(length):
    ppstring = ""
    for i in range(0, length):
        ppstring += '\t' + '%1.5e' %(1.0)
    ppstring += '\n' 
    return ppstring


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

    return out_file


def create_grp_file(design_matrix, grp_file_vector, output_dir, model_name):

    import os
    import numpy as np

    dimx = None
    dimy = None

    if len(design_matrix.shape) == 1:
        dimy = 1
        dimx = design_matrix.shape[0]
    else:
        dimx, dimy = design_matrix.shape

    filename = "grouping.grp"

    grp_file_vector = [int(x) for x in grp_file_vector]

    out_file = os.path.join(output_dir, model_name + ".grp")

    with open(out_file, "wt") as f:
        print('/NumWaves\t1', file=f)
        print('/NumPoints\t%d\n' %dimx, file=f)
        print('/Matrix', file=f)
        np.savetxt(f, grp_file_vector, fmt='%d', delimiter='\t')

    return out_file


def create_con_file(con_vecs, con_names, col_names, model_name,
                    current_output, out_dir):

    import os

    out_file = os.path.join(out_dir, model_name) + ".con"

    with open(out_file,'w+') as f:

        # write header
        num = 1

        for key in con_names:
            f.write("/ContrastName%s\t%s\n" %(num,key))
            num += 1

        f.write("/NumWaves\t%d\n" %len(col_names))
        f.write("/NumContrasts\t%d\n" %len(con_names))
        f.write("/PPheights%s" %create_dummy_string(len(con_vecs)))
        f.write("/RequiredEffect%s" %create_dummy_string(len(con_vecs)))
        f.write("\n\n")

        # print labels for the columns - mainly for double-checking your
        # model
        col_string = '\n'
        for col in col_names:
            col_string = col_string + col + '\t'
        print(col_string, '\n', file=f)

        # write data
        f.write("/Matrix\n")

        for vector in con_vecs:
            for v in vector:
                f.write("%1.5e\t" %v)
            f.write("\n")

    return out_file


def create_fts_file(ftest_list, con_names, model_name,
                    current_output, out_dir):

    import os
    import numpy as np

    try:

        print("\nFound f-tests in your model, writing f-tests file " \
              "(.fts)..\n")

        out_file = os.path.join(out_dir, model_name + '.fts')

        with open(out_file, 'w') as f:

            print('/NumWaves\t', len(con_names), file=f)
            print('/NumContrasts\t', len(ftest_list), file=f)

            # process each f-test
            ftst = []

            for ftest_string in ftest_list:

                ftest_vector = []
                
                cons_in_ftest = ftest_string.split(",")

                for con in con_names:
                    if con in cons_in_ftest:
                        ftest_vector.append(1)
                    else:
                        ftest_vector.append(0)

                ftst.append(ftest_vector)
        
            fts_n = np.array(ftst)

            # print labels for the columns - mainly for double-checking your
            # model
            col_string = '\n'
            for con in con_names:
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

    return out_file


def create_con_ftst_file(con_file, model_name, current_output, output_dir,
                         column_names, coding_scheme, group_sep):

    """
    Create the contrasts and fts file
    """

    import os
    import numpy as np

    column_names = [x for x in list(column_names) if 'participant_id' not in x]

    # Read the header of the contrasts file, which should contain the columns
    # of the design matrix and f-tests (if any)
    with open(con_file, "r") as f:
        evs = f.readline()
    evs = evs.rstrip('\r\n').split(',')
    
    if evs[0].strip().lower() != "contrasts":
        print("Error: first cell in contrasts file should contain " \
              "'Contrasts' ")
        raise Exception

    # remove "Contrasts" label and replace it with "Intercept"
    #evs[0] = "Intercept"

    # Count the number of f tests defined
    count_ftests = len([ev for ev in evs if "f_test" in ev ])

    # Whether any f tests are defined
    fTest = count_ftests > 0

    # Now read the actual contrasts
    try:
        contrasts_data = np.genfromtxt(con_file, names=True, delimiter=',',
                                       dtype=None)
    except:
        print("Error: Could not successfully read in contrast file: ",con_file)
        raise Exception

    lst = contrasts_data.tolist()
    # lst = list of rows of the contrast matrix (each row represents a
    # contrast, i.e. a name of the contrast, and then coefficients for each of
    # the design matrix columns, and finally coefficients for each of the
    # f tests specifying whether this contrast is part of that particular
    # f test).

    if isinstance(lst, tuple):
        lst = [lst]

    ftst = []
    fts_columns = []
    contrasts = []
    contrast_names = []

    length = len(list(lst[0]))

    for contr in lst:
        # tp = tuple in the format (contrast_name, 0, 0, 0, 0, ...)
        #      with the zeroes being the vector of contrasts for that contrast

        # extract the name of the contrast
        contrast_names.append(contr[0])

        # create a list of integers that is the vector for the contrast
        # ex. [0, 1, 1, 0, ..]
        con_vector = list(contr)[1:(length-count_ftests)]

        # fts_vector tells us which f-tests this contrast is a part of.
        fts_vector = list(contr)[(length-count_ftests):length]
        fts_columns.append(fts_vector)

        # add Intercept column
        if not group_sep:
            if False:
                # The following insertion gives an error further down the
                # line, because this suggests that there will be an intercept
                # column in the design matrix but such an intercept is never
                # actually added.
                if coding_scheme == "Treatment":
                    con_vector.insert(0, 0)
                elif coding_scheme == "Sum":
                    con_vector.insert(0, 1)

        contrasts.append(con_vector)

    # contrast_names = list of the names of the contrasts (not regressors)
    # contrasts = list of lists with the contrast vectors

    num_EVs_in_con_file = len(contrasts[0])

    contrasts = np.array(contrasts, dtype=np.float16)

    fts_columns = np.array(fts_columns)

    # if there are f-tests, create the array for them
    if fTest:
        ## TODO: Probably it would be more accurate to check that each
        ## f test itself contains enough contrasts, rather than whether
        ## there are in principle enough contrasts to form f tests.
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

        #raise Exception(err_string)
        print(err_string)
        return None, None

    for design_mat_col, con_csv_col in zip(column_names, evs[1:]):
        if con_csv_col not in design_mat_col:
            errmsg = "\n\n[!] CPAC says: The names of the EVs in your " \
                     "custom contrasts .csv file do not match the names or " \
                     "order of the EVs in the design matrix. Please make " \
                     "sure these are consistent.\nDesign matrix EV columns: "\
                     "%s\nYour contrasts matrix columns: %s\n\n" \
                     % (column_names, evs[1:])

            print(errmsg)
            return None, None        

    out_file = os.path.join(output_dir, model_name + '.con')

    with open(out_file,"wt") as f:

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
            if "contrast" not in ev and "Contrast" not in ev:
                col_string = col_string + ev + '\t'
        print(col_string, '\n', file=f)

        print('/Matrix', file=f)
   
        np.savetxt(f, contrasts, fmt='%1.5e', delimiter='\t')

    ftest_out_file = None
    if fTest:

        print("\nFound f-tests in your model, writing f-tests file (.fts)..\n")

        ftest_out_file = os.path.join(output_dir, model_name + '.fts')

        with open(ftest_out_file,"wt") as f:

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

    return out_file, ftest_out_file


def create_flame_model_files(design_matrix, col_names, contrasts_vectors,
                             contrast_names, custom_contrasts_csv, ftest_list,
                             group_sep, grouping_vector, coding_scheme,
                             model_name, output_measure, output_dir):

    if contrasts_vectors:
        con_file = create_con_file(contrasts_vectors, contrast_names,
                                   col_names, model_name, output_measure,
                                   output_dir)

        if len(ftest_list) > 0:
            fts_file = create_fts_file(ftest_list, contrast_names,
                                       model_name, output_measure, output_dir)
        else:
            fts_file = None

    elif custom_contrasts_csv:
        con_file, fts_file = create_con_ftst_file(custom_contrasts_csv,
            model_name, output_measure, output_dir, col_names, coding_scheme,
            group_sep)

    if not con_file:
        # don't write out the rest of the files if this didn't work out
        return None, None, None, None

    mat_file = write_mat_file(design_matrix, output_dir, model_name,
        col_names, output_measure)

    grp_file = create_grp_file(design_matrix, grouping_vector, output_dir,
        model_name)

    return mat_file, grp_file, con_file, fts_file

