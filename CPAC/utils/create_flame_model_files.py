
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

        print >>f, '/NumWaves\t%d' %dimy
        print >>f, '/NumPoints\t%d' %dimx
        print >>f, ppstring

        # print labels for the columns - mainly for double-checking your model
        col_string = '\n'

        for col in depatsified_EV_names:
            col_string = col_string + col + '\t'

        print >>f, col_string, '\n'

        print >>f, '/Matrix'

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

    out_file = os.path.join(output_dir, model_name + ".grp")

    with open(out_file, "wt") as f:

        print >>f, '/NumWaves\t1'
        print >>f, '/NumPoints\t%d\n' %dimx
        print >>f, '/Matrix'
        np.savetxt(f, grp_file_vector, fmt='%d', delimiter='\t')

    return out_file



def create_con_file(con_dict, col_names, model_name, current_output, out_dir):

    import os

    out_file = os.path.join(out_dir, model_name) + ".con"

    with open(out_file,'w+') as f:

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
        print >>f, col_string, '\n'

        # write data
        f.write("/Matrix\n")

        for key in con_dict:
            for v in con_dict[key]:
                f.write("%1.5e\t" %v)
            f.write("\n")

    return out_file



def create_fts_file(ftest_list, con_dict, model_name, current_output, \
                        out_dir):

    import os
    import numpy as np

    try:

        print "\nFound f-tests in your model, writing f-tests file " \
              "(.fts)..\n"

        out_file = os.path.join(out_dir, model_name + '.fts')

        with open(out_file, 'w') as f:

            print >>f, '/NumWaves\t', len(con_dict)
            print >>f, '/NumContrasts\t', len(ftest_list)

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
            print >>f, col_string, '\n'

            print >>f, '/Matrix'

            for i in range(fts_n.shape[0]):
                print >>f, ' '.join(fts_n[i].astype('str'))


    except Exception as e:

        filepath = os.path.join(out_dir, "model_files", current_output, \
                                    model_name + '.fts')

        errmsg = "\n\n[!] CPAC says: Could not create .fts file for " \
                    "FLAMEO or write it to disk.\nAttempted filepath: %s\n" \
                    "Error details: %s\n\n" % (filepath, e)

        raise Exception(errmsg)

    return out_file



def create_con_ftst_file(con_file, model_name, current_output, output_dir, \
                             column_names, coding_scheme, group_sep):

    """
    Create the contrasts and fts file
    """

    import os
    import numpy as np

    with open(con_file,"r") as f:
        evs = f.readline()

    evs = evs.rstrip('\r\n').split(',')
    count_ftests = 0

    # remove "Contrasts" label and replace it with "Intercept"
    evs[0] = "Intercept"

    fTest = False

    for ev in evs:
        if "f_test" in ev:
            count_ftests += 1

    if count_ftests > 0:
        fTest = True


    try:

        data = np.genfromtxt(con_file, names=True, delimiter=',', dtype=None)

    except:

        print "Error: Could not successfully read in contrast file: ",con_file
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

        # add Intercept column
        if group_sep == False:
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



    out_file = os.path.join(output_dir, model_name + '.con')

    with open(out_file,"wt") as f:

        idx = 1
        pp_str = '/PPheights'
        re_str = '/RequiredEffect'
        for name in contrast_names:

            print >>f, '/ContrastName%d' %idx, '\t', name
            pp_str += '\t%1.5e' %(1)
            re_str += '\t%1.5e' %(1)
            idx += 1

        print >>f, '/NumWaves\t', (contrasts.shape)[1]
        print >>f, '/NumContrasts\t', (contrasts.shape)[0]
        print >>f, pp_str
        print >>f, re_str + '\n'

        # print labels for the columns - mainly for double-checking your model
        col_string = '\n'
        for ev in evs:
            col_string = col_string + ev + '\t'
        print >>f, col_string, '\n'


        print >>f, '/Matrix'
   
        np.savetxt(f, contrasts, fmt='%1.5e', delimiter='\t')

    ftest_out_file = None
    if fTest:

        print "\nFound f-tests in your model, writing f-tests file (.fts)..\n"

        ftest_out_file = os.path.join(output_dir, model_name + '.fts')

        with open(ftest_out_file,"wt") as f:

            print >>f, '/NumWaves\t', (contrasts.shape)[0]
            print >>f, '/NumContrasts\t', count_ftests

            # print labels for the columns - mainly for double-checking your
            # model
            col_string = '\n'
            for con in contrast_names:
                col_string = col_string + con + '\t'
            print >>f, col_string, '\n'

            print >>f, '/Matrix'

            for i in range(fts_n.shape[0]):
                print >>f, ' '.join(fts_n[i].astype('str'))

    return out_file, ftest_out_file



def create_flame_model_files(design_matrix, col_names, contrasts_dict, \
	custom_contrasts_csv, ftest_list, group_sep, grouping_vector, \
	coding_scheme, model_name, output_measure, output_dir):

    mat_file = write_mat_file(design_matrix, output_dir, model_name, \
        col_names, output_measure)

    grp_file = create_grp_file(design_matrix, grouping_vector, output_dir, \
        model_name)

    if contrasts_dict:
        con_file = create_con_file(contrasts_dict, col_names, model_name, \
            output_measure, output_dir)

        if len(ftest_list) > 0:
	        fts_file = create_fts_file(ftest_list, contrasts_dict, \
                model_name, output_measure, output_dir)
        else:
            fts_file = None

    elif custom_contrasts_csv:
        con_file, fts_file = create_con_ftst_file(custom_contrasts_csv, \
            model_name, output_measure, output_dir, col_names, coding_scheme,\
            group_sep)

    return mat_file, grp_file, con_file, fts_file

