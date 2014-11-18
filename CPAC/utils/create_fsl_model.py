import os
import sys
import numpy as np
import csv
import yaml


def create_pheno_dict(gpa_fsl_yml):

    def read_phenotypic(pheno_file, ev_selections, subject_id_label):

        import csv
        import numpy as np

        ph = pheno_file

        # Read in the phenotypic CSV file into a dictionary named pheno_dict
        # while preserving the header fields as they correspond to the data
        p_reader = csv.DictReader(open(os.path.abspath(ph), 'rU'), skipinitialspace=True)

        #pheno_dict_list = []
        
        # dictionary to store the data in a format Patsy can use
        # i.e. a dictionary where each header is a key, and the value is a
        # list of all of that header's values
        pheno_data_dict = {}

        for line in p_reader:

            # here, each instance of 'line' is really a dictionary where the
            # keys are the pheno headers, and their values are the values of
            # each EV for that one subject - each iteration of this loop is
            # one subject

            for key in line.keys():

                if key not in pheno_data_dict.keys():
                    pheno_data_dict[key] = []

                # create a list within one of the dictionary values for that
                # EV if it is categorical; formats this list into a form
                # Patsy can understand regarding categoricals:
                #     example: { ADHD: ['adhd1', 'adhd1', 'adhd0', 'adhd1'] }
                #                instead of just [1, 1, 0, 1], etc.
                if 'categorical' in ev_selections.keys():
                    if key in ev_selections['categorical']:
                        pheno_data_dict[key].append(key + str(line[key]))

                    elif key == subject_id_label:
                        pheno_data_dict[key].append(line[key])

                    else:
                        pheno_data_dict[key].append(float(line[key]))

                elif key == subject_id_label:
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


    # pheno_data_dict gets loaded with the phenotypic data, in a dictionary
    # formatted for proper use with Patsy
    pheno_data_dict = read_phenotypic(gpa_fsl_yml.pheno_file, gpa_fsl_yml.ev_selections, gpa_fsl_yml.subject_id_label)

    return pheno_data_dict






def check_multicollinearity(matrix):

    U, s, V = np.linalg.svd(matrix)

    max_singular = np.max(s)
    min_singular = np.min(s)

    print max_singular, ' ~~~~~~ ', min_singular
    print 'RANK ~~~~~~ ', np.linalg.matrix_rank(matrix)
    if min_singular == 0:

        return 1

    else:

        condition_number = float(max_singular)/float(min_singular)
        print 'condition_number %f' % condition_number
        if condition_number > 30:

            return 1

    return 0





def create_mat_file(data, col_names, model_name, outputModelFilesDirectory):

    """
    create the .mat file
    """

    dimx = None
    dimy = None

    if len(data.shape) == 1:
        dimy = 1
        dimx = data.shape[0]
    else:
        dimx, dimy = data.shape


    ppstring = '/PPheights'

    for i in range(0, dimy):

        ppstring += '\t' + '%1.5e' %(1.0)

    ppstring += '\n'


    f = open(os.path.join(outputModelFilesDirectory, model_name + '.mat'), 'w')

    print >>f, '/NumWaves\t%d' %dimy
    print >>f, '/NumPoints\t%d' %dimx
    print >>f, ppstring

    # print labels for the columns - mainly for double-checking your model
    col_string = '\n'

    for col in col_names:
        col_string = col_string + col + '\t'

    print >>f, col_string, '\n'


    print >>f, '/Matrix'

    np.savetxt(f, data, fmt='%1.5e', delimiter='\t')

    f.close()



def create_grp_file(data, model_name, gp_var, outputModelFilesDirectory):

    """
    create the grp file
    """

    dimx = None
    dimy = None
    if len(data.shape) == 1:
        dimy = 1
        dimx = data.shape[0]
    else:
        dimx, dimy = data.shape
    data = np.ones(dimx)

    if not (gp_var == None):
        i = 1
        for key in sorted(gp_var.keys()):

            for index in gp_var[key]:
                data[index] = i

            i += 1


    f = open(os.path.join(outputModelFilesDirectory, model_name + '.grp'), 'w')

    print >>f, '/NumWaves\t1'
    print >>f, '/NumPoints\t%d\n' %dimx
    print >>f, '/Matrix'
    np.savetxt(f, data, fmt='%d', delimiter='\t')

    f.close()



def create_con_ftst_file(con_file, model_name, fTest, outputModelFilesDirectory):

    """
    Create the contrasts and fts file
    """

    evs = open(con_file, 'r').readline()
    evs = evs.rstrip('\r\n').split(',')
    count_ftests = 0

    for ev in evs:

        if 'f_test' in ev.lower():

            count_ftests += 1


    try:

        data = np.genfromtxt(con_file, names=True, delimiter=',', dtype=None)

    except:

        print "Error: Could not successfully read in contrast file: ", con_file
        raise Exception


    lst = data.tolist()


    ftst = []
    contrasts = []
    contrast_names = []

    length = None
    length = len(list(lst[0]))


    try:

        for tp in lst:

            contrast_names.append(tp[0])
            contrasts.append(list(tp)[1:length-count_ftests])

            if fTest:
                ftst.append(list(tp[length-count_ftests: length]))

        contrasts = np.array(contrasts, dtype=np.float16)
        
        if fTest:
            fts_n = np.array(ftst)

    except:
        print "\n\n" + "ERROR: Not enough contrasts for running f-tests." \
              "\n Tip: Do you have only one contrast in your contrasts file?" \
              " f-tests require more than one contrast." + "\n" + \
              "Either turn off f-tests or include more contrasts." + "\n" + \
              "Error name: create_fsl_model_0002" + "\n\n"
        raise Exception


    try:

        f = open(os.path.join(outputModelFilesDirectory, model_name + '.con'), 'w')

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
        print >>f, '/Matrix'
   
        np.savetxt(f, contrasts, fmt='%1.5e', delimiter='\t')

        f.close()

    except:
        print "Error: Could not create .con file."
        print ""
        raise Exception


    if fTest:

        try:

            fts_n = fts_n.T
            f = open(os.path.join(outputModelFilesDirectory, model_name + '.fts'), 'w')
            print >>f, '/NumWaves\t%d' % (contrasts.shape)[0]
            print >>f, '/NumContrasts\t%d\n' % count_ftests

            print >>f, '/Matrix'

            for i in range(fts_n.shape[0]):
                print >>f, ' '.join(fts_n[i].astype('str'))
            #np.savetxt(f, fts_n[None], fmt='%d', delimiter=' ')
            f.close()

        except:
            print "Error: Could not create .fts file."
            print ""
            raise Exception



"""
Class to set dictionary keys as map attributes
"""
class Configuration(object):
    def __init__(self, config_map):
        for key in config_map:
            if config_map[key] == 'None':
                config_map[key] = None
            setattr(self, key, config_map[key])



def pandas_alternate_organize_data(data, c):

    import pandas as ps
    import csv

    df = ps.DataFrame(data)

    categorical = []
    directional = []
    for i in range(0, len(c.columnsInModel)):

        if c.deMean[i]:
            col = c.columnsInModel[i]
            df[col] = df[col].astype('float32') - df[col].astype('float32').mean()

        if c.categoricalVsDirectional[i]:
            categorical.append(c.columnsInModel[i])
        else:
            directional.append(c.columnsInModel[i])

    #split on the grouping variable
    for name, group in df.groupby(c.groupingVariable):
        group[c.groupingVariable] = 1
        group[c.groupingVariable]= group[c.groupingVariable]
        df[c.groupingVariable + '__'+ name] = group[c.groupingVariable]
        df[c.groupingVariable + '__'+ name] = df[c.groupingVariable + '__'+ name].fillna(0)

        for col in directional:

            group[col] = group[col]
            df[col + '__'+ name] = group[col]
            df[col + '__'+ name] = df[col + '__'+ name].fillna(0)


    #split on (grouping variable and each of the (categoricals- grouping variable) )
    for col in categorical:

        if not (col == c.groupingVariable):
            for name, group in df.groupby([c.groupingVariable, col]):
                group[col] = 1
                df[col+'__'+'_'.join([str(e) for e in name])] = group[col]
                df[col+'__'+'_'.join([str(e) for e in name])] = df[col+'__'+'_'.join([str(e) for e in name])].fillna(0)

    for col in c.columnsInModel:
        del df[col]


    df.to_csv('./tempfile.csv', sep=',', index=False)

    print 'saved to abc.csv'
    sys.exit()



def split_directionals(gp, directional, data, c):

    for key in gp.keys():

        indices = gp[key]

        for col in directional:

            new_col = col + '__' + key

            for idx in range(0, len(data)):

                if idx in indices:
                    data[idx][new_col] = data[idx][col]
                else:
                    data[idx][new_col] = 0



def split_gp_var(gp, data, c):


    for key in gp.keys():

        indices = gp[key]

        new_col = c.groupingVariable + '__' + key

        for idx in range(0, len(data)):

            if idx in indices:
                data[idx][new_col] = 1
            else:
                data[idx][new_col] = 0



def group_by_gp_categorical(categorical, data, c):


    gp_cat_dict = {}
    for cat_col in categorical:

        if not (cat_col == c.groupingVariable):

            for idx in range(0, len(data)):
                new_col = '__'.join([str(cat_col), str(data[idx][cat_col]), str(c.groupingVariable), str(data[idx][c.groupingVariable])])
                if new_col in gp_cat_dict:
                    gp_cat_dict[new_col].append(idx)
                else:
                    gp_cat_dict[new_col] = [idx]



    for col in gp_cat_dict.keys():

        indices = gp_cat_dict[col]

        for idx in range(0, len(data)):

            if idx in indices:
                data[idx][col] = 1
            else:
                data[idx][col] = 0



def alternate_organize_data(data, c):

    categorical = []
    directional = []
    mean_cols = []

    groups_grouping_var = {}

    for i in range(0, len(c.columnsInModel)):

        if c.deMean[i]:
            col = c.columnsInModel[i]
            mean_cols.append(col)

        if c.categoricalVsDirectional[i]:
            categorical.append(c.columnsInModel[i])
        else:
            directional.append(c.columnsInModel[i])


    sum_ = {}

    idx = 0
    #take sum of each of columns to be demeaned
    for row in data:

        #on the side create groups for groups_grouping_var
        try:
            groups_grouping_var[str(row[c.groupingVariable])].append(idx)
        except:
            groups_grouping_var[str(row[c.groupingVariable])] = [idx]



        for col in mean_cols:
            try:
                sum_[col] += float(row[col])
            except:
                sum_[col] = float(row[col])

        idx += 1

    #take the mean
    for  row in data:

        for col in mean_cols:

            row[col] = float(row[col]) - float(sum_[col])/float(len(data))


    #split the directonal columns according to the groupingVariable
    split_directionals(groups_grouping_var, directional, data, c)
    #split the groupingVariable col
    split_gp_var(groups_grouping_var, data, c)
    #split categorical cols according to groupingVariable
    group_by_gp_categorical(categorical, data, c)

    #delete all original categorical and directional columns 
    for idx in range(0, len(data)):
        #del data[idx]['subject_id']
        for col in directional:
            del data[idx][col]

        for col in categorical:
            del data[idx][col]

    print '\t'.join(key for key in sorted(data[0].keys()) )
    for idx in range(0, len(data)):
        print '\t'.join(str(data[idx][key]) for key in sorted(data[idx].keys()))

    return data, data[0].keys(), groups_grouping_var




def run(config, fTest, param_file, pipeline_path, current_output, CPAC_run = False):

    # create_fsl_model.run()
    # this is called from cpac_group_analysis_pipeline.py

    # it collects the information the user provided for the FSL gpa model
    # which was saved in the group analysis FSL config .yaml file, and then
    # puts it into Patsy to create a design matrix
    #     it then also generates the contrast file from the contrasts the user
    #     provided
    #         ultimately this produces the .mat, .con and .grp files needed
    #         for FLAMEO for group analysis

    # see more info on Patsy:
    #     http://patsy.readthedocs.org/en/latest/overview.html


    # open the GROUP ANALYSIS FSL .YML CONFIG FILE, not the main pipeline
    # config .yml file!
    if CPAC_run:
        c = config
    else:
        try:
            c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
        except:
            raise Exception("Error in reading %s configuration file" % config)

    import csv
    import numpy as np

    # return the data from the phenotype file processed properly for Patsy
    # and load it into 'pheno_data_dict'
    #     format: dictionary, each key is the name of an EV, and its value is
    #             a LIST of values in order of the subjects
    #                 - categorical EVs are already renamed from '0,1,..' to
    #                   'EV0,EV1,..' with EV being the EV name
    #                 - EVs to be demeaned are already demeaned
    #                 - numerical EVs (non-categorical) are in a list which
    #                   have been converted into a NumPy array
    pheno_data_dict = create_pheno_dict(c)


    if param_file != None:

        ''' extract motion measures for insertion as EVs if selected '''
        # insert MeanFD or other measures into pheno_data_dict
        #     first, pull the measure values from the all_params .csv file written
        #     to the individual-level analysis output directory
        #     then, ensure the values are in the same order as the subject ids

        measures = ['MeanFD', 'MeanFD_Jenkinson', 'MeanDVARS']

        try:

            measure_dict = {}
            f = csv.DictReader(open(param_file,'rU'))

            for line in f:
                measure_map = {}
                for m in measures:
                    if line.get(m):
                        measure_map[m] = line[m]

                measure_dict[line['Subject']] = measure_map
           
        except:
            print '\n\n[!] CPAC says: Could not extract required information ' \
                  'from the parameters file.\n'
            print 'Path: ', param_file, '\n\n'
            raise Exception



        # function to demean measures the user included in the design formula
        # and then insert them in the right location in the pheno_data_dict
        def add_measure_to_pheno(measure_name):

            measure_list = []

            # create a blank list that is the proper length
            for sub in pheno_data_dict[c.subject_id_label]:
                measure_list.append(0)

            for subID in measure_dict.keys():

                # find matching subject IDs between the measure_dict and the
                # pheno_data_dict so we can insert measure values into the
                # pheno_data_dict
                for subject in pheno_data_dict[c.subject_id_label]:

                    if subject == subID:

                        # return the index (just an integer) of where in the
                        # pheno_data_dict list structure a subject ID is
                        idx = np.where(pheno_data_dict[c.subject_id_label]==subID)[0][0]

                        # insert Mean FD value in the proper point
                        measure_list[idx] = float(measure_dict[subID][measure_name])


            # time to demean the MeanFD values
            measure_sum = 0.0

            for measure in measure_list:
                measure_sum = measure_sum + measure

            measure_mean = measure_sum / len(measure_list)

            idx = 0

            for measure in measure_list:
                measure_list[idx] = measure - measure_mean
                idx += 1

            # add this new list to the pheno_data_dict
            pheno_data_dict[measure_name] = np.array(measure_list)

        ''' insert measures into pheno data '''
        # add measures selected in the design formula into pheno_data_dict
        # they are also demeaned prior
        for measure in measures:
            if measure in c.design_formula:
                add_measure_to_pheno(measure)



    if 'Measure_Mean' in c.design_formula:

        ''' extract the mean of derivative for each subject if selected '''
        # if the user has selected it to be part of their model, insert the mean
        # of the outputs included in group analysis (i.e. if running ReHo in
        # group-level analysis, have the mean of each subject's ReHo output
        # included as an EV in the phenotype - regress out the mean of measure
        #     pull the mean value from the output_means.csv file in the subject
        #     directory of the appropriate pipeline's output folder
        sub_means_dict = {}
        output_means_dict = {}

        for sub in pheno_data_dict[c.subject_id_label]:

            output_means_file = os.path.join(pipeline_path, sub, 'output_means_%s.csv' % sub)

            if os.path.exists(output_means_file):
            
                try:

                    output_means = csv.DictReader(open(output_means_file,'rU'))
                
                except:

                    print '\n\n[!] CPAC says: Could not open the output_means' \
                          '.csv file usually located in each subject\'s output ' \
                          'folder in the output directory.\n'
                    print 'Path: ', output_means_file, '\n\n'
                    raise Exception


                # pull in the output_means .csv as a dictionary
                for row in output_means:
                    sub_means_dict = row


                try:

                    output_means_dict[sub] = str(row[current_output])

                except:

                    print '\n\n[!] CPAC says: There is no mean value ' \
                          'stored for the output \'', current_output, \
                          '\' for subject \'', sub, '\'.\n'
                    print 'Path to means file: ', output_means_file, '\n'
                    print 'Possible situations:\n1. The output \'', \
                          current_output, '\' was not included in ' \
                          'individual-level analysis, but was included to ' \
                          'be run in group-level analysis.\n2. The means ' \
                          'file for this subject was not created properly.' \
                          '\n3. Individual-level analysis did not ' \
                          'complete properly.\n\n'
                    raise Exception


            else:
                print '\n\n[!] CPAC says: The output_means.csv file usually ' \
                      'located in each subject\'s output folder in the output ' \
                      'directory does not exist!\n'
                print 'Path not found: ', output_means_file, '\n\n'
                print 'Tip: Either check if individual-level analysis ' \
                      'completed successfully, or remove the measure mean ' \
                      'from your model design.\n\n'
                raise Exception

    
        # by the end of this for loop above, output_means_dict should look
        # something like this:
        #    {sub1: mean_val, sub2: mean_val, ..}
        #        as this code runs once per output, this dictionary contains
        #        the mean values of the one current output, right now




        ''' insert mean of derivatives into pheno data '''
        means_list = []

        # create a blank list that is the proper length
        for sub in pheno_data_dict[c.subject_id_label]:
            means_list.append(0)

        for subID in output_means_dict.keys():

            # find matching subject IDs between the output_means_dict and the
            # pheno_data_dict so we can insert mean values into the
            # pheno_data_dict
            for subject in pheno_data_dict[c.subject_id_label]:

                if subject == subID:

                    # return the index (just an integer) of where in the
                    # pheno_data_dict list structure a subject ID is
                    idx = np.where(pheno_data_dict[c.subject_id_label]==subID)[0][0]

                    # insert Mean FD value in the proper point
                    means_list[idx] = float(output_means_dict[subID])


        # time to demean the means!
        means_sum = 0.0

        for mean in means_list:

            means_sum = means_sum + mean

        measure_mean = means_sum / len(means_list)

        idx = 0

        for mean in means_list:

            means_list[idx] = mean - measure_mean
            idx += 1


        ''' insert means into pheno data if selected '''

        # add this new list to the pheno_data_dict
        pheno_data_dict['Measure_Mean'] = np.array(means_list)
   


    # make sure the group analysis output directory exists
    try:

        if not os.path.exists(c.output_dir):

            os.makedirs(c.output_dir)

    except:

        print '\n\n[!] CPAC says: Could not successfully create the group ' \
              'analysis output directory:\n', c.output_dir, '\n\nMake sure ' \
              'you have write access in this file structure.\n\n\n'
        raise Exception




    ''' create the Patsy design matrix '''
    # parse through ev_selections, find the categorical names within the
    # design formula and insert C(<name>, Sum) into the design formula
    #     this is required for Patsy to process the categorical EVs properly
    #     when generating the design matrix (this goes into the .mat file)
    formula = c.design_formula

    coding_scheme = c.coding_scheme[0]


    if 'categorical' in c.ev_selections.keys():
        for EV_name in c.ev_selections['categorical']:

            if coding_scheme == 'Treatment':
                formula = formula.replace(EV_name, 'C(' + EV_name + ')')
            elif coding_scheme == 'Sum':
                formula = formula.replace(EV_name, 'C(' + EV_name + ', Sum)')



    # create the actual design matrix using Patsy
    import patsy


    # drop pickles of the inputs meant for Patsy so you can manually test it
    # later if needed
    #import pickle
    #pickle.dump(formula, open(c.output_dir + '/' + "formula.p", "wb" ) )
    #pickle.dump(pheno_data_dict, open(c.output_dir + '/' + "data_dict.p", "wb" ) )


    try:
        dmatrix = patsy.dmatrix(formula, pheno_data_dict, NA_action='raise')
    except:
        print '\n\n[!] CPAC says: Design matrix creation wasn\'t ' \
                'successful - do the terms in your formula correctly ' \
                'correspond to the EVs listed in your phenotype file?\n'
        print 'Phenotype file provided: '
        print c.pheno_file, '\n\n'
        raise Exception




    ### CREATE CONTRAST FILE

    # parse in user-input contrast strings that were selected, and generate
    # the contrast file (.con)

    def greater_than(dmat, a, b, coding):
        c1 = positive(dmat, a, coding)
        c2 = positive(dmat, b, coding)
        return c1-c2

    def positive(dmat, a, coding):

        if coding == "Treatment":

            evs = dmat.design_info.column_name_indexes
            con = np.zeros(dmat.shape[1])

            print "a: ", a
            print "evs: ", evs

            if a in evs:
                print "a is in evs."
                con[evs[a]] = 1
            else:
                print "a is not in evs."
                #it is a dropped term so make all other terms in that category at -1
                term = a.split('[')[0]
                print "term: ", term
                for ev in evs:
                    if ev.startswith(term):
                        con[evs[ev]]= -1

            # make Intercept 0
            con[0] = 0

            return con

        elif coding == "Sum":

            evs = dmat.design_info.column_name_indexes
            con = np.zeros(dmat.shape[1])
            if a in evs:
                con[evs[a]] = 1
            else:
                #it is a dropped term so make all other terms in that category at -1
                term = a.split('[')[0]
                for ev in evs:
                    if ev.startswith(term):
                        con[evs[ev]]= -1

            # make Intercept 1
            con[0] = 1

            return con


    def negative(dmat, a, coding):
        con = 0-positive(dmat, a, coding)
        return con

    def create_dummy_string(length):
        ppstring = ""
        for i in range(0, length):
            ppstring += '\t' + '%1.5e' %(1.0)
        ppstring += '\n' 
        return ppstring

    def create_con_file(con_dict, col_names, file_name, out_dir):
        with open(os.path.join(out_dir, file_name)+".con",'w+') as f:
            #write header
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

            #write data
            f.write("/Matrix\n")

            for key in con_dict:
                for v in con_dict[key]:
                    f.write("%1.5e\t" %v)
                f.write("\n")

   

    contrasts = c.contrasts
    contrasts_list = []
    contrasts_dict = {}

    # collect the user-selected contrast strings into a list of strings
    for contrast in contrasts.keys():
        if contrasts[contrast] == True:

            contrasts_list.append(contrast)


    # take the contrast strings and process them appropriately
    #     extract the two separate contrasts (if there are two), and then
    #     identify which are categorical - adapting the string if so
    def process_contrast(operator):

        parsed_EVs_in_contrast = []

        EVs_in_contrast = parsed_contrast.split(operator)

        if '' in EVs_in_contrast:
            EVs_in_contrast.remove('')


        for EV in EVs_in_contrast:

            skip = 0

            # they need to be put back into Patsy formatted header titles
            # because the dmatrix gets passed into the function that writes
            # out the contrast matrix
            if 'categorical' in c.ev_selections.keys():
                for cat_EV in c.ev_selections['categorical']:

                    if cat_EV in EV:

                        # handle interactions
                        if ":" in EV:
                            temp_split_EV = EV.split(":")
                            for interaction_EV in temp_split_EV:
                                if cat_EV in interaction_EV:
                                    current_EV = interaction_EV
                        else:
                            current_EV = EV

                        if coding_scheme == 'Treatment':
                            cat_EV_contrast = EV.replace(EV, 'C(' + cat_EV + ')[T.' + current_EV + ']')

                        elif coding_scheme == 'Sum':
                            cat_EV_contrast = EV.replace(EV, 'C(' + cat_EV + ', Sum)[S.' + current_EV + ']')

                        parsed_EVs_in_contrast.append(cat_EV_contrast)
                        skip = 1
                    
            if skip == 0:

                parsed_EVs_in_contrast.append(EV)

            # handle interactions
            if ":" in EV and len(parsed_EVs_in_contrast) == 2:

                parsed_EVs_in_contrast = [parsed_EVs_in_contrast[0] + ":" + parsed_EVs_in_contrast[1]]

            if ":" in EV and len(parsed_EVs_in_contrast) == 3:

                parsed_EVs_in_contrast = [parsed_EVs_in_contrast[0], parsed_EVs_in_contrast[1] + ":" + parsed_EVs_in_contrast[2]]


        return parsed_EVs_in_contrast



    # parse the user-input contrast strings
    for contrast in contrasts_list:
        # each 'contrast' is a string the user input of the desired contrast

        # remove all spaces
        parsed_contrast = contrast.replace(' ', '')


        EVs_in_contrast = []
        parsed_EVs_in_contrast = []

        if '>' in parsed_contrast:

            print "parsed contrast: ", parsed_contrast

            parsed_EVs_in_contrast = process_contrast('>')

            print "parsed EVs in contrast: ", parsed_EVs_in_contrast

            contrasts_dict[parsed_contrast] = greater_than(dmatrix, parsed_EVs_in_contrast[0], parsed_EVs_in_contrast[1], coding_scheme)

            print "contrasts_dict[parsed_contrast]: ", contrasts_dict[parsed_contrast]


        elif '<' in parsed_contrast:

            parsed_EVs_in_contrast = process_contrast('<')

            contrasts_dict[parsed_contrast] = greater_than(dmatrix, parsed_EVs_in_contrast[1], parsed_EVs_in_contrast[0], coding_scheme)


        else:

            contrast_string = parsed_contrast.replace('+',',+,')
            contrast_string = contrast_string.replace('-',',-,')
            contrast_items = contrast_string.split(',')

            if '' in contrast_items:
                contrast_items.remove('')

            if '+' in contrast_items and len(contrast_items) == 2:

                parsed_EVs_in_contrast = process_contrast('+')

                contrasts_dict[parsed_contrast] = positive(dmatrix, parsed_EVs_in_contrast[0], coding_scheme)

            elif '-' in contrast_items and len(contrast_items) == 2:

                parsed_EVs_in_contrast = process_contrast('-')

                contrasts_dict[parsed_contrast] = negative(dmatrix, parsed_EVs_in_contrast[0], coding_scheme)


            if len(contrast_items) > 2:

                idx = 0
                for item in contrast_items:

                    # they need to be put back into Patsy formatted header titles
                    # because the dmatrix gets passed into the function that writes
                    # out the contrast matrix
                    if 'categorical' in c.ev_selections.keys():
                        for cat_EV in c.ev_selections['categorical']:

                            if cat_EV in item:

                                if coding_scheme == 'Treatment':
                                    item = item.replace(item, 'C(' + cat_EV + ')[T.' + item + ']')

                                elif coding_scheme == 'Sum':
                                    item = item.replace(item, 'C(' + cat_EV + ', Sum)[S.' + item + ']')


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

                                contrast_vector = positive(dmatrix, item)

                                if parsed_contrast not in contrasts_dict.keys():
                                    contrasts_dict[parsed_contrast] = contrast_vector
                                else:
                                    contrasts_dict[parsed_contrast] += contrast_vector


                            if contrast_items[idx-1] == '-':

                                contrast_vector = negative(dmatrix, item)

                                if parsed_contrast not in contrasts_dict.keys():
                                    contrasts_dict[parsed_contrast] = contrast_vector
                                else:
                                    contrasts_dict[parsed_contrast] += contrast_vector


                    idx += 1        



    ### CREATE the .mat, .con, and .grp files

    # convert the Patsy-generated design matrix into a NumPy array
    data = np.asarray((dmatrix))




    if check_multicollinearity(np.array(data)) == 1:

        print '\n\n[!] CPAC warns: Detected multicollinearity in the ' \
              'computed group-level analysis model. Please double-check ' \
              'your model design.\n\n'



    if c.group_sep == True and (c.grouping_var == None or (c.grouping_var not in c.design_formula)):
        print '\n\n[!] CPAC says: Model group variances separately is ' \
              'enabled, but the grouping variable set is either set to ' \
              'None, or was not included in the model as one of the EVs.\n'
        print 'Design formula: ', c.design_formula
        print 'Grouping variable: ', c.grouping_var, '\n\n'
        raise Exception



    # prep data and column names if user decides to model group variances
    # separately

    if c.group_sep == True:

        EV_options = []
        grouping_options = []
        new_options = []

        # take in what the grouping variable is. get the names of the options.
        for col_name in dmatrix.design_info.column_names:

            if col_name != 'Intercept':
                
                skip = 0

                if 'categorical' in c.ev_selections.keys():
                    for cat_EV in c.ev_selections['categorical']:
                        if cat_EV in col_name:

                            if coding_scheme == 'Treatment':
                                cat_EV_stripped = col_name.replace('C(' + cat_EV + ')[T.', '')
                            elif coding_scheme == 'Sum':
                                cat_EV_stripped = col_name.replace('C(' + cat_EV + ', Sum)[S.', '')

                            cat_EV_stripped = cat_EV_stripped.replace(']', '')
                            EV_options.append(cat_EV_stripped)
                            skip = 1

                if skip == 0:
                    EV_options.append(col_name)



        idx = 1

        for ev in EV_options:

            if c.grouping_var in ev:

                grouping_variable_info = []

                grouping_variable_info.append(ev)
                grouping_variable_info.append(idx)

                grouping_options.append(grouping_variable_info)

                # grouping_var_idx is the column numbers in the design matrix
                # which holds the grouping variable (and its possible levels)


            idx += 1


        # all the categorical values/levels of the grouping variable
        grouping_var_levels = []

        for gv_idx in grouping_options:
            for subject in dmatrix:

                level_label = '__' + gv_idx[0] + '_' + str(subject[gv_idx[1]])

                if level_label not in grouping_var_levels:
                    grouping_var_levels.append(level_label)


        # make the new header for the reorganized data
        for ev in EV_options:
            if c.grouping_var not in ev:
                for level in grouping_var_levels:
                    new_options.append(ev + level)
            elif c.grouping_var in ev:
                new_options.append(ev)



        grouped_data = []

        for subject in dmatrix:

            new_row = []

            # put in the Intercept first
            new_row.append(subject[0])

            for option in grouping_options:

                grouping_var_id = option[1]
                gp_var_value = subject[grouping_var_id]
                gp_var_label = '_' + str(gp_var_value)


                for orig_idx in range(1,len(subject)):

                    # if the current ev_value in the current subject line is the
                    # grouping variable
                    if orig_idx == grouping_var_id:
                        new_row.append(subject[orig_idx])

                    else:

                        for new_header in new_options:

                            if EV_options[orig_idx-1] in new_header:
                                if gp_var_label in new_header:
                                    new_row.append(subject[orig_idx])
                                else:
                                    new_row.append(0)

          

            grouped_data.append(new_row)


        data = np.array(grouped_data, dtype=np.float16)

        new_options.insert(0, 'Intercept')

        column_names = new_options

    else:

        data = np.array(data, dtype=np.float16)

        column_names = dmatrix.design_info.column_names



    '''
    still need:
        contrast handling?
    '''


    try:
        create_mat_file(data, column_names, c.model_name, c.output_dir)
    except:
        print '\n\n[!] CPAC says: Could not create .mat file during ' \
                  'group-level analysis model file generation.\n'
        print 'Attempted output directory: ', c.output_dir, '\n\n'
        raise Exception

    try:
        create_grp_file(data, c.model_name, c.grouping_var, c.output_dir)
    except:
        print '\n\n[!] CPAC says: Could not create .grp file during ' \
                  'group-level analysis model file generation.\n'
        print 'Attempted output directory: ', c.output_dir, '\n\n'
        raise Exception

    try:
        create_con_file(contrasts_dict, column_names, c.model_name, c.output_dir)
    except:
        print '\n\n[!] CPAC says: Could not create .con file during ' \
                  'group-level analysis model file generation.\n'
        print 'Attempted output directory: ', c.output_dir, '\n\n'
        raise Exception







