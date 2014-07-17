import os
import sys
import numpy as np
import csv
import yaml


# what this file now needs to do:

# - demean the selected EVs
# - create the pheno data dict like it does in the modelconfig

# - take the formula and put that and the data into dmatrix, with Sum encoding
# - parse the contrasts strings, use Aimi's code in that part
# - write the .mat and .con files, and .grp as well
# - figure out multicollinearity (if patsy doesn't do this)
# - figure out modeling group variances separately






def create_pheno_dict(gpa_fsl_yml):

    # \/ \/ DUPLICATED CODE! PUT THIS FUNCTION IN ONE PLACE,
    #       OR REMOVE THE ORIGINAL ONE
    def read_phenotypic(pheno_file, ev_selections):

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

            for key in line.keys():

                if key not in pheno_data_dict.keys():
                    pheno_data_dict[key] = []

                # create a list within one of the dictionary values for that
                # EV if it is categorical; formats this list into a form
                # Patsy can understand regarding categoricals:
                #     example: { ADHD: ['adhd1', 'adhd1', 'adhd0', 'adhd1'] }
                #                instead of just [1, 1, 0, 1], etc.
                if key in ev_selections['categorical']:
                    pheno_data_dict[key].append(key + str(line[key]))

                else:
                    pheno_data_dict[key].append(line[key])



        # this needs to run after each list in each key has been fully
        # populated above
        for key in pheno_data_dict.keys():


            # converts non-categorical EV lists into NumPy arrays
            # so that Patsy may read them in properly
            if key not in ev_selections['categorical']:
            
                pheno_data_dict[key] = np.array(pheno_data_dict[key])


            # demean the EVs marked for demeaning
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



        return pheno_data_dict


    # pheno_data_dict gets loaded with the phenotypic data, in a dictionary
    # formatted for proper use with Patsy
    pheno_data_dict = read_phenotypic(gpa_fsl_yml.pheno_file, gpa_fsl_yml.ev_selections)

    print pheno_data_dict
    raise





def filter_phenotypic(c):

    """
    The Main purpose of this function is to extract phenotypic data from the phenotypic file for the subjects in the subject list.
    The phenotypic data corresponding to column names specified in columnsInModel variable in the config_fsl.py is extracted.

    Parameters
    ----------

    c : The configuration file object containing all the variables specified in the configuration file. 

    Returns
    -------

    f_r: List of maps. Each map contains data corresponding to a row in the phenotypic file, but only for the columns specified
    in columnsInModel variable in the config_fsl.py

    """

    ph = c.phenotypicFile
    sublist = c.subjectListFile
 
    f = open(sublist, 'rU')

    subjects = f.readlines()

    sub_dict = {}


    # Read in the phenotypic CSV file into a dictionary named pheno_dict
    # while preserving the header fields as they correspond to the data
    p_reader = csv.DictReader(open(os.path.abspath(ph), 'rU'), skipinitialspace=True)

    pheno_dict_list = []
    for line in p_reader:
    
        pheno_dict_list.append(line)
        
        # pheno_dict_list is a list of dictionaries of phenotype header items
        # matched to their values, which also includes subject IDs
            
        # i.e. [{'header1': 'value', 'header2': 'value'}, {'header1': 'value', 'header2': 'value'}, ..]
            
        # these dictionaries are UNORDERED, i.e. header items ARE NOT ORDERED



    # Creates a dictionary sub_dict which stores the amount of
    # instances of each subject, per subject
  
    for subject in subjects:

        subject = subject.rstrip('\r\n')
        subject = subject.replace(',','_')


        if not subject in sub_dict:
            sub_dict[subject] = 1
        else:
            sub_dict[subject] += 1



    # Iterate over phenotypic CSV file removing any fields that are not listed
    # in the columnsInModel parameter in the FSL config file, and then write
    # the remaining fields into record_dict dictionary

    f_r = []


    try:

        record_dict = {}
      
        for record in pheno_dict_list:
            
            # "record" is a dictionary of phenotype header items matched
            # to their values and it is UNORDERED
            
            # if record[c.subjectColumn] (which is a subject ID number)
            # is in sub_dict or not
            if record[c.subjectColumn] in sub_dict:
                
                # record.keys() is a list of phenotype header items.
                # here, we remove any header items which are not included in the
                # 'columnsInModel' field in the group analysis FSL config file
                for rec in record.keys():
                    
                    if (not (rec in c.columnsInModel)) and not (c.subjectColumn == rec):
                        
                        del record[rec]
                            

                # record_dict is a dictionary of dictionaries, where the key is the
                # subject ID and the value is a dictionary of header items matched
                # to their values (WARNING: header items unordered)
                record_dict[record[c.subjectColumn]] = record
                

    except:
    
        print "Error processing phenotypic file data: ", ph
        print "Try checking your phenotype .csv file in a text editor " \
                "and ensuring the fields are separated by commas.\n\n"
        print "\n"
        raise Exception



    for subject in subjects:

        subject = subject.rstrip('\r\n')
        subject = subject.replace(',','_')

        try:

            # record_dict.keys() is a list of subject IDs within
            # the record_dict dictionary
            if subject in record_dict.keys():
                
                # f_r is a list of dictionaries of phenotype header items
                # matched to their values - like record_dict above, except
                # the dictionaries are not matched to a subject ID
                
                # HEADER ITEMS STILL UNORDERED
                
                # HOWEVER, the subject ID is present within each dictionary
                # in the f_r list, matched with a key named after the subject
                # column header item
                
                f_r.append(record_dict[subject])

        except:

            print "Exception: Could not read from record lookup table for subject #: ", subject
            raise Exception


    return f_r



def organize_data(filter_data, c):

    """

    The main purpose of this function is to identify the categorical and directional columns in the model,
    demean the categorical columns and organize the directional columns.

    Parameters
    ----------

    filter_data : List of maps. Each map contains data corresponding to a row in the phenotypic file, but only for the columns specified
    in columnsInModel variable in the config_fsl.py

    c : The configuration file object containing all the variables specified in the configuration file.

    Returns
    -------

    filter_data : List of maps. Each map contains data corresponding to a row in the phenotypic file, but only for the columns specified
    in columnsInModel variable in the config_fsl.py. The Directional columns get split according to the number of values they have. 

    field_names : The field names are the column names in the final model file

    """


    mean = {}
    mean_cols = []
    directional_cols = []
    directional_map = {}

    ### line up columns for the model
    for i in range(0, len(c.columnsInModel)):

        if c.deMean[i]:
            mean_cols.append(c.columnsInModel[i])

        # UNCOMMENT THIS AFTER ITS ENABLED
        #if c.categoricalVsDirectional[i]:
        #    directional_cols.append(c.columnsInModel[i])



    for data in filter_data:

        # filter_data is a LIST of dictionaries, one for each subject
        # so 'data' is each dictionary - each one has the phenotypic header
        # items (EVs) as its keys, and the values as its entries

        # there is one dictionary ('data') per subject, and the subject ID
        # is one of the entries in each dictionary

        for col in mean_cols:

            # mean_cols is a LIST of pheno header items (EVs) which are marked
            # to be demeaned (should only be continuous EVs)

            try:

                if not col in mean:
                    # if the current EV-to-be-demeaned does not already have
                    # an entry in the dictionary 'mean', then ensure the
                    # current subject's value for it is a float and put it in
                    mean[col] = float(data[col])

                else:
                    # otherwise, if an entry in 'mean' already exists for the
                    # current EV in the iteration, then simply ensure the
                    # current subject's value for it is a float and then ADD
                    # it to the current sum of this continuous EV's summation
                    mean[col] += float(data[col])

            except ValueError, e:
                print 'error ', e, ' for column: ', col, 'in data row: ', data
                raise


        for col in directional_cols:

            # directional_cols is a LIST of pheno header items (EVs) which are
            # marked as categorical/non-continuous
            # (i.e. are labels, like sex M/F)

            # here, each 'col' is an EV that has been marked as categorical

            # load the phenotype-file VALUE of the current subject's 'col'
            # (categorical EV) into 'val'
            val = data[col]

            # concatenate the name of the categorical EV with its value
            new_col = col + '__' + val


            if not col in directional_map:

                # add this 'categoricalEV__value' into this categorical
                # dictionary, if an entry for this categorical EV has not been
                # created in this 'directional_map' dictionary yet
                directional_map[col] = [new_col]

            else:

                # if categorical values for this EV have already been loaded
                # into 'directional_map', then load what is already in the
                # dictionary into 'val'
                val = directional_map[col]

                # if the current subject's value for this particular
                # categorical EV does not exist in what was already in
                # 'directional_map', then append the current value (the
                # 'categoricalEV__value') into 'val' (now a list?) and then
                # load this list 'val' back into 'directional_map' under the
                # key for this particular categorical EV
                if not new_col in val:

                    val.append(new_col)
                    directional_map[col] = val


    # the 'mean' dictionary now has a key for each EV (pheno header item),
    # with each key corresponding to a value which is the sum of all the
    # subjects' float values for that particular continuous EV

    # the 'directional_map' dictionary now has a key for each EV, with each
    # key corresponding to a LIST which contains all of the subjects'
    # categorical/non-continuous values for that particular categorical EV,
    # and these values are in the 'categoricalEV__value' format



    idx = 0
    for data in filter_data:

        # once again, each 'data' is each subject's collection of EVs and
        # their values

        for col in directional_map.keys():

            # make 'val' the current subject's current categorical EV's value
            val = data[col]

            # make 'vals' a list of all of the values of the current
            # categorical EV
            vals = directional_map[col]

            # remove the entry for the current subject's current categorical
            # EV from the dictionary 'data' (which is all of the EVs of the
            # current subject)
            del data[col]


            for value in vals:

                # 'value' is each categorical EV value (from directional_map)
                # in the 'categoricalEV__value' format

                column_name, v = value.rsplit('__',1)

                # iterating over each categorical EV value from
                # directional_map, check if the current subject's current
                # categorical EV value ('val', which holds steady) matches:
                if v == val:
                    # if it does, add an entry into 'data', the dictionary
                    # which has all of the current subject's EVs and their
                    # values, where the entry's key is the
                    # 'categoricalEV__value', and the value is '1'
                    data[value] = '1'

                else:
                    # if not, do the same as above, except make the value '0'
                    data[value] = '0'

            # now that 'for value in vals:' is complete for the current
            # subject's current categorical EV, the dictionary 'data' has been
            # cleared of all of its original categorical EV entries and new
            # entries have been added for each categorical EV found in
            # 'directional_map' which contain a '1' if this particular subject
            # had an EV with the same VALUE as the VALUE of any of the
            # categorical EVs in 'directional_map'

        filter_data[idx] = data

        idx += 1


    # by now, the dictionary 'filter_data' contains an entry for each subject,
    # with all of the categorical EVs 'broken out' into separate columns -
    # so if a categorical EV like sex can be either 'male' or 'female', it
    # creates a column for both options and places a 1 in the column that
    # corresponds for that subject



    for col in mean.keys():

        mean[col] = float(mean[col])/float(len(filter_data))


    idx = 0
    for data in filter_data:

        # this demeans the continuous EVs - only demeans the EVs that were
        # marked to be demeaned (only they were added to the 'mean'
        # dictionary)
        
        for col in mean.keys():

            val = 0.0
            val = float(data[col])
            val = val - mean[col]

            data[col] = str(val)
        filter_data[idx] = data
        idx += 1


    try:
        zeroth = filter_data[0]
    except:
        print "\n\n" + "ERROR: Subject information did not match properly" \
              " between the phenotypic file and group analysis subject " \
              "list.\n Tip: Double-check subject names." + "\n" + \
        "Error name: create_fsl_model_0001" + "\n\n"
        raise Exception


    field_names = [c.subjectColumn]

    keys = zeroth.keys()

    for ignore, directional_columns in directional_map.items():

        for column_name in sorted(directional_columns):

            field_names.append(column_name)



    for col in sorted(mean.keys()):

        field_names.append(col)

    for f_n in filter_data[0].keys():

        if not (f_n in field_names):
            field_names.append(f_n)



    return filter_data, field_names



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



def write_data(model_data, field_names, c):

    """
    The main purpose of this function is to populate the Model CSV File.
    This file becomes the basis for the FSL group Analysis inputs (.mat, .grp, .con, .fts etc)

    Parameters
    ----------

    model_data : List of maps. Each map contains data corresponding to a row in the phenotypic file, but only for the columns specified
    in columnsInModel variable in the config_fsl.py. The Directional columns get split according to the number of values they have. 

    field_names : The field names are the column names in the final model file

    Returns
    -------

    The populated Model File

    """

    evs = open(c.contrastFile, 'r').readline()
    evs = evs.rstrip('\n')
    evs = evs.rstrip('\n')
    evs = evs.split(',')
    evs = [ev.replace("\"", '') for ev in evs]

    new_evs = []
    

    for ev in evs:
        if (ev in field_names):
            new_evs.append(ev)

    # evs is now a list of contrast file header items
    evs = list(new_evs)
    del new_evs

    # new_field_names is a list of evs with the subjectID
    # column added at the beginning
    new_field_names = [c.subjectColumn] + evs


    try:

        csvPath = c.outputModelFilesDirectory + '/' + c.outputModelFile
        f = open(csvPath, 'wt')

    except:

        print "Could not open the output model file: ", csvPath
        print ""
        raise Exception


    try:

        # make fieldnames=new_field_names so that when the unordered
        # phenotype data is written to the model file, the header items
        # and their corresponding values will be in the correct order
        writer = csv.DictWriter(f, fieldnames=new_field_names)

        header = dict((n, n) for n in new_field_names)

        writer.writerow(header)

        dropped_columns_a = []
        dropped_columns_b = []
        dropped_columns_a = [name for name in field_names if not (name in list(set(field_names) & set(new_field_names))) ]
        dropped_columns_b = [name for name in new_field_names if not (name in list(set(field_names) & set(new_field_names))) ]
        dropped_columns = list(set(dropped_columns_a + dropped_columns_b))


        if not (len(dropped_columns) == 0):
            print 'dropping columns(not specified in contrasts) from the model ', dropped_columns


        new_data = []
        
        # model_data is a LIST of dictionaries of the phenotype
        # header items matched to their values
        for data in model_data:

            data_row = []
            for name in field_names:

                if not name in field_names:
                    del data[name]
                else:
                    if not (c.subjectColumn in name):
                        data_row.append(float(data[name]))

            new_data.append(list(data_row))
         
            writer.writerow(data)
                

        detect = 0

        detect = check_multicollinearity(np.array(new_data))

        if detect == 1:

            print 'Detected Multicollinearity in the computed Model. Please check %s ' % c.outputModelFile
            #sys.exit(0)

    finally:
        f.close()



def create_mat_file(data, model_name, outputModelFilesDirectory):

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



def run(config, fTest, CPAC_run = False):

    if CPAC_run:
        c = config
    else:
        try:
            c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
        except:
            raise Exception("Error in reading %s configuration file" % config)


    create_pheno_dict(c)

    ###This generates the model file

    ###parse the phenotypic file and pickup subjects and phenotypic
    ###columns that user requires to be in the model file with demeaning
    ###and splitting columns for categorical variables


    if c.modelGroupVariancesSeparately == 1 and (c.groupingVariable == None or (not c.groupingVariable in c.columnsInModel)):
        raise ValueError('modelGroupVariancesSeparately is set to 1 but groupingVariable not one of the columns in model')


    try:
        if not os.path.exists(c.outputModelFilesDirectory):

            os.makedirs(c.outputModelFilesDirectory)

    except OSError, e:

        print 'Error: ', e, ' while trying to create outputModelFilesDirectory'
        raise


    filter_data = filter_phenotypic(c)

    dumpFile = open(os.path.join(os.getcwd(), 'dumpFile.txt'), 'wt')
    print >>dumpFile, 'filter_data:\n\n', filter_data

    model_ready_data = None
    field_names = None
    gp_var = None


    if c.modelGroupVariancesSeparately == 0:

        model_ready_data, field_names = organize_data(filter_data, c)

        print >>dumpFile, '\n\nmodel_ready_data:\n\n', model_ready_data
        print >>dumpFile, '\n\nfield_names:\n\n', field_names
        dumpFile.close()
   
    else:

        model_ready_data, field_names, gp_var = alternate_organize_data(filter_data, c)


    write_data(model_ready_data, field_names, c)


    ###generate the final FSL .grp, .mat, .con, .fts files 
    model = c.outputModelFilesDirectory + '/' + c.outputModelFile

    con = c.contrastFile
    model_name = c.modelName

    rdr = csv.DictReader(open(model, "rb"))
    no_of_columns = len(rdr.fieldnames)

    tuple_columns = tuple([n for n in range(1, no_of_columns)])
    data = np.loadtxt(open(model, 'rb'), delimiter=',', skiprows=1, usecols=tuple_columns)
        

    data_lst = data.tolist()


    data = []

    for tp in data_lst:
        data.append(tp)

    print "Length of data list: ", len(data[:])
    print ""
    data = np.array(data, dtype=np.float16)



    try:
        create_mat_file(data, model_name, c.outputModelFilesDirectory)
    except:
        print "Error: Could not create .mat file."
        print ""
        raise Exception

    try:
        create_grp_file(data, model_name, gp_var, c.outputModelFilesDirectory)
    except:
        print "Error: Could not create .grp file."
        print ""
        raise Exception


    create_con_ftst_file(con, model_name, fTest, c.outputModelFilesDirectory)




