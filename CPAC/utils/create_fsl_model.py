import os
import sys
import argparse
import numpy as np
import csv

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

    p_reader = csv.DictReader(file(os.path.abspath(ph)), skipinitialspace=True)


    f = open(sublist, 'r')

    subjects = f.readlines()

    sub_dict = {}
    p_dict = {}


    for subject in subjects:

        subject = subject.rstrip('\r\n')
        if not subject in sub_dict:
            sub_dict[subject] = 1
        else:
            sub_dict[subject] += 1

    final_reader = []
    f_r = []

    record_dict = {}
    for record in p_reader:

        if record[c.subjectColumn] in sub_dict:


            for rec in record.keys():
                if (not (rec in c.columnsInModel)) and not (c.subjectColumn == rec):
                    del record[rec]
            record_dict[record[c.subjectColumn]] = record

    for subject in subjects:
        subject = subject.rstrip('\r\n')
        f_r.append(record_dict[subject])

    #print f_r

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
    lens_grouping_variable = {}
    gp_var = {}
    order = []

    ### line up columns for the model
    for i in range(0, len(c.columnsInModel)):

        if c.deMean[i]:
            mean_cols.append(c.columnsInModel[i])

        if c.categoricalVsDirectional[i]:
            directional_cols.append(c.columnsInModel[i])


    idx = 0
    for data in filter_data:

        if c.modelGroupVariancesSeparately == 1:
            if not (data[c.groupingVariable] in gp_var.keys()):
                gp_var[data[c.groupingVariable]] = [idx]
                order.append(data[c.groupingVariable])
            else:
                val = gp_var[data[c.groupingVariable]]
                val.append(idx)
                gp_var[data[c.groupingVariable]] = val
            idx += 1

        for col in mean_cols:

            try:

                if c.modelGroupVariancesSeparately == 0:
                    if not col in mean:

                        mean[col] = float(data[col])
                    else:
                        mean[col] += float(data[col])

                else:

                    key = col + '#' + str(data[c.groupingVariable])

                    if not key in mean:
                        mean[key] = float(data[col])
                        lens_grouping_variable[key] = 1
                    else:
                        mean[key] += float(data[col])
                        lens_grouping_variable[key] += 1

            except ValueError, e:
                print 'error ', e, ' for data row: ', data
                raise


        for col in directional_cols:

            val = data[col]

            new_col = col + '__' + val

            if not col in directional_map:

                directional_map[col] = [new_col]

            else:

                val = directional_map[col]

                if not new_col in val:

                    val.append(new_col)
                    directional_map[col] = val


    idx = 0
    for data in filter_data:

        for col in directional_map.keys():

            val = data[col]

            vals = directional_map[col]

            del data[col]
            for value in vals:

                column_name, v = value.split('__')

                if v == val:

                    data[value] = '1'

                else:
                    data[value] = '0'

        filter_data[idx] = data

        idx += 1


    for col in mean.keys():

        if c.modelGroupVariancesSeparately == 1:
            mean[col] = float(mean[col])/float(lens_grouping_variable[col])
        else:
            mean[col] = float(mean[col])/float(len(filter_data))
    idx = 0
    for data in filter_data:

        #print data
        for col in mean.keys():

            val = 0.0
            if c.modelGroupVariancesSeparately == 1:
                val = float(data[col.split('#')[0]])
            else:
                val = float(data[col])
            val = val - mean[col]

            if c.modelGroupVariancesSeparately == 1:
                if data[c.groupingVariable+ '__'+ col.split('#')[1]] == '1':
                    data[col.split('#')[0]] = str(val)
            else:
                data[col] = str(val)
        filter_data[idx] = data
        idx += 1


    zeroth = filter_data[0]

    field_names = [c.subjectColumn]

    keys = zeroth.keys()

    if c.groupingVariable in directional_map:
        gp_variable_vals = sorted(directional_map[(c.groupingVariable)])

        for val in gp_variable_vals:

            field_names.append(val)

        del directional_map[(c.groupingVariable)]

    for ignore, directional_columns in directional_map.items():

        for column_name in sorted(directional_columns):

            field_names.append(column_name)



    for col in sorted(mean.keys()):

        if c.modelGroupVariancesSeparately == 1:
            field_names.append(col.split('#')[0])
        else:
            field_names.append(col)

    for f_n in filter_data[0].keys():

        if not (f_n in field_names):
            field_names.append(f_n)


    return filter_data, field_names, gp_var, order


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
    #print model_data
    #print field_names

    evs = open(c.contrastFile, 'r').readline()
    evs = evs.rstrip('\r\n').split(',')

    idx = 0
    new_evs = []
    for ev in evs:
        if (ev in field_names):
            new_evs.append(ev)

    evs = list(new_evs)
    del new_evs

    new_field_names = [c.subjectColumn] + evs


    #print model_data[0]
    #print new_field_names

    f = open(c.outputModelFile, 'wt')


    try:
        writer = csv.DictWriter(f, fieldnames=new_field_names)

        header = dict((n, n) for n in new_field_names)

        writer.writerow(header)

        dropped_columns = []
        dropped_columns = [name for name in field_names if not (name in list(set(field_names) & set(new_field_names))) ]

        if not (len(dropped_columns) == 0):
            print 'dropping columns(not specified in contrasts) from the model ', dropped_columns

        new_data = []
        for data in model_data:

            data_row = []
            for name in field_names:

                if not name in new_field_names:
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

def create_grp_file(data, model_name, gp_var, order, outputModelFilesDirectory):

    """
    create the grp file
    """

    dimx, dimy = data.shape
    data = np.ones(dimx)


    if not (gp_var == []):
        i = 1
        for key in order:

            for index in gp_var[key]:
                data[index] = i

            i += 1


    f = open(os.path.join(outputModelFilesDirectory, model_name + '.grp'), 'w')

    print >>f, '/NumWaves\t1'
    print >>f, '/NumPoints\t%d\n' %dimx
    print >>f, '/Matrix'
    np.savetxt(f, data, fmt='%d', delimiter='\t')

    f.close()

def create_con_ftst_file(con_file, model_name, outputModelFilesDirectory):

    """
    Create the contrasts and fts file
    """
    evs = open(con_file, 'r').readline()
    evs = evs.rstrip('\r\n').split(',')
    count_ftests = 0

    for ev in evs:

        if 'f_test' in ev.lower():

            count_ftests += 1



    data = np.genfromtxt(con_file, names=True, delimiter=',', dtype=None)

    lst = data.tolist()

    ftst = []
    contrasts = []
    contrast_names = []

    length = None
    length = len(list(lst[0]))

    for tp in lst:

        contrast_names.append(tp[0])
        contrasts.append(list(tp)[1:length-count_ftests])

        ftst.append(list(tp[length-count_ftests: length]))


    contrasts = np.array(contrasts, dtype=np.float16)

    fts_n = np.array(ftst)
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

    fts_n = fts_n.T
    f = open(os.path.join(outputModelFilesDirectory, model_name + '.fts'), 'w')
    print >>f, '/NumWaves\t%d' % (contrasts.shape)[0]
    print >>f, '/NumContrasts\t%d\n' % count_ftests

    print >>f, '/Matrix'

    for i in range(fts_n.shape[0]):
        print >>f, ' '.join(fts_n[i].astype('str'))
    #np.savetxt(f, fts_n[None], fmt='%d', delimiter=' ')
    f.close()

"""
Class to set dictionary keys as map attributes
"""
class Configuration(object):
    def __init__(self, config_map):
        for key in config_map:
            if config_map[key] == 'None':
                config_map[key] = None
            setattr(self, key, config_map[key])

def run(config):

    try:
        c = Configuration(yaml.load(open(os.path.realpath(config), 'r')))
    except:
        raise Exception("Error in reading %s configuration file"%config) 

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
    model_ready_data, field_names, gp_var, order = organize_data(filter_data, c)
    write_data(model_ready_data, field_names, c)


    ###generate the final FSL .grp, .mat, .con, .fts files 
    import csv
    model = c.outputModelFile

    con = c.contrastFile
    model_name = c.modelName

    rdr = csv.DictReader(open(model, "rb"))
    no_of_columns = len(rdr.fieldnames)

    tuple_columns = tuple([n for n in range(1, no_of_columns)])
    data = np.loadtxt(open(model, 'rb'), delimiter=',', skiprows=1, usecols=tuple_columns)

    data_lst = data.tolist()

    data = []

    for tp in data_lst:
        data.append(tp[:])

    print len(data[:])
    data = np.array(data, dtype=np.float16)
    create_mat_file(data, model_name, c.outputModelFilesDirectory)
    create_grp_file(data, model_name, gp_var, order, c.outputModelFilesDirectory)
    create_con_ftst_file(con, model_name, c.outputModelFilesDirectory)