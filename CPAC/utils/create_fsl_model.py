import os
import sys
import argparse
import numpy as np
import csv
from collections import OrderedDict as od

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

    final_reader = []
    f_r = []

    for record in p_reader:

        if record['sub'] in sub_dict:

            for rec in record.keys():
                if (not (rec in c.columnsInModel)) and not ('sub' == rec):
                    del record[rec]
            f_r.append(record)


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

    for i in range(0, len(c.columnsInModel)):

        if c.deMean[i]:
            mean_cols.append(c.columnsInModel[i])

        if c.categoricalVsDirectional[i]:
            directional_cols.append(c.columnsInModel[i])


    for data in filter_data:

        for col in mean_cols:

            if not col in mean:

                mean[col] = float(data[col])

            else:
                mean[col] += float(data[col])

        for col in directional_cols:

            val = data[col]

            new_col = col + '#' + val

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

                column_name, v = value.split('#')

                if v == val:

                    data[value] = '1'

                else:
                    data[value] = '0'

        filter_data[idx] = data

        idx += 1


    for col in mean.keys():

        mean[col] = float(mean[col])/float(len(filter_data))

    idx = 0
    for data in filter_data:

        for col in mean_cols:

            val = float(data[col])
            val = val - mean[col]

            data[col] = str(val)
        filter_data[idx] = data
        idx += 1


    zeroth = filter_data[0]

    field_names = ['sub']

    keys = zeroth.keys()

    gp_variable_vals = sorted(directional_map[(c.groupingVariable)[0]])

    for val in gp_variable_vals:

        field_names.append(val)

    del directional_map[(c.groupingVariable)[0]]

    for ignore, directional_columns in directional_map.items():

        for column_name in sorted(directional_columns):

            field_names.append(column_name)



    for col in sorted(mean.keys()):

        field_names.append(col)

    for f_n in filter_data[0].keys():

        if not (f_n in field_names):
            field_names.append(f_n)


    return filter_data, field_names



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

    print model_data
    print field_names
    f = open(c.outputModelFile, 'wt')


    try:
        writer = csv.DictWriter(f, fieldnames=field_names)

        header = dict((n, n) for n in field_names)

        writer.writerow(header)


        for data in model_data:

            writer.writerow(data)

    finally:
        f.close()


def create_mat_file(data, model_name):

    """
    create the .mat file
    """

    dimx, dimy = data.shape

    ppstring = '/PPheights'

    for i in range(0, dimy):

        ppstring += '    ' + '%1.5e' %(1.0)

    ppstring += '\n'

    f = open(model_name + '.mat', 'a')

    print >>f, '/NumWaves    %d' %dimy
    print >>f, '/NumPoints    %d' %dimx
    print >>f, ppstring



    print >>f, '/Matrix'
    np.savetxt(f, data, fmt='%1.5e', delimiter='    ')

    f.close()

def create_grp_file(data, model_name):

    """
    create the grp file
    """

    dimx, dimy = data.shape
    data = np.ones(dimx)

    f = open(model_name + '.grp', 'a')

    print >>f, '/NumWaves    1'
    print >>f, '/NumPoints    %d\n' %dimx
    print >>f, '/Matrix'
    np.savetxt(f, data, fmt='%d', delimiter='    ')

    f.close()

def create_con_ftst_file(con_file, model_name):

    """
    Create the contrasts and fts file
    """
    data = np.genfromtxt(con_file, names=True, delimiter=',', dtype=None)

    lst = data.tolist()

    ftst = []
    contrasts = []
    contrast_names = []

    length = None
    length = len(list(lst[0]))

    for tp in lst:

        contrast_names.append(tp[0])
        contrasts.append(list(tp)[1:length-1])
        ftst.append((tp[length-1:length])[0])


    contrasts = np.array(contrasts, dtype=np.float16)

    fts_n = np.array(ftst, dtype=np.int)
    f = open(model_name + '.con', 'a')

    idx = 1
    pp_str = '/PPheights'
    re_str = '/RequiredEffect'
    for name in contrast_names:

        print >>f, '/ContrastName%d' %idx, '    ', name
        pp_str += '    %1.5e' %(1)
        re_str += '    %1.5e' %(1)
        idx += 1


    print >>f, '/NumWaves    ', (contrasts.shape)[1]
    print >>f, '/NumContrasts    ', (contrasts.shape)[0]
    print >>f, pp_str
    print >>f, re_str + '\n'
    print >>f, '/Matrix'


    np.savetxt(f, contrasts, fmt='%1.5e', delimiter='    ')

    f.close()


    f = open(model_name + '.fts', 'a')
    print >>f, '/NumWaves    ', (contrasts.shape)[0]
    print >>f, '/NumContrasts    1\n'

    print >>f, '/Matrix'

    np.savetxt(f, fts_n[None], fmt='%d', delimiter=' ')
    f.close()



def run(config, mode):

    path, fname = os.path.split(os.path.realpath(config))
    sys.path.append(path)
    c = __import__(fname.split('.')[0])


    if mode == 'generate_model_csv':
        ###This generates the model file

        ###parse the phenotypic file and pickup subjects and phenotypic
        ###columns that user requires to be in the model file with demeaning
        ###and splitting columns for categorical variables
        filter_data = filter_phenotypic(c)
        model_ready_data, field_names = organize_data(filter_data, c)
        write_data(model_ready_data, field_names, c)

    else:
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
        create_mat_file(data, model_name)
        create_grp_file(data, model_name)
        create_con_ftst_file(con, model_name)

