import os

def first_pass_organizing_files(qc_path):

    from CPAC.qc.utils import append_to_files_in_dict_way

    qc_files = os.listdir(qc_path)
    strat_dict = {}

    for file_ in sorted(qc_files, reverse=True):

        if not ('.txt' in file_):
            continue

        file_ = os.path.join(qc_path, file_)
        str_ = os.path.basename(file_)

        str_ = str_.replace('paths_file_', '')

        str_ = str_.replace('scan_', '')

#        str_ = str_.replace('_', '')
        str_ = str_.replace('.txt', '')


        str_ = str_.replace('____', '_')
        str_ = str_.replace('___', '_')
        str_ = str_.replace('__', '_')

        if '_hp_' in str_ and '_fwhm_' in str_ and not ('_bandpass_freqs_' in str_):

            str_, fwhm_val = str_.split('_fwhm_')

            fwhm_val = '_fwhm_' + fwhm_val

            str_, hp_lp_ = str_.split('_hp_')
            hp_lp_ = '_hp_'  + hp_lp_

            str_ = str_ + fwhm_val + hp_lp_

        if strat_dict.keys() == []:

            strat_dict[str_] = [file_]

        else:

            flag_ = 0
            for key_ in strat_dict.keys():

                if str_ in key_:
                    append_to_files_in_dict_way(strat_dict[key_], file_)
                    flag_ = 1

            if flag_ == 1:

                os.system('rm -f %s' % file_)

            else:

                strat_dict[str_] = [file_]


def second_pass_organizing_files(qc_path):

    from CPAC.qc.utils import append_to_files_in_dict_way

    qc_files = os.listdir(qc_path)

    strat_dict = {}

    for file_ in sorted(qc_files, reverse=True):

        if not ('.txt' in file_):
            continue
        str_ = file_
        file_ = os.path.join(qc_path, file_)

        str_ = str_.replace('paths_file_', '')
        str_ = str_.replace('scan_', '')
        str_ = str_.replace('.txt', '')
        str_ = str_.replace('____', '_')
        str_ = str_.replace('___', '_')
        str_ = str_.replace('__', '_')
        fwhm_val_ = ''

        #organize all derivatives excluding alff falff
        if '_bandpass_freqs_' in str_:

            if not str_ in strat_dict:
                strat_dict[str_] = [file_]
            else:
                print 'Error: duplicate keys for files in QC 2nd file_org pass: %s %s' % (strat_dict[str_], file_)
                raise

        #organize alff falff
        elif ('_hp_' in str_) and ('_lp_' in str_):


            key_ = ''
            key_1 = ''
            hp_lp_ = ''
            if '_fwhm_' in str_:
                key_1 = ''
                key_, hp_lp_ = str_.split('_hp_')
                ignore, fwhm_val_ = hp_lp_.split('_fwhm_')
                hp_lp_ = '_hp_' + ignore
                key_1 = '_fwhm_' + fwhm_val_
            else:
                key_, hp_lp_ = str_.split('_hp_')
                hp_lp_ = '_hp_' + hp_lp_


            flag_ = 0
            for key in strat_dict.keys():

                if (key_ in key) and (key_1 in key):

                    append_to_files_in_dict_way(strat_dict[key], file_)
                    str_ = strat_dict[key][0].replace('.txt', '')
                    new_fname = str_ + hp_lp_ + '.txt'

                    os.system('mv %s %s' %(strat_dict[key][0], new_fname))
                    del strat_dict[key]
                    flag_ = 1


                if flag_ == 1:
                    os.system('rm -f %s' % file_)

        else:
            if not str_ in strat_dict:
                strat_dict[str_] = [file_]
            else:
                print 'Error: duplicate keys for files in QC 2nd file_org pass: %s %s' % (strat_dict[str_], file_)
                raise



def populate_htmls(gp_html, sub_html, subj,
                subj_path,\
                meanFD, \
                meanDVARS, \
                mean_rms,\
                max_rms,\
                snr_val):

    f_ = None
    
        
    
    if not os.path.isfile(gp_html):

        f_ = open(gp_html, 'w')

        print >> f_, '<html>'
        print >> f_, '<head>'
        print >> f_, '<title>QC GROUP</title>'
        print >> f_, '</head>'
        print >> f_, '<h2><b>Group Comparison Page</b></h2>'
        print >> f_, '<form name="QC Group Form\" method=\"POST\">'
        print >> f_, "<table border=\"1\" cellspacing=\"1\" cellpadding=\"5\">"
        print >> f_, "    <tr>"
        print >> f_, "        <td>Subject</td>"
        print >> f_, "        <td>MEAN_FD</td>"
        print >> f_, "        <td>MEAN_DVARS</td>"
        print >> f_, "        <td>Mean_Relative_RMS_Displacement</td>"
        print >> f_, "        <td>Max_Relative_RMS_Displacement</td>"
        print >> f_, "        <td>Mean Functional SNR</td>"   ###
#         print >> f_, "        <td>Include in Group Analysis</td>"
#         print >> f_, "        <td>Comments</td>"
        print >> f_, "    </tr>"

    else:

        f_ = open(gp_html, 'a')

    print >> f_, "        <td><a href='%s'>%s</a></td>" % (sub_html, subj)
    print >> f_, "        <td>%s</td>" % (meanFD)
    print >> f_, "        <td>%s</td>" % (meanDVARS)
    print >> f_, "        <td>%s</td>" % (mean_rms)
    print >> f_, "        <td>%s</td>" % (max_rms)
    print >> f_, "        <td>%s</td>" % (snr_val)

#     print >> f_, "        <td><input type=\"checkbox\" name=\"gp_inc\" value=\"Y\">Y<br></td>"
# #    print >> f_, "        <td><textarea name=\"comments\" cols=\"3\" rows=\"1\"><br></td>"
#     print >> f_, "        <td>No Comments Yet!</td>"
    print >> f_, "    </tr>"

    f_.close()


def prep_resources(pip_path):

    #for pipeline in sorted(pipelines):

    #    pip_path = os.path.join(output_path, pipeline)


    subjects = os.listdir(pip_path)
    subjects = [subj for subj in subjects if os.path.isdir(os.path.join(pip_path, subj))]

    for subj in sorted(subjects):

        subj_path = os.path.join(pip_path, subj)

        if os.path.isdir(subj_path):
            qc_path = os.path.join(subj_path, 'qc_files_here')
            path_files_here = os.path.join(subj_path, 'path_files_here')
            path_copy = os.path.join(qc_path, 'path')

            if not os.path.exists(path_copy):
                os.makedirs(path_copy)

            os.system('cp %s/*.txt %s/' % (path_files_here, path_copy))
            first_pass_organizing_files(path_copy)
            second_pass_organizing_files(path_copy)

            #os.system('ls %s' %(path_copy))

            print 'mv %s/*.txt %s/' % (path_copy, qc_path)
            os.system('mv %s/*.txt %s/' % (path_copy, qc_path))
            os.system('rm -rf %s' %(path_copy))


def get_path_files_in_dict(qc_path, html_files):

    dict_ = {}
    files_ = os.listdir(qc_path)

    files_ = [file_ for file_ in files_ if 'paths_file' in file_]

    for html_ in html_files:

        str_html = html_
        str_html = str_html.split('.html')[0]
        str_html = str_html.replace('qc_', '')
        str_html = str_html.replace('__', '_')
        str_html = str_html.replace('___', '_')

        for file_ in files_:

            str_ = file_
            str_ = str_.split('paths_file_')[1]
            str_ = str_.replace('_.txt', '')
            str_ = str_.replace('.txt', '')
            str_ = str_.replace('__', '_')
            str_ = str_.replace('___', '_')


            if str_ in str_html:
                
                if not (str_html in dict_):
                    dict_[str_html] = [os.path.join(qc_path, file_)]
                else:

                    val = dict_[str_html]
                    val.append(os.path.join(qc_path, file_))
                    dict_[str_html]


    return dict(dict_)



def get_resources_for_strategy(files_, resources, special_res, dict_):


    things_to_replace = ['sca_tempreg_z_maps_', 'temp_reg_map_', \
                         'z_score_', '_warp', '_maths', \
                        'temp_reg_map_z_']

    for file_ in files_:

        f = open(file_, 'r')

        lines = f.readlines()

        lines = [line.rstrip('\r\n') for line in lines]

        for line in lines:

            for resource in resources:

                if resource in line:

                    if not resource in dict_:

                        dict_[resource] = [line]

                    else:
                        val = list(dict_[resource])
                        dict_[resource] = val.append(line)


            for res in special_res:

                if res in line:

                    mask_ = ''
                    map_ = ''

                    if not 'sca_seed_Z_to_standard':
                        mask_ = os.path.basename(os.path.dirname(line))
                    else:
                        mask_ = os.path.basename(os.path.dirname(os.path.dirname(line)))


                    mask_ = mask_.replace('_mask_', '')
                    mask_ = mask_.replace('_roi_', '')
                    mask_ = mask_.replace('_spatial_map_', '')
                    map_ = os.path.splitext(os.path.splitext(os.path.basename(line))[0])[0]

                    for thing in things_to_replace:

                        map_ = map_.replace(thing, '')



                    key_ = res + '--' + mask_ + '--' + map_


                    if not key_ in dict_:
                        dict_[key_] = [line]

                    else:
                        val = list(dict_[key_])
                        dict_[key_] = val.append(line)







# def organize_resources(output_path, pipelines):
# 
#     from CPAC.utils.create_all_qc import get_path_files_in_dict
#     from CPAC.utils.create_all_qc import get_resources_for_strategy
# 
#     resources = ['/alff_Z_to_standard/', '/falff_Z_to_standard/', \
#                 '/reho_Z_to_standard/', \
#                 '/vmhc_z_score_stat_map/']
# 
#     special_res = ['/sca_roi_Z_to_standard/', '/sca_seed_Z_to_standard/', \
#                    '/dr_tempreg_maps_z_files/', '/sca_tempreg_maps_z_files/', \
#                    '/centrality_outputs_zscore/']
# 
# 
#     resource_dict = {}
#     for pipeline in sorted(pipelines):
# 
#         pip_path = os.path.join(output_path, pipeline)
#         subjects = os.listdir(pip_path)
#         subjects = [subj for subj in subjects if os.path.isdir(os.path.join(pip_path, subj))]
# 
#         for subj in sorted(subjects):
# 
#             subj_path = os.path.join(pip_path, subj)
#             qc_path = os.path.join(subj_path, 'qc_files_here')
#             html_files = os.listdir(qc_path)
#             html_files = [html for html in html_files if not (html.endswith('_0.html') or html.endswith('_1.html') or html.endswith('.txt'))]
#             dict_ = get_path_files_in_dict(qc_path, html_files)
# 
#             for strat, resource_files in dict_.items():
# 
#                 if not strat in resource_dict:
# 
#                     per_strat_dict = {}
#                     per_strat_dict = get_resources_for_strategy(resource_files, resources, special_res, per_strat_dict)
#                     resource_dict[strat] = dict(per_strat_dict)
# 
#                 else:
# 
#                     per_strat_dict = dict(resource_dict[strat])
#                     per_strat_dict = get_resources_for_strategy(resource_files, resources, special_res, per_strat_dict)
#                     resource_dict[strat] = dict(per_strat_dict)
# 
# 
#     return resource_dict


def get_power_params(qc_path, file_):

    import csv

    subj_dir = os.path.dirname(qc_path)

    subj_dir = os.path.join(subj_dir, 'power_params')

    scans = os.listdir(subj_dir)

    meanFD = None
    meanDVARS = None

    for scan in scans:

        if scan in file_:

            subj_dir = os.path.join(subj_dir, scan)

            threshold = '_threshold_'

            if 'SCRUB_' in file_:
                threshold += (file_.split('SCRUB_')[1]).split('_')[0]
                subj_dir = os.path.join(subj_dir, threshold)
            else:

                subj_dir = os.path.join(subj_dir, os.listdir(subj_dir)[0])


            params_file = os.path.join(subj_dir, os.listdir(subj_dir)[0])
            csv_file = csv.DictReader(open(params_file, 'rb'), delimiter=',')
            print 'params_file: ', params_file
            print csv.list_dialects()

            line = None

            for line in csv_file:

                meanFD = line['MeanFD']
                meanDVARS = line['MeanDVARS']
                print meanFD, meanDVARS
                
                return meanFD, meanDVARS


def get_motion_params(qc_path, file_):

    import csv

    subj_dir = os.path.dirname(qc_path)

    subj_dir = os.path.join(subj_dir, 'motion_params')

    scans = os.listdir(subj_dir)

    mean_RMS = None
    max_RMS = None

    for scan in scans:

        if scan in file_:

            subj_dir = os.path.join(subj_dir, scan)


            params_file = os.path.join(subj_dir, os.listdir(subj_dir)[0])
            csv_file = csv.DictReader(open(params_file, 'rb'), delimiter=',')

            line = None

            for line in csv_file:

                mean_RMS = line['Mean_Relative_RMS_Displacement']
                max_RMS = line['Max_Relative_RMS_Displacement']

                return mean_RMS, max_RMS


def write_closing_tags(file_):

    f_ = open(file_, 'a')
    print >> f_, "</table>"
    print >> f_, "<br>"
    print >> f_, "<br>"
    #print >> f_, "<input type=\"submit\" value=\"create group subject list\" />"
    print >> f_, "</html>"

    f_.close()


def make_group_htmls(pip_path):

    import os

    files_ = []
    #pipelines = os.listdir(output_path)
    #pipelines = [pipeline for pipeline in pipelines if 'pipeline_' in pipeline]

    #pipelines = [pipeline for pipeline in pipelines if os.path.isdir(os.path.join(output_path, pipeline))]   ###
    prep_resources(pip_path)
#    organize_resources(output_path, pipelines)


    #for pipeline in sorted(pipelines):
    if os.path.exists(pip_path):
        #pip_path = os.path.join(output_path, pipeline)

        os.system('rm -f %s/qc_*.html'  % pip_path)
        subjects = os.listdir(pip_path)
        subjects = [subj for subj in subjects if os.path.isdir(os.path.join(pip_path, subj))]

        for subj in sorted(subjects):

            subj_path = os.path.join(pip_path, subj)
            

            print "subj_path: ", subj_path, '\n'
            
            snr_file = ''
            sval = ''
            
            
            ### read average snr value from the file
            for root, dirs, files in os.walk(subj_path + '/qc/snr_val'):
                if 'average_snr_file.txt' in files:
                    snr_file = os.path.join(root + '/average_snr_file.txt')


            if os.path.exists(snr_file):
                sval = open(snr_file, 'r').readline()
                sval = '%.2f' % float(sval)
                print 'snr value: ', sval


            if os.path.isdir(subj_path):
                qc_path = os.path.join(subj_path, 'qc_files_here')
                path_files_here = os.path.join(subj_path, 'path_files_here')
                path_copy = os.path.join(qc_path, 'path')
                os.makedirs(path_copy)
                os.system('cp %s/*.txt %s/' % (path_files_here, path_copy))
                first_pass_organizing_files(path_copy)
                second_pass_organizing_files(path_copy)
                ###os.system('ls %s' %(path_copy))
                print 'mv %s/*.txt %s/' % (path_copy, qc_path)
                os.system('mv %s/*.txt %s/' % (path_copy, qc_path))
                os.system('rm -rf %s' %(path_copy))
                html_files = os.listdir(qc_path)
                html_files = [html for html in html_files if not (html.endswith('_0.html') or html.endswith('_1.html') or html.endswith('.txt'))]

                for file_ in sorted(html_files):

                    html_ = os.path.join(pip_path, file_)

                    print 'qc_path: ', qc_path
                    print 'file_: ', file_

                    try:
                        meanFD, meanDvars = get_power_params(qc_path, file_)

                    except:
                        print "Error: some power params lost."
                        meanFD = 0.0
                        meanDvars = 0.0
                        pass

                    try:
                        mean_rms, max_rms = get_motion_params(qc_path, file_)
                    except:
                        print "Error: some motion params lost."
                        mean_rms = 0.0
                        max_rms = 0.0
                        pass                        


                    populate_htmls(html_, os.path.join(qc_path, file_), subj, \
                    subj_path, meanFD, meanDvars, mean_rms, max_rms, sval)
                    if not html_ in files_:
                        files_.append(html_)

        for file_ in files_:
            write_closing_tags(file_)





def run(pipeline_path):

    #from CPAC.utils.create_all_qc import make_group_htmls

    make_group_htmls(pipeline_path)

