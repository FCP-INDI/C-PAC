import commands
import numpy as np
import matplotlib
import pkg_resources as p
matplotlib.use('Agg')


def append_to_files_in_dict_way(list_files, file_):

    """
    Combine files so at each resource in file appears exactly once

    Parameters
    ----------

    list_files : list

    file_ : string

    Returns
    -------

    None

    Notes
    -----

    Writes contents of file_ into list_files, ensuring list_files finally has
    each resource appearing exactly once

    """

    f_1 = open(file_, 'r')

    lines = f_1.readlines()

    lines = [line.rstrip('\r\n') for line in lines]

    one_dict = {}

    for line in lines:

        if not line in one_dict:
            one_dict[line] = 1



    f_1.close()

    for f_ in list_files:

        two_dict = {}
        f_2 = open(f_, 'r')
        lines = f_2.readlines()
        f_2.close()
        f_2 = open(f_, 'w')
        lines = [line.rstrip('\r\n') for line in lines]
        for line in lines:

            if not line in one_dict:
                two_dict[line] = 1

        for key in one_dict:
                if not key in two_dict:
                    two_dict[key] = 1

        for key in two_dict:

            print >> f_2, key

        f_2.close



def first_pass_organizing_files(qc_path):


    """
    First Pass at organizing qc txt files

    Parameters
    ----------

    qc_path : string
        existing path of qc_files_here directory

    Returns
    -------

    None

    Notes
    -----

    Combines files with same strategy. First pass combines file names,
    where one file name is substring of the other.

    """

    import os
    from CPAC.qc.utils import append_to_files_in_dict_way

    qc_files = os.listdir(qc_path)
    strat_dict = {}

    for file_ in sorted(qc_files, reverse=True):

        if not ('.txt' in file_):
            continue

        file_ = os.path.join(qc_path, file_)
        str_ = os.path.basename(file_)

        str_ = str_.replace('qc_', '')

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

                print str_, ' ~~ ', key_
                if str_ in key_:
                    append_to_files_in_dict_way(strat_dict[key_], file_)
                    flag_ = 1

            if flag_ == 1:

                os.system('rm -f %s' % file_)

            else:

                strat_dict[str_] = [file_]



def second_pass_organizing_files(qc_path):


    """
    Second Pass at organizing qc txt files

    Parameters
    ----------

    qc_path : string
        existing path of qc_files_here directory

    Returns
    -------

    None

    Notes
    -----

    Combines files with same strategy. combines files for derivative 
    falff , alff with others

    """

    import os
    from CPAC.qc.utils import append_to_files_in_dict_way

    qc_files = os.listdir(qc_path)

    strat_dict = {}
    got_hp_lp = 0
    got_bp = 0
    for file_ in sorted(qc_files, reverse=True):

        if not ('.txt' in file_):
            continue
        str_ = file_
        file_ = os.path.join(qc_path, file_)

        str_ = str_.replace('qc_scan_', '')
        str_ = str_.replace('.txt', '')
        str_ = str_.replace('____', '_')
        str_ = str_.replace('___', '_')
        str_ = str_.replace('__', '_')
        print '~~>', str_
        fwhm_val_ = ''

        #organize all derivatives excluding alff falff
        if '_bandpass_freqs_' in str_:

            if not str_ in strat_dict:
                strat_dict[str_] = [file_]
            else:
                print 'Error: duplicate keys for files in QC 2nd file_org pass: %s %s' % (strat_dict[str_], file_)
                raise

            print strat_dict
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

                print key, '~~', key_, '~~', key_1
                if (key_ in key) and (key_1 in key):

                    append_to_files_in_dict_way(strat_dict[key], file_)
                    str_ = strat_dict[key][0].replace('.txt', '')
                    new_fname = str_ + hp_lp_ + '.txt'

                    print '~~>', new_fname
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



def organize(dict_, all_ids, png_, new_dict):

    """
    Organizes pngs according to their IDS in new_dict dictionary

    Parameters
    ----------

    dict_ : dictionary
        dict containing png id no and png type(montage/plot/hist)

    all_ids : list
        list of all png id numbers

    png_ : string
        path to png

    new_dict : dictionary
        dictionary containg ids and png lists


    Returns
    -------

    all_ids : list
        list of png id nos
    """

    for id_no, png_type in dict_.items():

        if png_type in png_:
            if not id_no in new_dict.keys():
                new_dict[id_no] = [png_]
            else:
                list_ = new_dict[id_no]
                list_.append(png_)
                new_dict[id_no] = list(list_)

            if not id_no in all_ids:
                all_ids.append(id_no)

    return  all_ids



def grp_pngs_by_id(pngs_, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id):

    """
    Groups pngs by their ids

    Parameters
    ----------

    pngs_ : list
        list of all pngs

    qc_montage_id_a : dictionary
        dictionary of axial montages key : id no
        value is list of png types 

    qc_montage_id_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of png types 

    qc_plot_id : dictionary
          dictionary of plot pngs key : id no
          value is list of png types

    qc_hist_id : dictionary
          dictionary of histogram pngs key : id no
          value is list of png types


    Returns
    -------
    dict_a : dictionary
        dictionary of axial montages key : id no
        value is list of paths to axial montages

    dict_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of paths to sagittal montages

    dict_plot : dictionary
          dictionary of plot pngs key : id no
          value is list of paths to plots

    dict_hist : dictionary
          dictionary of histogram pngs key : id no
          value is list of paths to histogram pngs

    all_ids : list
        list of png id nos
    """

    from CPAC.qc.utils import organize

    dict_a = {}
    dict_s = {}
    dict_hist = {}
    dict_plot = {}

    all_ids = []
    for png_ in pngs_:


        all_ids = organize(qc_montage_id_a, all_ids, png_, dict_a)
        all_ids = organize(qc_montage_id_s, all_ids, png_, dict_s)
        all_ids = organize(qc_plot_id, all_ids, png_, dict_plot)
        all_ids = organize(qc_hist_id, all_ids, png_, dict_hist)

    return dict(dict_a), dict(dict_s), dict(dict_hist), dict(dict_plot), list(all_ids)



def add_head(f_html_, f_html_0, f_html_1):

    """
    Write HTML Headers to various html files

    Parameters
    ----------

    f_html_ : string
        path to main html file

    f_html_0 : string
        path to navigation bar html file

    f_html_1 : string
        path to html file contaning pngs and plots


    Returns
    -------

    None

    """
    print >>f_html_, "<html>"
    print >>f_html_, "<head>"
    print >>f_html_, "<title>QC</title>"
    print >>f_html_, "</head>"
    print >>f_html_, ""
    print >>f_html_, "<frameset cols=\"20%,80%\">"
    print >>f_html_, ""
    print >>f_html_, "    <frame src=\"%s\" name=\"menu\"><frame src=\"%s\" name=\"content\"></frameset>" %(f_html_0.name, f_html_1.name)
    print >>f_html_, ""
    print >>f_html_, "</html>"

    print >>f_html_0, "<html>"
    print >>f_html_0, "<link href=\"%s\" rel=\"stylesheet\" media=\"screen\">"%(p.resource_filename('CPAC',"GUI/resources/html/_static/nature.css"))
    print >>f_html_0, "<link href=\"%s\" rel=\"stylesheet\" media=\"screen\">"%(p.resource_filename('CPAC',"GUI/resources/html/_static/pygments.css"))
    print >>f_html_0, "<head>"
    print >>f_html_0, "<base target=\"content\">"
    print >>f_html_0, "</head>"
    print >>f_html_0, "<body>"
    print >>f_html_0, "<div>"
    print >>f_html_0, "<div class=\"sphinxsidebarwrapper\">"
    print >>f_html_0, "<p class=\"logo\"><a href=\"#\">"
    print >>f_html_0, "<img class=\"logo\" src=\"%s\" alt=\"Logo\"/>"%(p.resource_filename('CPAC',"GUI/resources/html/_static/cpac_logo.jpg"))
    print >>f_html_0, "</a></p>"
    print >>f_html_0, "<h3><a href=\"#\">Table Of Contents</a></h3>"
    print >>f_html_0, "<ul>"

    print >>f_html_1, '<link href="default.css" rel="stylesheet" type="text/css" />'
    print >>f_html_1, "<html>"
    print >>f_html_1, "</style>"
    print >>f_html_1, "<body>"
    print >>f_html_1, "<a name='reverse'>"



def add_tail(f_html_, f_html_0, f_html_1):


    """
    Write HTML Tail Tags to various html files

    Parameters
    ----------

    f_html_ : string
        path to main html file

    f_html_0 : string
        path to navigation bar html file

    f_html_1 : string
        path to html file contaning pngs and plots


    Returns
    -------

    None

    """
#     print >>f_html_0, "</ul>"
#     print >>f_html_0, "<br>"
#     print >>f_html_0, "<br>"
#     print >>f_html_0, "<form class=\"navbar-form pull-left\" method=\"post\" action=\"\"> \
# <br> \
# Recommendation for QC <br>\
# <input type=\"radio\" name=\"group1\" value=\"1\"> 1<br> \
# <input type=\"radio\" name=\"group1\" value=\"2\"> 2<br> \
# <input type=\"radio\" name=\"group1\" value=\"3\"> 3<br> \
# <input type=\"radio\" name=\"group1\" value=\"4\"> 4<br> \
# <input type=\"radio\" name=\"group1\" value=\"5\" checked>5<br> \
#  <br> Notes\
# <br> \
# <textarea name=\"comments\" cols=\"25\" rows=\"5\"> \
# Enter QC Notes here... \
# </textarea><br> \
# <br> \
# <input type=\"submit\" value=\"Submit\" /> \
# <br> \
# <br> \
# <br> \
# </form>"
    print >>f_html_0, "</ul>"
    print >>f_html_0, "</div>"
    print >>f_html_0, "</div>"
    print >>f_html_0, "</body>"
    print >>f_html_0, "</html>"
    print >>f_html_1, "</body>"
    print >>f_html_1, "</html>"



def feed_line_nav(id_,
               image_name,
               anchor,
               f_html_0,
               f_html_1):

    """
    Write to navigation bar html file

    Parameters
    ----------

    id_ : string
        id of the image

    anchor : string
        anchor id of the image

    image_name : string
        name of image
    
    f_html_0 : string
        path to navigation bar html file

    f_html_1 : string
        path to html file contaning pngs and plots


    Returns
    -------

    None

    """

    ### add general user readable link names for QC navigation bar
    image_readable = image_name
    if image_name == 'skullstrip_vis':
        image_readable = 'Visual Result of Skull Strip'
    if image_name == 'csf_gm_wm':
        image_readable = 'Grey Matter, White Matter & CSF'
    if image_name == 'snr':
        image_readable = 'Signal to Noise Ratio'
    if image_name.find('snr_hist') > -1:
        image_readable = 'Histogram of Signal to Noise Ratio'
    if image_name.find('mni_normalized') > -1:
        image_readable = 'MNI Edge Overlapped on Normalized Anatomical'
    if image_name == 'mean_func_with_t1_edge':
        image_readable = 'T1 Edge Overlapped on Mean Functional Image'
    if image_name == 'mean_func_with_mni_edge':
        image_readable = 'MNI Edge Overlapped on Mean Functional Image'
    if image_name.find('movement_trans_plot') >-1:
        image_readable = 'Head Displacement Plot'
    if image_name.find('movement_rot_plot') >-1:
        image_readable = 'Head Rotation Plot'
    if image_name.find('fd_plot') > -1:
        image_readable = 'Framewise Displacement Plot'
    if image_name == 'sca_roi_smooth':
        image_readable = 'Seed-based Correlation Analysis'
    if image_name == 'sca_roi_smooth_hist':
        image_readable = 'Histogram of Seed-based Correlation Analysis'
    if image_name == 'centrality_smooth':
        image_readable = 'Network Centrality'
    if image_name == 'centrality_smooth_hist':
        image_readable = 'Histogram of Network Centrality'
    if image_name == 'temporal_dual_regression_smooth':
        image_readable = 'Temporal Dual Regression'
    if image_name == 'temporal_dual_regression_smooth_hist':
        image_readable = 'Histogram of Temporal Dual Regression'
    if image_name == 'vmhc_smooth':
        image_readable = 'Voxel-Mirrored Homotopic Connectivity'
    if image_name == 'vmhc_smooth_hist':
        image_readable = 'Histogram of Voxel-Mirrored Homotopic Connectivity'
    if image_name == 'reho_smooth':
        image_readable = 'Regional Homogeneity'
    if image_name == 'reho_smooth_hist':
        image_readable = 'Histogram of Regional Homogeneity'
    if image_name == 'alff_smooth':
        image_readable = 'Amplitude of Low-Frequency Fluctuation'
    if image_name == 'alff_smooth_hist':
        image_readable = 'Histogram of Amplitude of Low-Frequency Fluctuation'
    if image_name == 'falff_smooth':
        image_readable = 'fractional Amplitude of Low-Frequency Fluctuation'
    if image_name == 'falff_smooth_hist':
        image_readable = 'Histogram of fractional Amplitude of Low-Frequency Fluctuation' 

    print >>f_html_0, "<li><a href='%s#%s'> %s </a></li>" % (f_html_1.name, \
        anchor, image_readable)   ###



def feed_line_body(image_name, anchor, image, f_html_1):
    """
    Write to html file that has to contain images

    Parameters
    ----------

    image_name : string
        name of image

    anchor : string
        anchor id of the image

    image : string
        path to the image

    f_html_1 : string
        path to html file contaning pngs and plots


    Returns
    -------

    None

    """

    ### add general user readable link names for QC body
    image_readable = image_name
    if image_name == 'skullstrip_vis':
        image_readable = 'Visual Result of Skull Strip'
    if image_name == 'csf_gm_wm':
        image_readable = 'Grey Matter, White Matter & CSF'
    if image_name == 'snr':
        image_readable = 'Signal to Noise Ratio'
    if image_name.find('snr_hist') > -1:
        image_readable = 'Histogram of Signal to Noise Ratio'
    if image_name.find('mni_normalized') > -1:
        image_readable = 'MNI Edge Overlapped on Normalized Anatomical'
    if image_name == 'mean_func_with_t1_edge':
        image_readable = 'T1 Edge Overlapped on Mean Functional Image'
    if image_name == 'mean_func_with_mni_edge':
        image_readable = 'MNI Edge Overlapped on Mean Functional Image'
    if image_name.find('movement_trans_plot') >-1:
        image_readable = 'Head Displacement Plot'
    if image_name.find('movement_rot_plot') >-1:
        image_readable = 'Head Rotation Plot'
    if image_name.find('fd_plot') > -1:
        image_readable = 'Framewise Displacement Plot'
    if image_name == 'sca_roi_smooth':
        image_readable = 'Seed-based Correlation Analysis'
    if image_name == 'sca_roi_smooth_hist':
        image_readable = 'Histogram of Seed-based Correlation Analysis'
    if image_name == 'centrality_smooth':
        image_readable = 'Network Centrality'
    if image_name == 'centrality_smooth_hist':
        image_readable = 'Histogram of Network Centrality'
    if image_name == 'temporal_dual_regression_smooth':
        image_readable = 'Temporal Dual Regression'
    if image_name == 'temporal_dual_regression_smooth_hist':
        image_readable = 'Histogram of Temporal Dual Regression'
    if image_name == 'vmhc_smooth':
        image_readable = 'Voxel-Mirrored Homotopic Connectivity'
    if image_name == 'vmhc_smooth_hist':
        image_readable = 'Histogram of Voxel-Mirrored Homotopic Connectivity'
    if image_name == 'reho_smooth':
        image_readable = 'Regional Homogeneity'
    if image_name == 'reho_smooth_hist':
        image_readable = 'Histogram of Regional Homogeneity'
    if image_name == 'alff_smooth':
        image_readable = 'Amplitude of Low-Frequency Fluctuation'
    if image_name == 'alff_smooth_hist':
        image_readable = 'Histogram of Amplitude of Low-Frequency Fluctuation'
    if image_name == 'falff_smooth':
        image_readable = 'fractional Amplitude of Low-Frequency Fluctuation'
    if image_name == 'falff_smooth_hist':
        image_readable = 'Histogram of fractional Amplitude of Low-Frequency Fluctuation' 


    print >>f_html_1, "<h3><a name='%s'>%s</a> <a href='#reverse'>TOP</a></h3>" %(anchor, image_readable)   ###

    ###data_uri = open(image, 'rb').read().encode('base64').replace('\n', '')
    ###img_tag = '<br><img src="data:image/png;base64,{0}">'.format(data_uri)
    
    img_tag = "<br><img src='%s', alt='%s'>" %(image, image_readable)   ###
    print >>f_html_1, img_tag


def get_map_id(str_, id_):

    """
    Returns the proper map name given identifier for it

    Parameters
    ----------

    str_ : string
        string containing text for identifier

    id_ : string
        string for identifier

    Returns
    -------

    map_id : string
        proper name for a map
    """
    map_id = None

    if 'centrality' in id_:
        str_ = str_.split('_centrality_a.png')[0]
        type_, str_ = str_.split(id_)
        str_ = str_.split('_')[0]

        type_ = type_.replace('_', '')
        map_id = '_'.join([type_, id_, str_])
        return map_id

    else:
        str_ = str_.split(id_)[1]
        str_ = str_.split('_')[0]
        map_id = '_'.join([id_, str_])
        return map_id


def get_map_and_measure(png_a):

    """
    Extract Map name and Measure name from png

    Parameters
    ----------

    png_a : string
        name of png

    Returns
    -------

    map_name : string
        proper name for map

    measure_name : string
        proper name for measure    

    """

    import os
    from CPAC.qc.utils import get_map_id

    measure_name = None
    map_name = None
    if '_fwhm_' in png_a:

        measure_name = os.path.basename(os.path.dirname(os.path.dirname(png_a)))

    else:
        measure_name = os.path.basename(os.path.dirname((png_a)))

    str_ = os.path.basename(png_a)

    if 'sca_seeds' in png_a:
        map_name = 'seed'

    if 'sca_roi' in png_a:

        map_name = get_map_id(str_, 'ROI_number_')

    if 'temporal_regression_sca' in png_a:

        map_name = get_map_id(str_, 'roi_')

    if 'temporal_dual_regression' in png_a:

        map_name = get_map_id(str_, 'map_z_')

    if 'centrality' in png_a:

        try:
            map_name = get_map_id(str_, 'centrality_')
        except:
            pass

        try:
            map_name = get_map_id(str_, 'lFCD_')
        except:
            pass


    return map_name, measure_name


def feed_lines_html(id_,
                        dict_a,
                        dict_s,
                        dict_hist,
                        dict_plot,
                        qc_montage_id_a,
                        qc_montage_id_s,
                        qc_plot_id,
                        qc_hist_id,
                        f_html_0,
                        f_html_1):

    """
    Write HTML Tags to various html files and embeds images

    Parameters
    ----------

    dict_a : dictionary
        dictionary of axial montages key : id no
        value is list of paths to axial montages

    dict_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of paths to sagittal montages

    dict_plot : dictionary
          dictionary of plot pngs key : id no
          value is list of paths to plots

    dict_hist : dictionary
          dictionary of histogram pngs key : id no
          value is list of paths to histogram pngs

    qc_montage_id_a : dictionary
        dictionary of axial montages key : id no
        value is list of png types 

    qc_montage_id_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of png types 

    qc_plot_id : dictionary
          dictionary of plot pngs key : id no
          value is list of png types

    qc_hist_id : dictionary
          dictionary of histogram pngs key : id no
          value is list of png types

    f_html_0 : string
        path to navigation bar html file

    f_html_1 : string
        path to html file contaning pngs and plots


    Returns
    -------

    None

    """
    from CPAC.qc.utils import feed_line_nav
    from CPAC.qc.utils import feed_line_body
    from CPAC.qc.utils import get_map_and_measure

    #print 'id_ :', id_
    if id_ in dict_a:

        dict_a[id_] = sorted(dict_a[id_])
        dict_s[id_] = sorted(dict_s[id_])

        if id_ in dict_hist:

            dict_hist[id_] = sorted(dict_hist[id_])

        idxs = len(dict_a[id_])
        for idx in range(0, idxs):

            png_a = dict_a[id_][idx]
            png_s = dict_s[id_][idx]
            png_h = None

            if id_ in dict_hist:
                png_h = dict_hist[id_][idx]

            measure_name = None
            map_name = None

            if idxs > 1:
                map_name, measure_name =get_map_and_measure(png_a)


            id_a = str(id_)
            id_s = str(id_) + '_s'
            id_h = str(id_) + '_' + str(id_)

            image_name_a = None
            image_name_h = None

            image_name_a_nav = qc_montage_id_a[id_].replace('_a', '')
            if id_ in qc_hist_id:
                image_name_h_nav = qc_hist_id[id_]
            if map_name is not None:
                image_name_a = 'Measure: ' + qc_montage_id_a[id_].replace('_a', '') + '    Mask: '+ measure_name + '   Map: ' + map_name
                if id_ in qc_hist_id:
                    image_name_h = 'Measure: ' + qc_hist_id[id_]  + '    Mask:'+ measure_name + '    Map: ' + map_name
            else:
                image_name_a = qc_montage_id_a[id_].replace('_a', '')
                if id_ in qc_hist_id:
                    image_name_h = qc_hist_id[id_]

            if idx != 0:
                id_a = '_'.join([id_a, str(idx), 'a'])
                id_s = '_'.join([id_s, str(idx), 's'])
                id_h = '_'.join([id_h, str(idx), 'h' ])
            

            if idx == 0:

                # add general user readable link names for QC navigation bar
                if image_name_a_nav == 'skullstrip_vis':
                    image_readable = 'Visual Result of Skull Strip'
                if image_name_a_nav == 'csf_gm_wm':
                    image_readable = 'Grey Matter, White Matter & CSF'
                if image_name_a_nav == 'snr':
                    image_readable = 'Signal to Noise Ratio'
                if image_name_a_nav == 'snr_hist':
                    image_readable = 'Histogram of Signal to Noise Ratio'
                if image_name_a_nav == 'mean_func_with_t1_edge':
                    image_readable = 'T1 Edge Overlapped on Mean Functional Image'
                if image_name_a_nav == 'mean_func_with_mni_edge':
                    image_readable = 'MNI Edge Overlapped on Mean Functional Image'
                if image_name_a_nav == 'movement_trans_plot':
                    image_readable = 'Head Displacement Plot'
                if image_name_a_nav == 'movement_rot_plot':
                    image_readable = 'Head Rotation Plot'
                if image_name_a_nav == 'fd_plot':
                    image_readable = 'Framewise Displacement Plot'
                if image_name_a_nav == 'sca_roi_smooth':
                    image_readable = 'Seed-based Correlation Analysis'
                if image_name_a_nav == 'sca_roi_smooth_hist':
                    image_readable = 'Histogram of Seed-based Correlation Analysis'
                if image_name_a_nav == 'centrality_smooth':
                    image_readable = 'Network Centrality'
                if image_name_a_nav == 'centrality_smooth_hist':
                    image_readable = 'Histogram of Network Centrality'
                if image_name_a_nav == 'temporal_dual_regression_smooth':
                    image_readable = 'Temporal Dual Regression'
                if image_name_a_nav == 'temporal_dual_regression_smooth_hist':
                    image_readable = 'Histogram of Temporal Dual Regression'
                if image_name_a_nav == 'vmhc_smooth':
                    image_readable = 'Voxel-Mirrored Homotopic Connectivity'
                if image_name_a_nav == 'vmhc_smooth_hist':
                    image_readable = 'Histogram of Voxel-Mirrored Homotopic Connectivity'
                if image_name_a_nav == 'reho_smooth':
                    image_readable = 'Regional Homogeneity'
                if image_name_a_nav == 'reho_smooth_hist':
                    image_readable = 'Histogram of Regional Homogeneity'
                if image_name_a_nav == 'alff_smooth':
                    image_readable = 'Amplitude of Low-Frequency Fluctuation'
                if image_name_a_nav == 'alff_smooth_hist':
                    image_readable = 'Histogram of Amplitude of Low-Frequency Fluctuation'
                if image_name_a_nav == 'falff_smooth':
                    image_readable = 'fractional Amplitude of Low-Frequency Fluctuation'
                if image_name_a_nav == 'falff_smooth_hist':
                    image_readable = 'Histogram of fractional Amplitude of Low-Frequency Fluctuation'      

                feed_line_nav(id_, \
                            image_name_a_nav, \
                            #image_readable, \
                            id_a, \
                            f_html_0, \
                            f_html_1)


            feed_line_body(
                        image_name_a, \
                        id_a, \
                        png_a, \
                        f_html_1)

            feed_line_body(
                        '', \
                        id_s, \
                        png_s, \
                        f_html_1)


            if id_ in dict_hist.keys():

                if idx == 0:
                    feed_line_nav(id_, \
                              image_name_h_nav, \
                              id_h, \
                              f_html_0, \
                              f_html_1)

                feed_line_body(
                        image_name_h, \
                        id_h, \
                        png_h, \
                        f_html_1)

    if id_ in dict_plot:
        id_a = str(id_)
        image_name = qc_plot_id[id_]
        png_a = dict_plot[id_][0]
        print 'png_: ', png_a
        feed_line_nav(id_, \
                    image_name, \
                    id_a, \
                    f_html_0, \
                    f_html_1)

        feed_line_body(
                    image_name, \
                    id_a, \
                    png_a, \
                    f_html_1)


def make_page(file_, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id):

    """
    Make Page

    Parameters
    ----------

    file_ : string
        path to qc path file

    qc_montage_id_a : dictionary
        dictionary of axial montages key : id no
        value is list of png types 

    qc_montage_id_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of png types 

    qc_plot_id : dictionary
          dictionary of plot pngs key : id no
          value is list of png types

    qc_hist_id : dictionary
          dictionary of histogram pngs key : id no
          value is list of png types


    Returns
    -------

    None

    """
    import os
    from CPAC.qc.utils import grp_pngs_by_id, add_head, add_tail, \
                            feed_lines_html



    f_ = open(file_, 'r')
    pngs_ = [line.rstrip('\r\n') for line in f_.readlines()]
    f_.close()

    html_f_name = file_.replace('.txt', '')
 
    html_f_name = html_f_name.replace("'", "")   ###fixed in front, don't have to now

    html_f_name_0 = html_f_name + '_0.html'
    html_f_name_1 = html_f_name + '_1.html'
    html_f_name = html_f_name + '.html'

    
    try:

        f_html_ = open(html_f_name, 'wb')
        f_html_0 = open(html_f_name_0, 'wb')
        f_html_1 = open(html_f_name_1, 'wb')


        dict_a, dict_s, dict_hist, dict_plot, all_ids = grp_pngs_by_id(pngs_, qc_montage_id_a, \
                                            qc_montage_id_s, qc_plot_id, qc_hist_id)

        #for k, v in dict_plot.items():
        #        print '_a~~~> ', k, v

        add_head(f_html_, f_html_0, f_html_1)


        for id_ in sorted(all_ids):
            feed_lines_html(id_,
                              dict_a,
                              dict_s,
                              dict_hist,
                              dict_plot,
                              qc_montage_id_a,
                              qc_montage_id_s,
                              qc_plot_id,
                              qc_hist_id,
                              f_html_0,
                              f_html_1)

        add_tail(f_html_, f_html_0, f_html_1)

        f_html_.close()
        f_html_0.close()
        f_html_1.close()

    except:

        print '\n\nWarning: QC HTML pages could not be generated because' \
              ' the file names have become too long. Clear out the .html' \
              ' files present in the output directory and in the' \
              ' qc_files_here folders within the subject directories.\n\n'
        raise Exception

    
def make_qc_pages(qc_path, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id):


    """
    Calls make page

    Parameters
    ----------

    qc_path : string
        path to qc_files_here directory

    qc_montage_id_a : dictionary
        dictionary of axial montages key : id no
        value is list of png types 

    qc_montage_id_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of png types 

    qc_plot_id : dictionary
          dictionary of plot pngs key : id no
          value is list of png types

    qc_hist_id : dictionary
          dictionary of histogram pngs key : id no
          value is list of png types


    Returns
    -------

    None

    """
    import os
    from CPAC.qc.utils import make_page

    qc_files = os.listdir(qc_path)

    for file_ in qc_files:

        if not (file_.endswith('.txt')):
            continue
        #actually make the html page for the file_
        make_page(os.path.join(qc_path, file_), qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id)




def generateQCPages(qc_path, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id):

    """
    parafile = open('QC_input_para.txt', 'w')
    parafile.write(qc_path)
    parafile.write('qc_montage_id_a: ')
    for key, value in qc_montage_id_a.iteritems()
        parafile.write(key)
        parafile.write(value)
    parafile.write('qc_montage_id_s: ')
    for key, value in qc_montage_id_s.iteritems()
        parafile.write(key)
        parafile.write(value)
    parafile.write('qc_plot_id: ')
    for key, value in qc_plot_id.iteritems()
        parafile.write(key)
        parafile.write(value)
    parafile.write('qc_hist_id: ')
    for key, value in qc_hist_id.iteritems()
        parafile.write(key)
        parafile.write(value)
    parafile.close()
    """

    """
    Calls make_qc_page and organizes qc path files

    Parameters
    ----------

    qc_path : string
        path to qc_files_here directory

    qc_montage_id_a : dictionary
        dictionary of axial montages key : id no
        value is list of png types 

    qc_montage_id_s : dictionary
          dictionary of sagittal montages key : id no
          value is list of png types 

    qc_plot_id : dictionary
          dictionary of plot pngs key : id no
          value is list of png types

    qc_hist_id : dictionary
          dictionary of histogram pngs key : id no
          value is list of png types


    Returns
    -------

    None

    """
    import os
    from CPAC.qc.utils import first_pass_organizing_files, second_pass_organizing_files
    from CPAC.qc.utils import make_qc_pages

    #os.system('rm -rf %s/*.html')

    #according to preprocessing strategy combines the files
    first_pass_organizing_files(qc_path)
    #according to bandpass and hp_lp and smoothing iterables combines the files
    second_pass_organizing_files(qc_path)

    #generate pages from qc files
    make_qc_pages(qc_path, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id)



def make_edge(file_):

    """
    Make edge file from a scan image

    Parameters
    ----------

    file_ :    string
        path to the scan

    Returns
    -------

    new_fname : string
        path to edge file

    """

    import commands
    import os

    remainder, ext_ = os.path.splitext(file_)

    remainder, ext1_ = os.path.splitext(remainder)

    ext = ''.join([ext1_, ext_])

    new_fname = ''.join([remainder, '_edge', ext])
    new_fname = os.path.join(os.getcwd(), os.path.basename(new_fname))

    cmd = "3dedge3 -input %s -prefix %s" % (file_, new_fname)
    print cmd
    print commands.getoutput(cmd)

    return new_fname


def gen_func_anat_xfm(func_, ref_, xfm_, interp_):

    """
    Transform functional file (std dev) into anatomical space

    Parameters
    ----------

    func_ : string
        functional scan

    ref_ : string
        path to reference file

    xfm_ : string
        path to transformation mat file

    interp_ : string
        interpolation measure string

    Returns
    -------

    new_fname : string
        path to the transformed scan
    """

    import os
    import commands

    new_fname = os.path.join(os.getcwd(), 'std_dev_anat.nii.gz')

    cmd = 'applywarp --ref=%s --in=%s --out=%s --premat=%s --interp=%s' % (ref_, func_, new_fname, xfm_, interp_)
    print cmd

    print commands.getoutput(cmd)

    return new_fname


def gen_snr(std_dev, mean_func_anat):

    """
    Generate SNR file

    Parameters
    ----------

    std_dev : string
        path to std dev file in anat space

    mean_func_anat : string
        path to mean functional scan in anatomical space


    Returns
    -------

    new_fname : string
        path to the snr file
    """
    import os
    import commands

    new_fname = os.path.join(os.getcwd(), 'snr.nii.gz')

    cmd = '3dcalc -a %s -b %s -expr \"b/a\" -prefix %s' % (std_dev, mean_func_anat, new_fname)

    print cmd
    print commands.getoutput(cmd)

    return new_fname


###
def cal_snr_val(measure_file):

    """
    Calculate average snr value for snr image.

    Parameters
    ----------

    measure_file : string

        path to input nifti file


    Returns
    -------

    avg_snr_file : string

        a text file store average snr value

    """

    import numpy.ma as ma
    import nibabel as nb
    import os

    data = nb.load(measure_file).get_data()
    data_flat = data.flatten()
    data_no0 = data_flat[data_flat > 0]
    snr_val = ma.mean(data_no0)

    avg_snr_file = os.path.join(os.getcwd(), 'average_snr_file.txt')
    f = open(avg_snr_file, 'w')
    f.write(str(snr_val) + '\n')

    #f.write(str(measure_file) + '\n')
    #f.write(str(avg_snr_file) + '\n')

    f.close()

    return avg_snr_file





def gen_std_dev(mask_, func_):

    """
    Generate std dev file

    Parameters
    ----------
    
    mask_ : string
        path to whole brain mask file

    func_ : string
        path to functional scan

    Returns
    -------

    new_fname : string
        path to standard deviation file
    """

    import os
    import commands

    new_fname = os.path.join(os.getcwd(), 'std_dev.nii.gz')


    cmd = "3dTstat -stdev -mask %s -prefix %s %s" % (mask_, new_fname, func_)

    print commands.getoutput(cmd)

    return new_fname


def drange(min_, max_):

    """
    Generate list of float values in a specified range.

    Parameters
    ----------

    min_ : float
        Min value

    max_ : float
        Max value


    Returns
    -------

    range_ : list
        list of float values in the min_ max_ range
    
    """

    step = float(max_ - min_) /8.0
    range_ = []

    while min_ <= max_:

        range_.append(float('%.3f' % round(min_, 3)))
        min_ += step
    return range_


def gen_plot_png(arr, ex_vol, measure):

    """
    Generate Motion FD Plot. Shows which volumes were dropped.

    Parameters
    ----------

    arr : list
        Frame wise Displacements

    ex_vol : list
        Volumes excluded

    measure : string
        Label of the Measure


    Returns
    -------

    png_name : string
            path to the generated plot png
    """


    import matplotlib
    import commands
#    matplotlib.use('Agg')
    from matplotlib import pyplot
    matplotlib.rcParams.update({'font.size': 8})
    import matplotlib.cm as cm
    import numpy as np
    import os


    arr = np.loadtxt(arr)
    try:
        ex_vol = np.genfromtxt(ex_vol, delimiter=',', dtype=int)
        ex_vol = ex_vol[ex_vol > 0]
    except:
        ex_vol = []
    arr = arr[1:]
    del_el = [x for x in ex_vol if x < len(arr)]

    ex_vol = np.array(del_el)


    fig = pyplot.figure(figsize=(10, 6))
    pyplot.plot([i for i in xrange(len(arr))], arr, '-')
    fig.suptitle('%s plot with Mean %s = %0.4f' % (measure, measure, arr.mean()))
    if measure == 'FD' and len(ex_vol) > 0:

        pyplot.scatter(ex_vol, arr[ex_vol], c="red", zorder=2)

        for x in ex_vol:
            pyplot.annotate('( %d , %0.3f)' % (x, arr[x]), xy=(x, arr[x]), arrowprops=dict(facecolor='black', shrink=0.0))

    pyplot.xlabel('Volumes')
    pyplot.ylabel('%s' % measure)
    png_name = os.path.join(os.getcwd(), '%s_plot.png' % measure)
    fig.savefig(os.path.join(os.getcwd(), png_name))
    pyplot.close()
    matplotlib.rcdefaults()
    return png_name


def gen_motion_plt(motion_parameters):

    """
    Function to Generate Matplotlib plot for motion.
    Separate plots for Translation and Rotation are generated.

    Parameters
    ----------

    motion_parameters: string
                    Motion Parameters file

    Returns
    -------

    translation_plot : string
        path to translation plot

    rotation_plot : string
        path to rotation plot

    """

    import matplotlib
    import commands
#    matplotlib.use('Agg')
    import matplotlib.cm as cm
    import numpy as np
    from matplotlib import pyplot as plt
    import math
    import os

    png_name1 = 'motion_trans_plot.png'
    png_name2 = 'motion_rot_plot.png'
    data = np.loadtxt(motion_parameters)

    data_t = data.T

    translation_plot = None
    rotation_plot = None

    titles1 = ['x', 'y', 'z']
    titles2 = ['roll', 'pitch', 'yaw']

    plt.gca().set_color_cycle(['red', 'green', 'blue'])
    plt.plot(data_t[0])
    plt.plot(data_t[1])
    plt.plot(data_t[2])
    plt.legend(['x', 'y', 'z'], loc='upper right')
    plt.ylabel('Translation (mm)')
    plt.xlabel('Volume')

    plt.savefig(os.path.join(os.getcwd(), png_name1))

    plt.close()

    for i in range(3, 6):

        for j in range(len(data_t[i])):

            data_t[i][j] = math.degrees(data_t[i][j])

    plt.gca().set_color_cycle(['red', 'green', 'blue'])
    plt.plot(data_t[3])
    plt.plot(data_t[4])
    plt.plot(data_t[5])
    plt.legend(['roll', 'pitch', 'yaw'], loc='upper right')
    plt.ylabel('Rotation (degrees)')
    plt.xlabel('Volume')


    plt.savefig(os.path.join(os.getcwd(), png_name2))

    plt.close()

    translation_plot = os.path.join(os.getcwd(), png_name1)
    rotation_plot = os.path.join(os.getcwd(), png_name2)

    return translation_plot, rotation_plot


def gen_histogram(measure_file, measure):

    """
    Generates Histogram Image of intensities for a given input
    nifti file.

    Parameters
    ----------

    measure_file : string

                path to input nifti file

    measure : string

        Name of the measure label in the plot


    Returns
    -------

    hist_path : string

        Path to the generated histogram png

    """
    hist_path = None

    from CPAC.qc.utils import make_histogram
    import os
    m_ = measure
    if isinstance(measure_file, list):
        hist_path = []
        for file_ in measure_file:
            measure = m_
            if 'sca_roi' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                fname = fname.split('ROI_number_')[1]
                fname = 'SCA_ROI_Number_' + fname.split('_')[0]
                measure = fname
            if 'temporal_regression_sca' in measure.lower():

                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                fname = fname.split('z_maps_roi_')[1]
                fname = 'z_maps_roi_' + fname.split('_')[0]
                measure = fname

            if 'temporal_dual_regression' in measure.lower():

                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                fname = fname.split('map_z_')[1]
                fname = 'map_z_'+ fname.split('_')[0]
                measure = fname

            if 'centrality' in measure.lower():

                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])

                if 'centrality_' in fname:
                    type_, fname = fname.split('centrality_')
                    fname = type_ + 'centrality_' + fname.split('_')[0]
                    measure = fname
                elif 'lFCD_' in fname:
                    fname = 'lFCD_' + fname.split('_')[1]
                    measure = fname

            hist_path.append(make_histogram(file_, measure))

    else:

        print "measure_file: ", measure_file
        print "measure: ", measure

        hist_path = make_histogram(measure_file, measure)

    return hist_path



def make_histogram(measure_file, measure):

    """
    Generates Histogram Image of intensities for a given input
    nifti file.

    Parameters
    ----------

    measure_file : string

                path to input nifti file

    measure : string

        Name of the measure label in the plot


    Returns
    -------

    hist_path : string

        Path to the generated histogram png

    """

    import matplotlib
    import commands
#    matplotlib.use('Agg')
    from matplotlib import pyplot
    import matplotlib.cm as cm
    import numpy as np
    import nibabel as nb
    import os

    data = nb.load(measure_file).get_data()
    data_flat = data.flatten(order='F')
    y, binEdges = np.histogram(data_flat[data_flat != 0], bins=100)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

    fig = pyplot.figure()
    fig.suptitle('%s intensity plot' % measure)
    pyplot.plot(bincenters, y, '-')
    pyplot.xlabel('intensity')
    pyplot.ylabel('# of voxels')

    png_name = os.path.join(os.getcwd(), '%s_hist_plot.png' % measure)
    fig.savefig(os.path.join(os.getcwd(), png_name))

    pyplot.close()
    hist_path = os.path.join(os.getcwd(), png_name)

    """
    ###
    hist_file = os.path.join(os.getcwd(), '%s_hist_path_file.txt' % measure)
    fl = open(hist_file, 'w')
    fl.write(str(measure_file) + '\n')
    fl.write(str(hist_path) + '\n')

    fl.close()
    """

    return hist_path


def drop_percent_(measure_file, percent_):

    """
    Zeros out voxels in measure filewhose intensity doesnt fall
    in percent_ of voxel intensities

    Parameters
    ----------

    measure_file : string
                Input nifti file

    percent_ : percentage of the voxels to keep

    
    Returns
    -------

    modified_measure_file : string
                    measure_file with 1 - percent_ voxels zeroed out 
    

    """

    import nibabel as nb
    import numpy as np
    import os
    import commands

    img = nb.load(measure_file)

    data = img.get_data()

    x, y, z = data.shape


    max_val= float(commands.getoutput('fslstats %s -P %f' %(measure_file, percent_)))

    for i in range(x):

        for j in range(y):

            for k in range(z):
                if data[i][j][k] > 0.0:
                    if data[i][j][k] >= max_val:
                        data[i][j][k] = 0.0

    save_img = nb.Nifti1Image(data, header=img.get_header(), affine=img.get_affine())

    f_name = os.path.basename(os.path.splitext(os.path.splitext(measure_file)[0])[0])

    saved_name = None
    saved_name_correct_header = None
    ext = None

    if '.nii.gz' in measure_file:
        ext = '.nii.gz'
    else:
        ext = '.nii'

    saved_name = '%s_%d_%s' % (f_name, percent_, ext)
    saved_name_correct_header = '%s_%d%s' % (f_name, percent_, ext)
    save_img.to_filename(saved_name)

    commands.getoutput("3dcalc -a %s -expr 'a' -prefix %s" % (saved_name, saved_name_correct_header))

    modified_measure_file = os.path.join(os.getcwd(), saved_name_correct_header)
    return modified_measure_file


def get_spacing(across, down, dimension):

    """
    Get Spacing in slices to be selected for montage
    display varying in given dimension

    Parameters
    ----------

    across : integer
        # images placed horizontally in montage

    down : integer
        # images stacked vertically in montage


    Returns
    -------

    space : integer
        # of images to skip before displaying next one

    """

#    across = 6
#    down = 3
    space = 10


    prod = (across*down*space)

    if prod > dimension:
        while(across*down*space) > dimension:
            space -= 1
    else:
        while(across*down*space) < dimension:
            space += 1


    return space


def determine_start_and_end(data, direction, percent):

    """
    Determine start slice and end slice in data file in
    given direction with at least threshold percent of voxels
    at start and end slices.

    Parameters
    ----------

    data : string
        input nifti file

    direction : string
        axial or sagittal

    percent : float
        percent(from total) of non zero voxels at starting and ending slice


    Returns
    -------

    start : integer
            Index of starting slice

    end : integer
            Index of the last slice

    """

    x, y, z = data.shape

    xx1 = 0
    xx2 = x - 1
    zz1 = 0
    zz2 = z - 1
    total_non_zero_voxels = len(np.nonzero(data.flatten())[0])
    thresh = percent * float(total_non_zero_voxels)
    start = None
    end = None


    if 'axial' in direction:

        while(zz2 > 0):

            d = len(np.nonzero(data[:, :, zz2].flatten())[0])
            if float(d) > thresh:
                break

            zz2 -= 1

        while(zz1 < zz2):

            d = len(np.nonzero(data[:, :, zz1].flatten())[0])
            if float(d) > thresh:
                break

            zz1 += 1

        start =  zz1
        end = zz2

    else:
        while(xx2 > 0):
            d = len(np.nonzero(data[xx2, :, :].flatten())[0])
            if float(d) > thresh:
                break

            xx2 -= 1

        while(xx1 < xx2):

            d = len(np.nonzero(data[xx1, :, :].flatten())[0])
            if float(d) > thresh:
                break

            xx1 += 1

        start = xx1
        end = xx2

    return start, end



def montage_axial(overlay, underlay, png_name, cbar_name):

    """
    Draws Montage using overlay on Anatomical brain in Axial Direction
    calls make_montage_axial

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar 

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    import os
    from CPAC.qc.utils import make_montage_axial
    pngs = None
    if isinstance(overlay, list):
        pngs = []

        for ov in overlay:
            fname = os.path.basename(os.path.splitext(os.path.splitext(ov)[0])[0])
            pngs.append(make_montage_axial(ov, underlay, fname + '_' + png_name, cbar_name))
    else:

        pngs = make_montage_axial(overlay, underlay, png_name, cbar_name)

    png_name = pngs

    return png_name




def make_montage_axial(overlay, underlay, png_name, cbar_name):

    """
    Draws Montage using overlay on Anatomical brain in Axial Direction

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar 

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    import matplotlib
    import commands
#    matplotlib.use('Agg')
    import os
    matplotlib.rcParams.update({'font.size': 5})
    import matplotlib.cm as cm
    ###
    try:
        from mpl_toolkits.axes_grid1 import ImageGrid   
    except:
        from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib.colors as col
    import nibabel as nb
    import numpy as np
    from CPAC.qc.utils import determine_start_and_end, get_spacing

    Y = nb.load(underlay).get_data()
    X = nb.load(overlay).get_data()
    X = X.astype(np.float16)
    Y = Y.astype(np.float16)

    if  'skull_vis' in png_name:
        X[X < 20.0] = 0.0
    if 'skull_vis' in png_name or 't1_edge_on_mean_func_in_t1' in png_name or 'MNI_edge_on_mean_func_mni' in png_name:
        max_ = np.nanmax(np.abs(X.flatten()))
        X[X != 0.0] = max_
        print '^^', np.unique(X)
    z1, z2 = determine_start_and_end(Y, 'axial', 0.0001)
    spacing = get_spacing(6, 3, z2 - z1)
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    if ('snr' in png_name) or  ('reho' in png_name) or ('vmhc' in png_name) or ('sca_' in png_name) or ('alff' in png_name) or ('centrality' in png_name) or ('temporal_regression_sca' in png_name)  or ('temporal_dual_regression' in png_name):
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True, aspect=True, cbar_mode="single", cbar_pad=0.2, direction="row")
    else:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True, aspect=True, direction="row")

    zz = z1
    for i in range(6*3):

        if zz >= z2:
            break

        im = grid[i].imshow(np.rot90(Y[:, :, zz]), cmap=cm.Greys_r)
        zz += spacing

    x, y, z = X.shape
    X[X == 0.0] = np.nan
    max_ = np.nanmax(np.abs(X.flatten()))
    print '~~', max_


    zz = z1
    im = None
    print '~~~', z1, ' ', z2
    for i in range(6*3):


        if zz >= z2:
            break
        if cbar_name is 'red_to_blue':

            im = grid[i].imshow(np.rot90(X[:, :, zz]), cmap=cm.get_cmap(cbar_name), alpha=0.82, vmin=0, vmax=max_)   ###
        elif cbar_name is 'green':
            im = grid[i].imshow(np.rot90(X[:, :, zz]), cmap=cm.get_cmap(cbar_name), alpha=0.82, vmin=0, vmax=max_)
        else:
            im = grid[i].imshow(np.rot90(X[:, :, zz]), cmap=cm.get_cmap(cbar_name), alpha=0.82, vmin=- max_, vmax=max_)   

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

    if 'snr' in png_name:
        cbar.ax.set_yticks(drange(0, max_))

    elif  ('reho' in png_name) or ('vmhc' in png_name) or ('sca_' in png_name) or ('alff' in png_name) or ('centrality' in png_name) or ('temporal_regression_sca' in png_name) or ('temporal_dual_regression' in png_name):
        cbar.ax.set_yticks(drange(-max_, max_))


#    plt.show()
    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    matplotlib.rcdefaults()

    return png_name




def montage_sagittal(overlay, underlay, png_name, cbar_name):

    """
    Draws Montage using overlay on Anatomical brain in Sagittal Direction
    calls make_montage_sagittal

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar 

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    import os
    from CPAC.qc.utils import make_montage_sagittal
    pngs = None
    if isinstance(overlay, list):
        pngs = []

        for ov in overlay:
            fname = os.path.basename(os.path.splitext(os.path.splitext(ov)[0])[0])
            pngs.append(make_montage_sagittal(ov, underlay, fname + '_' + png_name, cbar_name))

    else:

        pngs = make_montage_sagittal(overlay, underlay, png_name, cbar_name)
    png_name = pngs

    return png_name


def make_montage_sagittal(overlay, underlay, png_name, cbar_name):

    """
    Draws Montage using overlay on Anatomical brain in Sagittal Direction

    Parameters
    ----------

    overlay : string
            Nifi file

    underlay : string
            Nifti for Anatomical Brain

    cbar_name : string
            name of the cbar 

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """
    from CPAC.qc.utils import determine_start_and_end, get_spacing
    import matplotlib
    import commands
#    matplotlib.use('Agg')
    import os
    import numpy as np
    matplotlib.rcParams.update({'font.size': 5})
    ###
    try:
        from mpl_toolkits.axes_grid1 import ImageGrid   
    except:
        from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import matplotlib.colors as col
    import nibabel as nb

    Y = nb.load(underlay).get_data()
    X = nb.load(overlay).get_data()
    X = X.astype(np.float16)
    Y = Y.astype(np.float16)


    if  'skull_vis' in png_name:
        X[X < 20.0] = 0.0
    if 'skull_vis' in png_name or 't1_edge_on_mean_func_in_t1' in png_name or 'MNI_edge_on_mean_func_mni' in png_name:
        max_ = np.nanmax(np.abs(X.flatten()))
        X[X != 0.0] = max_
        print '^^', np.unique(X)

    x1, x2 = determine_start_and_end(Y, 'sagittal', 0.0001)
    spacing = get_spacing(6, 3, x2 - x1)
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    if ('snr' in png_name) or  ('reho' in png_name) or ('vmhc' in png_name) or ('sca_' in png_name) or ('alff' in png_name) or ('centrality' in png_name) or ('temporal_regression_sca' in png_name)  or ('temporal_dual_regression' in png_name):
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True, aspect=True, cbar_mode="single", cbar_pad=0.5, direction="row")
    else:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True, aspect=True, cbar_mode="None", direction="row")

    xx = x1
    for i in range(6*3):

        if xx >= x2:
            break

        im = grid[i].imshow(np.rot90(Y[xx, :, :]), cmap=cm.Greys_r)
        grid[i].get_xaxis().set_visible(False)
        grid[i].get_yaxis().set_visible(False)
        xx += spacing

    x, y, z = X.shape
    X[X == 0.0] = np.nan
    max_ = np.nanmax(np.abs(X.flatten()))
    print '~~', max_
    xx = x1
    for i in range(6*3):


        if xx >= x2:
            break
        im = None
        if cbar_name is 'red_to_blue':

            im = grid[i].imshow(np.rot90(X[xx, :, :]), cmap=cm.get_cmap(cbar_name), alpha=0.82, vmin=0, vmax=max_)   ###
        elif cbar_name is 'green':
            im = grid[i].imshow(np.rot90(X[xx, :, :]), cmap=cm.get_cmap(cbar_name), alpha=0.82, vmin=0, vmax=max_)
        else:
            im = grid[i].imshow(np.rot90(X[xx, :, :]), cmap=cm.get_cmap(cbar_name), alpha=0.82, vmin=- max_, vmax=max_)   
        xx += spacing
    cbar = grid.cbar_axes[0].colorbar(im)

    if 'snr' in png_name:
        cbar.ax.set_yticks(drange(0, max_))

    elif  ('reho' in png_name) or ('vmhc' in png_name) or ('sca_' in png_name) or ('alff' in png_name) or ('centrality' in png_name) or ('temporal_regression_sca' in png_name)  or ('temporal_dual_regression' in png_name):
        cbar.ax.set_yticks(drange(-max_, max_))


#    plt.show()
    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()
    matplotlib.rcdefaults()

    return png_name


def montage_gm_wm_csf_axial(overlay_csf, overlay_wm, overlay_gm, underlay, png_name):

    """
    Draws Montage using GM WM and CSF overlays on Anatomical brain in Sagittal Direction

    Parameters
    ----------

    overlay_csf : string
            Nifi file CSF MAP

    overlay_wm : string
            Nifti file WM MAP

    overlay_gm : string
            Nifti file GM MAP

    underlay : string
            Nifti for Anatomical Brain

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """

    import os
    import matplotlib
    import commands
#    matplotlib.use('Agg')
    from CPAC.qc.utils import determine_start_and_end, get_spacing
    import numpy as np
    ###
    try:
        from mpl_toolkits.axes_grid1 import ImageGrid   
    except:
        from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib.colors as col
    import nibabel as nb
    import matplotlib.cm as cm

    Y = nb.load(underlay).get_data()
    z1, z2 = determine_start_and_end(Y, 'axial', 0.0001)
    spacing = get_spacing(6, 3, z2 - z1)
    X_csf = nb.load(overlay_csf).get_data()
    X_wm = nb.load(overlay_wm).get_data()
    X_gm = nb.load(overlay_gm).get_data()
    X_csf = X_csf.astype(np.float16)
    X_wm = X_wm.astype(np.float16)
    X_gm = X_gm.astype(np.float16)
    Y = Y.astype(np.float16)

    max_csf = np.nanmax(np.abs(X_csf.flatten()))
    X_csf[X_csf != 0.0] = max_csf
    max_wm = np.nanmax(np.abs(X_wm.flatten()))
    X_wm[X_wm != 0.0] = max_wm
    max_gm = np.nanmax(np.abs(X_gm.flatten()))
    X_gm[X_gm != 0.0] = max_gm
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))
    grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True, aspect=True, cbar_mode="None", direction="row")

    zz = z1
    for i in range(6*3):

        if zz >= z2:
            break

        im = grid[i].imshow(np.rot90(Y[:, :, zz]), cmap=cm.Greys_r)
        zz += spacing

    x, y, z = X_csf.shape
    X_csf[X_csf == 0.0] = np.nan
    X_wm[X_wm == 0.0] = np.nan
    X_gm[X_gm == 0.0] = np.nan
    print '~~', max_


    zz = z1
    im = None
    for i in range(6*3):


        if zz >= z2:
            break

        im = grid[i].imshow(np.rot90(X_csf[:, :, zz]), cmap=cm.get_cmap('green'), alpha=0.82, vmin=0, vmax=max_csf)   ###
        im = grid[i].imshow(np.rot90(X_wm[:, :, zz]), cmap=cm.get_cmap('blue'), alpha=0.82, vmin=0, vmax=max_wm)
        im = grid[i].imshow(np.rot90(X_gm[:, :, zz]), cmap=cm.get_cmap('red'), alpha=0.82, vmin=0, vmax=max_gm)   

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

#    plt.show()
    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    return png_name


def montage_gm_wm_csf_sagittal(overlay_csf, overlay_wm, overlay_gm, underlay, png_name):

    """
    Draws Montage using GM WM and CSF overlays on Anatomical brain in Sagittal Direction

    Parameters
    ----------

    overlay_csf : string
            Nifi file CSF MAP

    overlay_wm : string
            Nifti file WM MAP

    overlay_gm : string
            Nifti file GM MAP

    underlay : string
            Nifti for Anatomical Brain

    png_name : string
            Proposed name of the montage plot

    Returns
    -------

    png_name : Path to generated PNG

    """

    import os
    import matplotlib
    import commands
#    matplotlib.use('Agg')
    from CPAC.qc.utils import determine_start_and_end, get_spacing
    import numpy as np
    ###
    try:
        from mpl_toolkits.axes_grid1 import ImageGrid   
    except:
        from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib.colors as col
    import matplotlib.cm as cm
    import nibabel as nb

    Y = nb.load(underlay).get_data()
    x1, x2 = determine_start_and_end(Y, 'sagittal', 0.0001)
    spacing = get_spacing(6, 3, x2 - x1)
    X_csf = nb.load(overlay_csf).get_data()
    X_wm = nb.load(overlay_wm).get_data()
    X_gm = nb.load(overlay_gm).get_data()
    X_csf = X_csf.astype(np.float16)
    X_wm = X_wm.astype(np.float16)
    X_gm = X_gm.astype(np.float16)
    Y = Y.astype(np.float16)

    max_csf = np.nanmax(np.abs(X_csf.flatten()))
    X_csf[X_csf != 0.0] = max_csf
    max_wm = np.nanmax(np.abs(X_wm.flatten()))
    X_wm[X_wm != 0.0] = max_wm
    max_gm = np.nanmax(np.abs(X_gm.flatten()))
    X_gm[X_gm != 0.0] = max_gm
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))
    grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True, aspect=True, cbar_mode="None", direction="row")

    zz = x1
    for i in range(6*3):

        if zz >= x2:
            break

        im = grid[i].imshow(np.rot90(Y[zz, :, :]), cmap=cm.Greys_r)
        zz += spacing

    x, y, z = X_csf.shape
    X_csf[X_csf == 0.0] = np.nan
    X_wm[X_wm == 0.0] = np.nan
    X_gm[X_gm == 0.0] = np.nan
    print '~~', max_


    zz = x1
    im = None
    for i in range(6*3):


        if zz >= x2:
            break

        im = grid[i].imshow(np.rot90(X_csf[zz, :, :]), cmap=cm.get_cmap('green'), alpha=0.82, vmin=0, vmax=max_csf)   ###
        im = grid[i].imshow(np.rot90(X_wm[zz, :, :]), cmap=cm.get_cmap('blue'), alpha=0.82, vmin=0, vmax=max_wm)
        im = grid[i].imshow(np.rot90(X_gm[zz, :, :]), cmap=cm.get_cmap('red'), alpha=0.82, vmin=0, vmax=max_gm)   

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

#    plt.show()
    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    return png_name


def register_pallete(file_, cbar_name):

    """
    Registers color pallete to matplotlib

    Parameters
    ----------

    file_ : string
        file containing colors in hexadecimal formats in each line

    cbar_name : string
        Proposed name for the color bar


    Returns
    -------

    None

    """

    import matplotlib
    import commands
#    matplotlib.use('Agg')
    import matplotlib.colors as col
    import matplotlib.cm as cm
    f = open(file_, 'r')

    colors_ = f.readlines()

    colors = []

    for color in reversed(colors_):

        colors.append(color.rstrip('\r\n'))

    cmap3 = col.ListedColormap(colors, cbar_name)
    cm.register_cmap(cmap=cmap3)


def resample_1mm(file_):

    """
    Calls make_resample_1mm which resamples file to 1mm space

    Parameters
    ----------

    file_ : string
        path to the scan

    Returns
    -------

    new_fname : string
        path to 1mm resampled nifti file

    """

    from CPAC.qc.utils import make_resample_1mm

    new_fname = None
    if isinstance(file_, list):
        new_fname = []

        for f in file_:

            new_fname.append(make_resample_1mm(f))

    else:
        new_fname = make_resample_1mm(file_)

    return new_fname


def make_resample_1mm(file_):

    """
    Resamples input nifti file to 1mm space

    Parameters
    ----------

    file_ : string
        Input Nifti File

    Returns
    -------

    new_fname : string
            Input Nifti resampled to 1mm space
    """

    import os
    import commands

    remainder, ext_ = os.path.splitext(file_)

    remainder, ext1_ = os.path.splitext(remainder)

    ext = ''.join([ext1_, ext_])

    new_fname = ''.join([remainder, '_1mm', ext])
    new_fname = os.path.join(os.getcwd(), os.path.basename(new_fname))
    cmd = " 3dresample -dxyz 1.0 1.0 1.0 -prefix %s -inset %s " % (new_fname, file_)
    print cmd
    commands.getoutput(cmd)

    return new_fname
