import os
import re
import math
import base64
import commands
import pkg_resources as p

import numpy as np
import nibabel as nb
import numpy.ma as ma

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import matplotlib.cm as cm
from matplotlib import gridspec as mgs
from matplotlib.colors import ListedColormap

from nipype.interfaces import afni
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


derivative_descriptions = {
    'carpet': 'Carpet',

    'alff_smooth_hist': 'Histogram of Amplitude of Low-Frequency Fluctuation (smoothed)',
    'alff_smooth': 'Amplitude of Low-Frequency Fluctuation (smoothed)',
    'alff_to_standard': 'Amplitude of Low-Frequency Fluctuation',
    'alff_to_standard_hist': 'Histogram of Amplitude of Low-Frequency Fluctuation',
    'alff_to_standard_zstd': 'Amplitude of Low-Frequency Fluctuation (z-score standardized)',
    'alff_to_standard_zstd_hist': 'Histogram of Amplitude of Low-Frequency Fluctuation (z-score standardized)',
    'alff_to_standard_smooth': 'Amplitude of Low-Frequency Fluctuation (smoothed)',
    'alff_to_standard_smooth_hist': 'Histogram of Amplitude of Low-Frequency Fluctuation (smoothed)',
    'alff_to_standard_smooth_zstd': 'Amplitude of Low-Frequency Fluctuation (smoothed, z-score standardized)',
    'alff_to_standard_smooth_zstd_hist': 'Histogram of Amplitude of Low-Frequency Fluctuation (smoothed, z-score standardized)',

    'centrality_hist': 'Histogram of Network Centrality',
    'centrality_smooth_hist': 'Histogram of Network Centrality (smoothed)',
    'centrality_smooth_zstd_hist': 'Histogram of Network Centrality (smoothed, z-score standardized)',
    'centrality_smooth_zstd': 'Network Centrality (smoothed, z-score standardized)',
    'centrality_smooth': 'Network Centrality (smoothed)',
    'centrality_zstd_hist': 'Histogram of Network Centrality (z-score standardized)',
    'centrality_zstd_smooth_hist': 'Histogram of Network Centrality (z-score standardized, smoothed)',
    'centrality_zstd_smooth': 'Network Centrality (z-score standardized, smoothed)',
    'centrality_zstd': 'Network Centrality (z-score standardized)',
    'centrality': 'Network Centrality',
    
    'csf_gm_wm': 'Grey Matter, White Matter & CSF',

    'falff_smooth_hist': 'Histogram of Fractional Amplitude of Low-Frequency Fluctuation (smoothed)',
    'falff_smooth': 'Fractional Amplitude of Low-Frequency Fluctuation (smoothed)',
    'falff_to_standard': 'Fractional Amplitude of Low-Frequency Fluctuation',
    'falff_to_standard_hist': 'Histogram of Fractional Amplitude of Low-Frequency Fluctuation',
    'falff_to_standard_smooth': 'Fractional Amplitude of Low-Frequency Fluctuation (smoothed)',
    'falff_to_standard_smooth_hist': 'Histogram of Fractional Amplitude of Low-Frequency Fluctuation (smoothed)',
    'falff_to_standard_smooth_zstd': 'Fractional Amplitude of Low-Frequency Fluctuation (smoothed, z-score standardized)',
    'falff_to_standard_smooth_zstd_hist': 'Histogram of Fractional Amplitude of Low-Frequency Fluctuation (smoothed, z-score standardized)',
    'falff_to_standard_zstd': 'Fractional Amplitude of Low-Frequency Fluctuation (z-score standardized)',
    'falff_to_standard_zstd_hist': 'Histogram of Fractional Amplitude of Low-Frequency Fluctuation (z-score standardized)',

    'fd_plot': 'Framewise Displacement Plot',
    'mean_func_with_mni_edge': 'MNI Edge Overlapped on Mean Functional Image',
    'mean_func_with_t1_edge': 'T1 Edge Overlapped on Mean Functional Image',
    'mni_normalized_anatomical': 'MNI Edge Overlapped on Normalized Anatomical',
    'movement_rot_plot': 'Head Rotation Plot',
    'movement_trans_plot': 'Head Displacement Plot',

    'reho_smooth': 'Regional Homogeneity (smoothed)',
    'reho_smooth_hist': 'Histogram of Regional Homogeneity (smoothed)',
    'reho_to_standard': 'Regional Homogeneity',
    'reho_to_standard_hist': 'Histogram of Regional Homogeneity',
    'reho_to_standard_smooth': 'Regional Homogeneity (smoothed)',
    'reho_to_standard_smooth_hist': 'Histogram of Regional Homogeneity (smoothed)',
    'reho_to_standard_smooth_zstd': 'Regional Homogeneity (smoothed, z-score standardized)',
    'reho_to_standard_smooth_zstd_hist': 'Histogram of Regional Homogeneity (smoothed, z-score standardized)',
    'reho_to_standard_zstd': 'Regional Homogeneity (z-score standardized)',
    'reho_to_standard_zstd_hist': 'Histogram of Regional Homogeneity (z-score standardized)',

    'sca_roi_smooth': 'Seed-based Correlation Analysis (smoothed)',
    'sca_roi_smooth_hist': 'Histogram of Seed-based Correlation Analysis (smoothed)',
    'sca_roi_files_to_standard': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_fisher_zstd': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_fisher_zstd_hist': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_hist': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_smooth': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_smooth_fisher_zstd': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_smooth_fisher_zstd_hist': 'Seed-based Correlation Analysis',
    'sca_roi_files_to_standard_smooth_hist': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_files': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_files_hist': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_files_smooth': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_files_smooth_hist': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_zstat_files': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_zstat_files_hist': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_zstat_files_smooth': 'Seed-based Correlation Analysis',
    'sca_tempreg_maps_zstat_files_smooth_hist': 'Seed-based Correlation Analysis',

    
    'skullstrip_vis': 'Visual Result of Skull Strip',
    'snr_hist': 'Histogram of Signal to Noise Ratio',
    'snr': 'Signal to Noise Ratio',

    'temporal_dual_regression_smooth_hist': 'Histogram of Temporal Dual Regression',
    'temporal_dual_regression_smooth': 'Temporal Dual Regression',
    
    'vmhc_smooth': 'Voxel-Mirrored Homotopic Connectivity (smoothed)',
    'vmhc_smooth_hist': 'Histogram of Voxel-Mirrored Homotopic Connectivity (smoothed)',
    'vmhc_fisher_zstd': 'Fisher-Z transform map of Voxel-Mirrored Homotopic Connectivity (z-score standardized)',
    'vmhc_fisher_zstd_hist': 'Histogram of Fisher-Z transform map of Voxel-Mirrored Homotopic Connectivity (z-score standardized)',
    'vmhc_fisher_zstd_zstat_map': 'Z-Statistic map of Voxel-Mirrored Homotopic Connectivity (z-score standardized)',
    'vmhc_fisher_zstd_zstat_map_hist': 'Histogram of Z-Statistic map of Voxel-Mirrored Homotopic Connectivity (z-score standardized)',
    'vmhc_raw_score': 'Voxel-Mirrored Homotopic Connectivity',
    'vmhc_raw_score_hist': 'Histogram of Voxel-Mirrored Homotopic Connectivity',

    'dr_tempreg_maps_files_to_standard': 'Spatial Regression',
    'dr_tempreg_maps_files_to_standard_hist': 'Histogram of Spatial Regression',
    'dr_tempreg_maps_files_to_standard_smooth': 'Spatial Regression (smoothed)',
    'dr_tempreg_maps_files_to_standard_smooth_hist': 'Histogram of Spatial Regression (smoothed)',
    'dr_tempreg_maps_files_to_standard_smooth_zstd': 'Spatial Regression (smoothed, z-score standardized)',
    'dr_tempreg_maps_files_to_standard_smooth_zstd_hist': 'Histogram of Spatial Regression (smoothed, z-score standardized)',
    'dr_tempreg_maps_files_to_standard_zstd': 'Spatial Regression (z-score standardized)',
    'dr_tempreg_maps_files_to_standard_zstd_hist': 'Histogram of Spatial Regression (z-score standardized)',
    'dr_tempreg_maps_zstat_files_to_standard': 'Spatial Regression (z-score standardized)',
    'dr_tempreg_maps_zstat_files_to_standard_hist': 'Histogram of Spatial Regression (z-score standardized)',
    'dr_tempreg_maps_zstat_files_to_standard_smooth': 'Histogram of Spatial Regression (smoothed, z-score standardized)',
    'dr_tempreg_maps_zstat_files_to_standard_smooth_hist': 'Histogram of Spatial Regression (smoothed, z-score standardized)',
}


def append_to_files_in_dict_way(list_files, file_):
    """Combine files so at each resource in file appears exactly once.

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

    with open(file_, 'r') as f:
        lines = [line.rstrip('\r\n') for line in f.readlines()]
        one_dict = {line: 1 for line in lines}

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
    """First Pass at organizing qc txt files.

    Parameters
    ----------
    qc_path : string
        existing path of qc_html directory

    Returns
    -------
    None

    Notes
    -----
    Combines files with same strategy. First pass combines file names,
    where one file name is substring of the other.

    """

    if not os.path.exists(qc_path):
        os.makedirs(qc_path)

    qc_files = os.listdir(qc_path)
    strat_dict = {}

    for qc_file in sorted(qc_files, reverse=True):
        if not ('.txt' in qc_file):
            continue

        qc_file = os.path.join(qc_path, qc_file)
        qc_filename = os.path.basename(qc_file)

        qc_filename = qc_filename.replace('qc_', '')
        qc_filename = qc_filename.replace('scan_', '')
        qc_filename = qc_filename.replace('.txt', '')
        qc_filename = qc_filename.replace('____', '_')
        qc_filename = qc_filename.replace('___', '_')
        qc_filename = qc_filename.replace('__', '_')

        if '_hp_' in qc_filename and '_fwhm_' in qc_filename and \
                not ('_bandpass_freqs_' in qc_filename):
            qc_filename, fwhm_val = qc_filename.split('_fwhm_')

            fwhm_val = '_fwhm_' + fwhm_val

            qc_filename, hp_lp_ = qc_filename.split('_hp_')
            hp_lp_ = '_hp_' + hp_lp_

            qc_filename = qc_filename + fwhm_val + hp_lp_

        if strat_dict.keys() == []:
            strat_dict[qc_filename] = [qc_file]
        else:
            flag_ = 0
            for key_ in strat_dict.keys():
                if qc_filename in key_:
                    append_to_files_in_dict_way(strat_dict[key_], qc_file)
                    flag_ = 1

            if flag_ == 1:
                os.system('rm -f %s' % qc_file)
            else:
                strat_dict[qc_filename] = [qc_file]


def second_pass_organizing_files(qc_path):
    """Second Pass at organizing qc txt files.

    Parameters
    ----------
    qc_path : string
        existing path of qc_html directory

    Returns
    -------
    None

    Notes
    -----
    Combines files with same strategy. combines files for derivative 
    falff , alff with others

    """

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
        fwhm_val_ = ''

        # organize all derivatives excluding alff falff
        if '_bandpass_freqs_' in str_:
            if not str_ in strat_dict:
                strat_dict[str_] = [file_]
            else:
                print 'Error: duplicate keys for files in QC 2nd file_org ' \
                      'pass: %s %s' % (strat_dict[str_], file_)
                raise

        # organize alff falff
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
                print 'Error: duplicate keys for files in QC 2nd file_org ' \
                      'pass: %s %s' % (strat_dict[str_], file_)
                raise


def organize(dict_, all_ids, png_, new_dict):
    """Organizes pngs according to their IDS in new_dict dictionary

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

    return all_ids


def grp_pngs_by_id(pngs_, qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id):
    """Groups pngs by their ids.

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


def encode_to_url(f, type):
    with open(f, "rb") as image_file:
        b64 = str(base64.b64encode(image_file.read()).decode("utf-8"))
        return "data:" + type + ";" + "base64," + b64


def commonprefix(args, sep='/'):
	return os.path.commonprefix(args).rpartition(sep)[0]


def add_head(frameset_html_fd, menu_html_fd, content_html_fd, name):
    """Write HTML Headers to various html files.

    Parameters
    ----------
    frameset_html_fd : string
        path to main html file

    menu_html_fd : string
        path to navigation bar html file

    content_html_fd : string
        path to html file contaning pngs and plots

    Returns
    -------
    None

    """

    # Relativize files path to include on output
    html_menu_relative_name = os.path.join('qc_html', os.path.basename(menu_html_fd.name))
    html_content_relative_name = os.path.join('qc_html', os.path.basename(content_html_fd.name))

    frameset_html = """

<html>
    <head>
        <title>C-PAC QC</title>
    </head>
    <frameset cols="20%,80%">
        <frame src="{menu_file}" name="menu">
        <frame src="{content_file}" name="content">
    </frameset>
</html>

"""

    frameset_html_fd.write(frameset_html.format(
        menu_file=html_menu_relative_name,
        content_file=html_content_relative_name
    ))


    menu_html = """

<html>
    <head>
        <style>{css_nature}</style>
        <style>{css_pygments}</style>
        <base target="content">
    </head>

    <body bgcolor="#FFFF00">
        <div>
            <div class="sphinxsidebarwrapper">
                <p class="logo">
                    <a href="https://fcp-indi.github.io" target="website">
                        <img class="logo" src="{logo}" style="width:100%" alt="Logo"/>
                    </a>
                </p>
                <h3>Table Of Contents</h3>
                <ul>

"""

    with open(p.resource_filename('CPAC',"GUI/resources/html/_static/nature.css"), 'r') as content_file:
        css_nature_content = content_file.read()

    with open(p.resource_filename('CPAC',"GUI/resources/html/_static/pygments.css"), 'r') as content_file:
        css_pygments_content = content_file.read()

    menu_html_fd.write(menu_html.format(
        css_nature=css_nature_content,
        css_pygments=css_pygments_content,
        logo=encode_to_url(p.resource_filename('CPAC', "GUI/resources/html/_static/cpac_logo.jpg"), 'image/jpeg')
    ))


    content_html = """

<html>
    <body>
        <a name="reverse"></a>
        <h1>C-PAC Visual Data Quality Control Interface</h1>
        <h3>C-PAC Website: <a href=\"https://fcp-indi.github.io/\" target=\"website\">https://fcp-indi.github.io</a></h3>
        <h3>C-PAC Support Forum: <a href=\"https://groups.google.com/forum/#!forum/cpax_forum\" target=\"forum\">https://groups.google.com/forum/#!forum/cpax_forum</a></h3>
        <hr>
        <h3>Scan and strategy identifiers: {name}</h3>
    
"""

    content_html_fd.write(content_html.format(
        name=name
    ))


def add_tail(frameset_html_fd, menu_html_fd, content_html_fd):
    """Write HTML Tail Tags to various html files.

    Parameters
    ----------
    frameset_html_fd : string
        path to main html file

    menu_html_fd : string
        path to navigation bar html file

    content_html_fd : string
        path to html file contaning pngs and plots


    Returns
    -------
    None

    """

    menu_html_fd.write("""
    
                </ul>
            </div>
        </div>
    </body>
</html>

""")

    content_html_fd.write("""

    </body>
</html>
    
""")


def feed_line_nav(image_name, anchor, menu_html_fd, content_html_fd):
    """Write to navigation bar html file.

    Parameters
    ----------
    anchor : string
        anchor id of the image

    image_name : string
        name of image
    
    menu_html_fd : string
        path to navigation bar html file

    content_html_fd : string
        path to html file contaning pngs and plots

    Returns
    -------
    None

    """

    image_readable = derivative_descriptions[image_name]
        
    html_content_relative_name = os.path.join('qc_html', os.path.basename(content_html_fd.name))
    menu_html = """
                    <li><a href="{page}#{anchor}">{description}</a></li>
"""

    menu_html_fd.write(menu_html.format(
        page=html_content_relative_name,
        anchor=anchor,
        description=image_readable
    ))


def feed_line_body(image_name, anchor, image, content_html_fd):
    """Write to html file that has to contain images.

    Parameters
    ----------
    image_name : string
        name of image

    anchor : string
        anchor id of the image

    image : string
        path to the image

    content_html_fd : string
        path to html file contaning pngs and plots

    Returns
    -------
    None

    """

    folder = commonprefix([image, content_html_fd.name])

    html_rel = '/'.join(['..'] * content_html_fd.name.replace(folder + '/', '').count('/'))
    image_rel = image.replace(folder + '/', '')
    image_rel = '/'.join([html_rel, image_rel])

    description_html = """
        <h3><a name="{anchor}">{description}</a> <a href="#reverse">TOP</a></h3>
"""
    image_html = """
        <p><img src="{image}" alt="{description}"></p>
"""

    image_readable = image_name

    if image_name:
        image_readable = derivative_descriptions[image_name]

        content_html_fd.write(
            description_html.format(
                anchor=anchor,
                description=image_readable
            )
        )
        
    content_html_fd.write(
        image_html.format(
            image=image_rel,
            description=image_readable
        )
    )


def get_map_id(str_, id_):
    """Returns the proper map name given identifier for it.

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

    # so whatever goes into "type_" and then "map_id" becomes the "Map: "
    # Mask: should be the ROI nifti, but right now it's the nuisance strat...
    # Measure: should be eigenvector binarize etc., but it's just "centrality_outputs"

    if 'centrality' in id_ or 'lfcd' in id_:
        # TODO: way too reliant on a very specific string format
        # TODO: needs re-factoring
        str_ = str_.split('_a.png')[0]

        if 'lfcd' in str_:
            type_ = str_.rsplit('lfcd', 1)
        else:
            type_ = str_.rsplit(id_, 1)

        if len(type_) > 1:
            type_ = type_[1]

        if "_99_1mm_" in type_:
            type_ = type_.replace("_99_1mm_", "")
        map_id = type_

        '''
        str_ = str_.split('_')[0]
        type_ = type_.replace('_', '')
        map_id = '_'.join([type_, id_, str_])
        '''

        return map_id

    else:
        str_ = str_.split(id_)[1]
        str_ = str_.split('_')[0]
        map_id = '_'.join([id_, str_])

        return map_id


def get_map_and_measure(png_a):
    """Extract Map name and Measure name from png.

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

    measure_name = None
    map_name = None

    if '_fwhm_' in png_a:
        measure_name = os.path.basename(os.path.dirname(os.path.dirname(png_a)))
    else:
        measure_name = os.path.basename(os.path.dirname((png_a)))

    str_ = os.path.basename(png_a)

    if 'sca_tempreg' in png_a:
        map_name = get_map_id(str_, 'maps_')

    if 'sca_roi' in png_a:
        map_name = get_map_id(str_, 'roi_')

    if 'dr_tempreg' in png_a:
        map_name = get_map_id(str_, 'tempreg_maps_')

    if 'centrality' in png_a:
        map_name = get_map_id(str_, 'centrality_')

    return map_name, measure_name


def feed_lines_html(montage_id, montages_a, montages_s, histograms, dict_plot,
                    qc_montage_id_a, qc_montage_id_s, qc_plot_id, qc_hist_id,
                    menu_html_fd, content_html_fd):
    """Write HTML Tags to various html files and embeds images.

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

    if montage_id in montages_a:

        montages_a[montage_id] = sorted(montages_a[montage_id])
        montages_s[montage_id] = sorted(montages_s[montage_id])

        if montage_id in histograms:
            histograms[montage_id] = sorted(histograms[montage_id])

        idxs = len(montages_a[montage_id])

        for idx in range(0, idxs):
            png_a = montages_a[montage_id][idx]
            png_s = montages_s[montage_id][idx]
            png_h = None

            if montage_id in histograms:
                try:
                    png_h = histograms[montage_id][idx]
                except:
                    pass

            measure_name = None
            map_name = None

            if idxs > 1:
                map_name, measure_name = get_map_and_measure(png_a)

            id_a = str(montage_id)
            id_s = str(montage_id) + '_s'
            id_h = str(montage_id) + '_' + str(montage_id)

            image_name_a = None
            image_name_h = None

            image_name_a_nav = re.sub('_a$', '', qc_montage_id_a[montage_id])
            if montage_id in qc_hist_id:
                image_name_h_nav = qc_hist_id[montage_id]
                
            if map_name is not None:
                image_name_a = "Measure: {measure}; Mask: {mask}; Map: {map}".format(
                    measure=image_name_a_nav,
                    mask=measure_name,
                    map=map_name
                )

                if montage_id in qc_hist_id:
                    image_name_h = "Measure: {measure}; Mask: {mask}; Map: {map}".format(
                        measure=qc_hist_id[montage_id],
                        mask=measure_name,
                        map=map_name
                    )
            else:
                image_name_a = image_name_a_nav
                if montage_id in qc_hist_id:
                    image_name_h = qc_hist_id[montage_id]

            if idx != 0:
                id_a = '_'.join([id_a, str(idx), 'a'])
                id_s = '_'.join([id_s, str(idx), 's'])
                id_h = '_'.join([id_h, str(idx), 'h' ])

            if idx == 0:
                feed_line_nav(image_name_a_nav, id_a, menu_html_fd, content_html_fd)

            feed_line_body(image_name_a_nav, id_a, png_a, content_html_fd)
            feed_line_body(None, id_s, png_s, content_html_fd)

            if montage_id in histograms.keys():
                if idx == 0:
                    feed_line_nav(image_name_h_nav, id_h, menu_html_fd, content_html_fd)
                if png_h is not None:
                    feed_line_body(image_name_h_nav, id_h, png_h, content_html_fd)

    if montage_id in dict_plot:
        id_a = str(montage_id)
        image_name = qc_plot_id[montage_id]
        png_a = dict_plot[montage_id][0]
        feed_line_nav(image_name, id_a, menu_html_fd, content_html_fd)
        feed_line_body(image_name, id_a, png_a, content_html_fd)


def make_page(qc_file, sub_output_dir,
              qc_montage_id_a, qc_montage_id_s,
              qc_plot_id, qc_hist_id):
    """Convert a 'qc_html' text file in the CPAC output directory into
    a QC HTML page.

    Parameters
    ----------
    file_ : string
        path to qc path file

    sub_output_dir : string
        path to subject's output directory

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

    with open(qc_file, 'r') as f:
        qc_images = [line.rstrip('\r\n') for line in f.readlines()]

    frameset_html = qc_file.replace('.txt', '')
    frameset_html = frameset_html.replace("'", "")

    menu_html = frameset_html + '_navbar.html'
    content_html = frameset_html + '_page.html'

    frameset_html = "{0}.html".format(frameset_html.replace("qc_scan",
                                                            "QC-interface_scan"))
    log_dir = frameset_html.split('/qc_html')[0]
    frameset_html = frameset_html.replace("/qc_html", "")
    frameset_html = frameset_html.replace(log_dir, sub_output_dir)

    frameset_html_fd = open(frameset_html, 'wb')
    menu_html_fd = open(menu_html, 'wb')
    content_html_fd = open(content_html, 'wb')

    dict_a, dict_s, dict_hist, dict_plot, all_ids = \
        grp_pngs_by_id(qc_images,
                       qc_montage_id_a, qc_montage_id_s,
                       qc_plot_id, qc_hist_id)

    qc_path_file_id = os.path.basename(frameset_html).replace(".html", "")


    add_head(frameset_html_fd, menu_html_fd, content_html_fd, qc_path_file_id)

    for montage_id in sorted(all_ids):
        feed_lines_html(montage_id, dict_a, dict_s, dict_hist, dict_plot,
                        qc_montage_id_a, qc_montage_id_s, qc_plot_id,
                        qc_hist_id, menu_html_fd, content_html_fd)

    add_tail(frameset_html_fd, menu_html_fd, content_html_fd)


    frameset_html_fd.close()
    menu_html_fd.close()
    content_html_fd.close()

    
def make_qc_pages(qc_path, sub_output_dir, qc_montage_id_a, qc_montage_id_s,
                  qc_plot_id, qc_hist_id):
    """Generates a QC HTML file for each text file in the 'qc_html'
    folder in the CPAC output directory.

    Parameters
    ----------
    qc_path : string
        path to qc_html directory

    sub_output_dir : string
        path to subject's output directory

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
    qc_files = os.listdir(qc_path)

    for qc_file in qc_files:
        if not qc_file.endswith('.txt'):
            continue
        try:
            make_page(os.path.join(qc_path, qc_file), sub_output_dir,
                        qc_montage_id_a, qc_montage_id_s, qc_plot_id,
                        qc_hist_id)
        except IndexError as e:
            print('\n[!] Did not generate QC sub-page: {0}\n\nDetails:\n'
                  '{1}\n'.format(os.path.join(qc_path, qc_file), e))


def generate_qc_pages(qc_path, sub_output_dir,
                      qc_montage_id_a, qc_montage_id_s,
                      qc_plot_id, qc_hist_id):
    """Generates the QC HTML files populated with the QC images that were
    created during the CPAC pipeline run.

    This function runs after the pipeline is over.

    Parameters
    ----------
    qc_path : string
        path to qc_html directory

    sub_output_dir : string
        path to subject's output directory

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

    # according to preprocessing strategy combines the files
    first_pass_organizing_files(qc_path)

    # according to bandpass and hp_lp and smoothing iterables combines the
    # files
    second_pass_organizing_files(qc_path)

    make_qc_pages(qc_path, sub_output_dir, qc_montage_id_a, qc_montage_id_s,
                  qc_plot_id, qc_hist_id)


def cal_snr_val(measure_file):
    """Calculate average snr value for snr image.

    Parameters
    ----------
    measure_file : string
        path to input nifti file

    Returns
    -------
    avg_snr_file : string
        a text file store average snr value

    """

    data = nb.load(measure_file).get_data()
    data_flat = data.flatten()
    data_no0 = data_flat[data_flat > 0]
    snr_val = ma.mean(data_no0)

    avg_snr_file = os.path.join(os.getcwd(), 'average_snr_file.txt')
    with open(avg_snr_file, 'w') as f:
        f.write(str(snr_val) + '\n')

    return avg_snr_file



def drange(min_, max_):
    """Generate list of float values in a specified range.

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


def gen_plot_png(arr, measure, ex_vol=None):
    """Generate Motion FD Plot. Shows which volumes were dropped.

    Parameters
    ----------
    arr : list
        Frame wise Displacements

    measure : string
        Label of the Measure

    ex_vol : list
        Volumes excluded

    Returns
    -------
    png_name : string
            path to the generated plot png
    """

    matplotlib.rcParams.update({'font.size': 8})

    arr = np.loadtxt(arr)

    if ex_vol:
        try:
            ex_vol = np.genfromtxt(ex_vol, delimiter=',', dtype=int)
            ex_vol = ex_vol[ex_vol > 0]
        except:
            ex_vol = []
    else:
        ex_vol = []

    arr = arr[1:]
    del_el = [x for x in ex_vol if x < len(arr)]

    ex_vol = np.array(del_el)

    fig = plt.figure(figsize=(10, 6))
    plt.plot([i for i in xrange(len(arr))], arr, '-')
    fig.suptitle('%s plot with Mean %s = %0.4f' % (measure, measure,
                                                   arr.mean()))
    if measure == 'FD' and len(ex_vol) > 0:

        plt.scatter(ex_vol, arr[ex_vol], c="red", zorder=2)

        for x in ex_vol:
            plt.annotate('( %d , %0.3f)' % (x, arr[x]), xy=(x, arr[x]),
                            arrowprops=dict(facecolor='black', shrink=0.0))

    plt.xlabel('Volumes')
    plt.ylabel('%s' % measure)
    png_name = os.path.join(os.getcwd(), '%s_plot.png' % measure)
    fig.savefig(os.path.join(os.getcwd(), png_name))
    plt.close()
    matplotlib.rcdefaults()
    return png_name


def gen_carpet_plt(gm_mask, wm_mask, csf_mask, functional_to_standard, output):
    
    size = (950, 800)
    
    carpet_plot_path = os.path.join(os.getcwd(), output + '.png')

    func = nb.load(functional_to_standard).get_data()
    gm_voxels = func[nb.load(gm_mask).get_data().astype(bool)]
    wm_voxels = func[nb.load(wm_mask).get_data().astype(bool)]
    csf_voxels = func[nb.load(csf_mask).get_data().astype(bool)]
    del func

    data = np.concatenate((gm_voxels, wm_voxels, csf_voxels))
    seg = np.concatenate((
        np.ones(gm_voxels.shape[0]) * 1,
        np.ones(wm_voxels.shape[0]) * 2,
        np.ones(csf_voxels.shape[0]) * 3
    ))

    p_dec = 1 + data.shape[0] // size[0]
    if p_dec:
        data = data[::p_dec, :]
        seg = seg[::p_dec]

    t_dec = 1 + data.shape[1] // size[1]
    if t_dec:
        data = data[:, ::t_dec]
        
    interval = max((int(data.shape[-1] + 1) // 10, int(data.shape[-1] + 1) // 5, 1))
    xticks = list(range(0, data.shape[-1])[::interval])

    mycolors = ListedColormap(cm.get_cmap('tab10').colors[:4][::-1])

    gs = mgs.GridSpecFromSubplotSpec(1, 2, subplot_spec=mgs.GridSpec(1, 1)[0],
                                    width_ratios=[1, 100],
                                    wspace=0.0)
    ax0 = plt.subplot(gs[0])
    ax0.set_yticks([])
    ax0.set_xticks([])
    ax0.imshow(seg[:, np.newaxis], interpolation='none', aspect='auto',
            cmap=mycolors, vmin=1, vmax=4)
    ax0.grid(False)
    ax0.spines["left"].set_visible(False)
    ax0.spines["top"].set_visible(False)

    ax1 = plt.subplot(gs[1])
    ax1.imshow(data, interpolation='nearest', aspect='auto', cmap='gray')
    ax1.grid(False)
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_xticks(xticks)
    ax1.set_xlabel('time (frames)')
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')

    plt.savefig(carpet_plot_path, dpi=200, bbox_inches='tight')
    plt.close()

    return carpet_plot_path


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

    rotation_plot = os.path.join(os.getcwd(), 'motion_trans_plot.png')
    translation_plot = os.path.join(os.getcwd(), 'motion_rot_plot.png')

    data = np.loadtxt(motion_parameters).T

    plt.gca().set_color_cycle(['red', 'green', 'blue'])
    plt.plot(data[0])
    plt.plot(data[1])
    plt.plot(data[2])
    plt.legend(['roll', 'pitch', 'yaw'], loc='upper right')
    plt.ylabel('Rotation (degrees)')
    plt.xlabel('Volume')
    plt.savefig(rotation_plot)
    plt.close()

    plt.gca().set_color_cycle(['red', 'green', 'blue'])
    plt.plot(data[3])
    plt.plot(data[4])
    plt.plot(data[5])
    plt.legend(['x', 'y', 'z'], loc='upper right')
    plt.ylabel('Translation (mm)')
    plt.xlabel('Volume')
    plt.savefig(translation_plot)
    plt.close()

    return translation_plot, rotation_plot


def gen_histogram(measure_file, measure):
    """Generates Histogram Image of intensities for a given input nifti file.

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
                if 'ROI_' in fname:
                    fname = fname.rsplit('ROI_')[1]
                elif 'roi_' in fname:
                    fname = fname.rsplit('roi_')[1]
                fname = 'sca_roi_' + fname.split('_')[0]
                measure = fname
            if 'sca_tempreg' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                fname = fname.split('z_maps_roi_')[1]
                fname = 'sca_mult_regression_maps_roi_' + fname.split('_')[0]
                measure = fname
            if 'dr_tempreg' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                for i in ['temp_reg_map_', 'tempreg_map_', 'tempreg_maps_', 'temp_reg_maps_']:
                    if i in fname:
                        try:
                            fname = fname.rsplit(i)[1]
                            break
                        except IndexError:
                            continue
                fname = 'dual_regression_map_'+ fname.split('_')[0]
                measure = fname
            if 'centrality' in measure.lower():
                fname = os.path.basename(os.path.splitext(os.path.splitext(file_)[0])[0])
                type_, fname = fname.split('centrality_')
                fname = type_ + 'centrality_' + fname.split('_')[0]
                measure = fname

            hist_path.append(make_histogram(file_, measure))

    else:
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
    from matplotlib import pyplot as plt 
    import numpy as np
    import nibabel as nb
    import os

    data = nb.load(measure_file).get_data()
    data_flat = data.flatten(order='F')
    y, binEdges = np.histogram(data_flat[data_flat != 0], bins=100)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])

    fig = plt.figure()
    fig.suptitle('%s intensity plot' % measure)
    plt.plot(bincenters, y, '-')
    plt.xlabel('intensity')
    plt.ylabel('# of voxels')

    png_name = os.path.join(os.getcwd(), '%s_hist_plot.png' % measure)
    fig.savefig(os.path.join(os.getcwd(), png_name))

    plt.close()
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


def drop_percent(measure_file, percent):
    """
    Zeros out voxels in measure files whose intensity doesnt fall in percent
    of voxel intensities

    Parameters
    ----------

    measure_file : string
                Input nifti file

    percent : percentage of the voxels to keep

    
    Returns
    -------

    modified_measure_file : string
                    measure_file with 1 - percent voxels zeroed out
    """

    import os
    import nibabel as nb
    import numpy as np

    img = nb.load(measure_file)
    data = img.get_data()
    
    max_val = np.percentile(data[data != 0.0], percent)
    data[data >= max_val] = 0.0

    save_img = nb.Nifti1Image(data, header=img.get_header(), affine=img.get_affine())
   
    if '.nii.gz' in measure_file:
        ext = '.nii.gz'
    else:
        ext = '.nii'

    f_name = os.path.basename(os.path.splitext(os.path.splitext(measure_file)[0])[0])
    saved_name = '%s_%d_%s' % (f_name, percent, ext)
    save_img.to_filename(saved_name)

    modified_measure_file = os.path.join(os.getcwd(),
                                         saved_name)

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
    """Draws Montage using overlay on Anatomical brain in Axial Direction,
    calls make_montage_axial.

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

    pngs = None
    if isinstance(overlay, list):
        pngs = []
        for ov in overlay:
            fname = os.path.basename(os.path.splitext(os.path.splitext(ov)[0])[0])
            pngs.append(make_montage_axial(ov, underlay,
                                           fname + '_' + png_name, cbar_name))
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
    import os
    import matplotlib
    matplotlib.rcParams.update({'font.size': 5})
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.pyplot as plt
    import nibabel as nb
    import numpy as np

    Y = nb.load(underlay).get_data()
    X = nb.load(overlay).get_data()
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)

    if 'skull_vis' in png_name:
        X[X < 20.0] = 0.0
    if 'skull_vis' in png_name or \
            't1_edge_on_mean_func_in_t1' in png_name or \
                'MNI_edge_on_mean_func_mni' in png_name:
        max_ = np.nanmax(np.abs(X.flatten()))
        X[X != 0.0] = max_

    z1, z2 = determine_start_and_end(Y, 'axial', 0.0001)
    spacing = get_spacing(6, 3, z2 - z1)
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    if ('snr' in png_name) or ('reho' in png_name) or \
            ('vmhc' in png_name) or ('sca_' in png_name) or \
            ('alff' in png_name) or ('centrality' in png_name) or \
            ('dr_tempreg' in png_name):
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="single", cbar_pad=0.2,
                         direction="row")
    else:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, direction="row")

    zz = z1
    for i in range(6*3):
        if zz >= z2:
            break
        try:
            im = grid[i].imshow(np.rot90(Y[:, :, zz]), cmap=cm.Greys_r)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "      
                  "axial montage for {0}\n\nDetails:{1}. This error might occur because of a registration error encountered while using ANTs.\
                  Please refer to the png image located in your working directory for more insight."
                  "\n".format(png_name, e))
            pass
        zz += spacing

    x, y, z = X.shape
    X[X == 0.0] = np.nan
    max_ = np.nanmax(np.abs(X.flatten()))

    zz = z1
    im = None
    for i in range(6*3):
        if zz >= z2:
            break

        try:
            if cbar_name is 'red_to_blue':
                im = grid[i].imshow(np.rot90(X[:, :, zz]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            elif cbar_name is 'green':
                im = grid[i].imshow(np.rot90(X[:, :, zz]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            else:
                im = grid[i].imshow(np.rot90(X[:, :, zz]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=- max_, vmax=max_)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "axial montage for {0}\n\nDetails:{1}.This error might occur because of a registration error encountered while using ANTs.\
                   Please refer to the image located in your working directory for more insight"
                  "\n".format(png_name, e))
            pass

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

    if 'snr' in png_name:
        cbar.ax.set_yticks(np.linspace(0, max_, 8))

    elif ('reho' in png_name) or ('vmhc' in png_name) or \
            ('sca_' in png_name) or ('alff' in png_name) or \
            ('centrality' in png_name) or ('dr_tempreg' in png_name):
        cbar.ax.set_yticks(np.linspace(-max_, max_, 8))

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
    import os
    import numpy as np

    matplotlib.rcParams.update({'font.size': 5})

    try:
        from mpl_toolkits.axes_grid1 import ImageGrid   
    except:
        from mpl_toolkits.axes_grid import ImageGrid
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import nibabel as nb

    Y = nb.load(underlay).get_data()
    X = nb.load(overlay).get_data()
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)

    if 'skull_vis' in png_name:
        X[X < 20.0] = 0.0
    if 'skull_vis' in png_name or \
            't1_edge_on_mean_func_in_t1' in png_name or \
                'MNI_edge_on_mean_func_mni' in png_name:
        max_ = np.nanmax(np.abs(X.flatten()))
        X[X != 0.0] = max_

    x1, x2 = determine_start_and_end(Y, 'sagittal', 0.0001)
    spacing = get_spacing(6, 3, x2 - x1)
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    if ('snr' in png_name) or ('reho' in png_name) or \
            ('vmhc' in png_name) or ('sca_' in png_name) or \
            ('alff' in png_name) or ('centrality' in png_name) or \
            ('dr_tempreg' in png_name):
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="single", cbar_pad=0.5,
                         direction="row")
    else:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="None", direction="row")

    xx = x1
    for i in range(6*3):
        if xx >= x2:
            break

        try:
            im = grid[i].imshow(np.rot90(Y[xx, :, :]), cmap=cm.Greys_r)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "sagittal montage for {0}\n\nDetails:{1}.This error might occur because of a registration error encountered while using ANTs\
                   Please refer to the image located in your working directory for more insight"
                  "\n".format(png_name, e))
            pass

        grid[i].get_xaxis().set_visible(False)
        grid[i].get_yaxis().set_visible(False)
        xx += spacing

    x, y, z = X.shape
    X[X == 0.0] = np.nan
    max_ = np.nanmax(np.abs(X.flatten()))
    xx = x1
    for i in range(6*3):
        if xx >= x2:
            break
        im = None

        try:
            if cbar_name is 'red_to_blue':
                im = grid[i].imshow(np.rot90(X[xx, :, :]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            elif cbar_name is 'green':
                im = grid[i].imshow(np.rot90(X[xx, :, :]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=0, vmax=max_)
            else:
                im = grid[i].imshow(np.rot90(X[xx, :, :]),
                                    cmap=cm.get_cmap(cbar_name), alpha=0.82,
                                    vmin=- max_, vmax=max_)
        except IndexError as e:
            # TODO: send this to the logger instead
            print("\n[!] QC Interface: Had a problem with creating the "
                  "sagittal montage for {0}\n\nDetails:{1}.This error might occur because of a registration error encountered while using ANTs.\
                   Please refer to the image located in your working directory for more insight"
                  "\n".format(png_name, e))
            pass

        xx += spacing

    try:
        cbar = grid.cbar_axes[0].colorbar(im)

        if 'snr' in png_name:
            cbar.ax.set_yticks(np.linspace(0, max_, 8))
        elif ('reho' in png_name) or ('vmhc' in png_name) or \
                ('sca_' in png_name) or ('alff' in png_name) or \
                ('centrality' in png_name) or ('dr_tempreg' in png_name):
            cbar.ax.set_yticks(np.linspace(-max_, max_, 8))

    except AttributeError as e:
        # TODO: send this to the logger instead
        print("\n[!] QC Interface: Had a problem with creating the "
              "sagittal montage for {0}\n\nDetails:{1}"
              "\n".format(png_name, e))
        pass

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
    import numpy as np
    from mpl_toolkits.axes_grid import ImageGrid as ImageGrid
    import matplotlib.pyplot as plt
    import nibabel as nb
    import matplotlib.cm as cm

    Y = nb.load(underlay).get_data()
    z1, z2 = determine_start_and_end(Y, 'axial', 0.0001)
    spacing = get_spacing(6, 3, z2 - z1)
    X_csf = nb.load(overlay_csf).get_data()
    X_wm = nb.load(overlay_wm).get_data()
    X_gm = nb.load(overlay_gm).get_data()
    X_csf = X_csf.astype(np.float32)
    X_wm = X_wm.astype(np.float32)
    X_gm = X_gm.astype(np.float32)
    Y = Y.astype(np.float32)

    max_csf = np.nanmax(np.abs(X_csf.flatten()))
    X_csf[X_csf != 0.0] = max_csf
    max_wm = np.nanmax(np.abs(X_wm.flatten()))
    X_wm[X_wm != 0.0] = max_wm
    max_gm = np.nanmax(np.abs(X_gm.flatten()))
    X_gm[X_gm != 0.0] = max_gm
    fig = plt.figure(1)

    try:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                          aspect=True, cbar_mode="None", direction="row")
    except:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="None", direction="row")

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

    zz = z1
    im = None
    for i in range(6*3):
        if zz >= z2:
            break
        im = grid[i].imshow(np.rot90(X_csf[:, :, zz]), cmap=cm.get_cmap('green'), alpha=0.82, vmin=0, vmax=max_csf)
        im = grid[i].imshow(np.rot90(X_wm[:, :, zz]), cmap=cm.get_cmap('blue'), alpha=0.82, vmin=0, vmax=max_wm)
        im = grid[i].imshow(np.rot90(X_gm[:, :, zz]), cmap=cm.get_cmap('red'), alpha=0.82, vmin=0, vmax=max_gm)   

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

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

    import numpy as np
    from mpl_toolkits.axes_grid import ImageGrid as ImageGrid
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import nibabel as nb

    Y = nb.load(underlay).get_data()
    x1, x2 = determine_start_and_end(Y, 'sagittal', 0.0001)
    spacing = get_spacing(6, 3, x2 - x1)
    X_csf = nb.load(overlay_csf).get_data()
    X_wm = nb.load(overlay_wm).get_data()
    X_gm = nb.load(overlay_gm).get_data()
    X_csf = X_csf.astype(np.float32)
    X_wm = X_wm.astype(np.float32)
    X_gm = X_gm.astype(np.float32)
    Y = Y.astype(np.float32)

    max_csf = np.nanmax(np.abs(X_csf.flatten()))
    X_csf[X_csf != 0.0] = max_csf
    max_wm = np.nanmax(np.abs(X_wm.flatten()))
    X_wm[X_wm != 0.0] = max_wm
    max_gm = np.nanmax(np.abs(X_gm.flatten()))
    X_gm[X_gm != 0.0] = max_gm
    x, y, z = Y.shape
    fig = plt.figure(1)
    max_ = np.max(np.abs(Y))

    try:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                          aspect=True, cbar_mode="None", direction="row")
    except:
        grid = ImageGrid(fig, 111, nrows_ncols=(3, 6), share_all=True,
                         aspect=True, cbar_mode="None", direction="row")

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

    zz = x1
    im = None
    for i in range(6*3):
        if zz >= x2:
            break

        im = grid[i].imshow(np.rot90(X_csf[zz, :, :]),
                            cmap=cm.get_cmap('green'), alpha=0.82, vmin=0,
                            vmax=max_csf)
        im = grid[i].imshow(np.rot90(X_wm[zz, :, :]),
                            cmap=cm.get_cmap('blue'), alpha=0.82, vmin=0,
                            vmax=max_wm)
        im = grid[i].imshow(np.rot90(X_gm[zz, :, :]),
                            cmap=cm.get_cmap('red'), alpha=0.82, vmin=0,
                            vmax=max_gm)

        grid[i].axes.get_xaxis().set_visible(False)
        grid[i].axes.get_yaxis().set_visible(False)
        zz += spacing

    cbar = grid.cbar_axes[0].colorbar(im)

    plt.axis("off")
    png_name = os.path.join(os.getcwd(), png_name)
    plt.savefig(png_name, dpi=200, bbox_inches='tight')
    plt.close()

    return png_name


def register_pallete(colors_file, cbar_name):

    """
    Registers color pallete to matplotlib

    Parameters
    ----------

    colors_file : string
        file containing colors in hexadecimal formats in each line

    cbar_name : string
        Proposed name for the color bar


    Returns
    -------

    None

    """

    import matplotlib.colors as col
    import matplotlib.cm as cm
    
    with open(colors_file, 'r') as f:
        colors = [c.rstrip('\r\n') for c in reversed(f.readlines())]
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
    cmd = " 3dresample -dxyz 1.0 1.0 1.0 -prefix %s " \
          "-inset %s " % (new_fname, file_)
    commands.getoutput(cmd)

    return new_fname

