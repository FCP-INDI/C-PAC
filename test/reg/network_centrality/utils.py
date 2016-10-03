# test/reg/network_centrality/utils.py
#

'''
This module performs benchmarking and regression analysis tools for the
centrality calculations in CPAC and AFNI
'''


def download_inputs(img_list, sub_idx, inputs_dir):
    '''
    Function to download functional images for input to the centrality
    workflow from the fcp-indi S3 bucket on AWS

    Parameters
    ----------
    img_list : list
        a list of S3-relative file paths to download
    sub_idx : integer
        the index indicating which subject to download from the list
    inputs_dir : string
        filepath to the directory for storing local input files

    Returns
    -------
    local_path : string
        filepath to the locally-downloaded file
    '''

    # Import packages
    import os
    import urllib

    # Init variables
    sub_rel_path = img_list[sub_idx]
    url_path = 'https://s3.amazonaws.com/fcp-indi/' + sub_rel_path
    local_path = os.path.join(inputs_dir, 'sub_%d' % (sub_idx),
                              sub_rel_path.split('/')[-1])
    # Make directory for downloaded file
    if not os.path.exists(os.path.dirname(local_path)):
        os.makedirs(os.path.dirname(local_path))

    # Check to see if we should download
    if not os.path.exists(local_path):
        print 'Downloading %s to %s...\n' % (url_path, local_path)
        urllib.urlretrieve(url_path, local_path)

    # Return local path
    return local_path


def merge_img_paths(cpac_field, afni_field):
    '''
    Function for nipype JoinNode that will merge the lists of outputs
    produced by node iterable forking in the centrality workflow

    Parameters
    ----------
    cpac_field : list
        a list where each element is a list of binarized and weighted
        image files for a given subject - CPAC workflow outputs
    afni_field : list
        a list where each element is a list of binarized and weighted
        image files for a given subject - AFNI workflow outputs

    Returns
    -------
    map_yaml : string
        filepath to the mapping dictoinary yaml file between afni
        and cpac centrality outputs
    '''

    # Import packages
    import yaml
    import os

    # Path mapping dict
    path_map = {}

    # Iterate through the input lists to merge into a dictionary
    for idx, paths in enumerate(cpac_field):
        cpac_bin = paths[0]
        cpac_wght = paths[1]
        afni_paths = afni_field[idx]
        afni_bin = afni_paths[0]
        afni_wght = afni_paths[1]
        path_map[cpac_bin] = afni_bin
        path_map[cpac_wght] = afni_wght

    # Write dictionary to working dir
    with open('merged_paths.yml', 'w') as fout:
        fout.write(yaml.dump(path_map))

    # Return the mapping dictionary path
    map_yaml = os.path.abspath('merged_paths.yml')
    return map_yaml


def read_and_correlate(map_yaml):
    '''
    Read and correlate the paths from the mapping dictionary
    yaml

    Parameters
    ----------
    map_yaml : string
        filepath to the mapping dictoinary yaml file between afni
        and cpac centrality outputs

    Returns
    -------
    rho_dict : dictionary
        dictionary of pairwise concordances between the afni and
        cpac centrality implementations
    '''

    # Import packages
    import os
    import yaml
    import nibabel as nib

    from CPAC.utils.test_init import concordance

    # Init variables
    map_dict = yaml.load(open(map_yaml, 'r'))
    rho_dict = {}

    # Iteratae through mapping dict
    for cpac_nii, afni_nii in map_dict.items():
        cpac_arr = nib.load(cpac_nii).get_data()
        afni_arr = nib.load(afni_nii).get_data()
        rho = concordance(cpac_arr.flatten(), afni_arr.flatten())
        img_type = os.path.split(cpac_nii)[-1].split('.')[0].split('_')[-1]
        if rho_dict.has_key(img_type):
            rho_dict[img_type].append(rho)
        else:
            rho_dict[img_type] = [rho]

    # Return the concordance dictionary
    return rho_dict


def gen_boxplots(base_dir, rho_dict, img_desc):
    '''
    Function to generate a scatter plot of all of the images
    ran for a given centrality type

    Parameters
    ----------
    rho_dict : dictionary
        dictionary where keys are strings containing centrality
        output type and values are arrays of concordances
    img_desc : string
        a string describing the type of images being analyzed; this
        string will be used to title and name the plot png file

    Returns
    -------
    out_png : string
        filepath to the produced output png file
    '''

    # Import packages
    from collections import OrderedDict
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    # Set up plot
    rho_dict = OrderedDict(sorted(rho_dict.items()))
    plt.boxplot(rho_dict.values())
    plt.xticks(range(1, len(rho_dict)+1), rho_dict.keys(), rotation=45)
    plt.ylim([np.min(rho_dict.values())-0.1, 1.1])
    plt.title('CPAC-AFNI pairwise concordance boxplots: %s' % img_desc)
    plt.xlabel('Image type')
    plt.ylabel('Concordance')
    plt.grid()

    # Save figure
    fig = plt.gcf()
    fig.set_size_inches(16, 12)
    fig.tight_layout()

    # Output png
    out_png = os.path.join(base_dir, img_desc + '_boxplot.png')
    plt.savefig(out_png, dpi=200)

    # Clear and close
    plt.clf()
    plt.close()

    # Return png path
    return out_png


def merge_boxplots(base_dir):
    '''
    Function which merges the pre-computed box-plots via the
    merged_paths.yml files in the working directory
    '''

    # Import packages
    import os
    import yaml

    # Init variables
    yamls = []
    rho_dicts = {}
    merged_dict = {}

    # Collect yamls in base directory
    for root, dirs, files in os.walk(base_dir):
        yamls.extend([os.path.join(root, file) for file in files \
                      if file.endswith('merged_paths.yml')])

    # For each yaml
    for yaml in yamls:
        img_type = yaml.split(os.path.sep)[-3].rstrip('_test')
        rho_dicts[img_type] = read_and_correlate(yaml)

    # Expand dict
    for centrality, rho_dict in rho_dicts.items():
        for img_type, rhos in rho_dict.items():
            merged_dict['_'.join([centrality, img_type])] = rhos

    # Generate and return the output
    out_png = gen_boxplots(base_dir, merged_dict, 'merged')
    return out_png


def gen_scatterplot(base_dir, map_yaml, img_desc):
    '''
    Function to generate a scatter plot of all of the images
    ran for a given centrality type

    Parameters
    ----------
    base_dir : string
        filepath to the directory where the output should be written
    map_yaml : string
        filepath to the mapping dictoinary yaml file between afni
        and cpac centrality outputs
    img_desc : string
        a string describing the type of images being analyzed; this
        string will be used to title and name the plot png file

    Returns
    -------
    out_png : string
        filepath to the produced output png file
    '''

    # Import packages
    import os
    import yaml
    import numpy as np
    import matplotlib.pyplot as plt
    import nibabel as nib

    # Init variables
    cpac_bin = np.empty(0)
    afni_bin = np.empty(0)
    cpac_wght = np.empty(0)
    afni_wght = np.empty(0)
    map_dict = yaml.load(open(map_yaml, 'r'))

    # Extract and build pairwise arrays
    for cpac, afni in map_dict.items():
        cpac_arr = nib.load(cpac).get_data().flatten()
        afni_arr = nib.load(afni).get_data().flatten()
        if 'binarize' in cpac:
            cpac_bin = np.concatenate((cpac_bin, cpac_arr))
            afni_bin = np.concatenate((afni_bin, afni_arr))
        else:
            cpac_wght = np.concatenate((cpac_wght, cpac_arr))
            afni_wght = np.concatenate((afni_wght, afni_arr))

    # Get best fit lines and set up equation strs
    bin_fit = np.polyfit(cpac_bin, afni_bin, 1)
    wght_fit = np.polyfit(cpac_wght, afni_wght, 1)
    bin_eq_str = 'y = %.4fx + %.4f' % (bin_fit[0], bin_fit[1])
    wght_eq_str = 'y = %.4fx + %.4f' % (wght_fit[0], wght_fit[1])

    # Build plot
    bin_pts = plt.scatter(cpac_bin, afni_bin, color='b', alpha=0.4,
                          label='Binarized')
    wght_pts = plt.scatter(cpac_wght, afni_wght, color='r', alpha=0.4,
                           label='Weighted')
    plt.legend(handles=[bin_pts, wght_pts])
    plt.text(0.25*cpac_bin.max(), 0.75*afni_bin.max(), bin_eq_str, color='b')
    plt.text(0.75*cpac_bin.max(), 0.25*afni_bin.max(), wght_eq_str, color='r')
    plt.xlabel('C-PAC values')
    plt.ylabel('AFNI values')
    plt.title('CPAC-AFNI image intensities scatterplot: %s' % img_desc)
    plt.grid()

    # Save figure
    fig = plt.gcf()
    fig.set_size_inches(14, 9)

    # Output png
    out_png = os.path.join(base_dir, img_desc + '_scatter.png')
    plt.savefig(out_png, dpi=150)

    # Clear and close
    plt.clf()
    plt.close()

    # Return png path
    return out_png


def gen_runtime_plots(callback_log, img_desc):
    '''
    Function to generate box plots of the runtime and memory usage
    for the CPAC vs AFNI implementations

    Parameters
    ----------
    map_yaml : string
        filepath to the mapping dictoinary yaml file between afni
        and cpac centrality outputs
    img_desc : string
        a string describing the type of images being analyzed; this
        string will be used to title and name the plot png file

    Returns
    -------
    out_png : string
        filepath to the produced output png file
    '''

    # Import packages
    import matplotlib.pyplot as plt
    from nipype.utils.draw_gantt_chart import log_to_dict

    # Init variables
    nodes_list = log_to_dict(callback_log)

    # Isolate cpac and afni finish nodes
    cpac_nodes = [node for node in nodes_list if \
                  node['name'] == 'calculate_centrality' and node.has_key('finish')]
    afni_nodes = [node for node in nodes_list if \
                  node['name'] == 'afni_centrality' and node.has_key('finish')]

    # Gather cpac stats
    cpac_mems = [float(node['runtime_memory_gb']) for node in cpac_nodes]
    cpac_times = [node['duration'] for node in cpac_nodes]

    # Gather afni stats
    afni_mems = [float(node['runtime_memory_gb']) for node in afni_nodes]
    afni_times = [node['duration'] for node in afni_nodes]

    # Init pngs
    mem_png = callback_log.split('.log')[0] + '_memories.png'
    runtime_png = callback_log.split('.log')[0] + '_runtimes.png'

    # Set up plot - memory
    plt.boxplot([cpac_mems, afni_mems])
    plt.xticks(range(1, 3), ['C-PAC', 'AFNI'])
    plt.title('CPAC-AFNI memory usage for: %s' % img_desc)
    plt.xlabel('Implementation')
    plt.ylabel('Memory used (GB)')
    plt.grid()

    # Save figure
    fig = plt.gcf()
    fig.set_size_inches(16, 12)
    fig.tight_layout()

    plt.savefig(mem_png, dpi=200)

    # Clear and close
    plt.clf()
    plt.close()

    # Set up plot - runtime
    plt.boxplot([cpac_times, afni_times])
    plt.xticks(range(1, 3), ['C-PAC', 'AFNI'])
    plt.title('CPAC-AFNI runtime for: %s' % img_desc)
    plt.xlabel('Implementation')
    plt.ylabel('Runtime (seconds)')
    plt.grid()

    # Save figure
    fig = plt.gcf()
    fig.set_size_inches(16, 12)
    fig.tight_layout()

    plt.savefig(runtime_png, dpi=200)

    # Clear and close
    plt.clf()
    plt.close()
