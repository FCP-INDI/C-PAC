# collect_cent_stats.py
#
# Author: Daniel Clark, 2015

'''
'''
def collect_plot_correlations(base_dir, meth_type, thr_type):
    '''
    '''

    # Import packages
    import numpy as np
    import os

    # Init variables
    compare_dict = {}
    afni_niis = []
    folder_name = '_'.join([meth_type, thr_type])

    # Walk through to pick up the afni centrality maps
    for root, dirs, files in os.walk(base_dir):
        afni_niis.extend([os.path.join(root, fil) for fil in files \
                          if folder_name in root and \
                             fil.endswith('.nii.gz') and \
                             'centrality_test_afni' in root])

    # Get corresponding cpac maps
    cpac_niis = [nii.replace('centrality_test_afni', 'centrality_test_cpac') \
                 for nii in afni_niis]

    # Iterate through the afni maps and find correlations
    for idx, afni_nii in enumerate(afni_niis):
        cpac_nii = cpac_niis[idx]
        if os.path.exists(cpac_nii):




# Function to parse the logs in a given directory
def parse_logs(log_dir):
    '''
    '''

    # Import packages
    import copy
    import os
    import re

    # Init variables
    re_pat = '\d+.\d+'
    logs = []
    afni_run_times = {'deg-sparse' : {}, 'deg-rval' : {}, 'lfcd-rval' : {}}
    afni_memories = copy.deepcopy(afni_run_times)
    cpac_run_times = copy.deepcopy(afni_run_times)
    cpac_memories = copy.deepcopy(afni_run_times)

    # Collect logs
    for root, dirs, files in os.walk(log_dir):
        logs.extend([os.path.join(root, fil) for fil in files\
                     if fil.endswith('.log')])

    # Iterate through logs
    for log in logs:
        print 'Parsing log: %s...' % log

        # Read in log lines
        with open(log, 'r') as f:
            log_lines = f.readlines()

        # Filter out timestamps
        log_lines = [line.split(' : ')[-1] for line in log_lines]

        # Iterate through text stamps to pull out info
        for idx, line in enumerate(log_lines):
            if 'Utilizing AFNI centrality' in line:
                wflow_line_sp = log_lines[idx+1].split(' ')
                meth = wflow_line_sp[1]
                thr = wflow_line_sp[4]
                key = '-'.join([meth, thr])
                try:
                    afni_run_times[key][log] = float(re.findall(re_pat, log_lines[idx+2])[0])
                    afni_memories[key][log] = float(re.findall(re_pat, log_lines[idx+3])[0])
                except KeyError as exc:
                    err_msg = 'Key: %s not in dictionary - not analyzing!' % key
                    continue
                except IndexError as exc:
                    err_msg = 'Log either still be written to or has unfinished run, skipping...'
                    print err_msg
                    continue
            elif 'Utilizing C-PAC centrality' in line:
                wflow_line_sp = log_lines[idx+1].split(' ')
                meth = wflow_line_sp[1]
                thr = wflow_line_sp[4]
                key = '-'.join([meth, thr])
                try:
                    cpac_run_times[key][log] = float(re.findall(re_pat, log_lines[idx+2])[0])
                    cpac_memories[key][log] = float(re.findall(re_pat, log_lines[idx+3])[0])
                except KeyError as exc:
                    err_msg = 'Key: %s not in dictionary - not analyzing!' % key
                    print err_msg
                    continue
                except IndexError as exc:
                    err_msg = 'Log either still be written to or has unfinished run, skipping...'
                    print err_msg
                    continue
            else:
                continue

    return afni_run_times, afni_memories, cpac_run_times, cpac_memories


# Generate the histogram plot
def gen_hist(array, plt, title):
    '''
    '''

    # Import packages
    import numpy as np

    # Generate histogram
    hist, bins = np.histogram(array, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
    plt.title(title)

    # Return plot obj
    return plt


# Plot the histograms from ditionaries
def plot_histograms(run_times_dict, memories_dict):
    '''
    '''

    # Import packages
    import matplotlib.pyplot as plt
    import numpy as np

    # Init variables
    # Runtimes arrays
    lfcd_rval_runtimes = np.array([v for v in \
                                   run_times_dict['lfcd-rval'].values()],
                                   dtype='float32')
    deg_sparse_runtimes = np.array([v for v in \
                                   run_times_dict['deg-sparse'].values()],
                                   dtype='float32')
    deg_rval_runtimes = np.array([v for v in \
                                   run_times_dict['deg-rval'].values()],
                                   dtype='float32')

    # Memories arrays
    lfcd_rval_mems = np.array([v for v in \
                               memories_dict['lfcd-rval'].values()],
                               dtype='float32')
    deg_sparse_mems = np.array([v for v in \
                                memories_dict['deg-sparse'].values()],
                                dtype='float32')
    deg_rval_mems = np.array([v for v in \
                              memories_dict['deg-rval'].values()],
                              dtype='float32')

    # Histogram plots of runtimes
    plt.figure(1)

    plt.subplot(321)
    plt = gen_hist(deg_rval_runtimes, plt, 'deg-rval')
    plt.xlabel('Runtime (seconds)')
    plt.subplot(323)
    plt = gen_hist(deg_sparse_runtimes, plt, 'deg-sparse')
    plt.xlabel('Runtime (seconds)')

    plt.subplot(325)
    plt = gen_hist(lfcd_rval_runtimes, plt, 'lfcd-rval')
    plt.xlabel('Runtime (seconds)')

    # Memories
    plt.subplot(322)
    plt = gen_hist(deg_rval_mems, plt, 'deg-rval')
    plt.xlabel('Peak memory (GB)')
    plt.subplot(324)
    plt = gen_hist(deg_sparse_mems, plt, 'deg-sparse')
    plt.xlabel('Peak memory (GB)')

    plt.subplot(326)
    plt = gen_hist(lfcd_rval_mems, plt, 'lfcd-rval')
    plt.xlabel('Peak memory (GB)')


    # Return plot object
    return plt


# Main executable
if __name__ == '__main__':

    # Import packages
    import argparse
    import os

    from matplotlib.backends.backend_pdf import PdfPages

    # Init argparser
    parser = argparse.ArgumentParser(description=__doc__)
    # Required arguments
    parser.add_argument('-d', '--logs_dir', nargs=1, required=True,
                        type=str, help='Base directory of centrality logs')

    # Get arguments
    args = parser.parse_args()

    # Init variables
    logs_dir = args.logs_dir[0]

    # Call functions to plot
    afni_runs, afni_mems, cpac_runs, cpac_mems = parse_logs(logs_dir)

    # Generate and save figure
    afni_plt = plot_histograms(afni_runs, afni_mems)
    pdf_out = os.path.join(logs_dir, 'afni_cent_stats.pdf')
    pdf = PdfPages(pdf_out)
    pdf.savefig(afni_plt.figure(1))
    pdf.close()

    # Generate and save figure
    cpac_plt = plot_histograms(cpac_runs, cpac_mems)
    pdf_out = os.path.join(logs_dir, 'cpac_cent_stats.pdf')
    pdf = PdfPages(pdf_out)
    pdf.savefig(cpac_plt.figure(1))
    pdf.close()
