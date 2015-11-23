# collect_cent_stats.py
#
# Author: Daniel Clark, 2015

'''
'''

# Function to parse the logs in a given directory
def parse_logs(log_dir):
    '''
    '''

    # Import packages
    import copy
    import os
    import re

    # Init variables
    logs = []
    for root, dirs, files in os.walk(log_dir):
        logs.extend([os.path.join(root, fil) for fil in files\
                     if fil.endswith('.log')])
    run_times = {'deg-sparse' : {}, 'deg-rval' : {}, 'lfcd-rval' : {}}
    memories = copy.deepcopy(run_times)
    re_pat = '\d+.\d+'

    # Iterate through logs
    for log in logs:

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
                run_times[key][log] = re.findall(re_pat, log_lines[idx+2])
                memories[key][log] = re.findall(re_pat, log_lines[idx+3])
            elif 'Utilizing C-PAC centrality' in line:
                wflow_line_sp = log_lines[idx+1].split(' ')
                meth = wflow_line_sp[1]
                thr = wflow_line_sp[4]
                key = '-'.join([meth, thr])
                run_times[key][log] = re.findall(re_pat, log_lines[idx+2])
                memories[key][log] = re.findall(re_pat, log_lines[idx+3])
            else:
                continue

    return run_times, memories


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
    lfcd_pval_runtimes = np.array([v[0] for v in \
                                   run_times_dict['lfcd-pval'].values()],
                                   dtype='float32')
    lfcd_rval_runtimes = np.array([v[0] for v in \
                                   run_times_dict['lfcd-rval'].values()],
                                   dtype='float32')
    deg_pval_runtimes = np.array([v[0] for v in \
                                   run_times_dict['deg-pval'].values()],
                                   dtype='float32')
    deg_sparse_runtimes = np.array([v[0] for v in \
                                   run_times_dict['deg-sparse'].values()],
                                   dtype='float32')
    deg_rval_runtimes = np.array([v[0] for v in \
                                   run_times_dict['deg-rval'].values()],
                                   dtype='float32')

    # Memories arrays
    lfcd_pval_mems = np.array([v[0] for v in \
                               memories_dict['lfcd-pval'].values()],
                               dtype='float32')
    lfcd_rval_mems = np.array([v[0] for v in \
                               memories_dict['lfcd-rval'].values()],
                               dtype='float32')
    deg_pval_mems = np.array([v[0] for v in \
                              memories_dict['deg-pval'].values()],
                              dtype='float32')
    deg_sparse_mems = np.array([v[0] for v in \
                                memories_dict['deg-sparse'].values()],
                                dtype='float32')
    deg_rval_mems = np.array([v[0] for v in \
                              memories_dict['deg-rval'].values()],
                              dtype='float32')

    # Histogram plots of runtimes
    plt.figure(1)
    plt.subplot(221)
    plt = gen_hist(lfcd_pval_runtimes, plt, 'lfcd-pval')
    plt.xlabel('Runtime (seconds)')
    plt.subplot(222)
    plt = gen_hist(lfcd_rval_runtimes, plt, 'lfcd-rval')
    plt.xlabel('Runtime (seconds)')

    plt.subplot(223)
    plt = gen_hist(lfcd_pval_mems, plt, 'lfcd-pval')
    plt.xlabel('Peak memory (GB)')
    plt.subplot(224)
    plt = gen_hist(lfcd_rval_mems, plt, 'lfcd-rval')
    plt.xlabel('Peak memory (GB)')

    # Degree centrality
    plt.figure(2)
    plt.subplot(231)
    plt = gen_hist(deg_pval_runtimes, plt, 'deg-pval')
    plt.xlabel('Runtime (seconds)')
    plt.subplot(232)
    plt = gen_hist(deg_sparse_runtimes, plt, 'deg-sparse')
    plt.xlabel('Runtime (seconds)')
    plt.subplot(233)
    plt = gen_hist(deg_rval_runtimes, plt, 'deg-rval')
    plt.xlabel('Runtime (seconds)')

    plt.subplot(234)
    plt = gen_hist(deg_pval_mems, plt, 'deg-pval')
    plt.xlabel('Peak memory (GB)')
    plt.subplot(235)
    plt = gen_hist(deg_sparse_mems, plt, 'deg-sparse')
    plt.xlabel('Peak memory (GB)')
    plt.subplot(236)
    plt = gen_hist(deg_rval_mems, plt, 'deg-rval')
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
    runs, mems = parse_logs(logs_dir)
    plt = plot_histograms(runs, mems)

    # Save figure
    pdf_out = os.path.join(logs_dir, 'cent_stats.pdf')
    pdf = PdfPages(pdf_out)

    for fig in xrange(1, 3):
        pdf.savefig(plt.figure(fig))
    pdf.close()

