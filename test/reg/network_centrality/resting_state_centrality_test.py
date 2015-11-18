# test/reg/afni_network_centrality/resting_state_centrality_test.py
#
# Daniel Clark, 2015

'''
This module performs regression testing on the outputs from the network
centrality workflows in
CPAC/network_centrality/..
'''


# Run and record memory of function
def run_and_get_max_memory(func_tuple):
    '''
    Function to run and record memory usage of a function

    Parameters
    ----------
    func_tuple : tuple
        tuple contaning the function and any arguments in the form of
        (func, *args, **kwargs)

    Returns
    -------
    max_mem_gb : float
        the high watermark of memory usage by the function specified
    '''

    # Import packages
    import memory_profiler

    # Get memory
    max_mem = memory_profiler.memory_usage(func_tuple, max_usage=True)
    max_mem_gb = max_mem[0]/1024.0

    # Return memory watermark in GB
    return max_mem_gb


# Run and test centrality
def run_and_test_centrality(datafile, template, cent_imp, num_threads, memory_gb):
    '''
    Function to init, run, and test the outputs of the network
    centrality workflow

    Parameters
    ----------
    pass_thr : float
        the correlation threshold to be greater than to ensure the new
        centrality workflow outputs are accurate
    cent_imp : string
       either 'afni' or 'cpac' - indicating the type of centrality
       implementation to test
    '''

    # Import packages
    import datetime
    import logging
    import os

    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl

    from CPAC.utils import test_init
    from CPAC.network_centrality.afni_network_centrality import \
        create_afni_centrality_wf
    from CPAC.network_centrality.resting_state_centrality import \
        create_resting_state_graphs

    # Init variables
    ident_mat = '/usr/share/fsl/5.0/etc/flirtsch/ident.mat'
    meth_dict = {'deg' : 0.001,
                 'eig' : 0.001,
                 'lfcd' : 0.6}
    thr_types = ['pval', 'sparse', 'rval']

    # Workflow base directory
    test_dir = os.path.join(os.path.expanduser('~'), 'tests', 'centrality')
    wflow = pe.Workflow(name='centrality_test_%s' % cent_imp, base_dir=test_dir)

    # Init resample node
    resamp_wflow = pe.Node(fsl.FLIRT(), name='resamp_wf')
    resamp_wflow.inputs.interp = 'trilinear'
    resamp_wflow.inputs.apply_xfm = True
    resamp_wflow.inputs.in_matrix_file = ident_mat
    resamp_wflow.inputs.in_file = datafile
    resamp_wflow.inputs.reference = template

    # Init test log file
    log_path = os.path.join(os.getcwd(), '%s_centrality_test.log' % \
                                         os.path.basename(datafile))
    cent_test_log = test_init.setup_test_logger('cent_test_log', log_path,
                                                logging.INFO, to_screen=True)
    cent_test_log.info('Running centrality correlations tests. Storing log ' \
                       'in %s...' % log_path)

    # Log parameters
    cent_test_log.info('Centrality workflow parameters:\ninput img: %s\n' \
                       'template file: %s\nallocated memory (GB): %.3f\n' \
                       'thresholds: %s' % \
                       (datafile, template, memory_gb, str(meth_dict)))

    # For each threshold type
    for m_idx, meth_type in enumerate(meth_dict.keys()):
        # For each centrality method
        for t_idx, thr_type in enumerate(thr_types):
            threshold = meth_dict[meth_type]
            wf_name = meth_type + '_' + thr_type
            # Init afni implementation
            if cent_imp == 'afni':
                cent_test_log.info('Utilizing AFNI centrality...')
                # AFNI centrality takes in sparsity as percentage
                if thr_type == 'sparse':
                    threshold = threshold*100
                cent_wflow = create_afni_centrality_wf(wf_name,
                                                  meth_type, thr_type,
                                                  num_threads, memory_gb)
            # Init C-PAC python implementation
            elif cent_imp == 'cpac':
                cent_test_log.info('Utilizing C-PAC centrality...')
                cent_wflow = create_resting_state_graphs(wf_name, memory_gb)
                cent_wflow.inputs.inputspec.method_option = m_idx
                cent_wflow.inputs.inputspec.threshold_option = t_idx
                cent_wflow.inputs.inputspec.weight_options = [True, True]
            # Otherwise, raise error
            else:
                err_msg = 'Specify either \'afni\' or \'cpac\' for centrality '\
                          'implementation type!'
                raise Exception(err_msg)

            # Connect the resampled input to centrality workflow
            wflow.connect(resamp_wflow, 'out_file', cent_wflow, 'inputspec.datafile')
            # Add rest of parameters
            cent_wflow.inputs.inputspec.template = template
            cent_wflow.inputs.inputpsec.threshold = threshold

            # Log running status
            cent_test_log.info('Running %s workflow with %s thresholding ' \
                               'with a %.3f threshold...' % \
                               (meth_type, thr_type, threshold))
            # Run the workflow
            start = datetime.datetime.now()
            # Run workflow
            #wflow.run()
            max_used_mem_gb = run_and_get_max_memory((wflow.run,))
            # Get and log stats
            runtime = (datetime.datetime.now()-start).total_seconds()
            cent_test_log.info('Exection time: %.3f seconds\n'\
                               'Memory used: %.3f GB' % (runtime, max_used_mem_gb))


# Make module executable
if __name__ == '__main__':

    # Import packages
    import argparse

    # Init argparser
    parser = argparse.ArgumentParser(description=__doc__)

    # Required arguments
    parser.add_argument('-d', '--datafile', nargs=1, required=True,
                        type=str, help='Filepath to the functional mni')
    parser.add_argument('-t', '--template', nargs=1, required=True,
                        type=str, help='Filepath to the mask template')
    parser.add_argument('-c', '--cent_imp', nargs=1, required=True,
                        type=str, help='Centrality implementation: \'afni\' or \'cpac\'')
    parser.add_argument('-n', '--num_threads', nargs=1, required=True,
                        type=int, help='Number of threads to run')
    parser.add_argument('-m', '--memory_gb', nargs=1, required=True,
                        type=float, help='Alloted memory (GB)')

    # Parse arguments
    args = parser.parse_args()

    # Init variables
    datafile = args.datafile[0]
    template = args.template[0]
    cent_imp = args.cent_imp[0]
    num_threads = args.num_threads[0]
    memory_gb = args.memory_gb[0]

    # Run and test centrality
    run_and_test_centrality(datafile, template, cent_imp, num_threads, memory_gb)