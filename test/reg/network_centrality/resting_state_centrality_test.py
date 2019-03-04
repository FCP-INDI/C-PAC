# test/unit/network_centrality/resting_state_centrality_test.py
#

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Import packages
import unittest
import utils


# Test case for the run function
class RestingStateCentralityTestCase(unittest.TestCase):
    '''
    This is a test case for the CPAC/network_centrality subpackage
    '''

    # setUp method for the necessary arguments to run cpac_pipeline.run
    def setUp(self):
        '''
        Method to instantiate input arguments for the
        test case
        '''

        # Import packages
        import os
        import tempfile
        import urllib
        import yaml

        # Init variables
        self.rho_thresh = 0.99
        # Limit the amount of memory and threads for the test case
        self.mem_gb_limit = 50
        self.num_threads = 32
        # Number of subjects to run through centrality test/benchmark
        self.num_subs = 198
        # Identity matrix for resampling
        self.ident_mat = '/usr/share/fsl/5.0/etc/flirtsch/ident.mat'

        # Workflow base directory
        #self.base_dir = tempfile.mkdtemp()
        self.base_dir = '/home/ubuntu/centrality_results'
        # Make inputs directory
        self.inputs_dir = os.path.join(self.base_dir, 'inputs')
        if not os.path.exists(self.inputs_dir):
            os.makedirs(self.inputs_dir)

        # Init sublist yaml
        sublist_url = 'https://s3.amazonaws.com/fcp-indi/data/test_resources/'\
                      'cpac_resources/settings/resources/s3_subs_iba_trt.yml'
        sublist_filename = sublist_url.split('/')[-1]
        sublist_yaml = os.path.join(self.inputs_dir, sublist_filename)
        # Download sublist yaml
        if not os.path.exists(sublist_yaml):
            print('Downloading %s to %s\n' % (sublist_url, sublist_yaml))
            urllib.urlretrieve(sublist_url, sublist_yaml)

        # Init centrality mask file
        mask_url = 'https://s3.amazonaws.com/fcp-indi/data/test_resources/'\
                   'cpac_resources/settings/resources/benchmark_centrality_mask.nii.gz'
        mask_filename = mask_url.split('/')[-1]
        mask_path = os.path.join(self.inputs_dir, mask_filename)
        # Download centrality mask file
        if not os.path.exists(mask_path):
            print('Downloading %s to %s\n' % (mask_url, mask_path))
            urllib.urlretrieve(mask_url, mask_path)
        self.mask_path = mask_path

        # Read in yaml and set img_list attribute
        with open(sublist_yaml, 'r') as yml_in:
            self.img_list = yaml.load(yml_in)

    def _init_logger(self, centrality_method):
        '''
        Function to init logger filehandler path to write to method
        specific log name

        Parameters
        ----------
        centrality_method : string
            string describing the centrality method/threshold option

        Returns
        -------
        cb_logger : logging.Logger obj
            logger which writes 'callback' log results to file
        '''

        # Import packages
        import logging
        import os

        # Log path
        callback_logpath = os.path.join(self.base_dir, '%s.log' % centrality_method)
        # Add handler to callback log file
        cb_logger = logging.getLogger('callback')
        cb_logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler(callback_logpath)
        cb_logger.addHandler(handler)

        # Return logger
        return cb_logger

    def _init_centrality_wf(self, method, thresh_option, thresh):
        '''
        Create and return the AFNI/CPAC centrality run and merge
        workflow

        Parameters
        ----------
        method : string
            options are 'degree', 'eigenvector', 'lfcd'
        thresh_option : string
            options are 'sparsity', 'correlation'
        thresh : float
            threshold for simliarity matrix

        Returns
        -------
        wflow : nipype Workflow
            the complete workflow for running centrality
        '''

        # Import packages
        import nipype.interfaces.fsl as fsl
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util

        from CPAC.network_centrality.resting_state_centrality import \
            create_resting_state_graphs
        from CPAC.network_centrality.afni_network_centrality import \
            create_afni_centrality_wf
        import eigenvector_golden

        # Init workflow
        wflow = pe.Workflow(name='%s_%s_test' % (method, thresh_option))

        # Set up iterable input node
        input_node = pe.Node(util.Function(input_names=['img_list', 'sub_idx',
                                                        'inputs_dir'],
                                           output_names=['local_path'],
                                           function=utils.download_inputs),
                             name='inputspec')
        input_node.inputs.img_list = self.img_list
        input_node.inputs.inputs_dir = self.inputs_dir
        input_node.iterables = ('sub_idx', range(self.num_subs))

        # Set up resample node
        resamp_node = pe.Node(fsl.FLIRT(), name='resamp_wf')
        resamp_node.inputs.interp = 'trilinear'
        resamp_node.inputs.apply_xfm = True
        resamp_node.inputs.in_matrix_file = self.ident_mat
        resamp_node.inputs.reference = self.mask_path

        # Connect input node to resample
        wflow.connect(input_node, 'local_path', resamp_node, 'in_file')

        # Init the centrality workflows
        # Init variables
        wf_name = 'cpac_%s_%s' % (method, thresh_option)
        golden_node = pe.Node(util.Function(input_names=['nii_path',
                                                         'mask_path',
                                                         'thresh_type',
                                                         'thresh_val'],
                                            output_names=['centrality_outputs'],
                                            function=eigenvector_golden.eigen_centrality),
                              name=wf_name)
        golden_node.inputs.mask_path = self.mask_path
        golden_node.inputs.thresh_type = thresh_option
        golden_node.inputs.thresh_val = thresh
        wflow.connect(resamp_node, 'out_file', golden_node, 'nii_path')
        #cpac_wflow = create_resting_state_graphs(wf_name, 4.0)

        # Init workflow run parameters
        #cpac_wflow.inputs.inputspec.method_option = method
        #cpac_wflow.inputs.inputspec.threshold_option = thresh_option
        #cpac_wflow.inputs.inputspec.threshold = thresh

        # If it is sparsity thresholding, put into percentage for afni
        if thresh_option == 'sparsity':
            thresh = 100*thresh
        wf_name = 'afni_%s_%s' % (method, thresh_option)
        afni_wflow = create_afni_centrality_wf(wf_name, method, thresh_option,
                                               thresh, 1, 4.0)

        # Connect resampled functionalin to centrality workflow
        wflow.connect(resamp_node, 'out_file', afni_wflow, 'inputspec.in_file')
        #wflow.connect(resamp_node, 'out_file', cpac_wflow, 'inputspec.in_file')

        # Connect masks
        afni_wflow.inputs.inputspec.template = self.mask_path
        #cpac_wflow.inputs.inputspec.template = self.mask_path

        # Collect arrays MapNnode
        merge_outputs_node = pe.JoinNode(util.Function(input_names=['cpac_field',
                                                                    'afni_field'],
                                                       output_names=['map_yaml'],
                                                       function=utils.merge_img_paths),
                                         name='merge_img_paths',
                                         joinsource=input_node.name,
                                         joinfield=['cpac_field', 'afni_field'])

        # Connect the merge node from cpac/afni outputs
        #wflow.connect(cpac_wflow, 'outputspec.centrality_outputs',
        #              merge_outputs_node, 'cpac_field')
        wflow.connect(golden_node, 'centrality_outputs',
                      merge_outputs_node, 'cpac_field')
        wflow.connect(afni_wflow, 'outputspec.outfile_list',
                      merge_outputs_node, 'afni_field')

        # Return the complete workflow
        return wflow

    def _run_wf_and_map_outputs(self, method, thresh_option, thresh):
        '''
        Build and run the workflow for the desired centrality options
        and build the pairwise output mappings

        Parameters
        ----------
        method : string
            options are 'degree', 'eigenvector', 'lfcd'
        thresh_option : string
            options are 'sparsity', 'correlation'
        thresh : float
            threshold for simliarity matrix

        Returns
        -------
        map_yaml : string
            filepath to the mapping dictoinary yaml file between afni
            and cpac centrality outputs
        '''

        # Import packages
        import os
        try:
            from nipype.pipeline.plugins.callback_log import log_nodes_cb
        except:
            from nipype.utils.profiler impor log_nodes_cb

        # Init workflow
        centrality_wf = self._init_centrality_wf(method, thresh_option, thresh)
        centrality_wf.base_dir = self.base_dir

        centrality_wf.run(plugin='MultiProc',
                          plugin_args={'n_procs' : self.num_threads,
                                       'memory_gb' : self.mem_gb_limit,
                                       'status_callback' : log_nodes_cb})

        # Formulate mapping dictionary path
        map_yaml = os.path.join(self.base_dir, centrality_wf.name,
                                'merge_img_paths', 'merged_paths.yml')

        # Return the concordnace dictionary
        return map_yaml

    def test_degree_sparsity(self):
        '''
        Test AFNI and CPAC degree sparsity methods correlate
        '''

        # Init callback logger
        cb_logger = self._init_logger('degree_sparsity')

        # Run and correlate afni/cpac workflows
        deg_sparsity_map_yaml = \
            self._run_wf_and_map_outputs('degree', 'sparsity', 0.001)

        # De-init callback logger
        cb_logger.removeHandler(cb_logger.handlers[0])

        # Generate scatter plots
        out_png = utils.gen_scatterplot(self.base_dir, deg_sparsity_map_yaml,
                                        'degree_sparsity')

        # Pairwise correlate images
        degree_sparsity_results = utils.read_and_correlate(deg_sparsity_map_yaml)

        # Generate box plots
        out_png = utils.gen_boxplots(self.base_dir, degree_sparsity_results,
                                     'degree_sparsity')

        # Iterate through concordances and assert > 0.99
        for img_type, rho_list in degree_sparsity_results.items():
            err_msg = 'AFNI and C-PAC concordance: %.6f is too low for %s!'
            for rho in rho_list:
                self.assertGreater(rho, self.rho_thresh,
                                   msg=err_msg % (rho, img_type))

    def test_degree_correlation(self):
        '''
        Test AFNI and CPAC degree correlation methods correlate
        '''

        # Init callback logger
        cb_logger = self._init_logger('degree_correlation')

        # Run and correlate afni/cpac workflows
        degree_corr_map_yaml = \
            self._run_wf_and_map_outputs('degree', 'correlation', 0.6)

        # De-init callback logger
        cb_logger.removeHandler(cb_logger.handlers[0])

        # Generate scatter plots
        out_png = utils.gen_scatterplot(self.base_dir, degree_corr_map_yaml,
                                        'degree_correlation')

        # Pairwise correlate images
        degree_corr_results = utils.read_and_correlate(degree_corr_map_yaml)

        # Generate box plots
        out_png = utils.gen_boxplots(self.base_dir, degree_corr_results,
                                     'degree_correlation')

        # Iterate through concordances and assert > 0.99
        for img_type, rho_list in degree_corr_results.items():
            err_msg = 'AFNI and C-PAC concordance: %.6f is too low for %s!'
            for rho in rho_list:
                self.assertGreater(rho, self.rho_thresh,
                                   msg=err_msg % (rho, img_type))

    def test_eigen_sparsity(self):
        '''
        Test AFNI and CPAC eigenvector sparsity methods correlate
        '''

        # Init callback logger
        cb_logger = self._init_logger('eigen_sparsity')

        # Run and correlate afni/cpac workflows
        eigen_sparsity_map_yaml = \
            self._run_wf_and_map_outputs('eigenvector', 'sparsity', 0.001)

        # De-init callback logger
        cb_logger.removeHandler(cb_logger.handlers[0])

        # Generate scatter plots
        out_png = utils.gen_scatterplot(self.base_dir, eigen_sparsity_map_yaml,
                                        'eigen_sparsity')

        # Pairwise correlate images
        eigen_sparsity_results = utils.read_and_correlate(eigen_sparsity_map_yaml)

        # Generate box plots
        out_png = utils.gen_boxplots(self.base_dir, eigen_sparsity_results,
                                     'eigen_sparsity')

        # Iterate through concordances and assert > 0.99
        for img_type, rho_list in eigen_sparsity_results.items():
            err_msg = 'AFNI and C-PAC concordance: %.6f is too low for %s!'
            for rho in rho_list:
                self.assertGreater(rho, self.rho_thresh,
                                   msg=err_msg % (rho, img_type))

    def test_eigen_correlation(self):
        '''
        Test AFNI and CPAC eigenvector correlation methods correlate
        '''

        # Init callback logger
        cb_logger = self._init_logger('eigen_correlation')

        # Run and correlate afni/cpac workflows
        eigen_corr_map_yaml = \
            self._run_wf_and_map_outputs('eigenvector', 'correlation', 0.6)

        # De-init callback logger
        cb_logger.removeHandler(cb_logger.handlers[0])

        # Generate scatter plots
        out_png = utils.gen_scatterplot(self.base_dir, eigen_corr_map_yaml,
                                        'eigen_correlation')

        # Pairwise correlate images
        eigen_corr_results = utils.read_and_correlate(eigen_corr_map_yaml)

        # Generate box plots
        out_png = utils.gen_boxplots(self.base_dir, eigen_corr_results,
                                     'eigen_correlation')

        # Iterate through concordances and assert > 0.99
        for img_type, rho_list in eigen_corr_results.items():
            err_msg = 'AFNI and C-PAC concordance: %.6f is too low for %s!'
            for rho in rho_list:
                self.assertGreater(rho, self.rho_thresh,
                                   msg=err_msg % (rho, img_type))

    def test_lfcd_correlation(self):
        '''
        Test AFNI and CPAC lfcd correlation methods correlate
        '''

        # Init callback logger
        cb_logger = self._init_logger('lfcd_correlation')

        # Run and correlate afni/cpac workflows
        lfcd_corr_map_yaml = \
            self._run_wf_and_map_outputs('lfcd', 'correlation', 0.6)

        # De-init callback logger
        cb_logger.removeHandler(cb_logger.handlers[0])

        # Generate scatter plots
        out_png = utils.gen_scatterplot(self.base_dir, lfcd_corr_map_yaml,
                                        'lfcd_correlation')

        # Pairwise correlate images
        lfcd_corr_results = utils.read_and_correlate(lfcd_corr_map_yaml)

        # Generate box plots
        out_png = utils.gen_boxplots(self.base_dir, lfcd_corr_results, 'lfcd_correlation')

        # Iterate through concordances and assert > 0.99
        for img_type, rho_list in lfcd_corr_results.items():
            err_msg = 'AFNI and C-PAC concordance: %.6f is too low for %s!'
            for rho in rho_list:
                self.assertGreater(rho, self.rho_thresh,
                                   msg=err_msg % (rho, img_type))


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()
