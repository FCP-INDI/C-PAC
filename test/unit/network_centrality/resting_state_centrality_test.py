# test/unit/network_centrality/resting_state_centrality_test.py
#

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Import packages
import unittest


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
        import logging
        import os
        import tempfile
        import urllib
        import yaml

        # Init variables
        self.rho_thresh = 0.99
        # Limit the amount of memory and threads for the test case
        self.mem_gb_limit = 8.0
        self.num_threads = 1
        # Number of subjects to run through centrality test/benchmark
        self.num_subs = 3
        # Identity matrix for resampling
        self.ident_mat = '/usr/share/fsl/5.0/etc/flirtsch/ident.mat'

        # Workflow base directory
        #self.base_dir = tempfile.mkdtemp()
        self.base_dir = os.getcwd()
        # Make inputs directory
        self.inputs_dir = os.path.join(self.base_dir, 'inputs')
        if not os.path.exists(self.inputs_dir):
            os.makedirs(self.inputs_dir)
        # Log path
        self.callback_logpath = os.path.join(self.base_dir,
                                             'centrality_test_callback.log')
        # Add handler to callback log file
        cb_logger = logging.getLogger('callback')
        cb_logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler(self.callback_logpath)
        cb_logger.addHandler(handler)

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

        # Init workflow
        wflow = pe.Workflow(name='%s_%s_test' % (method, thresh_option))

        # Set up iterable input node
        input_node = pe.Node(util.Function(input_names=['img_list', 'sub_idx',
                                                        'inputs_dir'],
                                           output_names=['local_path'],
                                           function=download_inputs),
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
        resamp_node.interface.estimated_memory_gb = 2.0

        # Connect input node to resample
        wflow.connect(input_node, 'local_path', resamp_node, 'in_file')

        # Init the centrality workflows
        # Init variables
        wf_name = 'cpac_%s_%s' % (method, thresh_option)
        cpac_wflow = create_resting_state_graphs(wf_name, self.mem_gb_limit)

        # Init workflow run parameters
        cpac_wflow.inputs.inputspec.method_option = method
        cpac_wflow.inputs.inputspec.threshold_option = thresh_option
        cpac_wflow.inputs.inputspec.threshold = thresh

        # If it is sparsity thresholding, put into percentage for afni
        if thresh_option == 'sparsity':
            thresh = 100*thresh
        wf_name = 'afni_%s_%s' % (method, thresh_option)
        afni_wflow = create_afni_centrality_wf(wf_name, method, thresh_option,
                                            thresh, self.num_threads, self.mem_gb_limit)

        # Connect resampled functionalin to centrality workflow
        wflow.connect(resamp_node, 'out_file', afni_wflow, 'inputspec.in_file')
        wflow.connect(resamp_node, 'out_file', cpac_wflow, 'inputspec.in_file')

        # Connect masks
        afni_wflow.inputs.inputspec.template = self.mask_path
        cpac_wflow.inputs.inputspec.template = self.mask_path

        # Collect arrays MapNnode
        merge_outputs_node = pe.JoinNode(util.Function(input_names=['cpac_field',
                                                                    'afni_field'],
                                                       output_names=['map_yaml'],
                                                       function=merge_img_paths),
                                         name='merge_img_paths',
                                         joinsource=input_node.name,
                                         joinfield=['cpac_field', 'afni_field'])

        # Connect the merge node from cpac/afni outputs
        wflow.connect(cpac_wflow, 'outputspec.centrality_outputs',
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
        from nipype.pipeline.plugins.callback_log import log_nodes_cb

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

    def _gen_scatterplot(self, map_yaml, img_desc):
        '''
        Function to generate a scatter plot of all of the images
        ran for a given centrality type

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
        out_png = os.path.join(self.base_dir, img_desc + '_scatter.png')
        plt.savefig(out_png, dpi=150)

        # Clear and close
        plt.clf()
        plt.close()

        # Return png path
        return out_png

    def _read_and_correlate(self, map_yaml):
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

    def _gen_boxplots(self, rho_dict, img_desc):
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
        out_png = os.path.join(self.base_dir, img_desc + '_boxplot.png')
        plt.savefig(out_png, dpi=200)

        # Clear and close
        plt.clf()
        plt.close()

        # Return png path
        return out_png

    def _merge_boxplots(self):
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
        for root, dirs, files in os.walk(self.base_dir):
            yamls.extend([os.path.join(root, file) for file in files \
                          if file.endswith('merged_paths.yml')])
    
        # For each yaml
        for yaml in yamls:
            img_type = yaml.split(os.path.sep)[-3].rstrip('_test')
            rho_dicts[img_type] = self._read_and_correlate(yaml)

        # Expand dict
        for centrality, rho_dict in rho_dicts.items():
            for img_type, rhos in rho_dict.items():
                merged_dict['_'.join([centrality, img_type])] = rhos

        # Generate and return the output
        out_png = self._gen_boxplots(merged_dict, 'merged')
        return out_png

    def test_merge(self):
        '''
        '''

        # Import packages
        import os

        # Get output png
        out_png = self._merge_boxplots()
        merged_exists = os.path.exists(out_png)
        err_msg = 'Merge function could not write out merged box plots png!'
        self.assertTrue(merged_exists, msg=err_msg)

    def test_degree_sparsity(self):
        '''
        Test AFNI and CPAC degree sparsity methods correlate
        '''

        # Run and correlate afni/cpac workflows
        deg_sparsity_map_yaml = \
            self._run_wf_and_map_outputs('degree', 'sparsity', 0.001)

        # Generate scatter plots
        out_png = self._gen_scatterplot(deg_sparsity_map_yaml, 'degree_sparsity')

        # Pairwise correlate images
        degree_sparsity_results = self._read_and_correlate(deg_sparsity_map_yaml)

        # Generate box plots
        out_png = self._gen_boxplots(degree_sparsity_results, 'degree_sparsity')

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

        # Run and correlate afni/cpac workflows
        degree_corr_map_yaml = \
            self._run_wf_and_map_outputs('degree', 'correlation', 0.6)

        # Generate scatter plots
        out_png = self._gen_scatterplot(degree_corr_map_yaml, 'degree_correlation')

        # Pairwise correlate images
        degree_corr_results = self._read_and_correlate(degree_corr_map_yaml)

        # Generate box plots
        out_png = self._gen_boxplots(degree_corr_results, 'degree_correlation')

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

        # Run and correlate afni/cpac workflows
        eigen_sparsity_map_yaml = \
            self._run_wf_and_map_outputs('eigenvector', 'sparsity', 0.001)

        # Generate scatter plots
        out_png = self._gen_scatterplot(eigen_sparsity_map_yaml, 'eigen_sparsity')

        # Pairwise correlate images
        eigen_sparsity_results = self._read_and_correlate(eigen_sparsity_map_yaml)

        # Generate box plots
        out_png = self._gen_boxplots(eigen_sparsity_results, 'eigen_sparsity')

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

        # Run and correlate afni/cpac workflows
        eigen_corr_map_yaml = \
            self._run_wf_and_map_outputs('eigenvector', 'correlation', 0.6)

        # Generate scatter plots
        out_png = self._gen_scatterplot(eigen_corr_map_yaml, 'eigen_correlation')

        # Pairwise correlate images
        eigen_corr_results = self._read_and_correlate(eigen_corr_map_yaml)

        # Generate box plots
        out_png = self._gen_boxplots(eigen_corr_results, 'eigen_correlation')

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

        # Run and correlate afni/cpac workflows
        lfcd_corr_map_yaml = \
            self._run_wf_and_map_outputs('lfcd', 'correlation', 0.6)

        # Generate scatter plots
        out_png = self._gen_scatterplot(lfcd_corr_map_yaml, 'lfcd_correlation')

        # Pairwise correlate images
        lfcd_corr_results = self._read_and_correlate(lfcd_corr_map_yaml)

        # Generate box plots
        out_png = self._gen_boxplots(lfcd_corr_results, 'lfcd_correlation')

        # Iterate through concordances and assert > 0.99
        for img_type, rho_list in lfcd_corr_results.items():
            err_msg = 'AFNI and C-PAC concordance: %.6f is too low for %s!'
            for rho in rho_list:
                self.assertGreater(rho, self.rho_thresh,
                                   msg=err_msg % (rho, img_type))


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()