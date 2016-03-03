# test/unit/network_centrality/resting_state_centrality_test.py
#
#

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Import packages
import unittest


# Calculate concordance correlation coefficient
def _concordance(x, y):
    '''
    Return the concordance correlation coefficient as defined by
    Lin (1989)

    Parameters
    ----------
    x : list or array
        a list of array of length N of numbers
    y : list or array
        a list of array of length N of numbers

    Returns
    -------
    rho_c : numpy.float32
        the concordance value as a float
    '''

    # Import packages
    import numpy as np

    # Usage errors check
    x_shape = np.shape(x)
    y_shape = np.shape(y)
    if len(x_shape) != 1 or len(y_shape) != 1:
        err_msg = 'Inputs must be 1D lists or arrays.'
        raise ValueError(err_msg)
    elif x_shape != y_shape:
        err_msg = 'Length of the two inputs must be equal.\n'\
                'Length of x: %d\nLength of y: %d' % (len(x), len(y))
        raise ValueError(err_msg)

    # Init variables
    x_arr = np.array(x).astype('float64')
    y_arr = np.array(y).astype('float64')

    # Get pearson correlation
    rho = np.corrcoef(x_arr, y_arr)[0][1]

    # Get stdevs
    sigma_x = np.std(x_arr)
    sigma_y = np.std(y_arr)

    # Get means
    mu_x = np.mean(x_arr)
    mu_y = np.mean(y_arr)

    # Comput condordance
    rho_c = (2*rho*sigma_x*sigma_y) /\
            (sigma_x**2 + sigma_y**2 + (mu_x-mu_y)**2)

    # Return variables
    return rho_c


def _consolidate_results(arr_dict_list):
    '''
    '''

    # Import packages

    # Init variables
    results_dict = {}

    # Build dictionary
    for arr_dict in arr_dict_list:
        sub_id = arr_dict['sub_id']
        afni_arr = arr_dict['afni_arr']
        cpac_arr = arr_dict['cpac_arr']
        img_type = arr_dict['img_type']

        # Calculate concordance
        rho_c = _concordance(afni_arr, cpac_arr)

        # Populate dictionary
        if not results_dict.has_key(img_type):
            results_dict[img_type] = {sub_id : rho_c}
        else:
            results_dict[img_type][sub_id] = rho_c

    # Return results
    return results_dict


# Compare and get runtime stats
def _get_img_arrs(afni_output, cpac_output):
    '''
    '''

    # Import packages
    import os
    import nibabel as nb

    # Sub id
    sub_id = afni_output.split(os.path.sep)[-3]
    img_type = os.path.basename(afni_output).split('.')[0]

    # Read in images as arrays and get concordances
    afni_arr = nb.load(afni_output).get_data().flatten()
    cpac_arr = nb.load(cpac_output).get_data().flatten()

    # Create dictionary
    arr_dict = {'sub_id' : sub_id,
                'afni_arr' : afni_arr,
                'cpac_arr' : cpac_arr,
                'img_type' : img_type}

    # Return array dict
    return arr_dict


# Download input data
def _download_inputs(img_list, sub_idx, inputs_dir):
    '''
    '''

    # Import packages
    import os
    import urllib

    # Init variables
    sub_rel_path = img_list[sub_idx]
    url_path = 'https://s3.amazonaws.com/fcp-indi/' + sub_rel_path
    local_path = os.path.join(inputs_dir, 'sub_%d' % (sub_idx),
                              sub_rel_path.split('/')[-1])
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
    This class is a test case for the cpac_pipeline.run() function

    Inherits
    --------
    unittest.TestCase class

    Attributes (class):
    ------------------
    see unittest.TestCase documentation

    Attributes (instance):
    ----------------------
    '''

    # setUp method for the necessary arguments to run cpac_pipeline.run
    def setUp(self):
        '''
        Method to instantiate input arguments for the
        cpac_pipeline.run() method via instance attributes
        
        Parameters
        ----------
        self : CPACPippelineRunTestCase
            a unittest.TestCase-inherited class
        '''

        # Import packages
        import logging
        import os
        import tempfile
        import urllib
        import yaml

        # Init variables
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

    def _return_cpac_centrality_wf(self, method, thresh_option, thresh):
        '''
        '''

        # Import packages
        from CPAC.network_centrality.resting_state_centrality import \
            create_resting_state_graphs

        # Init variables
        wf_name = 'cpac_%s_%s' % (method, thresh_option)
        cpac_wf = create_resting_state_graphs(wf_name, self.mem_gb_limit)

        # Init workflow run parameters
        cpac_wf.inputs.inputspec.method_option = method
        cpac_wf.inputs.inputspec.threshold_option = thresh_option
        cpac_wf.inputs.inputspec.threshold = thresh

        # Return wf
        return cpac_wf

    def _return_afni_centrality_wf(self, method, thresh_option, thresh):
        '''
        '''

        # Import packages
        from CPAC.network_centrality.afni_network_centrality import \
            create_afni_centrality_wf

        # Init variables
        wf_name = 'afni_%s_%s' % (method, thresh_option)
        afni_wf = create_afni_centrality_wf(wf_name, method, thresh_option,
                                            thresh, self.num_threads, self.mem_gb_limit)

        # Return wf
        return afni_wf

    def _init_centrality_wf(self, method, thresh_option, thresh):
        '''
        '''

        # Import packages
        import nipype.interfaces.fsl as fsl
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util

        # Init workflow
        wflow = pe.Workflow(name='%s_%s_test' % (method, thresh_option))

        # Set up iterable input node
        input_node = pe.Node(util.Function(input_names=['img_list', 'sub_idx',
                                                        'inputs_dir'],
                                           output_names=['local_path'],
                                           function=_download_inputs),
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
        afni_wflow = self._return_afni_centrality_wf(method, thresh_option,
                                                     thresh*100)

        cpac_wflow = self._return_cpac_centrality_wf(method, thresh_option,
                                                     thresh)
        # Connect in centrality workflow
        wflow.connect(resamp_node, 'out_file', afni_wflow, 'inputspec.in_file')
        wflow.connect(resamp_node, 'out_file', cpac_wflow, 'inputspec.in_file')
        # Connect mask
        afni_wflow.inputs.inputspec.template = self.mask_path
        cpac_wflow.inputs.inputspec.template = self.mask_path

        # Collect arrays MapNnode
        collect_arrs_node = \
            pe.MapNode(util.Function(input_names=['afni_outputs',
                                                  'cpac_outputs'],
                                     output_names=['arr_dict'],
                                     function=_get_img_arrs),
                        name='extract_arrs',
                        iterfield=['afni_outputs', 'cpac_outputs'])

        # Connect arrays MapNode
        wflow.connect(afni_wflow, 'outputspec.outfile_list',
                      collect_arrs_node, 'afni_outputs')
        wflow.connect(cpac_wflow, 'outputspec.centrality_outputs',
                      collect_arrs_node, 'cpac_outputs')

        # Consolidate results node
        consolidate_results_node = \
            pe.Node(util.Function(input_names=['arr_dict_list'],
                                  output_names=['concord_dict'],
                                  function=_consolidate_results),
                    name='consolidate_results')

        # Connect arrays MapNode to consolidate results
        wflow.connect(collect_arrs_node, 'arr_dict',
                      consolidate_results_node, 'arr_dict_list')

        # Return the complete workflow
        return wflow

    def test_degree_sparsity(self):
        '''
        '''

        # Import packages
        from nipype.pipeline.plugins.callback_log import log_nodes_cb

        # Init test workflow
        deg_sparse_wf = self._init_centrality_wf('degree', 'sparsity', .001)
        deg_sparse_wf.base_dir = self.base_dir

        deg_sparse_wf.run(plugin='MultiProc',
                          plugin_args={'n_procs' : self.num_threads,
                                       'memory_gb' : self.mem_gb_limit,
                                       'status_callback' : log_nodes_cb})
        print 'hi'


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()