# test/unit/network_centrality/resting_state_centrality_test.py
#
#

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Import packages
import unittest
def _print_filenames(outfile_list):
    print '\n\n!!!!!!!', outfile_list

    return outfile_list
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
        import os
        import tempfile
        import urllib
        import yaml

        # Init variables
        # Limit the amount of memory and threads for the test case
        self.mem_gb_limit = 8.0
        self.num_threads = 1
        # Number of subjects to run through centrality test/benchmark
        self.num_subs = 1
        # Identity matrix for resampling
        self.ident_mat = '/usr/share/fsl/5.0/etc/flirtsch/ident.mat'

        # Workflow base directory
        self.base_dir = tempfile.mkdtemp()
        # Make inputs directory
        self.inputs_dir = os.path.join(tempfile.mkdtemp(), 'inputs')
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

    def _init_centrality_wf(self, implement, method, thresh_option, thresh):
        '''
        '''

        # Import packages
        import nipype.interfaces.fsl as fsl
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util

        # Init workflow
        wflow = pe.Workflow(name='%s_centrality_test' % implement)

        # Set up iterable input node
        input_node = pe.Node(util.Function(input_names=['img_list', 'sub_idx',
                                                        'inputs_dir', 'callback_f'],
                                           output_names=['local_path'],
                                           function=_download_inputs),
                             name='inputspec')
        input_node.inputs.img_list = self.img_list
        input_node.inputs.inputs_dir = self.inputs_dir
        input_node.iterables = ('sub_idx', range(self.num_subs))

        # Set up resample node
        resamp_node = pe.Node(fsl.FLIRT(), name='resamp_wf', iterfield=['in_file'])
        resamp_node.inputs.interp = 'trilinear'
        resamp_node.inputs.apply_xfm = True
        resamp_node.inputs.in_matrix_file = self.ident_mat
        resamp_node.inputs.reference = self.mask_path
        resamp_node.interface.estimated_memory_gb = 2.0

        # Connect input node to resample
        wflow.connect(input_node, 'local_path', resamp_node, 'in_file')

        # Init the centrality workflow
        if implement == 'afni':
            cent_wflow = self._return_afni_centrality_wf(method, thresh_option,
                                                         thresh)
        elif implement == 'cpac':
            cent_wflow = self._return_cpac_centrality_wf(method, thresh_option,
                                                         thresh)
        # Connect in centrality workflow
        wflow.connect(resamp_node, 'out_file', cent_wflow, 'inputspec.in_file')
        cent_wflow.inputs.inputspec.template = self.mask_path

        # Define outputs node
        output_node = pe.Node(util.IdentityInterface(fields=['outfile_list',
                                                             'oned_output']),
                              name='outputspec')
        if implement == 'afni':
            wflow.connect(cent_wflow, 'outputspec.outfile_list', output_node, 'outfile_list')
        else:
            wflow.connect(cent_wflow, 'outputspec.centrality_outputs', output_node, 'outfile_list')

        # Return the complete workflow
        return wflow



    def test_degree_sparsity(self):
        '''
        '''

        # Import packages
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util

        # Init test workflow
        test_wflow = pe.Workflow('face')
        afni_wf = self._init_centrality_wf('afni', 'degree', 'sparsity', .1)
        cpac_wf = self._init_centrality_wf('cpac', 'degree', 'sparsity', .001)
        print_node = pe.Node(util.Function(input_names=['outfile_list'],
                                           output_names=['outfile_list'],
                                           function=_print_filenames),
                             name='print_node')

        test_wflow.connect(afni_wf, 'outputspec.outfile_list', print_node, 'outfile_list')

        test_wflow.run(plugin='MultiProc',
                       plugin_args={'n_procs' : self.num_threads,
                                    'memory_gb' : self.mem_gb_limit})


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()