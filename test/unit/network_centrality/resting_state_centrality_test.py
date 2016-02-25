# test/unit/network_centrality/resting_state_centrality_test.py
#
#

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Import packages
import unittest

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

        # Limit the amount of memory for the test case
        self.mem_gb_limit = 8.0
        # Workflow base directory
        self.base_dir = tempfile.mkdtemp()
        # Make inputs directory
        self.inputs_dir = os.path.join(tempfile.mkdtemp(), 'inputs')
        if not self.inputs_dir:
            os.makedirs(self.inputs_dir)
        # Number of subjects to run through centrality test/benchmark
        self.num_subs = 1
        

        # Download sublist yaml
        sublist_url = 'https://s3.amazonaws.com/fcp-indi/data/test_resources/'\
                      'cpac_resources/settings/resources/s3_subs_iba_trt.yml'
        sublist_filename = self.sublist_url.split('/')[-1]
        sublist_yaml = os.path.join(self.inputs_dir, sublist_filename)

        # Download sublist yaml
        if not os.path.exists(sublist_yaml):
            print('Downloading %s to %s\n' % (sublist_url, sublist_yaml))
            urllib.urlretrieve(sublist_url, sublist_yaml,
                               reporthook=self._dl_progress)

        # Download centrality mask file
        mask_url = ''

        # Read in yaml and set img_list attribute
        with open(sublist_yaml, 'r') as yml_in:
            self.img_list = yaml.load(yml_in)

    def _dl_progress(self, count, block_size, total_size):
        '''
        '''

        # Import packages
        import sys

        # Write percentage to stdout
        percent = int(count*block_size*100/total_size)
        sys.stdout.write('\r' + '...%d%%' % percent)
        sys.stdout.flush()

    def _download_inputs(self, sub_idx):
        '''
        '''

        # Import packages
        import os
        import urllib

        # Init variables
        sub_rel_path = self.img_list[sub_idx]
        url_path = 'https://s3.amazonaws.com/fcp-indi/' + sub_rel_path
        local_path = os.path.join(self.inputs_dir, 'sub_%d' % (sub_idx),
                                  sub_rel_path.split('/')[-1])

        # Check to see if we should download
        if not os.path.exists(local_path):
            print 'Downloading %s to %s\n' % (url_path, local_path)
            urllib.urlretrieve(url_path, local_path,
                               reporthook=self._dl_progress)

        # Return local path
        return local_path

    def _return_cpac_centrality_wf(self, method_option, threshold_option,
                                 threshold, weight_options=[True, True]):
        '''
        '''

        # Import packages
        from CPAC.network_centrality.resting_state_centrality import \
            create_resting_state_graphs

        # Init variables
        wf_name = 'cpac_%s_%s_%f' % (method_option, threshold_option, threshold)
        cpac_wf = create_resting_state_graphs(wf_name, self.mem_gb_limit)

        # Init workflow run parameters
        cpac_wf.inputs.inputspec.method_option = method_option
        cpac_wf.inputs.inputspec.threshold_option = threshold_option
        cpac_wf.inputs.inputspec.threshold = threshold
        cpac_wf.inputs.inputspec.weight_options = weight_options

        # Return wf
        return cpac_wf

    def _init_centrality_wf(self, implement):
        '''
        '''

        # Import packages
        import nipype.interfaces.fsl as fsl
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util

        # Init workflow
        wflow = pe.Workflow(name='%s_centrality_test' % implement,
                            base_dir=self.base_dir)

        # Set up iterable input node
        input_node = pe.Node(util.Function(input_names=['sub_idx'],
                                           output_names=['local_path'],
                                           function=self._download_inputs),
                             name='inputspec')
        input_node.iterables = ('sub_idx', range(self.num_subs))

        # 
        resamp_wflow = pe.Node(fsl.FLIRT(), name='resamp_wf', iterfield=['in_file'])
        resamp_wflow.inputs.interp = 'trilinear'
        resamp_wflow.inputs.apply_xfm = True
        resamp_wflow.inputs.in_matrix_file = ident_mat
        resamp_wflow.inputs.reference = template
        resamp_wflow.interface.estimated_memory_gb = 2.0


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()