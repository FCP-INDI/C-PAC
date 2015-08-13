# test/unit/network_centrality/test_resting_state_centrality.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in
CPAC/network_centrality/resting_state_centrality.py
'''

# Import packages
import unittest

# Test case for the run function
class CentralityWorkflowTestCase(unittest.TestCase):
    '''
    This class is a test case for the resting_state_centrality workflow

    Inherits
    --------
    unittest.TestCase class

    Attributes (class):
    -------------------
    see unittest.TestCase documentation

    Attributes (instance):
    ----------------------
    wflow : nipype.pipeline.engine.Workflow instance
        the C-PAC resting_state_centrality workflow
    '''

    # Set up workflow
    def setUp(self):
        '''
        Function which inits a nipype workflow using network_centrality's
        create_resting_state_graphs() function to be used for unit testing

        Parameters
        ----------
        self : CentralityWorkflowTestCase
            a unittest.TestCase-inherited class

        Returns
        -------
        None
            this function does not return any values, but populates the
            instance attributes for:

            self.output_dirs : list
                a list of filepath strings of output base directories
            self.test_wflows : list
                a list of nipype.pipeline.engine.Workflow objects
        '''

        # Import packages
        import os
        from CPAC.network_centrality import resting_state_centrality
        from CPAC.utils import tests_init

        # Init variables
        test_wflows = {}
        output_dirs = tests_init.return_subj_measure_dirs('network_centrality')

        # Get template file path
        settings_dir = tests_init.return_resource_subfolder('settings')
        mask_template = os.path.join(settings_dir, 'resources', 'centrality',
                                     'benchmark_centrality_mask.nii.gz')

        # Init workflows for each centrality output directory
        for out_dir in output_dirs:
            # Init workflows
            wflow = resting_state_centrality.\
                    create_resting_state_graphs(allocated_memory=1.0)

            # Grab functional mni as subject input for that strategy
            func_mni_dir = out_dir.replace('network_centrality', 'functional_mni')
            func_mni = os.path.join(func_mni_dir, 'functional_mni_centrality.nii.gz')

            # Set up workflow parameters
            wflow.base_dir = out_dir.replace('output', 'tests')
            wflow.inputs.inputspec.subject = func_mni
            wflow.inputs.inputspec.template = mask_template

            # Make the key the strategy being used (last folder)
            wflow_strat = out_dir.split('/')[-1]
            test_wflows[wflow_strat] = wflow

        # Set centrality TestCase instance attributes
        self.output_dirs = output_dirs
        self.test_wflows = test_wflows

    # Test the ants registration strategy
    def test_ants_strategy(self):
        '''
        Function to run the centrality workflows for the ANTS
        registration strategy
        '''

        # Import packages
        import os
        from CPAC.utils import tests_init

        # Init variables
        ants_wflow = self.test_wflows['ants']
        smooth_dict = {}

        # Set up workflows and run each
        ants_wflow.base_dir = os.path.join(ants_wflow.base_dir, 'deg')
        ants_wflow.inputs.inputspec.method_option = 0
        ants_wflow.inputs.inputspec.weight_options = [True, True]
        ants_wflow.inputs.inputspec.threshold_option = 0
        ants_wflow.inputs.inputspec.threshold = 0.001
        print 'running degree centrality...'
        ants_wflow.run()

        # Get raw outputs for smoothing and z-score calculations
        deg_niis = tests_init.return_all_niis(ants_wflow.base_dir)
        

        # Set up workflows and run each
        ants_wflow.base_dir = ants_wflow.base_dir.replace('deg', 'eig')
        ants_wflow.inputs.inputspec.method_option = 1
        ants_wflow.inputs.inputspec.weight_options = [True, True]
        ants_wflow.inputs.inputspec.threshold_option = 0
        ants_wflow.inputs.inputspec.threshold = 0.001
        print 'running eigenvector centrality...'
        ants_wflow.run()

        # Set up workflows and run each
        ants_wflow.base_dir = ants_wflow.base_dir.replace('eig', 'lfcd')
        ants_wflow.inputs.inputspec.method_option = 2
        ants_wflow.inputs.inputspec.weight_options = [True, False]
        ants_wflow.inputs.inputspec.threshold_option = 2
        ants_wflow.inputs.inputspec.threshold = 0.6
        print 'running lfcd...'
        ants_wflow.run()

    # Smooth nifti file
    def smooth_nii_output(self, nii_file):
        '''
        '''

        # Import packages
        import nibabel as nib
        import numpy as np
        import scipy as sp

        # Init variables
        mask_file = self.test_wflows['ants'].inputs.inputspec.template
        smooth_arr = np.zeros(mask_file.shape, dtype=float)

        # Grab niftis as numpy arrays
        raw_arr = nib.load(nii_file).get_data()
        mask_arr = nib.load(mask_file).get_data()

        # Smooth input
        smooth_out = sp.ndimage.gaussian_filter(raw_arr, sig, order=0)

        # Get mask coordinates
        coords = np.argwhere(mask_arr)
        idx = 0
        for xyz in coords:
            x, y, z = xyz
            smooth_arr[x, y, z] = smooth_out[idx]
            idx += 1


    # Collect test outputs and compare
    def test_collect_and_compare(self):
        '''
        Function to collect the precomputed and test outputs and
        compare the images
        '''

        # Import packages
        import glob
        import os
        import nibabel as nb
        import numpy as np

        # Init variables
        outputs_to_test = {}

        # Grab precomputed outputs and corresponding test outputs
        # For each (strategy) precomputed output directory
        for out_dir in self.output_dirs:
            niis = glob.glob(os.path.join(out_dir, '*.nii.gz'))
            test_dir = out_dir.replace('output', 'tests')
            # For each precomputed output nifti
            for nii in niis:
                nii_file = os.path.basename(nii)
                f_list = []
                for root, dirs, files in os.walk(test_dir):
                    if files:
                        f_list.extend([os.path.join(root, file) for file in files \
                                  if file == nii_file])
                if len(f_list) > 1:
                    err_msg = 'More than one file found for %s in %s; '\
                              'please use only one' % (nii_file, str(f_list))
                    raise Exception(err_msg)
                elif len(f_list) == 0:
                    print 'No test outputs found for %s, skipping comparison' \
                          % nii_file
                else:
                    strat = os.path.basename(out_dir)
                    if not outputs_to_test.has_key(strat):
                        outputs_to_test[strat] = {nii : f_list[0]}
                    else:
                        outputs_to_test[strat][nii] = f_list[0]

        # Iterate through dictionary and assert correlations\
        pass_thr = 0.98
        err_msg = 'Test failed: correlation < %.3f' % pass_thr
        for strat, golden_vs_test in outputs_to_test.items():
            for golden, test in golden_vs_test.items():
                # Load in golden and test images for comparison
                img1 = nb.load(golden).get_data()
                img2 = nb.load(test).get_data()

                # Compute pearson correlation on flattened 4D images
                print 'Comparing %s outputs...' % \
                      (os.path.basename(test.rstrip('.nii.gz')))
                corr = np.corrcoef(img1.flatten(), img2.flatten())[0,1]
                print 'Correlation = %.3f' % corr

                # Assert the correlation is >= pass_threshold
                self.assertGreaterEqual(corr, pass_thr, err_msg)

# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()
