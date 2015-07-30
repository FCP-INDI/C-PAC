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
        '''

        # Import packages
        import copy
        import os

        # Init variables
        ants_deg_wflow = self.test_wflows['ants']
        ants_eig_wflow = copy.deepcopy(ants_deg_wflow)
        ants_lfcd_wflow = copy.deepcopy(ants_deg_wflow)

        # Set up workflows and run each
        ants_deg_wflow.base_dir = os.path.join(ants_deg_wflow.base_dir, 'deg')
        ants_deg_wflow.inputs.inputspec.method_option = 0
        ants_deg_wflow.inputs.inputspec.weight_options = [True, True]
        ants_deg_wflow.inputs.inputspec.threshold_option = 0
        ants_deg_wflow.inputs.inputspec.threshold = 0.001
        print 'running degree centrality...'
        ants_deg_wflow.run()

        # Set up workflows and run each
        ants_eig_wflow.base_dir = os.path.join(ants_eig_wflow.base_dir, 'eig')
        ants_eig_wflow.inputs.inputspec.method_option = 1
        ants_eig_wflow.inputs.inputspec.weight_options = [True, True]
        ants_eig_wflow.inputs.inputspec.threshold_option = 0
        ants_eig_wflow.inputs.inputspec.threshold = 0.001
        print 'running eigenvector centrality...'
        ants_eig_wflow.run()

        # Set up workflows and run each
        ants_lfcd_wflow.base_dir = os.path.join(ants_lfcd_wflow.base_dir, 'lfcd')
        ants_lfcd_wflow.inputs.inputspec.method_option = 2
        ants_lfcd_wflow.inputs.inputspec.weight_options = [True, False]
        ants_lfcd_wflow.inputs.inputspec.threshold_option = 0
        ants_lfcd_wflow.inputs.inputspec.threshold = 0.6
        print 'running lfcd...'
        ants_lfcd_wflow.run()


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()
