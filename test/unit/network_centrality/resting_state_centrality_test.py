# test/unit/network_centrality/test_resting_state_centrality.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in
CPAC/pipeline/cpac_pipeline.py
'''

# Import packages
import unittest
from test import RESOURCE_DIR

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
            self.wflow : nipype.pipeline.engine.Workflow object
            self.wflow.base_dir : string
        '''

        # Import packages
        import os
        from CPAC.network_centrality import resting_state_centrality

        # Init variables
        self.wflow = resting_state_centrality.create_resting_state_graphs(allocated_memory=1)
        self.wflow.base_dir = os.path.join(RESOURCE_DIR, 'network_centrality')


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()