# CPAC/tests/unit/pipeline/test_cpac_pipeline.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in
CPAC/pipeline/cpac_pipeline.py
'''

# Import packages
import unittest
from CPAC.pipeline.cpac_pipeline import run as cpac_pipeline_run
from CPAC.tests import RESOURCE_DIR

# Test case for the run function
class CentralityWorkflowTestCase(unittest.TestCase):
    '''
    This class is a test case for the resting_state_centrality workflow
    
    Inherits
    --------
    unittest.TestCase class
    '''
    
    # Set up workflow
    def setUp(self):
        '''
        Function which inits a nipype workflow using network_centrality's
        create_resting_state_graphs() function to be used for unit testing
    
        Parameters
        ----------
        meth_opt : integer
            0 - degree centrality, 1 - eigenvector centrality, 2 - lFCD
        thr_opt : integer
            0 - significance, 1 - sparsity, 2 - correlation
        thr : float
            the thresholding value at which to prune connections
        wght_opts : list (bool)
            a two-element list of booleans corresponding to the type of
            weighting centrality should perform when producing an output;
            e.g. [True, False] - perform binarized, not weighted centrality
        mask : string
            filepath to the location of the mask that will be used to
            isolate the portion of the input data, over which centrality
            will be analyzed
        data : string
            filepath to the location of the input functional data
        base_dir : string
            filepath to the base directory where nipype will work out of
            and save files to
        mem : float
            amount of allocated memory (in GB) that the centrality algorithm
            is permitted to use on the system
    
        Returns
        -------
           wflow : nipype.Workflow object
               the parameterized centrality nipype workflow to run
        '''
    
        # Import packages
        from CPAC.network_centrality import resting_state_centrality as graph
    
        # Init variables
        wflow = graph.create_resting_state_graphs(allocated_memory=mem)
    
        # Parameterize the workflow
        wflow.inputs.inputspec.method_option = meth_opt
        wflow.inputs.inputspec.threshold_option = thr_opt
        wflow.inputs.inputspec.threshold = thr
        wflow.inputs.inputspec.weight_options = wght_opts
        wflow.inputs.inputspec.template = mask
        wflow.inputs.inputspec.subject = data
        wflow.base_dir = base_dir
    
        # Return the workflow
        return wflow
    
    
    # Run the centrality workflow for all scenarios
    def main(data, mask, mem, base_dir):
        '''
        Function to calculate every kind of output for the network
        centrality workflow
    
        Parameters
        ----------
        data : string
            filepath to the location of the input functional data
        mask : string
            filepath to the location of the mask that will be used to
            isolate the portion of the input data, over which centrality
            will be analyzed
        mem : float
            amount of allocated memory (in GB) that the centrality algorithm
            is permitted to use on the system
        base_dir : string
            filepath to the base directory where nipype will work out of
            and save files to
    
        Returns
        -------
        None
            This function does not return any values; tt computes the 
            centrality outputs and writes them to disk
        '''
    
        # Import packages
    
        # Init variables
        meth_options = {0 : 'deg',
                        1 : 'eig',
                        2 : 'lfcd'}
        thr_options = {0 : 'sig',
                       1 : 'sp',
                       2 : 'r'}
        thr = 0.01
        wght_options = [True, True]
    
        # Iterate through the different methods
        for method, meth_dir in meth_options.items():
            # Run the workflow for each threshold type
            for thr_opt, thr_dir in thr_options.items():
                # Form the complete base directory
                wf_dir = '/'.join([base_dir, meth_dir, thr_dir])
                # Parameterize and return the workflow
                wflow = init_workflow(method,
                                      thr_opt,
                                      thr,
                                      wght_options,
                                      mask,
                                      data,
                                      wf_dir,
                                      mem)
                # And run the workflow
                wflow.run()


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()