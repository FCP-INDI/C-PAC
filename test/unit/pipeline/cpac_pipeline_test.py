# test/unit/pipeline/test_cpac_pipeline.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in
CPAC/pipeline/cpac_pipeline.py
'''

# Import packages
import unittest

# Test case for the run function
class CPACPipelineRunTestCase(unittest.TestCase):
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
    config_file : string
        filepath to the C-PAC pipeline configuration file
    sublist_file : string
        filepath to the C-PAC subject list file
    idx : integer
        integer index for the specific subject from the subject list to
        run
    config : CPAC.utils.configuration.Configuration object
        the pipeline configuration C-PAC internal data structure
    strategies : 
        the strategies to run for the pipeline
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


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()