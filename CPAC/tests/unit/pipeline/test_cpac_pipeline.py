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
class CPACPipelineRunTestCase(unittest.TestCase):
    '''
    This class is a test case for the cpac_pipeline.run function
    
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
        
        Returns
        -------
        None
            this function does not return any values, but populates the
            instance attributes for:
            self.config_file : string
            self.sublist_file : string
            self.idx : integer
            self.config : CPAC.utils.configuration.Configuration object
            self.strategies : list [dict]
        '''
        
        # Import packages
        import os
        import yaml
        from CPAC.utils.configuration import Configuration
        from CPAC.pipeline.cpac_runner import build_strategies
        
        # Init variables
        self.config_file = os.path.join(RESOURCE_DIR, 'test_pipeline_config.yml')
        self.sublist_file = os.path.join(RESOURCE_DIR, 'test_sublist.yml')
        self.idx = 0
        # Init Configuration class from config_file
        self.config = Configuration(yaml.load(open(self.config_file), 'r'))
        self.strategies = build_strategies(self.config)
    
    # Method to test AWS-S3 compatibility
    def test_s3_interaction(self):
        '''
        Method to test cpac_pipeline.run()'s ability to download data from
        '''
        
        # Import packages
        import os
        
        # AWS variables
        creds_path = os.path.join(RESOURCE_DIR, 'aws_creds.csv')
        bucket_name = 'fcp-indi'
        bucket_prefix = 'data/Projects/ABIDE_Initiative/'
        local_prefix = os.path.join(RESOURCE_DIR, 'inputs')
        
        # Run with AWS variables
        cpac_pipeline_run(self.config_file,
                          self.sublist_file,
                          self.idx,
                          self.strategies,
                          self.config.maskSpecificationFile,
                          self.config.roiSpecificationFile,
                          self.config.templateSpecificationFile,
                          p_name='test_s3',
                          creds_path=creds_path,
                          bucket_name=bucket_name,
                          bucket_prefix=bucket_prefix,
                          local_prefix=local_prefix)


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()