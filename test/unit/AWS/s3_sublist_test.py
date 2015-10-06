# test/unit/AWS/s3_sublist_test.py
#
# Author(s): Daniel Clark, 2015

'''
This module performs unit testing on functions in the s3_sublist
module in the CPAC/AWS subpackage
'''

# Import packages
import unittest

# Tets case for cpac datasink
class DataSinkTestCase(unittest.TestCase):
    '''
    This class is a test case for the s3_sublist module in
    CPAC/AWS

    Inherits
    --------
    unittest.TestCase class

    Attributes (class):
    ------------------
    see unittest.TestCase documentation

    Attributes (instance):
    '''

    # setUp method - init the input file
    def setUp(self):
        '''
        Method to instantiate input arguments for the
        AWS.fetch_creds() method via instance attributes

        Parameters
        ----------
        self : FetchCredsTestCase
            a unittest.TestCase-inherited class

        Returns
        -------
        None
            this function does not return any values, but populates the
            instance attributes for:

            self.creds_path : string
                filepath to AWS keys credentials file
        '''

