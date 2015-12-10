# test/unit/AWS/test_fetch_creds.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs testing on the functions in
CPAC/AWS/fetch_creds.py
'''

# Import packages
import unittest
from CPAC.AWS import fetch_creds

# Test case for the run function
class FetchCredsTestCase(unittest.TestCase):
    '''
    This class is a test case for the fetch_creds.py module

    Inherits
    --------
    unittest.TestCase class

    Attributes (class):
    ------------------
    see unittest.TestCase documentation

    Attributes (instance):
    aws_creds : string
        filepath to the csv file on disk with the AWS credentials
    db_creds : string
        filepath to the csv file on disk with the database login
        credentials
    '''

    # setUp method - init the creds_path
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
            self.creds_path
            self.db_creds
        '''

        # Init variables
        self.creds_path = AWS_CREDS
        self.db_creds = DB_CREDS
        self.bucket_name = BUCKET_NAME

    # Test getting AWS keys from creds file
    def test_return_aws_keys(self):
        '''
        Method to test the fetch_creds.return_aws_keys() function

        Parameters
        ----------
        self : FetchCredsTestCase
            a unittest.TestCase-inherited class

        Returns
        -------
        None
            this function does not return any values, but tests to make
            sure the fetch_creds.return_aws_keys() function returns
            a string for each of the two keys
        '''

        # Grab the keys
        keys = fetch_creds.return_aws_keys(self.creds_path)

        # Assert that there are two keys
        self.assertEqual(len(keys), 2)

        # Assert that they are each strings
        self.assertIsInstance(keys[0], str)
        self.assertIsInstance(keys[1], str)

    # Test getting S3 bucket
    def test_return_bucket(self):
        '''
        Method to test the fetch_creds.return_bucket() function

        Parameters
        ----------
        self : FetchCredsTestCase
            a unittest.TestCase-inherited class

        Returns
        -------
        None
            this function does not return any values, but tests to make
            sure the fetch_creds.return_bucket() function returns a
            bucket object
        '''

        # Import packages
        import boto3

        # Init variables
        bucket_name = 'fcp-indi'
        err_msg = 'Unable to get the S3 bucket because of faulty AWS '\
                  'credentials or boto package not found'

        # Grab the AWS bucket
        bucket = fetch_creds.return_bucket(None,
                                           'fcp-indi')

        # Assert that it is a boto bucket object
        self.assertIs(bucket.name, 'fcp-indi')


# Command-line run-able unittest module
if __name__ == '__main__':
    unittest.main()
