# test/unit/utils/datasource_test.py
#
# Contributing authors (please append):
# Daniel Clark

'''
This module performs unit testing on the functions and classes from
the CPAC/utils/datasource.py module
'''

# Import packages
import unittest

# Tets case for cpac datasink
class DataSinkTestCase(unittest.TestCase):
    '''
    This class is a test case for the DataSink class implemented in
    CPAC that inherits from nipype's DataSink class

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

            self.base_dir : string
                filepath to base directory of DataSink
            self.in_file : string
                filepath to nifti file to use as test input to datasink
            self.data_sink : nipype DataSink object
                object from nipype.interfaces.io.DataSink()
            self.cpac_ds : nipype Node object
                object from nipype.pipeline.engine.Node()
        '''

        # Import packages
        import os
        import nipype.pipeline.engine as pe
        from CPAC.utils import datasource, test_init

        # Init variables
        input_dir = test_init.return_resource_subfolder('input')
        subj_id = test_init.return_test_subj()

        # Datasink parameters
        base_dir = '/tmp/'

        # Get the input file as an anatomical scan from the CPAC_RESOURCES
        in_file = os.path.join(input_dir, 'site_1', subj_id, 'session_1',
                               'anat_1', 'mprage.nii.gz')

        # Init the datasinks
        data_sink = datasource.DataSink()
        ds_node = pe.Node(datasource.DataSink(), name='sinker_0')

        # Add instance attributes to TestCase
        self.base_dir = base_dir
        self.in_file = in_file
        self.data_sink = data_sink
        self.cpac_ds = cpac_ds

    # test datasink will write to base directory
    def run_datasink(self, data_sink):
        '''
        '''

        # Import packages
        import os

        # Init variables
        attr_folder = 'input_scan'
        container = 'subject'

        # Set up datasink
        data_sink.inputs.base_directory = self.base_dir
        data_sink.inputs.container = container
        # Feed input to input_file
        setattr(data_sink.inputs, attr_folder, self.in_file)

        # Run data_sink
        data_sink.run()

        # Check for output
        data_sink_outpath = os.path.join(self.base_dir, container,
                                         attr_folder,
                                         os.path.basename(self.in_file))

        # Check existence
        data_sink_out_exists = os.path.exists(data_sink_outpath)

        # Remove if it exists to keep clean
        if data_sink_out_exists:
            os.remove(data_sink_outpath)

        return data_sink_out_exists

    # test datasink node will write to base directory
    def test_datasinks(self):
        '''
        '''

        # Run the raw DataSink object
        ds_obj_out_exists = run_datasink(self, self.data_sink)
        # Assert the output was produced as expected
        self.assertTrue(ds_obj_out_exists)

        # Run the nipype node DataSink object
        ds_node_out_exists = run_datasink(self, self.data_sink_node)
        # Assert
        self.assertTrue(ds_node_out_exists)



# Make module executable
if __name__ == '__main__':
    unittest.main()

