# test/unit/nipype/datasource_test.py
#
# Author(s): Daniel Clark, 

'''
This module performs unit testing on the DataSink class from nipype
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
        Method to instantiate the TestCase

        Parameters
        ----------
        self : DatSinkTestCase
            a unittest.TestCase-inherited class

        Returns
        -------
        None
            this function does not return any values, but populates the
            instance attributes for:

            self.base_dir : string
                filepath to base directory of DataSink
            self.creds_path : string
                filepath to AWS keys credentials file
            self.data_sink : nipype DataSink object
                object from nipype.interfaces.io.DataSink()
            self.ds_node : nipype Node object
                object from nipype.pipeline.engine.Node()
            self.in_file : string
                filepath to nifti file to use as test input to datasink
        '''

        # Import packages
        import os
        import sys
        import tempfile

        import nipype.pipeline.engine as pe
        import nipype.interfaces.io as nio

        from CPAC.utils import test_init

        # Init variables
        try:
            input_dir = test_init.return_resource_subfolder('input')
            subj_id = test_init.return_test_subj()
            creds_path = test_init.return_aws_creds()
            bucket_name = test_init.default_bucket_name()
        except Exception as exc:
            print 'Unable to locate testing resources.\nError: %s' % exc
            sys.exit()

        # Datasink parameters
        base_dir = tempfile.mkdtemp()

        # Get the input file as an anatomical scan from th     e CPAC_RESOURCES
        in_file = os.path.join(input_dir, 'site_1', subj_id, 'session_1',
                               'anat_1', 'mprage.nii.gz')

        # Init the datasinks
        data_sink = nio.DataSink()
        ds_node = pe.Node(nio.DataSink(), name='sinker_0')

        # Add instance variables to TestCase
        self.base_dir = base_dir
        self.bucket_name = bucket_name
        self.creds_path = creds_path
        self.data_sink = data_sink
        self.ds_node = ds_node
        self.in_file = in_file

    # Collect the input file paths and form expected output files
    def _collect_inputs_outputs(self, data_sink):
        '''
        Method to collect input filepaths and return the expected
        output filepaths

        Parameters
        ----------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method
        data_sink : DataSink object
            data sink to test and to extract parameters from

        Returns
        -------
        in_out_dict : dictionary
            dictionary where the keys are tuples of (str path, str md5)
            and the values are corresponding string output paths
        '''

        # Import packages
        import os
        import hashlib

        # Init variables
        base_directory = data_sink.inputs.base_directory
        container = data_sink.inputs.container
        data_sink_src_d = data_sink.inputs._outputs
        in_out_dict = {}
        out_files = []
        s3_str = 's3://'

        # Explicitly lower-case the "s3"
        if base_directory.lower().startswith(s3_str):
            base_dir_sp = base_directory.split('/')
            base_dir_sp[0] = base_dir_sp[0].lower()
            base_directory = '/'.join(base_dir_sp)
            s3_strip = '/'.join(base_directory.split('/')[:3])

        # Get datasink out base
        out_prefix = os.path.join(base_directory, container)

        # If S3 in path, strip out s3_prefix (e.g. 's3://fcp-indi/')
        if base_directory.startswith(s3_str):
            out_prefix = out_prefix.replace(s3_strip, '')

        # Iterate over defined attribute keys and input sources
        for attr_folder, in_base in data_sink_src_d.items():
            in_files = []
            # Make the input base a list if it's not already
            if not type(in_base) is list:
                in_base = [in_base]
            # Iterate through each element of in_base
            for in_f in in_base:
                # If it's a directory, collect all files
                if os.path.isdir(in_f):
                    for root, dirs, files in os.walk(in_f):
                        in_files.extend([os.path.join(root, fil) for fil in files])
                # Otherwise, just append file
                else:
                    in_files.append(in_f)
                # Form temporary relative paths list
                rel_in_files = [i_f.replace(in_f, os.path.basename(in_f)) \
                                for i_f in in_files]
                # Extend the output files list
                out_files.extend([os.path.join(out_prefix, attr_folder, i_f) \
                                  for i_f in rel_in_files])

                # Finally add to dictionary of the MD5sums and out paths
                for idx, in_file in enumerate(in_files):
                    in_read = open(in_file, 'rb').read()
                    md5sum = hashlib.md5(in_read).hexdigest()
                    in_out_dict[(in_file, md5sum)] = out_files[idx]

        # Return the output paths
        return in_out_dict

    # Check to see if datasink output was produced
    def _check_output_exists(self, data_sink):
        '''
        Method to check for the existence of a datasink output

        Parameters
        ----------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method
        data_sink : DataSink object
            data sink to test and to extract parameters from

        Returns
        -------
        ds_out_exists : boolean
            flag indicating whether the expected output file exists
            from the datasink or not
        ds_out_msg : string
            message to forecast the status of the datasink output
        '''

        # Import packages
        import os
        import hashlib

        # Get in_files/md5sums and expected outputs
        in_out_dict = self._collect_inputs_outputs(data_sink)

        # Init outputs to positive
        ds_out_exists = True
        ds_out_msg = 'All files output successfully!'

        # Iterate through dictionary and set outputs to negative and break
        # if there is a problem
        for in_key, out_path in in_out_dict.items():
            in_path = in_key[0]
            in_md5sum = in_key[1]
            # Try and get the key info from S3
            try:
                in_read = open(out_path, 'rb').read()
            except IOError as exc:
                ds_out_msg = 'The specified file: %s does not exist.\nError: %s'\
                          % exc
                ds_out_exists = False
                break

            # Get md5sum from s3 key
            out_md5sum = hashlib.md5(in_read).hexdigest()
            if out_md5sum != in_md5sum:
                ds_out_msg = 'MD5 sums do not match for:\n%s and\n%s, there '\
                             'was an upload error. Try running datasink again'\
                             % (in_path, out_path)
                ds_out_exists = False
                break

        # Remove if it exists to free disk space
        if ds_out_exists:
            [os.remove(out_f) for out_f in in_out_dict.values()]

        # Return if it exists or not as boolean and output message
        return ds_out_exists, ds_out_msg

    # Check to see if datasink output was produced
    def _check_output_exists_s3(self, data_sink, bucket):
        '''
        Method to check for the existence of a datasink output

        Parameters
        ----------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        s3_out_exists : boolean
            flag indicating whether the expected output file exists
            from the datasink or not
        s3_out_msg : string
            message to forecast the status of the datasink output
        '''

        # Import packages
        from botocore.exceptions import ClientError

        # Init variables

        # Get in_files/md5sums and expected outputs
        in_out_dict = self._collect_inputs_outputs(data_sink)

        # Init outputs to positive
        s3_out_exists = True
        s3_out_msg = 'All files uploaded successfully!'

        # Iterate through dictionary and set outputs to negative and break
        # if there is a problem
        for in_key, out_path in in_out_dict.items():
            in_path = in_key[0]
            in_md5sum = in_key[1]
            if out_path.startswith('/'):
                out_path = out_path.lstrip('/')
            s3_key = bucket.Object(key=out_path)
            # Try and get the key info from S3
            try:
                s3_exists = s3_key.get()
            except ClientError as exc:
                s3_out_msg = 'The specified key: %s does not exist.\nError: %s'\
                          % (out_path, exc)
                s3_out_exists = False
                break

            # Get md5sum from s3 key
            out_md5sum = s3_key.e_tag.strip('"')
            if out_md5sum != in_md5sum:
                s3_out_msg = 'MD5 sums do not match for:\n%s and\n%s, there '\
                             'was an upload error. Try running datasink again'\
                             % (in_path, out_path)
                s3_out_exists = False
                break

        # Remove if it exists to free disk space
        if s3_out_exists:
            for out_path in in_out_dict.values():
                if out_path.startswith('/'):
                    out_path = out_path.lstrip('/')
                s3_obj = bucket.Object(key=out_path)
                s3_obj.delete()

        # Return outputs
        return s3_out_exists, s3_out_msg

    # Test datasink node will write to base directory
    def test_datasink_node(self):
        '''
        Method to test the datasink node for local output file

        Paramters
        ---------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        None
            this method does not return any value, but performs an
            assertion on whether the expected output exists or not
        '''

        # Import packages

        # Init variables
        attr_folder = 'input_scan'
        container = 'test_datasink'
        data_sink = self.data_sink

        # Set up datasink
        data_sink.inputs.base_directory = self.base_dir
        data_sink.inputs.container = container
        # Feed input to input_file
        setattr(data_sink.inputs, attr_folder, self.in_file)

        # Run data_sink
        data_sink.run()

        # Check if output exists
        ds_out_exists, out_path = \
            self._check_output_exists(data_sink)

        # Assert the output was produced as expected
        err_msg = 'Could not find output of datasink on disk: %s' % out_path
        self.assertTrue(ds_out_exists, msg=err_msg)

    # Test datasink node will write to base directory
    def test_datasink_node_folder(self):
        '''
        Method to test the datasink node for local output folder

        Paramters
        ---------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        None
            this method does not return any value, but performs an
            assertion on whether the expected output exists or not
        '''

        # Import packages
        import os

        # Init variables
        attr_folder = 'input_scan'
        container = 'test_datasink'
        data_sink = self.data_sink

        # Set up datasink
        data_sink.inputs.base_directory = self.base_dir
        data_sink.inputs.container = container
        # Feed input to input_file
        setattr(data_sink.inputs, attr_folder, os.path.dirname(self.in_file))

        # Run data_sink
        data_sink.run()

        # Check if output exists
        ds_out_exists, out_path = \
            self._check_output_exists(data_sink)

        # Assert the output was produced as expected
        err_msg = 'Could not find output of datasink on disk: %s' % out_path
        self.assertTrue(ds_out_exists, msg=err_msg)

    # Test resource pool node -> datsink node workflow will write output
    def test_datasink_res_pool_wf(self):
        '''
        Method to test the datasink running in a nipype workflow; this
        mimics the behavior implemented in the C-PAC pipeline and uses
        a resource pool node to feed the datasink node to produce the
        output file

        Paramters
        ---------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        None
            this method does not return any value, but performs an
            assertion on whether the expected output exists or not
        '''

        # Import packages
        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as niu

        # Init variables
        attr_folder = 'input_scan'
        container = 'test_datasink_res_pool_wf'
        data_sink = self.ds_node

        # Init workflow
        pipeline_wf = pe.Workflow(name='datasink_tester')
        pipeline_wf.base_dir = self.base_dir

        # Set up the dummy resource pool node
        res_pool_node = pe.Node(niu.IdentityInterface(fields=['out_file']),
                                name='res_pool')
        res_pool_node.inputs.out_file = self.in_file

        # Set up data sink node
        data_sink.inputs.base_directory = self.base_dir
        data_sink.inputs.container = container

        # Connect resource pool's out file to datasink node
        pipeline_wf.connect(res_pool_node, 'out_file', data_sink, attr_folder)

        # Run workflow
        pipeline_wf.run()

        # Check if output exists
        ds_out_exists, out_path = \
            self._check_output_exists(data_sink)

        # Assert the output was produced as expected
        err_msg = 'Could not find output of datasink on disk: %s' % out_path
        self.assertTrue(ds_out_exists, msg=err_msg)

    # Test datasink node writes to S3
    def test_s3_datasink_node(self):
        '''
        Method to test the datasink node for s3 output file

        Paramters
        ---------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        None
            this method does not return any value, but performs an
            assertion on whether the expected output exists or not on
            AWS S3
        '''

        # Import packages
        import os

        # Init variables
        attr_folder = 'input_scan'
        container = 'test_datasink'
        data_sink = self.data_sink

        # Init AWS variables
        s3_output_dir = os.path.join('S3://', self.bucket_name, 'data/unittest')

        # Set up datasink
        data_sink.inputs.base_directory = s3_output_dir
        data_sink.inputs.container = container
        data_sink.inputs.creds_path = self.creds_path

        # Feed input to input_file
        setattr(data_sink.inputs, attr_folder, self.in_file)

        # Run data_sink
        data_sink.run()

        # Grab datasink bucket
        bucket = data_sink._fetch_bucket(self.bucket_name)

        # Check if output exists
        s3_out_exists, s3_out_msg = self._check_output_exists_s3(data_sink, bucket)

        # Assert the output was produced as expected
        self.assertTrue(s3_out_exists, msg=s3_out_msg)

    # Test datasink node writes to S3
    def test_s3_datasink_node_folder(self):
        '''
        Method to test the datasink node for s3 output folder

        Paramters
        ---------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        None
            this method does not return any value, but performs an
            assertion on whether the expected output exists or not on
            AWS S3
        '''

        # Import packages
        import os

        # Init variables
        attr_folder = 'input_scan'
        container = 'test_datasink'
        data_sink = self.data_sink

        # Init AWS variables
        s3_output_dir = os.path.join('s3://', self.bucket_name, 'data/unittest')

        # Set up datasink
        data_sink.inputs.base_directory = s3_output_dir
        data_sink.inputs.container = container
        data_sink.inputs.creds_path = self.creds_path

        # Feed input to input_file
        setattr(data_sink.inputs, attr_folder, os.path.dirname(self.in_file))

        # Run data_sink
        data_sink.run()

        # Grab datasink bucket
        bucket = data_sink._fetch_bucket(self.bucket_name)

        # Check if output exists
        s3_out_exists, s3_out_msg = \
            self._check_output_exists_s3(data_sink, bucket)

        # Assert the output was produced as expected
        self.assertTrue(s3_out_exists, msg=s3_out_msg)

    # Test datasink node writes to S3
    def test_s3_datasink_rejects_anon(self):
        '''
        Method to test the datasink raises the proper error when trying
        to write to a non-public-writeable bucket with no credentials

        Paramters
        ---------
        self : unittest.TestCase
            instance method inherits self automatically; instance
            variables are used in the method

        Returns
        -------
        None
            this method does not return any value, but performs an
            assertion on whether the datasink will raise a ClientError
            if one tries to access an S3 bucket to write without
            credentials
        '''

        # Import packages
        import os
        from botocore.exceptions import ClientError

        # Init variables
        attr_folder = 'input_scan'
        container = 'test_datasink'
        data_sink = self.data_sink

        # Init AWS variables
        s3_output_dir = os.path.join('s3://', self.bucket_name, 'data/unittest')

        # Set up datasink
        data_sink.inputs.base_directory = s3_output_dir
        data_sink.inputs.container = container

        # Feed input to input_file
        setattr(data_sink.inputs, attr_folder, self.in_file)

        # Run data_sink
        self.assertRaises(ClientError, data_sink.run)


# Make module executable
if __name__ == '__main__':
    unittest.main()
