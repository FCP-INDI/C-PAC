# CPAC/utils/datasource.py
#
#

'''
This module contains classes and functions used to interface with data
access
'''

# Import packages
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio

# Class to track percentage of S3 file upload
class ProgressPercentage(object):
    '''
    Call-able class instsance (via __call__ method) that displays
    upload percentage of a file to S3
    '''

    def __init__(self, filename):
        '''
        '''

        # Import packages
        import threading
        import os

        # Initialize data attributes
        self._filename = filename
        self._size = float(os.path.getsize(filename))
        self._seen_so_far = 0
        self._lock = threading.Lock()

    def __call__(self, bytes_amount):
        '''
        '''

        # Import packages
        import sys

        # With the lock on, print upload status
        with self._lock:
            self._seen_so_far += bytes_amount
            percentage = (self._seen_so_far / self._size) * 100
            progress_str = '%d / %d (%.2f%%)\r'\
                           % (self._seen_so_far, self._size, percentage)

            # Write to stdout
            sys.stdout.write(progress_str)
            sys.stdout.flush()


# Custom DataSinkInputSpec class
class DataSinkInputSpec(nio.DataSinkInputSpec):
    '''
    '''

    # Import packages
    import traits.api as tap
    import nipype.interfaces.traits_extension as nit
    import nipype.interfaces.base as nib

    # Init inputspec data attributes
    base_directory = nit.Directory(
        desc='Path to the base directory for storing data.')
    container = tap.Str(
        desc='Folder within base directory in which to store output')
    parameterization = tap.Bool(True, usedefault=True,
                                   desc='store output in parametrized structure')
    strip_dir = nit.Directory(desc='path to strip out of filename')
    substitutions = nib.InputMultiPath(tap.Tuple(tap.Str, tap.Str),
                                   desc=('List of 2-tuples reflecting string '
                                         'to substitute and string to replace '
                                         'it with'))
    regexp_substitutions = \
        nib.InputMultiPath(tap.Tuple(tap.Str, tap.Str),
                           desc=('List of 2-tuples reflecting a pair of a '\
                                 'Python regexp pattern and a replacement '\
                                 'string. Invoked after string `substitutions`'))

    _outputs = tap.Dict(tap.Str, value={}, usedefault=True)
    remove_dest_dir = tap.Bool(False, usedefault=True,
                                  desc='remove dest directory when copying dirs')

    # AWS S3 data attributes
    creds_path = tap.Str(desc='Filepath to AWS credentials file for S3 bucket '\
                              'access')
    encrypt_bucket_keys = tap.Bool(desc='Flag indicating whether to use S3 '\
                                        'server-side AES-256 encryption')

    # Set call-able inputs attributes
    def __setattr__(self, key, value):
        import nipype.interfaces.traits_extension as nit

        if key not in self.copyable_trait_names():
            if not nit.isdefined(value):
                super(DataSinkInputSpec, self).__setattr__(key, value)
            self._outputs[key] = value
        else:
            if key in self._outputs:
                self._outputs[key] = value
            super(DataSinkInputSpec, self).__setattr__(key, value)


# Custom DataSink class
class DataSink(nio.DataSink):
    '''
    '''

    # Give obj .inputs and .outputs
    input_spec = DataSinkInputSpec
    output_spec = nio.DataSinkOutputSpec

    # Check for s3 in base directory
    def _check_s3_base_dir(self):
        '''
        '''

        # Import packages
        import os
        import sys

        # Init variables
        s3_str = 's3://'
        sep = os.path.sep
        base_directory = self.inputs.base_directory

        # Check if 's3://' in base dir
        if base_directory.startswith(s3_str):
            try:
                # Expects bucket name to be 's3://bucket_name/base_dir/..'
                bucket_name = base_directory.split(s3_str)[1].split(sep)[0]
                # Get the actual bucket object
                self.bucket = self._fetch_bucket(bucket_name)
            except Exception as exc:
                err_msg = 'Unable to access S3 bucket. Error:\n%s. Exiting...'\
                          % exc
                print err_msg
                sys.exit()
            # Bucket access was a success, set flag
            s3_flag = True
        # Otherwise it's just a normal datasink
        else:
            s3_flag = False

        # Return s3_flag
        return s3_flag

    # Function to return AWS secure environment variables
    def _return_aws_keys(self, creds_path):
        '''
        Method to return AWS access key id and secret access key using
        credentials found in a local file.

        Parameters
        ----------
        creds_path : string (filepath)
            path to the csv file with 'AWSAccessKeyId=' followed by access
            key in the first row and 'AWSSecretAccessKey=' followed by
            secret access key in the second row

        Returns
        -------
        aws_access_key_id : string
            string of the AWS access key ID
        aws_secret_access_key : string
            string of the AWS secret access key
        '''

        # Import packages
        import csv

        # Init variables
        csv_reader = csv.reader(open(creds_path, 'r'))

        # Grab csv rows
        row1 = csv_reader.next()[0]
        row2 = csv_reader.next()[0]

        # And split out for keys
        aws_access_key_id = row1.split('=')[1]
        aws_secret_access_key = row2.split('=')[1]

        # Return keys
        return aws_access_key_id, aws_secret_access_key

    # Fetch bucket object
    def _fetch_bucket(self, bucket_name):
        '''
        Method to a return a bucket object which can be used to interact
        with an AWS S3 bucket using credentials found in a local file.

        Parameters
        ----------
        bucket_name : string
            string corresponding to the name of the bucket on S3

        Returns
        -------
        bucket : boto3.resources.factory.s3.Bucket
            boto3 s3 Bucket object which is used to interact with files
            in an S3 bucket on AWS
        '''

        # Import packages
        import logging

        try:
            import boto3
            import botocore
        except ImportError as exc:
            err_msg = 'Boto3 package is not installed - install boto3 and '\
                      'try again.'
            raise Exception(err_msg)

        # Init variables
        creds_path = self.inputs.creds_path
        iflogger = logging.getLogger('interface')

        # Try and get AWS credentials if a creds_path is specified
        if creds_path:
            try:
                aws_access_key_id, aws_secret_access_key = \
                    self._return_aws_keys(creds_path)
            except Exception as exc:
                err_msg = 'There was a problem extracting the AWS credentials '\
                          'from the credentials file provided: %s. Error:\n%s'\
                          % (creds_path, exc)
                raise Exception(err_msg)
            # Init connection
            iflogger.info('Connecting to S3 bucket: %s with credentials from '\
                          '%s ...' % (bucket_name, creds_path))
            # Use individual session for each instance of DataSink
            # Better when datasinks are being used in multi-threading, see:
            # http://boto3.readthedocs.org/en/latest/guide/resources.html#multithreading
            session = boto3.session.Session(aws_access_key_id=aws_access_key_id,
                                            aws_secret_access_key=aws_secret_access_key)
        # Otherwise, connect anonymously
        else:
            iflogger.info('Connecting to AWS: %s anonymously...'\
                          % bucket_name)
            session = boto3.session.Session()

        # Explicitly declare a secure SSL connection for bucket object
        s3_resource = session.resource('s3', use_ssl=True)
        bucket = s3_resource.Bucket(bucket_name)

        # And try fetch the bucket with the name argument
        try:
            s3_resource.meta.client.head_bucket(Bucket=bucket_name)
        except botocore.exceptions.ClientError as exc:
            error_code = int(exc.response['Error']['Code'])
            if error_code == 403:
                err_msg = 'Access to bucket: %s is denied; check credentials'\
                          % bucket_name
                raise Exception(err_msg)
            elif error_code == 404:
                err_msg = 'Bucket: %s does not exist; check spelling and try '\
                          'again' % bucket_name
                raise Exception(err_msg)
            else:
                err_msg = 'Unable to connect to bucket: %s. Error message:\n%s'\
                          % (bucket_name, exc)
        except Exception as exc:
            err_msg = 'Unable to connect to bucket: %s. Error message:\n%s'\
                      % (bucket_name, exc)
            raise Exception(err_msg)

        # Return the bucket
        return bucket


    # Send up to S3 method
    def _upload_to_s3(self, src, dst):
        '''
        Method to upload outputs to S3 bucket instead of on local disk
        '''

        # Import packages
        import logging
        import os

        # Init variables
        bucket = self.bucket
        iflogger = logging.getLogger('interface')
        s3_str = 's3://'
        s3_prefix = os.path.join(s3_str, bucket.name)

        # If src is a directory, collect files (this assumes dst is a dir too)
        if os.path.isdir(src):
            src_files = []
            for root, dirs, files in os.walk(src):
                src_files.extend([os.path.join(root, fil) for fil in files])
            # Make the dst files have the dst folder as base dir
            dst_files = [os.path.join(dst, src_f.split(src)[1]) \
                         for src_f in src_files]
        else:
            src_files = [src]
            dst_files = [dst]

        # Iterate over src and copy to dst
        for src_idx, src_f in enumerate(src_files):
            # Get destination filename/keyname
            dst_f = dst_files[src_idx]
            dst_k = dst_f.replace(s3_prefix, '').lstrip('/')

            # Copy file up to S3 (either encrypted or not)
            iflogger.info('Copying %s to S3 bucket, %s, as %s...'\
                          % (src_f, bucket.name, dst_f))
            if self.inputs.encrypt_bucket_keys:
                extra_args = {'ServerSideEncryption' : 'AES256'}
                bucket.upload_file(src_f, dst_k, ExtraArgs=extra_args,
                                   Callback=ProgressPercentage(src_f))
            else:
                bucket.upload_file(src_f, dst_k,
                                   Callback=ProgressPercentage(src_f))

    # List outputs, main run routine
    def _list_outputs(self):
        """Execute this module.
        """

        # Import packages
        import logging
        import os
        import shutil

        import nipype.utils.filemanip as nuf
        import nipype.utils.misc as num
        from nipype import config

        # Init variables
        iflogger = logging.getLogger('interface')
        outputs = self.output_spec().get()
        out_files = []
        outdir = self.inputs.base_directory
        use_hardlink = num.str2bool(config.get('execution',
                                               'try_hard_link_datasink'))

        # If base directory isn't given, assume current directory
        if not nuf.isdefined(outdir):
            outdir = '.'

        # Check if base directory reflects S3-bucket upload
        s3_flag = self._check_s3_base_dir()
        if not s3_flag:
            outdir = os.path.abspath(outdir)

        # If container input is given, append that to outdir
        if nuf.isdefined(self.inputs.container):
            outdir = os.path.join(outdir, self.inputs.container)
        # Create the directory if it doesn't exist
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except OSError, inst:
                if 'File exists' in inst:
                    pass
                else:
                    raise(inst)

        # Iterate through outputs attributes {key : path(s)}
        for key, files in self.inputs._outputs.items():
            if not nuf.isdefined(files):
                continue
            iflogger.debug("key: %s files: %s" % (key, str(files)))
            files = nuf.filename_to_list(files)
            tempoutdir = outdir
            for d in key.split('.'):
                if d[0] == '@':
                    continue
                tempoutdir = os.path.join(tempoutdir, d)

            # flattening list
            if isinstance(files, list):
                if isinstance(files[0], list):
                    files = [item for sublist in files for item in sublist]

            # Iterate through passed-in source files
            for src in nuf.filename_to_list(files):
                # Format src and dst files
                src = os.path.abspath(src)
                if not os.path.isfile(src):
                    src = os.path.join(src, '')
                dst = self._get_dst(src)
                dst = os.path.join(tempoutdir, dst)
                dst = self._substitute(dst)
                path, _ = os.path.split(dst)

                # Create output directory if it doesnt exist
                if not os.path.exists(path):
                    try:
                        os.makedirs(path)
                    except OSError, inst:
                        if 'File exists' in inst:
                            pass
                        else:
                            raise(inst)

                # If we're uploading to S3
                if s3_flag:
                    self._upload_to_s3(src, dst)
                    out_files.append(dst)
                # Otherwise, copy locally src -> dst
                else:
                    # If src is a file, copy it to dst
                    if os.path.isfile(src):
                        iflogger.debug('copyfile: %s %s' % (src, dst))
                        nuf.copyfile(src, dst, copy=True, hashmethod='content',
                                     use_hardlink=use_hardlink)
                        out_files.append(dst)
                    # If src is a directory, copy entire contents to dst dir
                    elif os.path.isdir(src):
                        if os.path.exists(dst) and self.inputs.remove_dest_dir:
                            iflogger.debug('removing: %s' % dst)
                            shutil.rmtree(dst)
                        iflogger.debug('copydir: %s %s' % (src, dst))
                        nio.copytree(src, dst)
                        out_files.append(dst)

        # Return outputs dictionary
        outputs['out_file'] = out_files

        return outputs


def create_func_datasource(rest_dict, wf_name='func_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)


    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'scan'],
                                mandatory_inputs=True),
                        name='inputnode')
    inputnode.iterables = [('scan', rest_dict.keys())]

    selectrest = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['rest'],
                        function=get_rest),
                         name='selectrest')
    selectrest.inputs.rest_dict = rest_dict

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'rest',
                                                     'scan' ]),
                         name='outputspec')

    wf.connect(inputnode, 'scan', selectrest, 'scan')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(selectrest, 'rest', outputnode, 'rest')
    wf.connect(inputnode, 'scan', outputnode, 'scan')

    return wf


def get_rest(scan, rest_dict):
    return rest_dict[scan]



def create_anat_datasource(wf_name='anat_datasource'):

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    wf = pe.Workflow(name=wf_name)

    inputnode = pe.Node(util.IdentityInterface(
                                fields=['subject', 'anat'],
                                mandatory_inputs=True),
                        name='inputnode')

    outputnode = pe.Node(util.IdentityInterface(fields=['subject',
                                                     'anat' ]),
                         name='outputspec')

    wf.connect(inputnode, 'subject', outputnode, 'subject')
    wf.connect(inputnode, 'anat', outputnode, 'anat')

    return wf



def create_roi_mask_dataflow(dir_path, mask_type, wf_name='datasource_roi_mask'):

    import nipype.interfaces.io as nio
    import os

    wf = pe.Workflow(name=wf_name)

    if mask_type == 'roi':
        tab = 'ROI Average TSE'
    elif mask_type == 'voxel':
        tab = 'ROI Voxelwise TSE'
    elif mask_type == 'centrality':
        tab = 'Network Centrality'


    if '.nii' in dir_path:

        masks = []
        masks.append(dir_path)

    elif '.txt' in dir_path:
        
        masks = open(dir_path, 'r').readlines()

    else:

        print '\n\n[!] CPAC says: Your ROI/mask specification file (under ' \
              '%s options) either needs to be a NIFTI file (.nii or ' \
              '.nii.gz) of an ROI/mask or a text file (.txt) containing a ' \
              'list of NIFTI files of ROI/mask files.\nPlease change this ' \
              'in your pipeline configuration file and try again.\n\n' % tab
        raise Exception


    mask_dict = {}

    for mask_file in masks:

        mask_file = mask_file.rstrip('\r\n')

        if not os.path.exists(mask_file):
            err = '\n\n[!] CPAC says: One of your ROI/mask specification ' \
                  'files (under %s options) does not have a correct path ' \
                  'or does not exist.\nTip: If all the paths are okay, ' \
                  'then ensure there are no whitespaces or blank lines in ' \
                  'your ROI specification file.\n\n' % mask_type
            raise Exception(err)

        if mask_file.strip() == '' or mask_file.startswith('#'):
            continue

        base_file = os.path.basename(mask_file)
        base_name = ''
        if base_file.endswith('.nii'):
            base_name = os.path.splitext(base_file)[0]
        elif(base_file.endswith('.nii.gz')):
            base_name = os.path.splitext(os.path.splitext(base_file)[0])[0]
        else:
            err = "\n\n[!] CPAC says: One of your ROI/mask specification " \
                  "files (under %s options) does not have '.nii' or " \
                  "'.nii.gz' as an extension.\n\nMask file: %s\n\n" \
                  % (tab, mask_file)
            raise Exception(err)

        if not (base_name in mask_dict):
            mask_dict[base_name] = mask_file
        else:
            err = "\n\n[!] CPAC says: You have two or more ROI/mask files " \
            "with the same name - please make sure these files are named " \
            "differently.\n\nDuplicate name: %s\n\nNote: This can be " \
            "changed in the ROI/mask file you specified under the %s " \
            "options.\n\n" % (mask_file, tab)
            raise Exception(err)


    inputnode = pe.Node(util.IdentityInterface(
                            fields=['mask'],
                            mandatory_inputs=True),
                    name='inputspec')

    inputnode.iterables = [('mask', mask_dict.keys())]

    selectmask = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['out_file'],
                                       function=get_rest),
                         name='select_mask')
    selectmask.inputs.rest_dict = mask_dict

    outputnode = pe.Node(util.IdentityInterface(fields=['out_file']),
                         name='outputspec')

    wf.connect(inputnode, 'mask',
               selectmask, 'scan')

    wf.connect(selectmask, 'out_file',
               outputnode, 'out_file')
    return wf




def create_spatial_map_dataflow(dirPath, wf_name='datasource_maps'):

    import nipype.interfaces.io as nio
    import os


    wf = pe.Workflow(name=wf_name)
    spatial_maps = open(dirPath, 'r').readlines()

    spatial_map_dict = {}
    
    for spatial_map_file in spatial_maps:

        spatial_map_file = spatial_map_file.rstrip('\r\n')

        if not os.path.exists(spatial_map_file):
            print "\n\n" + "ERROR: One of your spatial map files (under Spatial" + \
            " Regression options) does not have a correct path or does not exist." + \
            "\n" + "Tip: If all the paths are okay, then ensure there are no" + \
            " whitespaces or blank lines in your spatial map specification file." + \
            "\n\n" + "Error name: datasource_0001" + "\n\n"
            raise Exception

        base_file = os.path.basename(spatial_map_file)
        base_name = ''
        try:
            if base_file.endswith('.nii'):
                base_name = os.path.splitext(base_file)[0]
            elif(base_file.endswith('.nii.gz')):
                base_name = os.path.splitext(os.path.splitext(base_file)[0])[0]
            else:
                raise Exception("File extension not in  .nii and .nii.gz File: %s" % spatial_map_file)
        except Exception, e:
            print('error in spatial_map_dataflow: ', e)

        if not (base_name in spatial_map_dict):
            spatial_map_dict[base_name] = spatial_map_file
        else:
            raise ValueError('Files with same name not allowed %s %s' % (spatial_map_file, spatial_map_dict[base_name]))

    inputnode = pe.Node(util.IdentityInterface(
                            fields=['spatial_map'],
                            mandatory_inputs=True),
                    name='inputspec')

    inputnode.iterables = [('spatial_map', spatial_map_dict.keys())]

    select_spatial_map = pe.Node(util.Function(input_names=['scan', 'rest_dict'],
                                       output_names=['out_file'],
                                       function=get_rest),
                         name='select_spatial_map')
    select_spatial_map.inputs.rest_dict = spatial_map_dict

    wf.connect(inputnode, 'spatial_map',
               select_spatial_map, 'scan')
    return wf



def create_grp_analysis_dataflow(wf_name='gp_dataflow'):

        import nipype.pipeline.engine as pe
        import nipype.interfaces.utility as util
        from CPAC.utils import select_model_files

        wf = pe.Workflow(name=wf_name)

        inputnode = pe.Node(util.IdentityInterface(fields=['ftest',
                                                           'grp_model',
                                                           'model_name'],
                                                   mandatory_inputs=True),
                            name='inputspec')

        selectmodel = pe.Node(util.Function(input_names=['model',
                                                         'ftest',
                                                         'model_name'],
                                           output_names=['fts_file',
                                                         'con_file',
                                                         'grp_file',
                                                         'mat_file'],
                                           function=select_model_files),
                             name='selectnode')

        wf.connect(inputnode, 'ftest',
                   selectmodel, 'ftest')
        wf.connect(inputnode, 'grp_model',
                   selectmodel, 'model')
        wf.connect(inputnode, 'model_name', selectmodel, 'model_name')



        outputnode = pe.Node(util.IdentityInterface(fields=['fts',
                                                            'grp',
                                                            'mat',
                                                            'con'],
                                mandatory_inputs=True),
                    name='outputspec')

        wf.connect(selectmodel, 'mat_file',
                   outputnode, 'mat')
        wf.connect(selectmodel, 'grp_file',
                   outputnode, 'grp') 
        wf.connect(selectmodel, 'fts_file',
                   outputnode, 'fts')
        wf.connect(selectmodel, 'con_file',
                   outputnode, 'con')


        return wf
