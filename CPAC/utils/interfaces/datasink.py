from builtins import object, zip, filter, range, open, str

import time
import glob
import fnmatch
import string
import json
import os
import os.path as op
import shutil
import re
import copy
import tempfile
from os.path import join, dirname
from shutil import SameFileError
from warnings import warn
from nipype.interfaces.io import IOBase, DataSinkInputSpec, DataSinkOutputSpec, ProgressPercentage, copytree

from nipype import config, logging
from nipype.utils.misc import human_order_sorted, str2bool

from nipype.utils.filemanip import (
    copyfile, simplify_list, ensure_list,
    get_related_files)

from nipype.interfaces.base import (CommandLineInputSpec, CommandLine, Directory, TraitedSpec,
                    traits, isdefined, File, InputMultiObject, InputMultiPath,
                    Undefined, Str)


iflogger = logging.getLogger('nipype.interface')


RETRY = 5
RETRY_WAIT = 5


def _get_head_bucket(s3_resource, bucket_name):
    """ Try to get the header info of a bucket, in order to
    check if it exists and its permissions
    """

    import botocore

    # Try fetch the bucket with the name argument
    err_msg = None
    for _ in range(RETRY):
        try:
            s3_resource.meta.client.head_bucket(Bucket=bucket_name)
            return

        except botocore.exceptions.ClientError as exc:
            error_code = int(exc.response['Error']['Code'])
            if error_code == 403:
                err_msg = 'Access to bucket: %s is denied; check credentials'\
                            % bucket_name
                break
            elif error_code == 404:
                err_msg = 'Bucket: %s does not exist; check spelling and try '\
                            'again' % bucket_name
                break
            else:
                err_msg = 'Unable to connect to bucket: %s. Error message:\n%s'\
                            % (bucket_name, exc)

        except Exception as exc:
            err_msg = 'Unable to connect to bucket: %s. Error message:\n%s'\
                        % (bucket_name, exc)

        time.sleep(RETRY_WAIT)

    if err_msg is not None:
        raise Exception(err_msg)


class DataSink(IOBase):
    """ Generic datasink module to store structured outputs
        Primarily for use within a workflow. This interface allows arbitrary
        creation of input attributes. The names of these attributes define the
        directory structure to create for storage of the files or directories.
        The attributes take the following form:
        string[[.[@]]string[[.[@]]string]] ...
        where parts between [] are optional.
        An attribute such as contrasts.@con will create a 'contrasts' directory
        to store the results linked to the attribute. If the @ is left out, such
        as in 'contrasts.con', a subdirectory 'con' will be created under
        'contrasts'.
        the general form of the output is::
           'base_directory/container/parameterization/destloc/filename'
           destloc = string[[.[@]]string[[.[@]]string]] and
           filename comesfrom the input to the connect statement.
        .. warning::
            This is not a thread-safe node because it can write to a common
            shared location. It will not complain when it overwrites a file.
        .. note::
            If both substitutions and regexp_substitutions are used, then
            substitutions are applied first followed by regexp_substitutions.
            This interface **cannot** be used in a MapNode as the inputs are
            defined only when the connect statement is executed.
        Examples
        --------
        >>> ds = DataSink()
        >>> ds.inputs.base_directory = 'results_dir'
        >>> ds.inputs.container = 'subject'
        >>> ds.inputs.structural = 'structural.nii'
        >>> setattr(ds.inputs, 'contrasts.@con', ['cont1.nii', 'cont2.nii'])
        >>> setattr(ds.inputs, 'contrasts.alt', ['cont1a.nii', 'cont2a.nii'])
        >>> ds.run()  # doctest: +SKIP
        To use DataSink in a MapNode, its inputs have to be defined at the
        time the interface is created.
        >>> ds = DataSink(infields=['contasts.@con'])
        >>> ds.inputs.base_directory = 'results_dir'
        >>> ds.inputs.container = 'subject'
        >>> ds.inputs.structural = 'structural.nii'
        >>> setattr(ds.inputs, 'contrasts.@con', ['cont1.nii', 'cont2.nii'])
        >>> setattr(ds.inputs, 'contrasts.alt', ['cont1a.nii', 'cont2a.nii'])
        >>> ds.run()  # doctest: +SKIP
    """

    # Give obj .inputs and .outputs
    input_spec = DataSinkInputSpec
    output_spec = DataSinkOutputSpec

    # Initialization method to set up datasink
    def __init__(self, infields=None, force_run=True, **kwargs):
        """
        Parameters
        ----------
        infields : list of str
            Indicates the input fields to be dynamically created
        """

        super(DataSink, self).__init__(**kwargs)
        undefined_traits = {}
        # used for mandatory inputs check
        self._infields = infields
        if infields:
            for key in infields:
                self.inputs.add_trait(key, traits.Any)
                self.inputs._outputs[key] = Undefined
                undefined_traits[key] = Undefined
        self.inputs.trait_set(trait_change_notify=False, **undefined_traits)
        if force_run:
            self._always_run = True

    # Get destination paths
    def _get_dst(self, src):
        # If path is directory with trailing os.path.sep,
        # then remove that for a more robust behavior
        src = src.rstrip(os.path.sep)
        path, fname = os.path.split(src)
        if self.inputs.parameterization:
            dst = path
            if isdefined(self.inputs.strip_dir):
                dst = dst.replace(self.inputs.strip_dir, '')
            folders = [
                folder for folder in dst.split(os.path.sep)
                if folder.startswith('_')
            ]
            dst = os.path.sep.join(folders)
            if fname:
                dst = os.path.join(dst, fname)
        else:
            if fname:
                dst = fname
            else:
                dst = path.split(os.path.sep)[-1]
        if dst[0] == os.path.sep:
            dst = dst[1:]
        return dst

    # Substitute paths in substitutions dictionary parameter
    def _substitute(self, pathstr):
        pathstr_ = pathstr
        if isdefined(self.inputs.substitutions):
            for key, val in self.inputs.substitutions:
                oldpathstr = pathstr
                pathstr = pathstr.replace(key, val)
                if pathstr != oldpathstr:
                    iflogger.debug('sub.str: %s -> %s using %r -> %r',
                                   oldpathstr, pathstr, key, val)
        if isdefined(self.inputs.regexp_substitutions):
            for key, val in self.inputs.regexp_substitutions:
                oldpathstr = pathstr
                pathstr, _ = re.subn(key, val, pathstr)
                if pathstr != oldpathstr:
                    iflogger.debug('sub.regexp: %s -> %s using %r -> %r',
                                   oldpathstr, pathstr, key, val)
        if pathstr_ != pathstr:
            iflogger.info('sub: %s -> %s', pathstr_, pathstr)
        return pathstr

    # Check for s3 in base directory
    def _check_s3_base_dir(self):
        '''
        Method to see if the datasink's base directory specifies an
        S3 bucket path; if it does, it parses the path for the bucket
        name in the form 's3://bucket_name/...' and returns it
        Parameters
        ----------
        Returns
        -------
        s3_flag : boolean
            flag indicating whether the base_directory contained an
            S3 bucket path
        bucket_name : string
            name of the S3 bucket to connect to; if the base directory
            is not a valid S3 path, defaults to '<N/A>'
        '''

        # Init variables
        s3_str = 's3://'
        bucket_name = '<N/A>'
        base_directory = self.inputs.base_directory

        if not isdefined(base_directory):
            s3_flag = False
            return s3_flag, bucket_name

        # Explicitly lower-case the "s3"
        if base_directory.lower().startswith(s3_str):
            base_dir_sp = base_directory.split('/')
            base_dir_sp[0] = base_dir_sp[0].lower()
            base_directory = '/'.join(base_dir_sp)

        # Check if 's3://' in base dir
        if base_directory.startswith(s3_str):
            # Expects bucket name to be 's3://bucket_name/base_dir/..'
            bucket_name = base_directory.split(s3_str)[1].split('/')[0]
            s3_flag = True
        # Otherwise it's just a normal datasink
        else:
            s3_flag = False

        # Return s3_flag
        return s3_flag, bucket_name

    # Function to return AWS secure environment variables
    def _return_aws_keys(self):
        '''
        Method to return AWS access key id and secret access key using
        credentials found in a local file.
        Parameters
        ----------
        self : nipype.interfaces.io.DataSink
            self for instance method
        Returns
        -------
        aws_access_key_id : string
            string of the AWS access key ID
        aws_secret_access_key : string
            string of the AWS secret access key
        '''

        # Import packages
        import os

        # Init variables
        creds_path = self.inputs.creds_path

        # Check if creds exist
        if creds_path and os.path.exists(creds_path):
            with open(creds_path, 'r') as creds_in:
                # Grab csv rows
                row1 = creds_in.readline()
                row2 = creds_in.readline()

            # Are they root or user keys
            if 'User Name' in row1:
                # And split out for keys
                aws_access_key_id = row2.split(',')[1]
                aws_secret_access_key = row2.split(',')[2]
            elif 'AWSAccessKeyId' in row1:
                # And split out for keys
                aws_access_key_id = row1.split('=')[1]
                aws_secret_access_key = row2.split('=')[1]
            else:
                err_msg = 'Credentials file not recognized, check file is correct'
                raise Exception(err_msg)

            # Strip any carriage return/line feeds
            aws_access_key_id = aws_access_key_id.replace('\r', '').replace(
                '\n', '')
            aws_secret_access_key = aws_secret_access_key.replace('\r',
                                                                  '').replace(
                                                                      '\n', '')
        else:
            aws_access_key_id = os.getenv('AWS_ACCESS_KEY_ID')
            aws_secret_access_key = os.getenv('AWS_SECRET_ACCESS_KEY')

        # Return keys
        return aws_access_key_id, aws_secret_access_key

    # Fetch bucket object
    def _fetch_bucket(self, bucket_name):
        '''
        Method to return a bucket object which can be used to interact
        with an AWS S3 bucket using credentials found in a local file.
        Parameters
        ----------
        self : nipype.interfaces.io.DataSink
            self for instance method
        bucket_name : string
            string corresponding to the name of the bucket on S3
        Returns
        -------
        bucket : boto3.resources.factory.s3.Bucket
            boto3 s3 Bucket object which is used to interact with files
            in an S3 bucket on AWS
        '''

        # Import packages
        try:
            import boto3
            import botocore
        except ImportError as exc:
            err_msg = 'Boto3 package is not installed - install boto3 and '\
                      'try again.'
            raise Exception(err_msg)

        # Init variables
        creds_path = self.inputs.creds_path

        # Get AWS credentials
        try:
            aws_access_key_id, aws_secret_access_key = \
                self._return_aws_keys()
        except Exception as exc:
            err_msg = 'There was a problem extracting the AWS credentials '\
                      'from the credentials file provided: %s. Error:\n%s'\
                      % (creds_path, exc)
            raise Exception(err_msg)

        # Try and get AWS credentials if a creds_path is specified
        if aws_access_key_id and aws_secret_access_key:
            # Init connection
            iflogger.info('Connecting to S3 bucket: %s with credentials...',
                          bucket_name)
            # Use individual session for each instance of DataSink
            # Better when datasinks are being used in multi-threading, see:
            # http://boto3.readthedocs.org/en/latest/guide/resources.html#multithreading
            session = boto3.session.Session(
                aws_access_key_id=aws_access_key_id,
                aws_secret_access_key=aws_secret_access_key)

        else:
            iflogger.info('Connecting to S3 bucket: %s with IAM role...',
                          bucket_name)

            # Lean on AWS environment / IAM role authentication and authorization
            session = boto3.session.Session()

        s3_resource = session.resource('s3', use_ssl=True)

        # And try fetch the bucket with the name argument
        try:
            _get_head_bucket(s3_resource, bucket_name)
        except Exception as exc:

            # Try to connect anonymously
            s3_resource.meta.client.meta.events.register(
                'choose-signer.s3.*', botocore.handlers.disable_signing)

            iflogger.info('Connecting to AWS: %s anonymously...', bucket_name)
            _get_head_bucket(s3_resource, bucket_name)

        # Explicitly declare a secure SSL connection for bucket object
        bucket = s3_resource.Bucket(bucket_name)

        # Return the bucket
        return bucket


    # Send up to S3 method
    def _upload_to_s3(self, bucket, src, dst):
        '''
        Method to upload outputs to S3 bucket instead of on local disk
        '''

        # Import packages
        import hashlib
        import os

        from botocore.exceptions import ClientError

        # Init variables
        s3_str = 's3://'
        s3_prefix = s3_str + bucket.name

        # Explicitly lower-case the "s3"
        if dst[:len(s3_str)].lower() == s3_str:
            dst = s3_str + dst[len(s3_str):]

        # If src is a directory, collect files (this assumes dst is a dir too)
        if os.path.isdir(src):
            src_files = []
            for root, dirs, files in os.walk(src):
                src_files.extend([os.path.join(root, fil) for fil in files])
            # Make the dst files have the dst folder as base dir
            dst_files = [
                os.path.join(dst,
                             src_f.split(src)[1]) for src_f in src_files
            ]
        else:
            src_files = [src]
            dst_files = [dst]

        # Iterate over src and copy to dst
        for src_idx, src_f in enumerate(src_files):
            # Get destination filename/keyname
            dst_f = dst_files[src_idx]
            dst_k = dst_f.replace(s3_prefix, '').lstrip('/')

            # See if same file is already up there
            try:
                dst_obj = bucket.Object(key=dst_k)
                dst_md5 = dst_obj.e_tag.strip('"')

                # See if same file is already there
                src_read = open(src_f, 'rb').read()
                src_md5 = hashlib.md5(src_read).hexdigest()
                # Move to next loop iteration
                if dst_md5 == src_md5:
                    iflogger.info('File %s already exists on S3, skipping...',
                                  dst_f)
                    continue
                else:
                    iflogger.info('Overwriting previous S3 file...')

            except ClientError:
                iflogger.info('New file to S3')

            # Copy file up to S3 (either encrypted or not)
            iflogger.info('Uploading %s to S3 bucket, %s, as %s...', src_f,
                          bucket.name, dst_f)
            if self.inputs.encrypt_bucket_keys:
                extra_args = {'ServerSideEncryption': 'AES256'}
            else:
                extra_args = {}

            retry_exc = None
            for _ in range(RETRY):
                try:
                    bucket.upload_file(
                        src_f,
                        dst_k,
                        ExtraArgs=extra_args,
                        Callback=ProgressPercentage(src_f)
                    )
                    break
                except Exception as exc:
                    time.sleep(RETRY_WAIT)
                    retry_exc = exc

            if retry_exc is not None:
                raise retry_exc

    # List outputs, main run routine
    def _list_outputs(self):
        """Execute this module.
        """

        # Init variables
        outputs = self.output_spec().get()
        out_files = []
        # Use hardlink
        use_hardlink = str2bool(
            config.get('execution', 'try_hard_link_datasink'))

        # Set local output directory if specified
        if isdefined(self.inputs.local_copy):
            outdir = self.inputs.local_copy
        else:
            outdir = self.inputs.base_directory
            # If base directory isn't given, assume current directory
            if not isdefined(outdir):
                outdir = '.'

        # Check if base directory reflects S3 bucket upload
        s3_flag, bucket_name = self._check_s3_base_dir()
        if s3_flag:
            s3dir = self.inputs.base_directory
            # If user overrides bucket object, use that
            if self.inputs.bucket:
                bucket = self.inputs.bucket
            # Otherwise fetch bucket object using name
            else:
                try:
                    bucket = self._fetch_bucket(bucket_name)
                # If encountering an exception during bucket access, set output
                # base directory to a local folder
                except Exception as exc:
                    s3dir = '<N/A>'
                    if not isdefined(self.inputs.local_copy):
                        local_out_exception = os.path.join(
                            os.path.expanduser('~'),
                            's3_datasink_' + bucket_name)
                        outdir = local_out_exception
                    # Log local copying directory
                    iflogger.info(
                        'Access to S3 failed! Storing outputs locally at: '
                        '%s\nError: %s', outdir, exc)
        else:
            s3dir = '<N/A>'

        # If container input is given, append that to outdir
        if isdefined(self.inputs.container):
            outdir = os.path.join(outdir, self.inputs.container)
            s3dir = os.path.join(s3dir, self.inputs.container)

        # If sinking to local folder
        if outdir != s3dir:
            outdir = os.path.abspath(outdir)
            # Create the directory if it doesn't exist
            if not os.path.exists(outdir):
                try:
                    os.makedirs(outdir)
                except OSError as inst:
                    if 'File exists' in inst.strerror:
                        pass
                    else:
                        raise (inst)

        # Iterate through outputs attributes {key : path(s)}
        for key, files in list(self.inputs._outputs.items()):
            if not isdefined(files):
                continue
            iflogger.debug("key: %s files: %s", key, str(files))
            files = ensure_list(files if files else [])
            tempoutdir = outdir
            if s3_flag:
                s3tempoutdir = s3dir
            for d in key.split('.'):
                if d[0] == '@':
                    continue
                tempoutdir = os.path.join(tempoutdir, d)
                if s3_flag:
                    s3tempoutdir = os.path.join(s3tempoutdir, d)

            # flattening list
            if files and isinstance(files, list):
                if isinstance(files[0], list):
                    files = [item for sublist in files for item in sublist]

            # Iterate through passed-in source files
            for src in ensure_list(files):
                # Format src and dst files
                src = os.path.abspath(src)
                if not os.path.isfile(src):
                    src = os.path.join(src, '')
                dst = self._get_dst(src)
                if s3_flag:
                    s3dst = os.path.join(s3tempoutdir, dst)
                    s3dst = self._substitute(s3dst)
                dst = os.path.join(tempoutdir, dst)
                dst = self._substitute(dst)
                path, _ = os.path.split(dst)

                # If we're uploading to S3
                if s3_flag:
                    self._upload_to_s3(bucket, src, s3dst)
                    out_files.append(s3dst)
                # Otherwise, copy locally src -> dst
                if not s3_flag or isdefined(self.inputs.local_copy):
                    # Create output directory if it doesn't exist
                    if not os.path.exists(path):
                        try:
                            os.makedirs(path)
                        except OSError as inst:
                            if 'File exists' in inst.strerror:
                                pass
                            else:
                                raise (inst)
                    try:
                        # If src == dst, it's already home
                        if (not os.path.exists(dst)) or (
                            os.stat(src) != os.stat(dst)
                        ):
                            # If src is a file, copy it to dst
                            if os.path.isfile(src):
                                iflogger.debug(f'copyfile: {src} {dst}')
                                copyfile(
                                    src,
                                    dst,
                                    copy=True,
                                    hashmethod='content',
                                    use_hardlink=use_hardlink)
                            # If src is a directory, copy
                            # entire contents to dst dir
                            elif os.path.isdir(src):
                                if (
                                    os.path.exists(dst) and
                                    self.inputs.remove_dest_dir
                                ):
                                    iflogger.debug('removing: %s', dst)
                                    shutil.rmtree(dst)
                                iflogger.debug('copydir: %s %s', src, dst)
                                copytree(src, dst)
                                out_files.append(dst)
                    except SameFileError:
                        iflogger.debug(f'copyfile (same file): {src} {dst}')

        # Return outputs dictionary
        outputs['out_file'] = out_files

        return outputs
