# CPAC/AWS/aws_utils.py
#
# Contributing authors:
# Daniel Clark

'''
This module contains functions which assist in interacting with AWS
services, including uploading/downloading data and file checking.
'''

# Class to track percentage of S3 file upload
class ProgressPercentage(object):
    '''
    Callable class instsance (via __call__ method) that displays
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
        if not hasattr(filename, 'content_length'):
            self._size = float(os.path.getsize(filename))
        else:
            self._size = float(filename.content_length)
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
            if self._size != 0:
                percentage = (self._seen_so_far / self._size) * 100
            else:
                percentage = 0
            progress_str = '%d / %d (%.2f%%)\r'\
                           % (self._seen_so_far, self._size, percentage)

            # Write to stdout
            sys.stdout.write(progress_str)
            sys.stdout.flush()


# Get the MD5 sums of files on S3
def md5_sum(bucket, prefix='', filt_str=''):
    '''
    Function to get the filenames and MD5 checksums of files stored in
    an S3 bucket and return this as a dictionary.

    Parameters
    ----------
    bucket : boto3 Bucket instance
        an instance of the boto3 S3 bucket class to download from
    prefix : string (optional), default=''
        the bucket prefix where all of the file keys are located
    filt_str : string (optional), defualt=''
        a string to filter the filekeys of interest;
        e.g. 'matrix_data' will only return filekeys with the string
        'matrix_data' in their filepath name

    Returns
    -------
    md5_dict : dictionary {str : str}
        a dictionary where the keys are the S3 filename and the values
        are the MD5 checksum values
    '''

    # Init variables
    blist = bucket.objects.filter(Prefix=prefix)
    md5_dict = {}

    # And iterate over keys to copy over new ones
    for bkey in blist:
        filename = str(bkey.key)
        if filt_str in filename:
            md5_sum = str(bkey.etag).strip('"')
            md5_dict[filename] = md5_sum
            print 'filename: %s' % filename
            print 'md5_sum: %s' % md5_sum

    # Return the dictionary
    return md5_dict


# Rename s3 keys from src_list to dst_list
def s3_rename(bucket, src_list, dst_list,
              keep_old=False, make_public=False):
    '''
    Function to rename files from an AWS S3 bucket via a copy and delete
    process. Uses all keys in src_list as the original names and renames
    the them to the corresponding keys in the dst_list.
    (e.g. src_list[9] --> dst_list[9])

    Parameters
    ----------
    bucket : boto3 Bucket instance
        an instance of the boto3 S3 bucket class to download from
    src_list : list (str)
        a list of relative paths of the files to delete from the bucket
    dst_list : list (str)
        a list of relative paths of the files to delete from the bucket
    keep_old : boolean (optional), default=False
        flag indicating whether to keep the src_list files
    make_public : boolean (optional), default=False
        set to True if files should be publically available on S3

    Returns
    -------
    None
        The function doesn't return any value, it deletes data from
        S3 and prints its progress and a 'done' message upon completion
    '''

    # Import packages
    from botocore.exceptions import ClientError

    # Check list lengths are equal
    if len(src_list) != len(dst_list):
        raise ValueError('src_list and dst_list are different lengths!')

    # Init variables
    num_files = len(src_list)

    # Check if the list lengths match 
    if num_files != len(dst_list):
        raise RuntimeError, 'src_list and dst_list must be the same length!'

    # And iterate over keys to copy over new ones
    for idx, src_f in enumerate(src_list):
        src_key = bucket.Object(key=src_f)
        try:
            src_key_exists = src_key.get()
        except ClientError as exc:
            print 'source file %s doesnt exist, skipping...' % src_f
            continue

        # Get corresponding destination file
        dst_key = dst_list[idx]
        dst_obj = bucket.Object(key=dst_key)
        try:
            dst_key_exists = dst_obj.get()
            print 'Destination key %s exists, skipping...' % dst_key
            continue
        except ClientError as exc:
            print 'copying source: ', str(src_f)
            print 'to destination: ', dst_key
            if make_public:
                print 'making public...'
                dst_obj.copy_from(CopySource=bucket.name + '/' + str(src_f), ACL='public-read')
            else:
                dst_obj.copy_from(CopySource=bucket.name + '/' + str(src_f))
            if not keep_old:
                src_key.delete()

        # Print status
        per = 100*(float(idx+1)/num_files)
        print 'Done renaming %d/%d\n%f%% complete' % (idx+1, num_files, per)

    # Done iterating through list
    print 'done!'


# Delete s3 keys based on input list
def s3_delete(bucket, bucket_keys):
    '''
    Method to delete files from an AWS S3 bucket that have the same
    names as those of an input list to a local directory.

    Parameters
    ----------
    bucket : boto3 Bucket instance
        an instance of the boto3 S3 bucket class to delete from
    bucket_keys : list
        a list of relative paths of the files to delete from the bucket

    Returns
    -------
    None
        The function doesn't return any value, it deletes data from
        S3 and prints its progress and a 'done' message upon completion
    '''

    # Init variables
    num_files = len(bucket_keys)

    # Iterate over list and delete S3 items
    for idx, bkey in enumerate(bucket_keys):
        try:
            print 'attempting to delete %s from %s...' % (bkey, bucket.name)
            bobj = bucket.Object(bkey)
            bobj.delete()
            per = 100*(float(idx+1)/num_files)
            print 'Done deleting %d/%d\n%f%% complete' % (idx+1, num_files, per)
        except Exception as exc:
            print 'Unable to delete bucket key %s because of: %s' % (bkey, exc)

    # Done iterating through list
    print 'done!'


# Download files from AWS S3 to local machine
def s3_download(bucket, bucket_keys, download_dir):
    '''
    Method to download files from an AWS S3 bucket that have the same
    names as those of an input list to a local directory.

    Parameters
    ----------
    bucket : boto3 Bucket instance
        an instance of the boto3 S3 bucket class to download from
    bucket_keys : list
        a list of relative paths of the files to download from the bucket
    downoad_dir : string
        local directory prefix to store the downloaded data

    Returns
    -------
    None
        The function doesn't return any value, it downloads data from
        S3 and prints its progress and a 'done' message upon completion
    '''

    # Impor packages
    import hashlib
    import os

    from botocore.exceptions import ClientError

    # Init variables
    num_files = len(bucket_keys)

    # Get filepaths from S3 with prefix
    for idx, bkey in enumerate(bucket_keys):
        # Create a new key from the bucket and set its contents
        bobj = bucket.Object(key=bkey)

        # See if need to upload
        try:
            # If it exists, compare md5sums
            bkey_exists = bobj.get()
        except ClientError as exc:
            '%s does not exist in bucket on S3! Skipping...' % bkey
            continue
        s3_md5 = bobj.e_tag.strip('"')
        # Get local path
        local_path = os.path.join(download_dir, bkey)
        # Create subdirs if necessary
        dirname = os.path.dirname(local_path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        # If it exists, check its md5 before skipping
        if os.path.exists(local_path):
            if os.path.isdir(local_path):
                continue
            in_read = open(local_path, 'rb').read()
            local_md5 = hashlib.md5(in_read).hexdigest()
            if local_md5 == s3_md5:
                print 'Skipping %s, already downloaded...' % bkey
            else:
                try:
                    print 'Overwriting %s...' % local_path
                    bucket.download_file(bkey, local_path,
                                         Callback=ProgressPercentage(bobj))
                except Exception as exc:
                    print 'Could not download file %s because of: %s, skipping..' \
                          % (bkey, exc)
        else:
            print 'Downloading %s to %s' % (bkey, local_path)
            bucket.download_file(bkey, local_path,
                                 Callback=ProgressPercentage(bobj))

        # Print status
        per = 100*(float(idx+1)/num_files)
        print 'finished file %d/%d\n%f%% complete\n' % (idx+1, num_files, per)

    # Done iterating through list
    print 'done!'


# Upload files to AWS S3
def s3_upload(bucket, src_list, dst_list, make_public=False, encrypt=False):
    '''
    Function to upload a list of data to an S3 bucket

    Parameters
    ----------
    bucket : boto3 Bucket instance
        an instance of the boto3 S3 bucket class to upload to
    src_list : list (str)
        list of filepaths as strings to upload to S3
    dst_list : list (str)
        list of filepaths as strings coinciding with src_list, such
        that src_list[1] gets uploaded to S3 with the S3 path given in
        dst_list[1]
    make_public : boolean (optional), default=False
        set to True if files should be publically read-able on S3
    encrypt : boolean (optional), default=False
        set to True if the uploaded files should overwrite what is
        already there

    Returns
    -------
    None
        The function doesn't return any value, it uploads data to S3
        and prints its progress and a 'done' message upon completion
    '''

    # Import packages
    import hashlib
    from botocore.exceptions import ClientError

    # Init variables
    num_files = len(src_list)
    s3_str = 's3://'
    extra_args = {}

    # If make public, pass to extra args
    if make_public:
        extra_args['ACL'] = 'public-read'

    # If encryption is desired init extra_args
    if encrypt:
        extra_args['ServerSideEncryption'] = 'AES256'

    # Check if the list lengths match 
    if num_files != len(dst_list):
        raise RuntimeError, 'src_list and dst_list must be the same length!'

    # For each source file, upload
    for idx, src_file in enumerate(src_list):
        # Get destination path
        dst_file = dst_list[idx]

        # Check for s3_prefix
        if src_file.startswith(s3_str):
            bucket_name = src_file.split('/')[2]
            src_file = src_file.replace(s3_str+bucket_name, '').lstrip('/')
        elif dst_file.startswith(s3_str):
            bucket_name = dst_file.split('/')[2]
            dst_file = dst_file.replace(s3_str+bucket_name, '').lstrip('/')

        # Print status
        print 'Uploading %s to S3 bucket %s as %s' % \
        (src_file, bucket.name, dst_file)

        # Create a new key from the bucket and set its contents
        dst_key = bucket.Object(key=dst_file)

        # See if need to upload
        try:
            # If it exists, compare md5sums
            dst_key_exists = dst_key.get()
            dst_md5 = str(dst_key.e_tag.strip('"'))
            src_read = open(src_file, 'rb').read()
            src_md5 = hashlib.md5(src_read).hexdigest()
            # If md5sums dont match, re-upload via except ClientError
            if src_md5 != dst_md5:
                bucket.upload_file(src_file, dst_file, ExtraArgs=extra_args,
                                   Callback=ProgressPercentage(src_file))
        except ClientError as exc:
            bucket.upload_file(src_file, dst_file, ExtraArgs=extra_args,
                               Callback=ProgressPercentage(src_file))

        per = 100*(float(idx+1)/num_files)
        print 'finished file %d/%d\n\n%f%% complete\n' % (idx+1, num_files, per)

    # Print when finished
    print 'Done!'


# Test write-access to bucket first
def test_bucket_access(creds_path, output_directory, subject_id):
    '''
    '''

    # Import packages
    import os
    import botocore.exceptions as bexc
    from CPAC.AWS import fetch_creds

    # Init variables
    s3_str = 's3://'
    test_file = '/tmp/test-output.txt'

    # Explicitly lower-case the "s3"
    if output_directory.lower().startswith(s3_str):
        out_dir_sp = output_directory.split('/')
        out_dir_sp[0] = out_dir_sp[0].lower()
        output_directory = '/'.join(out_dir_sp)

    # Get bucket name
    bucket_name = output_directory.replace(s3_str, '').split('/')[0]

    # Get bucket
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    # Create local file
    with open(test_file, 'w') as f:
        f.write('test123')
    f.close()

    # Formulate test ouput key in bucket path output directory
    rel_key_path = output_directory.replace(\
                   os.path.join(s3_str, bucket_name), '').lstrip('/')
    write_test_key = os.path.join(rel_key_path, 'test-output_%s.txt' % subject_id)

    # Attempt a write to bucket
    try:
        bucket.upload_file(test_file, write_test_key)
        print 'Confirmed S3 write access for CPAC output!'
        test_key = bucket.Object(key=write_test_key)
        test_key.delete()
        s3_write_access = True
    # Otherwise we set the access flag to false
    except bexc.ClientError:
        s3_write_access = False

    # Return the access flag
    return s3_write_access
