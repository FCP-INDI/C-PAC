# CPAC/AWS/aws_utils.py
#
# Contributing authors:
# Daniel Clark

'''
This module contains functions which assist in interacting with AWS
services, including uploading/downloading data and file checking.
'''

# Build and download a subject list
def build_download_sublist(bucket, bucket_prefix, local_prefix, sub_list):
    '''
    Function to build and download a subject list from S3
    
    Parameters
    ----------
    bucket : boto.s3.bucket.Bucket instance
        an instance of the boto S3 bucket class to download from
    bucket_prefix : string
        the bucket prefix where all of the file keys are located
    local_prefix : string
        the local disk prefix where all of the files should be downloaded to
    sub_list : list (dict)
        the C-PAC subject list that has inputs in local_prefix
    
    Returns
    -------
    None
        this function does not return a value, it downloads the subjects from
        the C-PAC subject list to disk from S3
    '''
    
    # Import packages
    import os
    
    # Init variables
    local_list = []
    for sub_dict in sub_list:
        local_list.append(sub_dict['anat'])
        local_list.extend([v for v in sub_dict['rest'].values()])

    # Substitute the prefixes to build S3 list to download from
    s3_list = [l.replace(local_prefix, bucket_prefix) for l in local_list]

    # Check already-existing files and remove from download lists
    local_rm = []
    s3_rm = []
    # Build remove-lists
    for i in range(len(local_list)):
        l = local_list[i]
        s = s3_list[i]
        if os.path.exists(l):
            local_rm.append(l)
            s3_rm.append(s)
    # Go through remove lists and remove files
    for l, s in zip(local_rm, s3_rm):
        local_list.remove(l)
        s3_list.remove(s)

    # Download the data to the local prefix
    s3_download(bucket, s3_list, local_prefix, bucket_prefix=bucket_prefix)
    
    # Check to see they all downloaded successfully
    for l in local_list:
        l_path = os.path.abspath(l)
        if not os.path.exists(l_path):
            raise IOError('S3 files were not all downloaded.\n'\
                          'Could not find: %s' % l_path)


# Collect all files in directory as the source list
def collect_subject_files(prefix_star, sub_id):
    '''
    Function to collect all of the files in a directory into a list of
    full paths
    
    Parameters
    ----------
    prefix_star : string
        filepath to the folder, in which, all of the sub-files are
        collected; this filepath should have a wildcard character of
        '*' so that glob can collect the files via the pattern given
    sub_id : string
        the subject id to look for in the output folder
    
    Returns
    -------
    src_list : list [str]
        a list of filepaths (as strings)
    '''

    # Import packages
    import glob
    import os

    # Init variables
    bases = glob.glob(prefix_star)
    src_list = []
    
    # For each pipeline
    for base in bases:
        # Iterate through directory
        for root, dirs, files in os.walk(base):
            # If it's in the subject's folder and there are files
            if sub_id in root and files:
                src_list.extend([os.path.join(root, f) for f in files])

    # Return the list
    return src_list


# Get the MD5 sums of files on S3
def md5_sum(bucket, prefix='', filt_str=''):
    '''
    Function to get the filenames and MD5 checksums of files stored in
    an S3 bucket and return this as a dictionary.

    Parameters
    ----------
    bucket : boto.s3.bucket.Bucket instance
        an instance of the boto S3 bucket class to download from
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
    blist = bucket.list(prefix)
    md5_dict = {}

    # And iterate over keys to copy over new ones
    for fkey in blist:
        filename = str(fkey.key)
        if filt_str in filename:
            md5_sum = str(fkey.etag).strip('"')
            md5_dict[filename] = md5_sum
            print 'filename: %s' % filename
            print 'md5_sum: %s' % md5_sum

    # Return the dictionary
    return md5_dict


# Rename s3 keys from src_list to dst_list
def s3_rename(bucket, src_list, dst_list,
              keep_old=False, overwrite=False, make_public=False):
    '''
    Function to rename files from an AWS S3 bucket via a copy and delete
    process. Uses all keys in src_list as the original names and renames
    the them to the corresponding keys in the dst_list.
    (e.g. src_list[9] --> dst_list[9])

    Parameters
    ----------
    bucket : boto.s3.bucket.Bucket instance
        an instance of the boto S3 bucket class to download from
    src_list : list (str)
        a list of relative paths of the files to delete from the bucket
    dst_list : list (str)
        a list of relative paths of the files to delete from the bucket
    keep_old : boolean (optional), default=False
        flag indicating whether to keep the src_list files
    overwrite : boolean (optional), default=False
        flag indicated whether to overwrite the files in dst_list
    make_public : boolean (optional), default=False
        set to True if files should be publically available on S3
    Returns
    -------
    None
        The function doesn't return any value, it deletes data from
        S3 and prints its progress and a 'done' message upon completion
    '''

    # Check list lengths are equal
    if len(src_list) != len(dst_list):
        raise ValueError('src_list and dst_list are different lengths!')

    # Init variables
    i = 0
    no_files = len(src_list)

    # And iterate over keys to copy over new ones
    for f in src_list:
        src_key = bucket.get_key(f)
        if not src_key:
            print 'source file %s doesnt exist, skipping...' % f
            continue
        dst_key = dst_list[i]
        dst_exists = bucket.get_key(dst_key)
        if not dst_exists or overwrite:
            print 'copying source: ', str(src_key.key)
            print 'to destination: ', dst_key
            src_key.copy(bucket.name, dst_key)
            if make_public:
                print 'making public...'
                dk = bucket.get_key(dst_key)
                dk.make_public()
            if not keep_old:
                src_key.delete()
        else:
            print '%s exists and not overwriting' % dst_key
        i += 1
        per = 100*(float(i)/no_files)
        print 'Done renaming %d/%d\n%f%% complete' % (i, no_files, per)


# Delete s3 keys based on input list
def s3_delete(bucket, in_list):
    '''
    Method to delete files from an AWS S3 bucket that have the same
    names as those of an input list to a local directory.

    Parameters
    ----------
    bucket : boto.s3.bucket.Bucket instance
        an instance of the boto S3 bucket class to download from
    in_list : list (str)
        a list of relative paths of the files to delete from the bucket

    Returns
    -------
    None
        The function doesn't return any value, it deletes data from
        S3 and prints its progress and a 'done' message upon completion
    '''

    # Init variables
    no_files = len(in_list)
    i = 0
    # Iterate over list and delete S3 items
    for f in in_list:
        i += 1
        try:
            print 'attempting to delete %s from %s...' % (f, bucket.name)
            k = bucket.get_key(f)
            k.delete()
            per = 100*(float(i)/no_files)
            print 'Done deleting %d/%d\n%f%% complete' % (i, no_files, per)
        except AttributeError:
            print 'No key found for %s on bucket %s' % (f, bucket.name)

        # Done iterating through list
        print 'done!'


# Download files from AWS S3 to local machine
def s3_download(bucket, in_list, local_prefix, bucket_prefix=''):
    '''
    Method to download files from an AWS S3 bucket that have the same
    names as those of an input list to a local directory.

    Parameters
    ----------
    bucket : boto.s3.bucket.Bucket instance
        an instance of the boto S3 bucket class to download from
    in_list : list (str)
        a list of relative paths of the files to download from the bucket
    local_prefix : string
        local directory prefix to store the downloaded data
    bucket_prefix : string (optional)
        bucket_prefix, if specified, will be substituted with
        local_prefix; otherwise, the local_prefix will only prepend the
        downloaded files

    Returns
    -------
    None
        The function doesn't return any value, it downloads data from
        S3 and prints its progress and a 'done' message upon completion
    '''

    # Impor packages
    import os

    # Init variables
    no_files = len(in_list)
    i = 0
    # Check for trailing '/'
    if not local_prefix.endswith('/'):
        local_prefix = local_prefix + '/'
    if bucket_prefix and not bucket_prefix.endswith('/'):
        bucket_prefix = bucket_prefix + '/'
    # For each item in the list, try to download it
    for f in in_list:
        i += 1
        remote_filename = bucket.name + ': ' + f
        if bucket_prefix:
            local_filename = f.replace(bucket_prefix, local_prefix)
        else:
            local_filename = os.path.join(local_prefix, f)
        # Check to see if the local folder setup exists or not
        local_folders = os.path.dirname(local_filename)
        if not os.path.isdir(local_folders):
            print 'creating %s on local machine' % local_folders
            os.makedirs(local_folders)
        # Attempt to download the file
        print 'attempting to download %s to %s...' % (remote_filename,
                                                      local_filename)
        try:
            if not os.path.exists(local_filename):
                k = bucket.get_key(f)
                k.get_contents_to_filename(local_filename)
                per = 100*(float(i)/no_files)
                print 'Done downloading %d/%d\n%f%% complete' % (i, no_files, per)
            else:
                print 'File %s already exists, skipping...' % local_filename
        except AttributeError:
            print 'No key found for %s on bucket %s' % (f, bucket.name)

    # Done iterating through list
    print 'done!'


# Upload files to AWS S3
def s3_upload(bucket, src_list, dst_list, make_public=False, overwrite=False):
    '''
    Function to upload a list of data to an S3 bucket

    Parameters
    ----------
    bucket : boto.s3.bucket.Bucket instance
        an instance of the boto S3 bucket class to download from
    src_list : list (str)
        list of filepaths as strings to upload to S3
    dst_list : list (str)
        list of filepaths as strings coinciding with src_list, such
        that src_list[1] gets uploaded to S3 with the S3 path given in
        dst_list[1]
    make_public : boolean (optional), default=False
        set to True if files should be publically available on S3
    overwrite : boolean (optional), default=False
        set to True if the uploaded files should overwrite what is
        already there

    Returns
    -------
    None
        The function doesn't return any value, it uploads data to S3
        and prints its progress and a 'done' message upon completion
    '''

    # Callback function for upload progress update
    def callback(complete, total):
        '''
        Method to illustrate file uploading and progress updates
        '''

        # Import packages
        import sys

        # Write ...'s to the output for loading progress
        sys.stdout.write('.')
        sys.stdout.flush()

    # Init variables
    no_files = len(src_list)
    i = 0

    # Check if the list lengths match 
    if no_files != len(dst_list):
        raise RuntimeError, 'src_list and dst_list must be the same length!'

    # For each source file, upload
    for src_file in src_list:
        # Get destination path
        dst_file = dst_list[i]
        # Print status
        print 'Uploading %s to S3 bucket %s as %s' % \
        (src_file, bucket.name, dst_file)

        # Create a new key from the bucket and set its contents
        k = bucket.new_key(dst_file)
        if k.exists() and not overwrite:
            print 'key %s already exists, skipping...' % dst_file
        else:
            k.set_contents_from_filename(src_file, cb=callback, replace=True)
        # Make file public if set to True
        if make_public:
            print 'make public()'
            k.make_public()
        i += 1
        per = 100*(float(i)/no_files)
        print 'finished file %d/%d\n%f%% complete\n' % \
        (i, no_files, per)

    # Print when finished
    print 'Done!'
