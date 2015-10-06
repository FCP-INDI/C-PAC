# CPAC/AWS/build_s3_sublist.py
#
# Author(s): Daniel Clark, 2015

'''
This module has functions to build a subject list from S3 filepaths
'''

# Return matching filepaths
def return_s3_filepaths(bucket, file_pattern):
    '''
    Function to return the filepaths from an S3 bucket given a file
    pattern

    Parameters
    ----------
    bucket : boto3.resource.Bucket object
        the bucket of interest to parse
    file_pattern : string
        file pattern with two '%s's indicating site and subject-level
        directories, respectively

    Returns
    -------
    s3_filepaths : list
        a list of strings of the filepaths from the S3 bucket
    '''

    # Import packages
    import fnmatch

    # Check for errors
    if file_pattern.count('%s') != 2:
        err_msg = 'Please provide \'%s\' in file pattern where site and '\
                  'subject folders are present'
        raise Exception(err_msg)

    # Init variables
    prefix = file_pattern.split('%s')[0].lstrip('/')

    # Get filepaths from S3 with prefix
    print 'Gathering files from S3 to parse...'
    s3_filepaths = []
    for s3_obj in bucket.objects.filter(Prefix=prefix):
        s3_filepaths.append(str(s3_obj.key))

    # File pattern filter
    s3_filepaths = fnmatch.filter(s3_filepaths, file_pattern.replace('%s', '*'))

    # Print how many found
    num_s3_files = len(s3_filepaths)
    print 'Found %d files!' % num_s3_files

    # Return the filepaths as a list
    return s3_filepaths


# Build the S3 C-PAC subject list
def build_s3_sublist(creds_path, bucket_name, anat_fp, func_fp):
    '''
    Function to build a C-PAC-compatible subject list from S3, given
    anatomical and functional file paths and AWS credentials

    Parameters
    ----------
    creds_path : string
        filepath to the AWS credentials file
    bucket_name : string
        the S3 bucket name to parse
    anat_fp : string
        the anatomical data file pattern
    func_fp : string
        the anatomical data file pattern

    Returns
    -------
    sublist : list
        C-PAC subject list of subject dictionaries
    '''

    # Import packages
    import os
    from CPAC.AWS import fetch_creds

    # Init variables
    tmp_dict = {}
    s3_str = 's3://'
    s3_prefix = s3_str + bucket_name

    # Attempt to get bucket
    try:
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)
    except Exception as exc:
        print 'There was an error in retrieving the bucket: %s' % bucket_name
        print 'Error message:\n%s' % exc.message
        return

    # Get anatomical filepaths from s3
    print 'Fetching anatomical files...'
    anat_s3_files = return_s3_filepaths(bucket, anat_fp)
    # Get functional filepaths from s3
    print 'Fetching functional files...'
    func_s3_files = return_s3_filepaths(bucket, func_fp)

    # Get folder level indices of site and subject
    fp_split = anat_fp.split('/')
    site_idx = fp_split.index('%s')
    subj_idx = fp_split[site_idx+1:].index('%s') + site_idx+1

    # Iterate through file paths and build subject list
    for anat in anat_s3_files:
        anat_sp = anat.split('/')
        site = anat_sp[site_idx]
        subj = anat_sp[subj_idx]
        sess = anat_sp[subj_idx+1]
        subj_d = {'anat' : os.path.join(s3_prefix, anat), 'rest' : {},
                  'subject_id' : subj, 'unique_id' : sess}
        tmp_key = '_'.join([subj, site, sess])
        tmp_dict[tmp_key] = subj_d

    # Now go through and populate functional scans dictionaries
    for func in func_s3_files:
        # Extract info from filepath
        func_sp = func.split('/')
        site = func_sp[site_idx]
        subj = func_sp[subj_idx]
        sess = func_sp[subj_idx+1]
        scan = func_sp[-2]
        # Build tmp key and get subject dictionary from tmp dictionary
        tmp_key = '_'.join([subj, site, sess])
        # Try and find the associated anat scan
        try:
            subj_d = tmp_dict[tmp_key]
        except KeyError as exc:
            print 'Unable to find anatomical image for %s. Skipping...' \
                  % tmp_key
            continue
        # Set the rest dictionary with the scan
        subj_d['rest'][scan] = os.path.join(s3_prefix, func)
        # And replace it back in the dictionary
        tmp_dict[tmp_key] = subj_d

    # Build a subject list from dictionary values
    sublist = [v for v in tmp_dict.values()]

    # Return the subject list
    return sublist


# Read in data config file
def read_data_config(data_config_path):
    '''
    Function to read in the data config yaml file C-PAC generates and
    extract the relevant paramters to build a subject list from an S3
    bucket

    Parameters
    ----------
    data_config_path : string
        filepath to the data_config.yml file

    Returns
    -------
    creds_path : string
        filepath to the csv file containing AWS keys for bucket access
    bucket_name : string
        name of the S3 bucket to download subject data from
    anat_template : string
        the file pattern template for the subject anatomical data
    func_template : string
        the file pattern template for the subject functional data
    '''

    # Import packages
    import yaml

    # Init variables
    s3_prefix = 's3://'
    data_config = yaml.load(open(data_config_path, 'r'))

    # Check anat/func filepath templates
    anat_fp_template = data_config['anatomicalTemplate']
    func_fp_template = data_config['functionalTemplate']

    # If they are S3 resources, expect 's3://bucket_name/...' first
    if anat_fp_template.startswith(s3_prefix) and \
       func_fp_template.startswith(s3_prefix):

        # Get relative filepath templates
        anat_template = anat_fp_template.lstrip(s3_prefix)
        func_template = func_fp_template.lstrip(s3_prefix)

        # Get bucket names from each template path
        anat_bucket_name = anat_template.split('/')[0]
        func_bucket_name = func_template.split('/')[0]

        # In case bucket names dont match, raise Exception
        if anat_bucket_name != func_bucket_name:
            err_msg = 'Bucket names do not match. Currently downloading '\
                      'data from separate buckets is not supported.\n'\
                      'Anat bucket name: %s\nFunc bucket name: %s'\
                      % (anat_bucket_name, func_bucket_name)
            raise Exception(err_msg)
        else:
            bucket_name = anat_bucket_name

        # Strip bucket name from path templates as well
        anat_template = anat_template.lstrip(bucket_name)
        func_template = func_template.lstrip(bucket_name)

    # See if a credentials path was specified
    try:
        creds_path = data_config['credentialsPath']
    except KeyError:
        creds_path = None
        print 'No AWS credentials path specified, accessing data anonymously...'

    # Return arguments
    return creds_path, bucket_name, anat_template, func_template
