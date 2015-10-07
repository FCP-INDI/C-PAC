# CPAC/AWS/s3_sublist.py
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
    prefix = file_pattern.split('%s')[0]

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


# Read in data config file and return s3 parameters
def return_s3_params(data_config_dict):
    '''
    Function to read in the data config dictionary C-PAC generates and
    extract the relevant paramters to build a subject list from an S3
    bucket

    Parameters
    ----------
    data_config_dict : dictionary
        data config dictionary loaded in from yaml file

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

    # Init variables
    s3_str = 's3://'

    # Check anat/func filepath templates
    anat_fp_template = data_config_dict['anatomicalTemplate']
    func_fp_template = data_config_dict['functionalTemplate']

    # If they are S3 resources, expect 's3://bucket_name/...' first
    if anat_fp_template.startswith(s3_str) and \
       func_fp_template.startswith(s3_str):

        # Get relative filepath templates
        anat_template = anat_fp_template.lstrip(s3_str)
        func_template = func_fp_template.lstrip(s3_str)

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
        anat_template = anat_template.lstrip(bucket_name).lstrip('/')
        func_template = func_template.lstrip(bucket_name).lstrip('/')

    # See if a credentials path was specified
    try:
        creds_path = data_config_dict['credentialsPath']
    except KeyError:
        creds_path = None
        print 'No AWS credentials path specified, accessing data anonymously...'

    # Return arguments
    return creds_path, bucket_name, anat_template, func_template


# Filter out unwanted sites/subjects
def filter_sub_paths(sub_paths, include_sites, include_subs,
                     exclude_subs, site_idx, subj_idx):
    '''
    Function to filter out unwanted sites and subjects from the
    collected filepaths

    Parameters
    ----------
    sub_paths : list
        a list of the subject data files
    include_sites : list or string
        indicates which sites to keep filepaths from
    include_subs : list or string
        indicates which subjects to keep filepaths from
    exclude_subs : list or string
        indicates which subjects to remove filepaths from
    site_idx : integer
        index of filepath, split by '/', where site directory is
    subj_idx : integer
        index of filepath, split by '/', where subj directory is

    Returns
    -------
    keep_subj_paths : list
        a list of the filepaths to use in the filtered subject list
    '''

    # Filter out sites that are not included
    if include_sites is not None:
        keep_site_paths = []
        if type(include_sites) is not list:
            include_sites = [include_sites]
        print 'Only including sites: %s' % include_sites
        for site in include_sites:
            site_matches = filter(lambda sp: sp.split('/')[site_idx] == site,
                                  sub_paths)
            keep_site_paths.extend(site_matches)
            #sub_paths = [sp for sp in sub_paths if sp not in site_matches]
    else:
        keep_site_paths = sub_paths

    # Only keep subjects in inclusion list or remove those in exclusion list
    if include_subs is not None and exclude_subs is not None:
        err_msg = 'Please only populate subjects to include or exclude '\
                  '- not both!'
        raise Exception(err_msg)
    # Include only
    elif include_subs is not None:
        keep_subj_paths = []
        if type(include_subs) is not list:
            include_subs = [include_subs]
        print 'Only including subjects: %s' % include_subs
        for inc_sub in include_subs:
            sub_matches = filter(lambda sp: \
                                 sp.split('/')[subj_idx] == inc_sub,
                                 keep_site_paths)
            keep_subj_paths.extend(sub_matches)
    # Or exclude only
    elif exclude_subs is not None:
        keep_subj_paths = []
        if type(exclude_subs) is not list:
            exclude_subs = [exclude_subs]
        print 'Including all subjects but: %s' % exclude_subs
        sub_excludes = []
        for exc_sub in exclude_subs:
            sub_excludes.extend(filter(lambda sp: \
                                sp.split('/')[subj_idx] == exc_sub,
                                keep_site_paths))
        # Prune out exclusions
        keep_subj_paths = [sp for sp in keep_site_paths \
                           if sp not in sub_excludes]
    else:
        keep_subj_paths = keep_site_paths

    # Return kept paths
    return keep_subj_paths


# Extract site-based scan parameters
def extract_scan_params(scan_params_csv):
    '''
    '''

    # Import packages
    import csv

    # Init variables
    csv_open = open(scan_params_csv, 'r')
    site_dict = {}

    # Init csv dictionary reader
    reader = csv.DictReader(csv_open)

    # Iterate through the csv and pull in parameters
    for dict_row in reader:
        site = dict_row['Site']
        site_dict[site] = {key.lower() : val for key, val in dict_row.items()\
                           if key != 'Site'}
        # Assumes all other fields are formatted properly, but TR might not
        site_dict[site]['tr'] = site_dict[site].pop('tr (seconds)')

    # Return site dictionary
    return site_dict


# Build the S3 C-PAC subject list
def build_s3_sublist(data_config_yml):
    '''
    Function to build a C-PAC-compatible subject list from S3, given
    anatomical and functional file paths and AWS credentials

    Parameters
    ----------
    data_config_yml : string
        filepath to the data config yaml file

    Returns
    -------
    sublist : list
        C-PAC subject list of subject dictionaries
    '''

    # Import packages
    import os
    import yaml
    from CPAC.AWS import fetch_creds

    # Init variables
    tmp_dict = {}

    # Load in data config from yaml to dictionary
    data_config_dict = yaml.load(open(data_config_yml, 'r'))

    # Extract s3_prefix
    s3_prefix = '/'.join(data_config_dict['anatomicalTemplate'].split('/')\
                         [:3])

    # Get s3-specific paramters from data config
    creds_path, bucket_name, anat_fp, func_fp = \
        return_s3_params(data_config_dict)

    # Attempt to get bucket
    try:
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)
    except Exception as exc:
        err_msg = 'There was an error in retrieving S3 bucket: %s.\nError: %s'\
                  %(bucket_name, exc)
        raise Exception(err_msg)

    # Get anatomical filepaths from s3
    print 'Fetching anatomical files...'
    anat_s3_paths = return_s3_filepaths(bucket, anat_fp)
    # Get functional filepaths from s3
    print 'Fetching functional files...'
    func_s3_paths = return_s3_filepaths(bucket, func_fp)

    # Get folder level indices of site and subject - anat
    fp_split = anat_fp.split('/')
    anat_site_idx = fp_split.index('%s')
    anat_subj_idx = fp_split[anat_site_idx+1:].index('%s') + anat_site_idx+1

    # Get folder level indices of site and subject - func
    fp_split = func_fp.split('/')
    func_site_idx = fp_split.index('%s')
    func_subj_idx = fp_split[func_site_idx+1:].index('%s') + func_site_idx+1

    # Get inclusion, exclusion, scan csv parameters, and output location
    include_subs = data_config_dict['subjectList']
    exclude_subs = data_config_dict['exclusionSubjectList']
    include_sites = data_config_dict['siteList']
    scan_params_csv = data_config_dict['scanParametersCSV']
    sublist_outdir = data_config_dict['outputSubjectListLocation']
    sublist_name = data_config_dict['subjectListName']
    sublist_out_yml = os.path.join(sublist_outdir,
                                   sublist_name + '_subject_list.yml')

    # Filter out unwanted anat and func filepaths
    anat_s3_paths = filter_sub_paths(anat_s3_paths, include_sites,
                                     include_subs, exclude_subs,
                                     anat_site_idx, anat_subj_idx)
    print 'Filtered down to %d anatomical files' % len(anat_s3_paths)

    func_s3_paths = filter_sub_paths(func_s3_paths, include_sites,
                                     include_subs, exclude_subs,
                                     func_site_idx, func_subj_idx)
    print 'Filtered down to %d functional files' % len(func_s3_paths)

    # Read in scan parameters and return site-based dictionary
    if scan_params_csv is not None:
        site_scan_params = extract_scan_params(scan_params_csv)

    # Iterate through file paths and build subject list
    for anat in anat_s3_paths:
        anat_sp = anat.split('/')
        site = anat_sp[anat_site_idx]
        subj = anat_sp[anat_subj_idx]
        sess = anat_sp[anat_subj_idx+1]
        subj_d = {'anat' : os.path.join(s3_prefix, anat), 'rest' : {},
                  'subject_id' : subj, 'unique_id' : sess}
        if scan_params_csv is not None:
            try:
                subj_d['scan_parameters'] = site_scan_params[site]
            except KeyError as exc:
                print 'Site %s missing from scan parameters csv, skipping...'\
                      % site
        tmp_key = '_'.join([subj, site, sess])
        tmp_dict[tmp_key] = subj_d

    # Now go through and populate functional scans dictionaries
    for func in func_s3_paths:
        # Extract info from filepath
        func_sp = func.split('/')
        site = func_sp[func_site_idx]
        subj = func_sp[func_subj_idx]
        sess = func_sp[func_subj_idx+1]
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

    # Write subject list out to out location
    with open(sublist_out_yml, 'w') as out_sublist:
        out_sublist.write(yaml.dump(sublist))

    # Return the subject list
    return sublist

