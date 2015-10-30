# CPAC/AWS/s3_sublist.py
#
# Author(s): Daniel Clark, 2015

'''
This module has functions to build a subject list from S3 and
local filepaths
'''

# Extract site-based scan parameters
def extract_scan_params(scan_params_csv):
    '''
    Function to extract the site-based scan parameters from a csv file
    and return a dictionary of their values

    Parameters
    ----------
    scan_params_csv : string
        filepath to the scan parameters csv file

    Returns
    -------
    site_dict : dictionary
        a dictionary where site names are the keys and the scan
        parameters for that site are the values stored as a dictionary
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


# Get site, ppant, session-level directory indicies
def return_dir_indices(path_template):
    '''
    Function to return the site, particpant, and session-level
    directory indicies based on splitting the path template by
    directory seperation '/'

    Parameters
    ----------
    path_template : string
        filepath template in the form of:
        's3://bucket_name/base_dir/{site}/{participant}/{session}/..
        ../file.nii.gz'

    Returns
    -------
    site_idx : integer
        the directory level of site folders
    ppant_idx : integer
        the directory level of participant folders
    sess_idx : integer
        the directory level of the session folders
    '''

    # Get folder level indices of site and subject - anat
    fp_split = path_template.split('/')
    site_idx = fp_split.index('{site}')
    ppant_idx = fp_split.index('{participant}')
    # Session level isn't required, but recommended
    try:
        sess_idx = fp_split.index('{session}')
    except ValueError as exc:
        sess_idx = ppant_idx+1

    # Return indices
    return site_idx, ppant_idx, sess_idx


# Return matching filepaths
def return_local_filepaths(path_template, bids_flag=False):
    '''
    Function to return the filepaths from local directories given a
    file pattern template

    Parameters
    ----------
    path_template : string
        filepath template in the form of:
        '/base_dir/{site}/{participant}/{session}/../file.nii.gz'
    bids_flag : boolean (optional); default=False
        flag to indicate if the dataset to gather is organized to the
        BIDS standard

    Returns
    -------
    local_filepaths : list
        a list of strings of the local filepaths
    '''

    # Import packages
    import glob
    import os

    # Gather local files
    if bids_flag:
        file_pattern = path_template
    else:
        file_pattern = path_template.replace('{site}', '*').\
                       replace('{participant}', '*').replace('{session}', '*')

    local_filepaths = glob.glob(file_pattern)

    # Get absolute paths
    local_filepaths = [os.path.abspath(fp) for fp in local_filepaths]

    # Print how many found
    num_local_files = len(local_filepaths)
    print 'Found %d files!' % num_local_files

    # Return the filepaths as a list
    return local_filepaths


# Return matching filepaths
def return_s3_filepaths(path_template, creds_path=None, bids_flag=False):
    '''
    Function to return the filepaths from an S3 bucket given a file
    pattern template and, optionally, credentials

    Parameters
    ----------
    path_template : string
        filepath template in the form of:
        's3://bucket_name/base_dir/{site}/{participant}/{session}/..
        ../file.nii.gz'; if bids_flag is set, path_template is just the
        base directory of the BIDS data set
    creds_path : string (optional); default=None
        filepath to a credentials file containing the AWS credentials
        to access the S3 bucket objects
    bids_flag : boolean (optional); default=False
        flag to indicate if the dataset to gather is organized to the
        BIDS standard

    Returns
    -------
    s3_filepaths : list
        a list of strings of the filepaths from the S3 bucket
    '''

    # Import packages
    import fnmatch
    import os
    from CPAC.AWS import fetch_creds

    # Check for errors
    if not bids_flag:
        if not ('{site}' in path_template and '{participant}' in path_template):
            err_msg = 'Please provide \'{site}\' and \'{particpant}\' in '\
                      'filepath template where site and participant-level '\
                      'directories are present'
            raise Exception(err_msg)

    # Init variables
    bucket_name = path_template.split('/')[2]
    s3_prefix = '/'.join(path_template.split('/')[:3])

    # Extract base prefix to search through in S3
    if bids_flag:
        prefix = path_template.split('*')[0].replace(s3_prefix, '').lstrip('/')
    else:
        prefix = path_template.split('{site}')[0].replace(s3_prefix, '').lstrip('/')

    # Attempt to get bucket
    try:
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)
    except Exception as exc:
        err_msg = 'There was an error in retrieving S3 bucket: %s.\nError: %s'\
                  %(bucket_name, exc)
        raise Exception(err_msg)

    # Get filepaths from S3 with prefix
    print 'Gathering files from S3 to parse...'
    s3_filepaths = []
    for s3_obj in bucket.objects.filter(Prefix=prefix):
        s3_filepaths.append(str(s3_obj.key))

    # Prepend 's3://bucket_name/' on found paths
    s3_filepaths = [os.path.join(s3_prefix, s3_fp) for s3_fp in s3_filepaths]

    # File pattern filter
    if bids_flag:
        file_pattern = path_template
    else:
        file_pattern = path_template.replace('{site}', '*').\
                       replace('{participant}', '*').replace('{session}', '*')

    # Get only matching s3 paths
    s3_filepaths = fnmatch.filter(s3_filepaths, file_pattern)

    # Print how many found
    num_s3_files = len(s3_filepaths)
    print 'Found %d files!' % num_s3_files

    # Return the filepaths as a list
    return s3_filepaths


# Get the pattern that captures scan type
def return_bids_template(base_dir, scan_type, creds_path=None):
    '''
    Function that returns the path template of the desired scan type
    from a BIDS dataset

    Parameters
    ----------
    base_dir : string
        base directory of the BIDS dataset
    scan_type : string
        type of scan; e.g. 'anat', 'func', etc.
    creds_path : string (optional); default=None
        filepath to a set of AWS credentials to access a BIDS dataset
        stored on S3 that isn't public

    Returns
    -------
    file_template : string
        regular expression-compatible file template indicating data
        path organization
    '''

    # Import packages
    import os
    from CPAC.AWS import fetch_creds

    # Init variables
    s3_str = 's3://'
    file_path = None

    # If base directory is in S3
    if base_dir.startswith(s3_str):
        bucket_name = base_dir.split('/')[2]
        s3_prefix = '/'.join(base_dir.split('/')[:3])

        # Extract base prefix to search through in S3
        prefix = base_dir.split('*')[0].replace(s3_prefix, '').lstrip('/')

        # Attempt to get bucket
        try:
            bucket = fetch_creds.return_bucket(creds_path, bucket_name)
        except Exception as exc:
            err_msg = 'There was an error in retrieving S3 bucket: %s.\nError: %s'\
                      %(bucket_name, exc)
            raise Exception(err_msg)

        # Get filepaths from S3 with prefix
        print 'Gathering files from S3 to parse...'
        for s3_obj in bucket.objects.filter(Prefix=prefix):
            file_path = s3_obj.key
            scan_dir = file_path.split('/')[-2]
            if scan_dir == scan_type:
                break
    # Else, the base directory is locally stored
    else:
        for root, dirs, files in os.walk(base_dir):
            if file_path:
                break
            for fil in files:
                file_path = os.path.join(root, fil)
                scan_dir = file_path.split('/')[-2]
                if fil.endswith('.nii.gz') and scan_dir == scan_type:
                    break
                else:
                    file_path = None

    # Now replace file_path intermediate dirs with *
    if file_path:
        rel_path = file_path.replace(base_dir, '').lstrip('/')
        interm_dirs = rel_path.split('/')[:-2]
        for imd in interm_dirs:
            file_path = file_path.replace(imd, '*')
    else:
        err_msg = 'Could not find any files in directory, check files!'
        raise Exception(err_msg)

    # Set template as any file *
    file_template = os.path.join(os.path.dirname(file_path), '*.nii.gz')

    # Return file pattern template
    return file_template


# Build the C-PAC subject list
def build_sublist(data_config_yml):
    '''
    Function to build a C-PAC-compatible subject list, given
    anatomical and functional file paths

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

    from CPAC.utils import bids_metadata

    # Init variables
    tmp_dict = {}
    s3_str = 's3://'

    # Load in data config from yaml to dictionary
    data_config_dict = yaml.load(open(data_config_yml, 'r'))

    # Get inclusion, exclusion, scan csv parameters, and output location
    data_format = data_config_dict['dataFormat'][0]
    bids_base_dir = data_config_dict['bidsBaseDir']
    anat_template = data_config_dict['anatomicalTemplate']
    func_template = data_config_dict['functionalTemplate']
    include_subs = data_config_dict['subjectList']
    exclude_subs = data_config_dict['exclusionSubjectList']
    include_sites = data_config_dict['siteList']
    scan_params_csv = data_config_dict['scanParametersCSV']
    sublist_outdir = data_config_dict['outputSubjectListLocation']
    sublist_name = data_config_dict['subjectListName']
    # Older data configs won't have this field
    try:
        creds_path = data_config_dict['awsCredentialsFile']
    except KeyError:
        creds_path = None

    # Change any 'None' to None of optional arguments
    if bids_base_dir == 'None':
        bids_base_dir = None
    if include_subs == 'None':
        include_subs = None
    if exclude_subs == 'None':
        exclude_subs = None
    if include_sites == 'None':
        include_sites = None
    if scan_params_csv == 'None':
        scan_params_csv = None
    if creds_path == 'None':
        creds_path = None

    # Get dataformat if BIDS or not
    if data_format == 'BIDS':
        bids_flag = True
        anat_template = return_bids_template(bids_base_dir, 'anat')
        func_template = return_bids_template(bids_base_dir, 'func')
    else:
        bids_flag = False

    # See if the templates are s3 files
    if s3_str in anat_template and s3_str in func_template or \
       s3_str in bids_base_dir:
        # Get anatomical filepaths from s3
        print 'Fetching anatomical files...'
        anat_paths = return_s3_filepaths(anat_template, creds_path, bids_flag)
        # Get functional filepaths from s3
        print 'Fetching functional files...'
        func_paths = return_s3_filepaths(func_template, creds_path, bids_flag)

    # If one is in S3 and the other is not, raise error - not supported
    elif (s3_str in anat_template and s3_str not in func_template) or \
         (s3_str not in anat_template and s3_str in func_template):
        err_msg = 'Both anatomical and functional files should either be '\
                  'on S3 or local. Separating the files is currently not '\
                  'supported.'
        raise Exception(err_msg)

    # Otherwise, it's local files
    else:
        # Get anatomical filepaths
        print 'Gathering anatomical files...'
        anat_paths = return_local_filepaths(anat_template, bids_flag)
        # Get functional filepaths
        print 'Gathering functional files...'
        func_paths = return_local_filepaths(func_template, bids_flag)

    # Get directory indicies
    # If data is BIDS
    if bids_flag:
        anat_site_idx = func_site_idx = None
        anat_path_dirs = anat_template.split('/')
        anat_ppant_idx = func_ppant_idx = anat_path_dirs.index('*')
        anat_sess_idx = func_sess_idx = len(anat_path_dirs)-3
        # If session index is ppant index, then no session dir
        if anat_ppant_idx == anat_sess_idx:
            anat_sess_idx = func_sess_idx = None
    # Get indices from {site} {ppant} {session} identifiers
    else:
        anat_site_idx, anat_ppant_idx, anat_sess_idx = \
            return_dir_indices(anat_template)
        func_site_idx, func_ppant_idx, func_sess_idx = \
            return_dir_indices(func_template)

    # Filter out unwanted anat and func filepaths
    anat_paths = filter_sub_paths(anat_paths, include_sites,
                                  include_subs, exclude_subs,
                                  anat_site_idx, anat_ppant_idx)
    print 'Filtered down to %d anatomical files' % len(anat_paths)
    func_paths = filter_sub_paths(func_paths, include_sites,
                                  include_subs, exclude_subs,
                                  func_site_idx, func_ppant_idx)
    print 'Filtered down to %d functional files' % len(func_paths)

    if len(anat_paths) == 0 or len(func_paths) == 0:
        err_msg = 'Unable to find any files after filtering sites and '\
                  'subjects! Check site and subject inclusion/exclusion fields!'
        raise Exception(err_msg)

    # Read in scan parameters and return site-based dictionary
    if scan_params_csv is not None:
        site_scan_params = extract_scan_params(scan_params_csv)

    # Iterate through file paths and build subject list
    for anat in anat_paths:
        anat_sp = anat.split('/')
        subj = anat_sp[anat_ppant_idx]
        sess = anat_sp[anat_sess_idx]
        if bids_flag:
            site = ''
        else:
            site = anat_sp[anat_site_idx]
        subj_d = {'anat' : anat, 'creds_path' : creds_path, 'rest' : {},
                  'subject_id' : subj, 'unique_id' : sess,
                  'scan_parameters': {}}
        tmp_key = '_'.join([subj, site, sess])
        tmp_dict[tmp_key] = subj_d

    # Now go through and populate functional scans dictionaries
    for func in func_paths:
        # Extract info from filepath
        func_sp = func.split('/')
        subj = func_sp[func_ppant_idx]
        sess = func_sp[func_sess_idx]
        if bids_flag:
            site = ''
            scan_params = bids_metadata.get_metadata_for_nifti(bids_base_dir,
                                                               func)
        else:
            site = func_sp[func_site_idx]
            if scan_params_csv is not None:
                try:
                    scan_params = site_scan_params[site]
                except KeyError as exc:
                    print 'Site %s missing from scan parameters csv, skipping...'\
                          % site
                    scan_params = None

        # If there is no scan sub-folder under session, make scan
        # the name of the image itself without extension
        if func_sess_idx == len(func_sp)-2:
            scan = func_sp[-1].split('.nii')[0]
        # Othwerwise, there use scan sub folder
        else:
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
        subj_d['rest'][scan] = func
        # And replace it back in the dictionary
        tmp_dict[tmp_key] = subj_d

    # Build a subject list from dictionary values
    sublist = [v for v in tmp_dict.values()]

    # Write subject list out to out location
    sublist_out_yml = os.path.join(sublist_outdir,
                                   'CPAC_subject_list_%s.yml' % sublist_name)
    with open(sublist_out_yml, 'w') as out_sublist:
        # Make sure YAML doesn't dump aliases (so it's more human read-able)
        noalias_dumper = yaml.dumper.SafeDumper
        noalias_dumper.ignore_aliases = lambda self, data: True
        out_sublist.write(yaml.dump(sublist, default_flow_style=False,
                                    Dumper=noalias_dumper))

    # Return the subject list
    return sublist
