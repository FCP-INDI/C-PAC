# CPAC/AWS/s3_sublist.py
#
# Author(s): Daniel Clark, 2015

'''
This module has functions to build a subject list from S3 and
local filepaths
'''

# Init global variables
site_kw = '{site}'
ppant_kw = '{participant}'
sess_kw = '{session}'
ser_kw = '{series}'
kw_strs = [site_kw, ppant_kw, sess_kw, ser_kw]

# Check format of filepath templates
def check_template_format(file_template, site_kw, ppant_kw, sess_kw, ser_kw):
    '''
    Function to validate the file templalte contains all required
    keywords

    Parameters
    ----------
    file_template : string
        the file pattern template to check
    site_kw : string
        the keyword indicating the site-level strings in the filepaths
    ppant_kw : string
        the keyword indicating the ppant-level strings in the filepaths
    sess_kw : string
        the keyword indicating the session-level strings in the filepaths
    ser_kw : string
        the keyword indicating the series-level strings in the filepaths

    Returns
    -------
    None

    Raises
    ------
    Exception
        raises an exception if the filepaths don't contain the required
        keywords or if the filepaths have repeated keywords
    '''

    # Check for ppant, series-level directories
    if not (ppant_kw in file_template and ser_kw in file_template):
        err_msg = 'Please provide \'%s\' and \'%s\' level directories in '\
                  'filepath template where participant and series-level '\
                  'directories are present in file template: %s' \
                  % (ppant_kw, ser_kw, file_template)
        raise Exception(err_msg)

    # Check to make sure all keywords are only used once in template
    if file_template.count(site_kw) > 1 or file_template.count(ppant_kw) > 1 or \
       file_template.count(sess_kw) > 1 or file_template.count(ser_kw) > 1:
        err_msg = 'There are multiple instances of key words in the provided '\
                  'file template: %s. Fix this and try again.' % file_template
        raise Exception(err_msg)


# Extract keyword from filepath
def extract_keyword_from_path(filepath, keyword, template):
    '''
    Function to extract the key string from a filepath, given a
    keyword and a file pattern template

    Parameters
    ----------
    filepath : string
        filepath to the file of interest
    keyword : string
        string of the keyword
    template : string
        file pattern template containing the keyword

    Returns
    -------
    key_str : string
        extracted string where the keyword was located in the filepath 
    '''

    # Init variables
    temp_split = template.split('/')
    fp_split = filepath.split('/')

    # Extract directory  name of the 
    kw_dirname = [dir for dir in temp_split if keyword in dir]

    # If the keyword is in the template, extract string from filepath
    if len(kw_dirname) == 1:
        # Get the directory fullname from template, as well as any surrounding
        kw_dirname = kw_dirname[0]
        kw_idx = temp_split.index(kw_dirname)
        # Extract directory with key string in it from filepath
        key_str = fp_split[kw_idx]

        # Get the prefix and suffix surrounding keyword
        kw_prefix = kw_dirname.split(keyword)[0]
        kw_suffix = kw_dirname.split(keyword)[-1]

        # Replace other keywords in prefix/suffix with wildcards '*'
        for kw in kw_strs:
            # If keyword was found, grab any text in that position to split out
            if kw in kw_prefix:
                kw_prefix = kw_prefix.replace(kw, '*')
            if kw in kw_suffix:
                kw_suffix = kw_suffix.replace(kw, '*')

        # Strip out prefix patterns
        # If it ends with a '*', get everything before it
        while kw_prefix.endswith('*'):
            kw_prefix = kw_prefix[:-1]

        # Make sure what is left is more than ''
        if kw_prefix != '':
            # Find the previous '*' from the right
            prev_star_in_prefix = kw_prefix.rfind('*')
            # If there is '*', grab from it to end of prefix as delim
            if prev_star_in_prefix > 0:
                prefix_delim = kw_prefix[prev_star_in_prefix+1:]
            # Otherwise, just use the whole prefix as delim
            else:
                prefix_delim = kw_prefix

            # Split the filepath by prefix delim
            prefix_list = key_str.split(prefix_delim)

            # If it was found and split-able, take everything past the first
            # ocurrence of the delim (because delim could also be in suffix)
            if len(prefix_list) > 1:
                key_str = prefix_delim.join(prefix_list[1:])
            # Othwerwise, it wasn't found or split-able, use what we had
            else:
                key_str = prefix_list[0]

        # Strip out suffix patterns
        # If it starts with a '*', get everything after it
        while kw_suffix.startswith('*'):
            kw_suffix = kw_suffix[1:]

        # Make sure what is left is more than ''
        if kw_suffix != '':
            # Find the next '*' from the left
            next_star_in_suffix = kw_suffix.find('*')
            # If there is another '*', grab non-wildcards up until '*' as delim
            if next_star_in_suffix > 0:
                suffix_delim = kw_suffix[:next_star_in_suffix]
            # Otherwise, just use the whole prefix as delim
            else:
                suffix_delim = kw_suffix

            # Split the filepath by suffix delim
            suffix_list = key_str.split(suffix_delim)

            # If it was found and split-able, take everything up to the last
            # ocurrence of the delim (because delim could also be in prefix)
            if len(suffix_delim) > 1:
                key_str = suffix_delim.join(suffix_list[:-1])
            # Otherwise, it wasn't found or split-able, use what we had
            else:
                key_str = suffix_list[0]

        # Check to see if we split out everything, if so, just grab
        # whole directory
        if key_str == '':
            print 'Could not distinguish %s from filepath %s using the file ' \
                  'pattern template %s.\nInstead, using entire directory: %s ' \
                  'for keyword %s.\nCheck data organization and file ' \
                  'pattern template' % (keyword, filepath, template, keyword)
            key_str = fp_split[kw_idx]
    else:
        print 'Keyword %s not found in template %s' % (keyword, template)
        key_str = ''

    # Remove any nifti extensions
    key_str = key_str.rstrip('.gz').rstrip('nii')

    # Return the key string
    return key_str


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
def filter_sub_paths(sub_paths, include_sites, include_subs, exclude_subs,
                     site_kw, ppant_kw, path_template):
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
    path_template : string
        filepath template in the form of:
        '.../base_dir/{site}/{participant}/{session}/..
        {series}/file.nii.gz'

    Returns
    -------
    keep_subj_paths : list
        a list of the filepaths to use in the filtered subject list
    '''

    # Check if {site} was specified
    if site_kw in path_template and include_sites is not None:
        # Filter out sites that are not included
        keep_site_paths = []
        if type(include_sites) is not list:
            include_sites = [include_sites]
        print 'Only including sites: %s' % include_sites
        site_matches = filter(lambda sp: \
                              extract_keyword_from_path(sp,
                                                        site_kw,
                                                        path_template) in \
                              include_sites, sub_paths)
        keep_site_paths.extend(site_matches)
    else:
        print 'Not filtering out any potential sites...'
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
        subj_matches = filter(lambda sp: \
                              extract_keyword_from_path(sp,
                                                        ppant_kw,
                                                        path_template) in \
                              include_subs, sub_paths)
        keep_subj_paths.extend(subj_matches)
    # Or exclude only
    elif exclude_subs is not None:
        keep_subj_paths = []
        if type(exclude_subs) is not list:
            exclude_subs = [exclude_subs]
        print 'Including all subjects but: %s' % exclude_subs
        subj_matches = filter(lambda sp: \
                              extract_keyword_from_path(sp,
                                                        ppant_kw,
                                                        path_template) not in \
                              exclude_subs, sub_paths)
        keep_subj_paths.extend(subj_matches)

    else:
        keep_subj_paths = keep_site_paths

    # Filter out any duplicates
    keep_subj_paths = list(set(keep_subj_paths))
    print 'Filtered down to %d files' % len(keep_subj_paths)

    # Return kept paths
    return keep_subj_paths


# Return matching filepaths
def return_local_filepaths(file_pattern):
    '''
    Function to return the filepaths from local directories given a
    file pattern template

    Parameters
    ----------
    file_pattern : string
        regexp and glob compatible file pattern for gathering local
        files

    Returns
    -------
    local_filepaths : list
        a list of strings of the local filepaths
    '''

    # Import packages
    import glob
    import os

    # Gather local files
    local_filepaths = glob.glob(file_pattern)

    # Get absolute paths
    local_filepaths = [os.path.abspath(fp) for fp in local_filepaths]

    # Print how many found
    num_local_files = len(local_filepaths)
    print 'Found %d files!' % num_local_files

    # Return the filepaths as a list
    return local_filepaths


# Return matching filepaths
def return_s3_filepaths(file_pattern, creds_path=None):
    '''
    Function to return the filepaths from an S3 bucket given a file
    pattern template and, optionally, credentials

    Parameters
    ----------
    file_pattern : string
        regexp and glob compatible file pattern for gathering local
        files
    creds_path : string (optional); default=None
        filepath to a credentials file containing the AWS credentials
        to access the S3 bucket objects

    Returns
    -------
    s3_filepaths : list
        a list of strings of the filepaths from the S3 bucket
    '''

    # Import packages
    import fnmatch
    import os
    import re

    from CPAC.AWS import fetch_creds

    # Init variables
    bucket_name = file_pattern.split('/')[2]
    s3_prefix = '/'.join(file_pattern.split('/')[:3])

    # Find non regular expression patterns to get prefix
    s3_rel_path = file_pattern.replace(s3_prefix, '').lstrip('/')
    fp_split = s3_rel_path.split('/')
    for idx, dir in enumerate(fp_split):
        # Search for characters that we're permitting to be in folder names
        ok_chars = [char in dir for char in ['_', '-', '^']]
        if any(ok_chars):
            continue

        # Put escape '\' character in front of any non alphanumeric character
        reg_dir = re.escape(dir)

        # If we find a non alpha-numeric character in string, break
        if reg_dir != dir:
            break

    # Use all dirs up to first directory with non alpha-numeric character
    # as prefix
    prefix = '/'.join(fp_split[:idx])

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
    s3_filepaths = fnmatch.filter(s3_filepaths, file_pattern)

    # Print how many found
    num_s3_files = len(s3_filepaths)
    print 'Found %d files!' % num_s3_files

    # Return the filepaths as a list
    return s3_filepaths


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

    # Init variables
#     site_kw = '{site}'
#     ppant_kw = '{participant}'
#     sess_kw = '{session}'
#     ser_kw = '{series}'
#     kw_strs = [site_kw, ppant_kw, sess_kw, ser_kw]

    tmp_dict = {}
    s3_str = 's3://'

    # Load in data config from yaml to dictionary
    data_config_dict = yaml.load(open(data_config_yml, 'r'))

    # Get inclusion, exclusion, scan csv parameters, and output location
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

    # Check templates for proper formatting
    check_template_format(anat_template, site_kw, ppant_kw, sess_kw, ser_kw)
    check_template_format(func_template, site_kw, ppant_kw, sess_kw, ser_kw)

    # Replace keywords with glob regex wildcards
    anat_pattern = anat_template.replace(site_kw, '*').replace(ppant_kw, '*').\
                   replace(sess_kw, '*').replace(ser_kw, '*')
    func_pattern = func_template.replace(site_kw, '*').replace(ppant_kw, '*').\
                   replace(sess_kw, '*').replace(ser_kw, '*')

    # See if the templates are s3 files
    if s3_str in anat_template.lower() and s3_str in func_template.lower():
        # Get anatomical filepaths from s3
        print 'Fetching anatomical files...'
        anat_paths = return_s3_filepaths(anat_pattern, creds_path)
        # Get functional filepaths from s3
        print 'Fetching functional files...'
        func_paths = return_s3_filepaths(func_pattern, creds_path)

    # If one is in S3 and the other is not, raise error - not supported
    elif (s3_str in anat_template.lower() and s3_str not in func_template.lower()) or \
         (s3_str not in anat_template.lower() and s3_str in func_template.lower()):
        err_msg = 'Both anatomical and functional files should either be '\
                  'on S3 or local. Separating the files is currently not '\
                  'supported.'
        raise Exception(err_msg)

    # Otherwise, it's local files
    else:
        # Get anatomical filepaths
        print 'Gathering anatomical files...'
        anat_paths = return_local_filepaths(anat_pattern)
        # Get functional filepaths
        print 'Gathering functional files...'
        func_paths = return_local_filepaths(func_pattern)

    # Filter out unwanted anat and func filepaths
    print 'Filtering anatomical files...'
    anat_paths = filter_sub_paths(anat_paths, include_sites,
                                  include_subs, exclude_subs,
                                  site_kw, ppant_kw, anat_template)

    print 'Filtering functional files...'
    func_paths = filter_sub_paths(func_paths, include_sites,
                                  include_subs, exclude_subs,
                                  site_kw, ppant_kw, func_template)

    # If all data is filtered out, raise exception
    if len(anat_paths) == 0 or len(func_paths) == 0:
        err_msg = 'Unable to find any files after filtering sites and ' \
                  'subjects! Check site and subject inclusion fields as well ' \
                  'as filepaths and template pattern!'
        raise Exception(err_msg)

    # Read in scan parameters and return site-based dictionary
    if scan_params_csv is not None:
        site_scan_params = extract_scan_params(scan_params_csv)
    else:
        site_scan_params = {}

    # Iterate through file paths and build subject list
    for anat in anat_paths:
        site = extract_keyword_from_path(anat, site_kw, anat_template)
        ppant = extract_keyword_from_path(anat, ppant_kw, anat_template)
        session = extract_keyword_from_path(anat, sess_kw, anat_template)
        series = extract_keyword_from_path(anat, ser_kw, anat_template)

        # Init dictionary
        subj_d = {'anat' : anat, 'creds_path' : creds_path, 'rest' : {},
                  'subject_id' : ppant, 'unique_id' : '_'.join([site, session])}

        # Check for scan parameters
        if scan_params_csv is not None:
            try:
                subj_d['scan_parameters'] = site_scan_params[site]
            except KeyError as exc:
                print 'Site %s missing from scan parameters csv, skipping...'\
                      % site

        # Test to make sure subject key isn't present already
        # Should be a unique entry for every anatomical image
        tmp_key = '_'.join([site, ppant, session])
        if tmp_dict.has_key(tmp_key):
            err_msg = 'Key for anatomical file already exists: %s\n'\
                      'Either duplicate scan or data needs re-organization to '\
                      'differentiate between subjects from different sites/sessions'
            raise Exception(err_msg)

        tmp_dict[tmp_key] = subj_d

    # Now go through and populate functional scans dictionaries
    for func in func_paths:
        site = extract_keyword_from_path(func, site_kw, func_template)
        ppant = extract_keyword_from_path(func, ppant_kw, func_template)
        session = extract_keyword_from_path(func, sess_kw, func_template)
        series = extract_keyword_from_path(func, ser_kw, func_template)

        # Build tmp key and get subject dictionary from tmp dictionary
        tmp_key = '_'.join([site, ppant, session])

        # Try and find the associated anat scan
        try:
            subj_d = tmp_dict[tmp_key]
        except KeyError as exc:
            print 'Unable to find anatomical image for %s. Skipping...' \
                  % tmp_key
            continue

        # Set the rest dictionary with the scan
        subj_d['rest'][series] = func
        # And replace it back in the dictionary
        tmp_dict[tmp_key] = subj_d

    # Build a subject list from dictionary values
    sublist = []
    for data_bundle in tmp_dict.values():
        if data_bundle['anat'] != '' and data_bundle['rest'] != {}:
            sublist.append(data_bundle)
    # Check to make sure subject list has at least one valid data bundle
    if len(sublist) == 0:
        err_msg = 'Unable to find both anatomical and functional data for any '\
                  'subjects. Check data organization and file patterns!'
        raise Exception(err_msg)

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
