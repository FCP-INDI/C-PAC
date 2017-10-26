# utils/build_sublist.py
#

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


def read_subj_txtfile(filepath):
    with open(filepath,"r") as f:
        sub_ids = f.read().splitlines()
        sub_ids = filter(None,sub_ids)
    return sub_ids


# Check for glob-style patterns
def check_for_glob_patterns(delim, filepath, suffix_flag=False):
    '''
    Function to check the delimeter for any glob patterns (?, []);
    if some are found, they are replaced by the filepath's characters
    in those locations - this allows for a more accurate keyword
    extraction

    Parameters
    ----------
    delim : string
        the non-wildcard characters before/after the desired keyword
    filepath : string
        filepath to the file of interest
    suffix_flag : boolean
        flag to indicate if delimeter is prefix or suffix

    Returns
    -------
    delim : string
        the glob-matched delimeter with the filepath
    '''

    # If it's a suffix, traverse backwards
    if suffix_flag:
        delim = delim[::-1]
        filepath = filepath[::-1]
        first_brak = ']'
        last_brak = '['
    else:
        first_brak = '['
        last_brak = ']'

    # Init variables
    delim_list = list(delim)
    filepath_test = filepath
    in_brak_flg = False

    # Iterate through prefix delimeter
    for idx, char in enumerate(delim):
        if in_brak_flg:
            if char == last_brak:
                in_brak_flg = False
            delim_list[idx] = ''
            continue
        if char == '?' or char == first_brak:
            if char == first_brak:
                in_brak_flg = True
            const_str = ''.join(delim[:idx])
            fp_idx = filepath.find(const_str)
            filepath_test = filepath[fp_idx+len(const_str):]
            char_match = filepath_test[0]
            delim_list[idx] = char_match
            delim = ''.join(delim_list)

    # Join list
    delim = ''.join(delim_list)

    # If it's a suffix, reverse it
    if suffix_flag:
        delim = delim[::-1]

    # Return the glob-matches prefix delimiter
    return delim


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

    # Import packages
    import logging

    # Get logger
    logger = logging.getLogger('sublist_builder')

    # Check for ppant, series-level directories
    if not ppant_kw in file_template:
        err_msg = 'Please provide \'%s\' level directories in '\
                  'filepath template where participant-level '\
                  'directories are present in file template: %s' \
                  % (ppant_kw, file_template)
        logger.error(err_msg)
        raise Exception(err_msg)

    # Check to make sure all keywords are only used once in template
    if file_template.count(site_kw) > 1 or file_template.count(ppant_kw) > 1 or \
       file_template.count(sess_kw) > 1 or file_template.count(ser_kw) > 1:
        err_msg = 'There are multiple instances of key words in the provided '\
                  'file template: %s. Fix this and try again.' % file_template
        logger.error(err_msg)
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

    # Import packages
    import logging

    # Init variables
    # creates a list of directory levels of the filepath, and of the template
    temp_split = template.split('/')
    fp_split = filepath.split('/')

    # Get logger
    # logger = logging.getLogger('sublist_builder')

    # Extract directory name of the keyword, from the template
    kw_dirname = [dir for dir in temp_split if keyword in dir]

    # If the keyword is in the template, extract string from filepath
    if len(kw_dirname) > 0:
        # Get the directory fullname from template, as well as any surrounding
        kw_dirname = kw_dirname[0]
        # kw_idx is the index of where in the template path the keyword is
        kw_idx = temp_split.index(kw_dirname)
        # Extract directory with key string in it from filepath, i.e. if this
        # is for {participant}, key_str will be the participant ID string from
        # the filepath
        key_str = fp_split[kw_idx]
        # Get the prefix and suffix surrounding keyword
        kw_prefix = kw_dirname.split(keyword)[0]
        kw_suffix = kw_dirname.split(keyword)[1]

        # Replace other keywords in prefix/suffix with wildcards '*'
        for kw in kw_strs:
            # If keyword was found, grab any text in that position to split
            # out
            if kw in kw_prefix:
                kw_prefix = kw_prefix.replace(kw, '*')
            if kw in kw_suffix:
                kw_suffix = kw_suffix.replace(kw, '*')

        # Strip out prefix patterns
        # If it ends with a '*', get everything before it
        while kw_prefix.endswith('*'):
            kw_prefix = kw_prefix[:-1]

        # Make sure what is left is more than ''
        # This will not run if the keyword was the only thing in the directory
        # level
        if kw_prefix != '':
            # Find the previous '*' from the right
            prev_star_in_prefix = kw_prefix.rfind('*')
            # If there is '*', grab from it to end of prefix as delim
            if prev_star_in_prefix >= 0:
                prefix_delim = kw_prefix[prev_star_in_prefix+1:]
            # Otherwise, just use the whole prefix as delim
            else:
                prefix_delim = kw_prefix
            # Check for glob-style characters in delimeter
            prefix_delim = check_for_glob_patterns(prefix_delim, key_str)

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
        # This will not run if the keyword was the only thing in the directory
        # level
        if kw_suffix != '':
            # Find the next '*' from the left
            next_star_in_suffix = kw_suffix.find('*')
            # If there is another '*', grab non-wildcards up until '*' as
            # delim
            if next_star_in_suffix >= 0:
                suffix_delim = kw_suffix[:next_star_in_suffix]
            # Otherwise, just use the whole prefix as delim
            else:
                suffix_delim = kw_suffix
            # Check for glob-style characters in delimeter
            suffix_delim = check_for_glob_patterns(suffix_delim, key_str,
                                                   suffix_flag=True)

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
            msg = 'Could not distinguish %s from filepath %s using the ' \
                  'file pattern template %s.\nInstead, using entire ' \
                  'directory: %s for %s.\nCheck data organization and ' \
                  'file pattern template' \
                  % (keyword, filepath, template, key_str, keyword)
            logger.info(msg)
            key_str = fp_split[kw_idx]
    else:
        # logger.info('Keyword %s not found in template %s'
        #             % (keyword, template))
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
        site_dict[site] = {key.lower(): val for key, val in dict_row.items()\
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

    # Import packages
    import os
    import logging

    # Init variables
    logger = logging.getLogger('sublist_builder')

    # Check if {site} was specified
    keep_site_paths = []
    if site_kw in path_template and include_sites is not None:
        # Filter out sites that are not included
        if type(include_sites) is not list:
            include_sites = [include_sites]
        logger.info('Only including sites: %s' % include_sites)
        site_matches = filter(lambda sp: \
                              extract_keyword_from_path(sp,
                                                        site_kw,
                                                        path_template) in \
                              include_sites, sub_paths)
        # note which site IDs in "subjects to include" are missing
        missing = list(include_sites)
        for site_id in site_matches:
            for include in include_sites:
                if include in site_id:
                    if include in missing:
                        missing.remove(include)
        if len(missing) > 0:
            logger.info("Site IDs marked in 'Sites to Include' not found:" \
                        "\n%s" % str(missing))
        keep_site_paths.extend(site_matches)
    else:
        logger.info('Not filtering out any potential sites...')
        keep_site_paths = sub_paths

    # Only keep subjects in inclusion list or remove those in exclusion list
    if include_subs is not None and exclude_subs is not None:
        err_msg = 'Please only populate subjects to include or exclude '\
                  '- not both!'
        logger.error(err_msg)
        raise Exception(err_msg)
    # Include only
    elif include_subs is not None:
        keep_subj_paths = []
        if ".txt" in include_subs:
            if os.path.exists(include_subs):
                include_subs = read_subj_txtfile(include_subs)
            else:
                err = "\n\n[!] The filepath to the subject inclusion " \
                      "list does not exist!\nFilepath: %s\n\n" \
                      % include_subs
                raise Exception(err)
        if type(include_subs) is not list:
            include_subs = [include_subs]
        logger.info('Only including subjects: %s' % include_subs)
        subj_matches = filter(lambda sp: \
                              extract_keyword_from_path(sp,
                                                        ppant_kw,
                                                        path_template) in \
                              include_subs, sub_paths)
        # note which participant IDs in "subjects to include" are missing
        missing = list(include_subs)
        for subj_path in subj_matches:
            for include in include_subs:
                if include in subj_path:
                    if include in missing:
                        missing.remove(include)
        if len(missing) > 0:
            logger.info("Participant IDs marked in 'Subjects to Include' "\
                        "not found:\n%s" % str(missing))
        keep_subj_paths.extend(subj_matches)

        if keep_site_paths:
            keep_subj_paths = list(set(keep_subj_paths) & set(keep_site_paths))

    # Or exclude only
    elif exclude_subs is not None:
        keep_subj_paths = []
        if ".txt" in exclude_subs:
            if os.path.exists(exclude_subs):
                exclude_subs = read_subj_txtfile(exclude_subs)
            else:
                err = "\n\n[!] The filepath to the subject exclusion " \
                      "list does not exist!\nFilepath: %s\n\n" \
                      % exclude_subs
                raise Exception(err)
        if type(exclude_subs) is not list:
            exclude_subs = [exclude_subs]
        logger.info('Including all subjects but: %s' % exclude_subs)
        subj_matches = filter(lambda sp: \
                              extract_keyword_from_path(sp,
                                                        ppant_kw,
                                                        path_template) not in \
                              exclude_subs, sub_paths)
        keep_subj_paths.extend(subj_matches)

        if keep_site_paths:
            keep_subj_paths = list(set(keep_subj_paths) & set(keep_site_paths))

    else:
        keep_subj_paths = keep_site_paths

    # Filter out any duplicates
    keep_subj_paths = list(set(keep_subj_paths))
    logger.info('Filtered down to %d files' % len(keep_subj_paths))

    # Return kept paths
    return keep_subj_paths


# Get site, ppant, session-level directory indicies
def return_dir_indices(path_template):
    '''
    Function to return the site, participant, and session-level
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
    # Session level isn't required, but recommended
    sess_idx = None
    sess_extra = None
    for dir_level in fp_split:
        if "{site}" in dir_level:
            site_idx = fp_split.index(dir_level)
            site_extra = filter(bool,dir_level.split("{site}"))
        if "{participant}" in dir_level:
            ppant_idx = fp_split.index(dir_level)
            ppant_extra = filter(bool,dir_level.split("{participant}"))
        if "{session}" in dir_level:
            sess_idx = fp_split.index(dir_level)
            sess_extra = filter(bool,dir_level.split("{session}"))

    # Return extra characters
    return site_idx, ppant_idx, sess_idx, site_extra, ppant_extra, sess_extra


# Return matching filepaths
def return_local_matched_paths(path_template, bids_flag=False):
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
    matched_paths : list
        a list of strings of the local filepaths
    '''

    # Import packages
    import glob
    import logging
    import os

    # Get logger
    logger = logging.getLogger('sublist_builder')

    # Gather local files
    # if bids_flag:
    #     file_pattern = path_template
    # else:
    file_pattern = path_template.replace('{site}', '*').\
                   replace('{participant}', '*').\
                   replace('{session}', '*').\
                   replace('{series}', '*')

    # in case the user doesn't end their glob path with ".nii" or ".nii.gz"
    if ".nii" not in file_pattern:
        "".join([file_pattern, ".nii*"])

    local_filepaths = glob.glob(file_pattern)

    if len(local_filepaths) == 0:
        err = "\n\n[!] No data filepaths were found with the path template " \
              "given:\n%s\n\n" % path_template
        raise Exception(err)

    # Restrict filepaths and pattern to be of same directory depth
    # as fnmatch will expand /*/ recursively to .../*/*/...
    matched_paths = []
    for lfp in local_filepaths:
        s3_split = lfp.split('/')
        fp_split = file_pattern.split('/')
        if len(s3_split) == len(fp_split):
            matched_paths.append(lfp)

    # Get absolute paths
    matched_paths = [os.path.abspath(fp) for fp in matched_paths]

    # Print how many found
    num_local_files = len(matched_paths)
    logger.info('Found %d files!' % num_local_files)

    # Return the filepaths as a list
    return matched_paths


def return_local_filepaths(base_dir):

    import os

    print('Gathering file paths from directory:\n{0}\n'.format(base_dir))

    local_filepaths = []
    for root, dirs, files in os.walk(base_dir):
        for filename in files:
            fullpath = os.path.join(root, filename)
            local_filepaths.append(fullpath)

    print('{0} file paths found.\n'.format(str(len(local_filepaths))))

    return local_filepaths


def return_s3_filepaths(base_dir, creds_path=None):
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
    matched_s3_paths : list
        a list of strings of the filepaths from the S3 bucket
    '''

    # Import packages
    import logging
    import os

    from indi_aws import fetch_creds

    # # Check for errors
    # if not bids_base_dir:
    #     if not ('{site}' in path_template and '{participant}' in path_template):
    #         err_msg = 'Please provide \'{site}\' and \'{particpant}\' in '\
    #                   'filepath template where site and participant-level '\
    #                   'directories are present'
    #         raise Exception(err_msg)

    # if running this with "Custom" (non-BIDS) file templates
    if '{site}' in base_dir:
        base_dir = base_dir.split('{site}')[0]
    elif '{participant}' in base_dir:
        base_dir = base_dir.split('{participant}')[0]

    # Init variables
    bucket_name = base_dir.split('/')[2]
    s3_prefix = '/'.join(base_dir.split('/')[:3])

    # Get logger
    logger = logging.getLogger('sublist_builder')

    # Extract base prefix to search through in S3
    prefix = base_dir.replace(s3_prefix, '').lstrip('/')

    # Attempt to get bucket
    try:
        bucket = fetch_creds.return_bucket(creds_path, bucket_name)
    except Exception as exc:
        err_msg = 'There was an error in retrieving S3 bucket: %s.\n' \
                  'Error: %s' % (bucket_name, exc)
        logger.error(err_msg)
        raise Exception(err_msg)

    # Get filepaths from S3 with prefix
    logger.info('Gathering files from S3 to parse...')
    s3_filepaths = []
    for s3_obj in bucket.objects.filter(Prefix=prefix):
        s3_filepaths.append(str(s3_obj.key))

    # Prepend 's3://bucket_name/' on found paths
    s3_filepaths = [os.path.join(s3_prefix, s3_fp) for s3_fp in s3_filepaths]

    return s3_filepaths


def process_s3_paths(s3_filepaths, path_template):

    import logging
    import fnmatch

    # Get logger
    logger = logging.getLogger('sublist_builder')

    # File pattern filter
    file_pattern = path_template.replace('{site}', '*').replace('{participant}', '*').replace('{session}', '*').replace('{series}', '*')

    # Get only matching s3 paths
    s3_filepaths = fnmatch.filter(s3_filepaths, file_pattern)

    # Restrict filepaths and pattern to be of same directory depth
    # as fnmatch will expand /*/ recursively to .../*/*/...
    matched_s3_paths = []
    for s3fp in s3_filepaths:
        s3_split = s3fp.split('/')
        fp_split = file_pattern.split('/')
        if len(s3_split) == len(fp_split):
            matched_s3_paths.append(s3fp)

    # Print how many found
    num_s3_files = len(matched_s3_paths)
    logger.info('Found %d files!' % num_s3_files)

    # Return the filepaths as a list
    return matched_s3_paths


def parse_BIDS_filepaths(bids_filepaths, bids_base_dir, creds_path=None):

    sub_dict = {}
    json_dict = {}

    for path in bids_filepaths:
        # get the stuff after the base dir by itself
        dirs_path = path.split(bids_base_dir)[1].lstrip('/')
        dirs = dirs_path.split('/')
        filename = dirs[-1]

        if '.nii' in path:
            if '.nii' not in filename or '_' not in filename:
                continue
            # 'anat' or 'func'
            type = dirs[-2]

            if type != 'anat' and type != 'func':
                print '\nFilepath not in BIDS format! Skipping..\n' \
                      'Filepath: {0}\n'.format(path)
                continue

            sub_id = None
            ses_id = None
            for dir in dirs:
                if 'sub-' in dir and '.nii' not in dir:
                    sub_id = dir
                    sub_idx = dirs.index(dir)
                if 'ses-' in dir and '.nii' not in dir:
                    ses_id = dir

            if sub_idx and sub_idx != 0:
                site_id = dirs[sub_idx - 1]
            else:
                site_id = None

            # parse filename and get scan ID
            if type == 'func':
                scan_id = None
                run_id = None
                for file_part in filename.split('_'):
                    if 'task-' in file_part:
                        scan_id = file_part.split('task-')[1]
                    if 'run-' in file_part:
                        run_id = file_part
                if not scan_id:
                    continue
                if run_id:
                    scan_id = '_'.join([scan_id, run_id])

            # assign defaults if needed
            if not site_id:
                site_id = 'site-1'
            if not ses_id:
                ses_id = 'ses-1'

            key = ('{0}_{1}_{2}'.format(sub_id, site_id, ses_id))

            # populate sub_dict
            if key not in sub_dict.keys():
                if type == 'anat':
                    sub_dict[key] = {'anat': path,
                                     'subject_id': sub_id,
                                     'unique_id': ses_id,
                                     'site': site_id,
                                     'creds_path': creds_path}
                if type == 'func':
                    sub_dict[key] = {'func': {scan_id: path},
                                     'creds_path': creds_path}
            else:
                tmp_dict = sub_dict[key]

                if type == 'anat':
                    if 'anat' not in tmp_dict.keys():
                        tmp_dict.update({'anat': path})
                elif type == 'func':
                    if 'func' not in tmp_dict.keys():
                        tmp_dict['func'] = {scan_id: path}
                    else:
                        if scan_id not in tmp_dict['func']:
                            tmp_dict['func'].update({scan_id: path})

                sub_dict[key] = tmp_dict

        elif '.json' in path:
            site_id = None
            sub_id = None
            ses_id = None
            scan_id = None

            if len(dirs) > 2:
                # more specific json
                for dir in dirs:
                    if 'sub-' in dir and '.json' not in dir:
                        sub_id = dir
                        sub_idx = dirs.index(dir)
                    if 'ses-' in dir and '.json' not in dir:
                        ses_id = dir
                    if '.json' in dir:
                        scan_id = dir.replace('.json', '')
            else:
                # site-level json
                if len(dirs) > 1:
                    if '.json' in dirs[1]:
                        site_id = dirs[0]
                        if 'task-' in dirs[1]:
                            scan_id = dirs[1].replace('.json', '')
                else:
                    site_id = 'site-1'
                    if '.json' in dirs[0]:
                        scan_id = dirs[0].replace('.json', '')

            json_dict[(site_id, sub_id, ses_id, scan_id)] = path

    # grab the scan parameter JSONs, if they are there
    for id_tuple in json_dict.keys():
        json_path = json_dict[id_tuple]
        site_id = id_tuple[0]
        sub_id = id_tuple[1]
        ses_id = id_tuple[2]
        scan_id = id_tuple[3]

        if scan_id:
            if 'task-' not in scan_id:
                continue
        else:
            continue

        if site_id and not sub_id and not ses_id:
            # site-level json
            for id_string in sub_dict.keys():
                if site_id in id_string:
                    tmp_dict = sub_dict[id_string]
                    if 'func' not in tmp_dict.keys():
                        tmp_dict['func'] = {'scan_parameters': json_path}
                    else:
                        if scan_id not in tmp_dict['func']:
                            tmp_dict['func'].update({'scan_parameters': json_path})
                    sub_dict[id_string] = tmp_dict

        if sub_id:
            # more specific
            for id_string in sub_dict.keys():
                if sub_id in id_string:
                    if ses_id:
                        if ses_id in id_string:
                            tmp_dict = sub_dict[id_string]
                    else:
                        tmp_dict = sub_dict[id_string]
                    if 'func' not in tmp_dict.keys():
                        tmp_dict['func'] = {'scan_parameters': json_path}
                    else:
                        if scan_id not in tmp_dict['func']:
                            tmp_dict['func'].update({'scan_parameters': json_path})
                    sub_dict[id_string] = tmp_dict

    return sub_dict


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
    import logging
    import os
    import yaml

    from CPAC.utils import bids_metadata
    from CPAC.utils.utils import setup_logger

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

    # Get C-PAC logger if it's available
    log_path = os.path.join(sublist_outdir,
                            'sublist_build_%s.log' % \
                            os.path.basename(data_config_yml).split('.')[0])
    logger = setup_logger('sublist_builder', log_path, logging.INFO, to_screen=True)

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

    '''
    # Get dataformat if BIDS or not
    if data_format == 'BIDS':
        bids_flag = True
        anat_template = return_bids_template(bids_base_dir, 'anat')
        func_template = return_bids_template(bids_base_dir, 'func')
    else:
        bids_flag = False

    ### WE NEED A TEMPLATE FOR RETURN S3 PATHS???

    # Get directory indices
    # If data is BIDS
    if bids_flag:
        anat_site_idx = func_site_idx = None
        anat_path_dirs = anat_template.split('/')

        anat_ppant_idx = func_ppant_idx = \
            anat_path_dirs.index('{participant}')

        if '{site}' in anat_path_dirs:
            anat_site_idx = func_site_idx = anat_path_dirs.index('{site}')
        else:
            anat_site_idx = func_site_idx = None

        if '{session}' in anat_path_dirs:
            anat_sess_idx = func_sess_idx = anat_path_dirs.index('{session}')
        else:
            anat_sess_idx = func_sess_idx = None
    '''

    # Read in scan parameters and return site-based dictionary
    if scan_params_csv is not None:
        site_scan_params = extract_scan_params(scan_params_csv)
    else:
        site_scan_params = {}

    # Iterate through file paths and build subject list
    if data_format == 'BIDS':
        if 's3://' in bids_base_dir:
            # on the AWS S3 bucket
            bids_filepaths = return_s3_filepaths(bids_base_dir, creds_path)
        else:
            # local
            bids_filepaths = return_local_filepaths(bids_base_dir)

        sub_dict = parse_BIDS_filepaths(bids_filepaths, bids_base_dir,
                                        creds_path)

    elif data_format == 'Custom':
        # See if the templates are s3 files
        if s3_str in anat_template and s3_str in func_template:
            # Get anatomical filepaths from s3
            print 'Fetching anatomical files...'
            anat_s3_paths = return_s3_filepaths(anat_template, creds_path)
            anat_paths = process_s3_paths(anat_s3_paths, anat_template)
            # Get functional filepaths from s3
            print 'Fetching functional files...'
            func_s3_paths = return_s3_filepaths(func_template, creds_path)
            func_paths = process_s3_paths(func_s3_paths, func_template)

        # If one is in S3 and the other is not, raise error - not supported
        elif (
                s3_str in anat_template.lower() and s3_str not in func_template.lower()) or \
                (
                        s3_str not in anat_template.lower() and s3_str in func_template.lower()):
            err_msg = 'Both anatomical and functional files should either be ' \
                      'on S3 or local. Separating the files is currently not ' \
                      'supported.'
            raise Exception(err_msg)

        # Otherwise, it's local files
        else:
            # Get anatomical filepaths
            print 'Gathering anatomical files...'
            anat_paths = return_local_matched_paths(anat_template)
            # Get functional filepaths
            print 'Gathering functional files...'
            func_paths = return_local_matched_paths(func_template)

        #TODO: this might be in the wrong spot!!!
        # Filter out unwanted anat and func filepaths
        logger.info('Filtering anatomical files...')
        anat_paths = filter_sub_paths(anat_paths, include_sites,
                                      include_subs, exclude_subs,
                                      site_kw, ppant_kw, anat_template)
        print 'Filtered down to %d anatomical files' % len(anat_paths)

        func_paths = filter_sub_paths(func_paths, include_sites,
                                      include_subs, exclude_subs,
                                      site_kw, ppant_kw, func_template)

        # If all data is filtered out, raise exception
        if len(anat_paths) == 0 or len(func_paths) == 0:
            err_msg = 'Unable to find any files after filtering sites and ' \
                      'subjects! Check site and subject inclusion fields as ' \
                      'well as filepaths and template pattern!'
            logger.error(err_msg)
            raise Exception(err_msg)

        for anat in anat_paths:
            subj = extract_keyword_from_path(anat, "{participant}",
                                             anat_template)
            try:
                sess = extract_keyword_from_path(anat, "{session}",
                                                 anat_template)
            except TypeError:
                sess = "ses-1"
            site = extract_keyword_from_path(anat, "{site}", anat_template)

            subj_d = {'anat': anat, 'creds_path': creds_path, 'func': {},
                      'subject_id': subj, 'unique_id': sess}
            tmp_key = '_'.join([subj, site, sess])
            tmp_dict[tmp_key] = subj_d

        # Now go through and populate functional scans dictionaries
        for func in func_paths:
            func_sp = func.split('/')
            subj = extract_keyword_from_path(func, "{participant}",
                                             func_template)
            try:
                sess = extract_keyword_from_path(func, "{session}",
                                                 func_template)
            except TypeError:
                sess = "ses-1"
            site = extract_keyword_from_path(func, "{site}", func_template)
            scan_params = None
            if scan_params_csv is not None:
                try:
                    scan_params = site_scan_params[site]
                except KeyError as exc:
                    print 'Site %s missing from scan parameters csv, ' \
                          'skipping...' % site

            # Build tmp key and get subject dictionary from tmp dictionary
            tmp_key = '_'.join([subj, site, sess])
            # Try and find the associated anat scan
            try:
                subj_d = tmp_dict[tmp_key]
            except KeyError as exc:
                logger.info('Unable to find anatomical image for %s. '
                            'Skipping...' % tmp_key)
                continue

            # If there is no scan sub-folder, make scan
            # the name of the image itself without extension
            if "{series}" in func_template:
                scan = extract_keyword_from_path(func, "{series}",
                                                 func_template)
            else:
                scan = func_sp[-1].split('.nii')[0]

            # Set the rest dictionary with the scan
            subj_d['func'][scan] = func
            subj_d['func']["scan_parameters"] = scan_params
            # And replace it back in the dictionary
            tmp_dict[tmp_key] = subj_d

        sub_dict = tmp_dict

    # all done
    # Build a subject list from dictionary values
    sublist = []
    for data_bundle in sub_dict.values():
        if 'anat' not in data_bundle.keys() or \
                        'func' not in data_bundle.keys():
            continue
        if data_bundle['anat'] != '' and data_bundle['func'] != {}:
            sublist.append(data_bundle)
    # Check to make sure subject list has at least one valid data bundle
    if len(sublist) == 0:
        err_msg = 'Unable to find both anatomical and functional data for ' \
                  'any subjects. Check data organization and file patterns!'
        logger.error(err_msg)
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

    print('Participant list created - {0}'.format(sublist_out_yml))

    # Return the subject list
    return sublist
