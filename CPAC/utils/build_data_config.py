


# this will go into core tools
def gather_file_paths(base_directory, verbose=False):

    # ideas: return number of paths, optionally instead
    #        or write out, optionally, a text file with all paths, for easy
    #        searching

    # test cases:
    #   - proper data directory
    #   - empty directory

    import os

    path_list = []

    for root, dirs, files in os.walk(base_directory):
        for path in files:
            fullpath = os.path.join(root, path)
            path_list.append(fullpath)

    if verbose:
        print "Number of paths: {0}".format(len(path_list))

    return path_list


def extract_scan_params(scan_params_csv):
    """
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
    """

    # Import packages
    import csv

    # Init variables
    csv_open = open(scan_params_csv, 'r')
    site_dict = {}

    # Init csv dictionary reader
    reader = csv.DictReader(csv_open)

    placeholders = ['None', 'NONE', 'none', 'All', 'ALL', 'all', '', ' ']

    # Iterate through the csv and pull in parameters
    for dict_row in reader:

        if dict_row['Site'] in placeholders:
            site = 'All'
        else:
            site = dict_row['Site']

        ses = 'All'
        if 'Session' in dict_row.keys():
            if dict_row['Session'] not in placeholders:
                ses = dict_row['Session']

        if ses != 'All':
            # for scan-specific scan parameters (less common)
            if site in site_dict.keys():
                site_dict[site][ses] = {key.lower(): val
                                        for key, val in dict_row.items()
                                        if key != 'Site' and key != 'Session'}
            else:
                site_dict[site] = {ses: {key.lower(): val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and key != 'Session'}}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            site_dict[site][ses]['tr'] = \
                site_dict[site][ses].pop('tr (seconds)')

        else:
            # site-specific scan parameters only (more common)
            if site not in site_dict.keys():
                site_dict[site] = \
                    {ses: {key.lower(): val for key, val in dict_row.items()
                           if key != 'Site' and key != 'Session'}}
            else:
                site_dict[site][ses] = \
                {key.lower(): val for key, val in dict_row.items()
                 if key != 'Site' and key != 'Session'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            site_dict[site][ses]['tr'] = \
                site_dict[site][ses].pop('tr (seconds)')

    # Return site dictionary
    return site_dict


def format_incl_excl_dct(site_incl_list=None, participant_incl_list=None,
                         session_incl_list=None, scan_incl_list=None):

    incl_dct = {}

    if isinstance(site_incl_list, str):
        if '.txt' in site_incl_list:
            with open(site_incl_list, 'r') as f:
                incl_dct['sites'] = [x for x in f.readlines() if x != '']
    elif isinstance(site_incl_list, list):
        incl_dct['sites'] = site_incl_list

    if isinstance(participant_incl_list, str):
        if '.txt' in participant_incl_list:
            with open(participant_incl_list, 'r') as f:
                incl_dct['participants'] = \
                    [x.rstrip('\n') for x in f.readlines() if x != '']
    elif isinstance(participant_incl_list, list):
        incl_dct['participants'] = participant_incl_list

    if isinstance(session_incl_list, str):
        if '.txt' in session_incl_list:
            with open(session_incl_list, 'r') as f:
                incl_dct['sessions'] = [x for x in f.readlines() if x != '']
    elif isinstance(session_incl_list, list):
        incl_dct['sessions'] = session_incl_list

    if isinstance(scan_incl_list, str):
        if '.txt' in scan_incl_list:
            with open(scan_incl_list, 'r') as f:
                incl_dct['scans'] = [x for x in f.readlines() if x != '']
    elif isinstance(scan_incl_list, list):
        incl_dct['scans'] = scan_incl_list

    return incl_dct


def get_nonBIDS_data(anat_template, func_template, scan_params_dct=None,
                     aws_creds_path=None, inclusion_dct=None,
                     exclusion_dct=None):

    # go over the file paths, validate for nifti's?
    # work with the template

    # TODO: handle both json and csv scan parameters

    # test cases:
    #   - no anats, no funcs
    #   - combination of .nii and .nii.gz?
    #   - throw error/warning if anat and func templates are identical
    #   - all permutations of scan parameters json/csv's at different levels

    import glob

    # should have the {participant} label at the very least
    if '{participant}' not in anat_template or \
            '{participant}' not in func_template:
        print "ERROR"

    keywords = ['{site}', '{participant}', '{session}', '{scan}']

    # make globby templates, to use them to filter down the path_list into
    # only paths that will work with the templates
    anat_glob = anat_template
    func_glob = func_template

    for keyword in keywords:
        if keyword in anat_glob:
            anat_glob = anat_glob.replace(keyword, '*')
        if keyword in func_glob:
            func_glob = func_glob.replace(keyword, '*')

    # presumably, the paths contained in each of these pools should be anat
    # and func files only, respectively, if the templates were set up properly
    anat_pool = glob.glob(anat_glob)
    func_pool = glob.glob(func_glob)

    # okay, a few cases here:
    #   - one participant label alone in directory
    #   - one participant label in a string with other things
    #   - multiple participant labels

    data_dct = {}

    # pull out the site/participant/etc. IDs from each path and connect them
    for anat_path in anat_pool:

        # NIFTIs only
        if '.nii' not in anat_path:
            continue

        path_dct = {}

        site_parts = anat_template.split('{site}')

        partic_parts = []
        for part in site_parts:
            partic_parts = partic_parts + part.split('{participant}')
        ses_parts = []
        for part in partic_parts:
            ses_parts = ses_parts + part.split('{session}')

        new_template = anat_template
        new_path = anat_path

        for idx in range(0, len(ses_parts)):
            part1 = ses_parts[idx]
            try:
                part2 = ses_parts[idx+1]
            except IndexError:
                break

            label = new_template.split(part1, 1)[1]
            label = label.split(part2, 1)[0]

            id = new_path.split(part1, 1)[1]
            id = id.split(part2, 1)[0]

            if label not in path_dct.keys():
                path_dct[label] = id
            else:
                if path_dct[label] != id:
                    print "warning"

            new_template = new_template.replace(part1, '', 1)
            new_template = new_template.replace(label, '', 1)

            new_path = new_path.replace(part1, '', 1)
            new_path = new_path.replace(id, '', 1)

        sub_id = path_dct['{participant}']

        if '{site}' in path_dct.keys():
            site_id = path_dct['{site}']
        else:
            site_id = 'site-1'

        if '{session}' in path_dct.keys():
            ses_id = path_dct['{session}']
        else:
            ses_id = 'ses-1'
            
        if inclusion_dct:
            if 'sites' in inclusion_dct.keys():
                if site_id not in inclusion_dct['sites']:
                    continue
            if 'sessions' in inclusion_dct.keys():
                if ses_id not in inclusion_dct['sessions']:
                    continue
            if 'participants' in inclusion_dct.keys():
                if sub_id not in inclusion_dct['participants']:
                    continue
                    
        if exclusion_dct:
            if 'sites' in exclusion_dct.keys():
                if site_id in exclusion_dct['sites']:
                    continue
            if 'sessions' in exclusion_dct.keys():
                if ses_id in exclusion_dct['sessions']:
                    continue
            if 'participants' in exclusion_dct.keys():
                if sub_id in exclusion_dct['participants']:
                    continue

        temp_sub_dct = {'subject_id': sub_id,
                        'unique_id': ses_id,
                        'site': site_id,
                        'anat': anat_path,
                        'creds_path': aws_creds_path}

        if site_id not in data_dct.keys():
            data_dct[site_id] = {}
        if sub_id not in data_dct[site_id].keys():
            data_dct[site_id][sub_id] = {}
        if ses_id not in data_dct[site_id][sub_id].keys():
            data_dct[site_id][sub_id][ses_id] = temp_sub_dct
        else:
            #TODO: THROW ERROR/WARNING!!!
            pass

    # functional time
    for func_path in func_pool:

        # NIFTIs only
        if '.nii' not in func_path:
            continue

        path_dct = {}

        site_parts = func_template.split('{site}')

        partic_parts = []
        for part in site_parts:
            partic_parts = partic_parts + part.split('{participant}')
        ses_parts = []
        for part in partic_parts:
            ses_parts = ses_parts + part.split('{session}')
        scan_parts = []
        for part in ses_parts:
            scan_parts = scan_parts + part.split('{scan}')

        new_template = func_template
        new_path = func_path

        for idx in range(0, len(scan_parts)):
            part1 = scan_parts[idx]
            try:
                part2 = scan_parts[idx+1]
            except IndexError:
                break

            label = new_template.split(part1, 1)[1]
            label = label.split(part2, 1)[0]

            id = new_path.split(part1, 1)[1]
            id = id.split(part2, 1)[0]

            if label not in path_dct.keys():
                path_dct[label] = id
            else:
                if path_dct[label] != id:
                    print "warning"

            new_template = new_template.replace(part1, '', 1)
            new_template = new_template.replace(label, '', 1)

            new_path = new_path.replace(part1, '', 1)
            new_path = new_path.replace(id, '', 1)

        sub_id = path_dct['{participant}']

        if '{site}' in path_dct.keys():
            site_id = path_dct['{site}']
        else:
            site_id = 'site-1'

        if '{session}' in path_dct.keys():
            ses_id = path_dct['{session}']
        else:
            ses_id = 'ses-1'

        if '{scan}' in path_dct.keys():
            scan_id = path_dct['{scan}']
        else:
            scan_id = 'func-1'

        if inclusion_dct:
            if 'sites' in inclusion_dct.keys():
                if site_id not in inclusion_dct['sites']:
                    continue
            if 'sessions' in inclusion_dct.keys():
                if ses_id not in inclusion_dct['sessions']:
                    continue
            if 'participants' in inclusion_dct.keys():
                if sub_id not in inclusion_dct['participants']:
                    continue
            if 'scans' in inclusion_dct.keys():
                if scan_id not in inclusion_dct['scans']:
                    continue

        if exclusion_dct:
            if 'sites' in exclusion_dct.keys():
                if site_id in exclusion_dct['sites']:
                    continue
            if 'sessions' in exclusion_dct.keys():
                if ses_id in exclusion_dct['sessions']:
                    continue
            if 'participants' in exclusion_dct.keys():
                if sub_id in exclusion_dct['participants']:
                    continue
            if 'scans' in exclusion_dct.keys():
                if scan_id in exclusion_dct['scans']:
                    continue

        temp_func_dct = {scan_id: func_path}

        scan_params = None
        if scan_params_dct:
            if site_id in scan_params_dct.keys():
                if ses_id in scan_params_dct[site_id].keys():
                    # site-specific and session-specific scan params
                    scan_params = scan_params_dct[site_id][ses_id]
                elif 'All' in scan_params_dct[site_id].keys():
                    # site-specific scan params, same across all sessions
                    scan_params = scan_params_dct[site_id]['All']
            elif 'All' in scan_params_dct.keys():
                if ses_id in scan_params_dct['All'].keys():
                    # session-specific scan params, same across all sites
                    scan_params = scan_params_dct['All'][ses_id]
                elif 'All' in scan_params_dct['All'].keys():
                    # same scan params across all sites and sessions
                    scan_params = scan_params_dct['All']['All']

        if scan_params:
            temp_func_dct.update({'scan_parameters': scan_params})

        if site_id not in data_dct.keys():
            print "error"
        if sub_id not in data_dct[site_id].keys():
            print "error"
        if ses_id not in data_dct[site_id][sub_id].keys():
            print "error"

        if 'func' not in data_dct[site_id][sub_id][ses_id].keys():
            data_dct[site_id][sub_id][ses_id]['func'] = temp_func_dct
        else:
            data_dct[site_id][sub_id][ses_id]['func'].update(temp_func_dct)

    return data_dct


def run(data_settings_yml):

    import os
    import yaml

    print "Generating data configuration file.."

    with open(data_settings_yml, "r") as f:
        settings_dct = yaml.load(f)

    print "settings: ", settings_dct

    # local (not on AWS S3 bucket), non-BIDS files
    if 'Custom' in settings_dct['dataFormat'] or 'custom' in settings_dct['dataFormat']:

        params_dct = None
        if settings_dct['scanParametersCSV']:
            if '.csv' in settings_dct['scanParametersCSV']:
                params_dct = \
                    extract_scan_params(settings_dct['scanParametersCSV'])

        print params_dct

        incl_dct = format_incl_excl_dct(settings_dct['siteList'],
                                        settings_dct['subjectList'],
                                        settings_dct['sessionList'],
                                        settings_dct['scanList'])
        print incl_dct
        excl_dct = format_incl_excl_dct(settings_dct['exclusionSiteList'],
                                        settings_dct['exclusionSubjectList'],
                                        settings_dct['exclusionSessionList'],
                                        settings_dct['exclusionScanList'])
        print excl_dct
        data_dct = get_nonBIDS_data(settings_dct['anatomicalTemplate'],
                                    settings_dct['functionalTemplate'],
                                    params_dct,
                                    settings_dct['awsCredentialsFile'],
                                    incl_dct, excl_dct)
        print data_dct
        # get some data
        num_sites = len(data_dct.keys())
        num_subs = num_sess = num_scan = 0
        for site in data_dct.keys():
            num_subs += len(data_dct[site])
            for sub in data_dct[site].keys():
                num_sess += len(data_dct[site][sub])
                for session in data_dct[site][sub].keys():
                    for scan in data_dct[site][sub][session]['func'].keys():
                        if scan != "scan_parameters":
                            num_scan += 1

        if len(data_dct) > 0:
            data_config_outfile = \
                os.path.join(settings_dct['outputSubjectListLocation'],
                             "data_config_{0}.yml".format(settings_dct['subjectListName']))

            #TODO: WRITE IT OUT IN CPAC SUBLIST FORMAT!!!!
            #TODO: also what's with the &id001 thing in scan params

            with open(data_config_outfile, "wt") as f:
                yaml.dump(data_dct, f)

            if os.path.exists(data_config_outfile):
                print "\nCPAC DATA SETTINGS file entered:" \
                      "\n{0}".format(data_settings_yml)
                print "\nCPAC DATA CONFIGURATION file created:" \
                      "\n{0}\n".format(data_config_outfile)
                print "Number of:"
                print "...sites: {0}".format(num_sites)
                print "...participants: {0}".format(num_subs)
                print "...participant-sessions: {0}".format(num_sess)
                print "...functional scans: {0}".format(num_scan)




