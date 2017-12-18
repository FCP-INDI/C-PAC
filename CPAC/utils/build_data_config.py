

def gather_file_paths(base_directory, verbose=False):

    # this will go into core tools eventually

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


def pull_s3_sublist(data_folder, creds_path=None):

    import os
    from indi_aws import fetch_creds

    if creds_path:
        creds_path = os.path.abspath(creds_path)

    s3_path = data_folder.split("s3://")[1]
    bucket_name = s3_path.split("/")[0]
    bucket_prefix = s3_path.split(bucket_name + "/")[1]

    s3_list = []
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    # ensure slash at end of bucket_prefix, so that if the final
    # directory name is a substring in other directory names, these
    # other directories will not be pulled into the file list
    if "/" not in bucket_prefix[-1]:
        bucket_prefix += "/"

    # Build S3-subjects to download
    for bk in bucket.objects.filter(Prefix=bucket_prefix):
        s3_list.append(str(bk.key).replace(bucket_prefix,""))

    return s3_list


def extract_scan_params_csv(scan_params_csv):
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

        sub = "All"
        if "Participant" in dict_row.keys():
            if dict_row["Participant"] not in placeholders:
                sub = dict_row["Participant"]

        ses = 'All'
        if 'Session' in dict_row.keys():
            if dict_row['Session'] not in placeholders:
                ses = dict_row['Session']

        if ses != 'All':
            # for session-specific scan parameters
            if site not in site_dict.keys():
                site_dict[site] = {}
            if sub not in site_dict[site].keys():
                site_dict[site][sub] = {}

            site_dict[site][sub][ses] = {key.lower(): val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and
                                         key != 'Participant' and
                                         key != 'Session'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            site_dict[site][sub][ses]['tr'] = \
                site_dict[site][sub][ses].pop('tr (seconds)')

        elif sub != "All":
            # participant-specific scan parameters
            if site not in site_dict.keys():
                site_dict[site] = {}
            if sub not in site_dict[site].keys():
                site_dict[site][sub] = {}

            site_dict[site][sub][ses] = {key.lower(): val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and
                                         key != 'Participant' and
                                         key != 'Session'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            site_dict[site][sub][ses]['tr'] = site_dict[site][sub][ses].pop('tr (seconds)')

        else:
            # site-specific scan parameters only
            if site not in site_dict.keys():
                site_dict[site] = {}
            if sub not in site_dict[site].keys():
                site_dict[site][sub] = {}

            site_dict[site][sub][ses] = {key.lower(): val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and
                                         key != 'Participant' and
                                         key != 'Session'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            site_dict[site][sub][ses]['tr'] = \
                site_dict[site][sub][ses].pop('tr (seconds)')

    # Return site dictionary
    return site_dict


def extract_scan_params_json(scan_params_json):

    import json

    with open(scan_params_json, "r") as f:
        params_dct = json.load(f)

    return params_dct


def format_incl_excl_dct(site_incl_list=None, participant_incl_list=None,
                         session_incl_list=None, scan_incl_list=None):

    incl_dct = {}

    if isinstance(site_incl_list, str):
        if '.txt' in site_incl_list:
            with open(site_incl_list, 'r') as f:
                incl_dct['sites'] = [x.rstrip("\n") for x in f.readlines() if x != '']
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
                incl_dct['sessions'] = [x.rstrip("\n") for x in f.readlines() if x != '']
    elif isinstance(session_incl_list, list):
        incl_dct['sessions'] = session_incl_list

    if isinstance(scan_incl_list, str):
        if '.txt' in scan_incl_list:
            with open(scan_incl_list, 'r') as f:
                incl_dct['scans'] = [x.rstrip("\n") for x in f.readlines() if x != '']
    elif isinstance(scan_incl_list, list):
        incl_dct['scans'] = scan_incl_list

    return incl_dct


def get_BIDS_data_dct(bids_base_dir, aws_creds_path=None, inclusion_dct=None,
                      exclusion_dct=None):

    import os
    import glob

    anat_sess = os.path.join(bids_base_dir,
                             "sub-{participant}/ses-{session}/anat/sub-"
                             "{participant}_ses-{session}_T1w.nii.gz")
    anat = os.path.join(bids_base_dir,
                        "sub-{participant}/anat/sub-{participant}_ses-"
                        "{session}_T1w.nii.gz")

    func_sess = os.path.join(bids_base_dir,
                             "sub-{participant}/ses-{session}/func/sub-"
                             "{participant}_ses-{session}_task-{scan}_"
                             "bold.nii.gz")
    func = os.path.join(bids_base_dir,
                        "sub-{participant}/func/sub-{participant}_task-"
                        "{scan}_bold.nii.gz")

    fmap_phase_sess = os.path.join(bids_base_dir,
                                   "sub-{participant}/ses-{session}/fmap/"
                                   "sub-{participant}_ses-{session}_phase"
                                   "diff.nii.gz")
    fmap_phase = os.path.join(bids_base_dir,
                              "sub-{participant}/fmap/sub-{participant}"
                              "_phasediff.nii.gz")

    fmap_mag_sess = os.path.join(bids_base_dir,
                                 "sub-{participant}/ses-{session}/fmap/"
                                 "sub-{participant}_ses-{session}_"
                                 "magnitude1.nii.gz")
    fmap_mag = os.path.join(bids_base_dir,
                            "sub-{participant}/fmap/sub-{participant}"
                            "_magnitude1.nii.gz")

    sess_glob = os.path.join(bids_base_dir, "sub-*/ses-*/*")

    site_jsons = glob.glob(os.path.join(bids_base_dir, "*.json"))
    sub_jsons = glob.glob(os.path.join(bids_base_dir, "sub-*/*.json"))
    ses_jsons = glob.glob(os.path.join(bids_base_dir, "sub-*/ses-*/*.json"))

    # setting site level to "All" because in BIDS datasets, it is either
    # assumed that the dataset is from one site, or that site information is
    # encoded in a separate file
    # TODO: if site information is ever encoded into the JSON scan parameter
    # TODO: file name, this must be implemented here

    scan_params_dct = {}

    if len(ses_jsons) > 0:
        for json_file in ses_jsons:
            ids = os.path.basename(json_file).rstrip(".json").split("_")
            for id in ids:
                if "sub-" in id:
                    sub_id = id.replace("sub-", "")
                if "ses-" in id:
                    ses_id = id.replace("ses-", "")

            if "All" not in scan_params_dct.keys():
                scan_params_dct["All"] = {}
            if sub_id not in scan_params_dct["All"].keys():
                scan_params_dct["All"][sub_id] = {}

            scan_params_dct["All"][sub_id][ses_id] = json_file

    if len(sub_jsons) > 0:
        for json_file in sub_jsons:
            ids = os.path.basename(json_file).rstrip(".json").split("_")
            for id in ids:
                if "sub-" in id:
                    sub_id = id.replace("sub-", "")

            if "All" not in scan_params_dct.keys():
                scan_params_dct["All"] = {}
            if sub_id not in scan_params_dct["All"].keys():
                scan_params_dct["All"][sub_id] = {}

            scan_params_dct["All"][sub_id]["All"] = json_file

    if len(site_jsons) > 0:
        if len(site_jsons) > 1:
            # TODO: find out site-level naming convention for scan param JSONs
            # TODO: improve this message (if it is necessary)
            print "More than one dataset-level JSON file!!!"
        else:
            if "All" not in scan_params_dct.keys():
                scan_params_dct["All"] = {}
            if "All" not in scan_params_dct["All"].keys():
                scan_params_dct["All"]["All"] = {}

            for json_file in site_jsons:
                scan_params_dct["All"]["All"]["All"] = json_file

    if len(glob.glob(sess_glob)) > 0:
        # if there is a session level in the BIDS dataset
        data_dct = get_nonBIDS_data(anat_sess, func_sess, scan_params_dct,
                                    fmap_phase_template=fmap_phase_sess,
                                    fmap_mag_template=fmap_mag_sess,
                                    aws_creds_path=aws_creds_path,
                                    inclusion_dct=inclusion_dct,
                                    exclusion_dct=exclusion_dct)
    else:
        # no session level
        data_dct = get_nonBIDS_data(anat, func, scan_params_dct,
                                    fmap_phase_template=fmap_phase,
                                    fmap_mag_template=fmap_mag,
                                    aws_creds_path=aws_creds_path,
                                    inclusion_dct=inclusion_dct,
                                    exclusion_dct=exclusion_dct)

    return data_dct


def get_nonBIDS_data(anat_template, func_template, scan_params_dct=None,
                     fmap_phase_template=None, fmap_mag_template=None,
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

    import os
    import glob

    # should have the {participant} label at the very least
    if '{participant}' not in anat_template:
        err = "\n[!] The {participant} keyword is missing from your " \
              "anatomical path template.\n"
        raise Exception(err)
    if '{participant}' not in func_template:
        err = "\n[!] The {participant} keyword is missing from your " \
              "functional path template.\n"
        raise Exception(err)

    keywords = ['{site}', '{participant}', '{session}', '{scan}']

    # make globby templates, to use them to filter down the path_list into
    # only paths that will work with the templates
    anat_glob = anat_template
    func_glob = func_template

    # backwards compatibility
    if '{series}' in anat_glob:
        anat_glob = anat_glob.replace('{series}', '{scan}')
    if '{series}' in func_glob:
        func_glob = func_glob.replace('{series}', '{scan}')

    if '{scan}' in anat_glob:
        err = "\n[!] CPAC does not support multiple anatomical scans at " \
              "this time. You are seeing this message because you have a " \
              "{scan} or {series} keyword in your anatomical path template.\n"
        raise Exception(err)

    for keyword in keywords:
        if keyword in anat_glob:
            anat_glob = anat_glob.replace(keyword, '*')
        if keyword in func_glob:
            func_glob = func_glob.replace(keyword, '*')

    # presumably, the paths contained in each of these pools should be anat
    # and func files only, respectively, if the templates were set up properly
    anat_pool = glob.glob(anat_glob)
    func_pool = glob.glob(func_glob)

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
                skip = False
            else:
                if path_dct[label] != id:
                    warn = "\n\n[!] WARNING: While parsing your input data " \
                           "files, a file path was found with conflicting " \
                           "IDs for the same data level.\n\n" \
                           "File path: {0}\n" \
                           "Level: {1}\n" \
                           "Conflicting IDs: {2}, {3}\n\n" \
                           "This file has not been added to the data " \
                           "configuration.".format(anat_path, label,
                                                   path_dct[label], id)
                    print warn
                    skip = True
                    break

            new_template = new_template.replace(part1, '', 1)
            new_template = new_template.replace(label, '', 1)

            new_path = new_path.replace(part1, '', 1)
            new_path = new_path.replace(id, '', 1)

        if skip:
            continue

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
            #TODO: finish this
            warn = "\n\n[!] WARNING: Duplicate site-participant-session " \
                   "entry found in your input data directory!\n\n"
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
                skip = False
            else:
                if path_dct[label] != id:
                    warn = "\n\n[!] WARNING: While parsing your input data " \
                           "files, a file path was found with conflicting " \
                           "IDs for the same data level.\n\n" \
                           "File path: {0}\n" \
                           "Level: {1}\n" \
                           "Conflicting IDs: {2}, {3}\n\n" \
                           "This file has not been added to the data " \
                           "configuration.".format(func_path, label,
                                                   path_dct[label], id)
                    print warn
                    skip = True
                    break

            new_template = new_template.replace(part1, '', 1)
            new_template = new_template.replace(label, '', 1)

            new_path = new_path.replace(part1, '', 1)
            new_path = new_path.replace(id, '', 1)

        if skip:
            continue

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

        # scan parameters time
        scan_params = None

        if scan_params_dct:

            #TODO: implement scan-specific level once scan-level nesting is
            #TODO: available

            if site_id in scan_params_dct.keys():
                if sub_id in scan_params_dct[site_id].keys():
                    if ses_id in scan_params_dct[site_id][sub_id].keys():
                        # site, sub, session specific scan params
                        scan_params = scan_params_dct[site_id][sub_id][ses_id]
                    elif 'All' in scan_params_dct[site_id][sub_id].keys():
                        # site and sub specific scan params, same across all
                        # sessions
                        scan_params = scan_params_dct[site_id][sub_id]['All']
                    elif "All" in scan_params_dct[site_id].keys():
                        if ses_id in scan_params_dct[site_id]["All"].keys():
                            # site and session specific scan params, same
                            # across all subs
                            scan_params = scan_params_dct[site_id]["All"][ses_id]
                        elif "All" in scan_params_dct[site_id]["All"].keys():
                            # site specific scan params, same across all subs
                            # and sessions
                            scan_params = scan_params_dct[site_id]["All"]["All"]

            elif "All" in scan_params_dct.keys():
                if sub_id in scan_params_dct["All"].keys():
                    if ses_id in scan_params_dct["All"][sub_id].keys():
                        # sub and session specific scan params, same across
                        # all sites
                        scan_params = scan_params_dct["All"][sub_id][ses_id]
                    elif "All" in scan_params_dct["All"][sub_id].keys():
                        # sub specific scan params, same across all sites and
                        # sessions
                        scan_params = scan_params_dct["All"][sub_id]["All"]
                elif "All" in scan_params_dct["All"].keys():
                    if ses_id in scan_params_dct["All"]["All"].keys():
                        # session-specific scan params, same across all sites
                        # and subs
                        scan_params = scan_params_dct["All"]["All"][ses_id]
                    elif "All" in scan_params_dct["All"]["All"].keys():
                        # same scan params across all sites and sessions
                        scan_params = scan_params_dct["All"]["All"]["All"]

        if scan_params:
            temp_func_dct.update({'scan_parameters': scan_params})

        #TODO: fill these
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

    if fmap_phase_template and fmap_mag_template:
        # if we're doing the whole field map distortion correction thing

        # make globby templates, to use them to filter down the path_list into
        # only paths that will work with the templates
        fmap_phase_glob = fmap_phase_template
        fmap_mag_glob = fmap_mag_template

        # backwards compatibility
        if '{series}' in fmap_phase_glob:
            fmap_phase_glob = fmap_phase_glob.replace('{series}', '{scan}')
        if '{series}' in fmap_mag_glob:
            fmap_mag_glob = fmap_mag_glob.replace('{series}', '{scan}')

        for keyword in keywords:
            if keyword in fmap_phase_glob:
                fmap_phase_glob = fmap_phase_glob.replace(keyword, '*')
            if keyword in fmap_mag_glob:
                fmap_mag_glob = fmap_mag_glob.replace(keyword, '*')

        # presumably, the paths contained in each of these pools should be
        # field map files only, if the templates were set up properly
        fmap_phase_pool = glob.glob(fmap_phase_glob)
        fmap_mag_pool = glob.glob(fmap_mag_glob)

        for fmap_phase in fmap_phase_pool:

            path_dct = {}

            site_parts = fmap_phase_template.split('{site}')

            partic_parts = []
            for part in site_parts:
                partic_parts = partic_parts + part.split('{participant}')
            ses_parts = []
            for part in partic_parts:
                ses_parts = ses_parts + part.split('{session}')
            scan_parts = []
            for part in ses_parts:
                scan_parts = scan_parts + part.split('{scan}')

            new_template = fmap_phase_template
            new_path = fmap_phase

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
                    skip = False
                else:
                    if path_dct[label] != id:
                        warn = "\n\n[!] WARNING: While parsing your input data " \
                               "files, a file path was found with conflicting " \
                               "IDs for the same data level.\n\n" \
                               "File path: {0}\n" \
                               "Level: {1}\n" \
                               "Conflicting IDs: {2}, {3}\n\n" \
                               "This file has not been added to the data " \
                               "configuration.".format(fmap_phase, label,
                                                       path_dct[label], id)
                        print warn
                        skip = True
                        break

                new_template = new_template.replace(part1, '', 1)
                new_template = new_template.replace(label, '', 1)

                new_path = new_path.replace(part1, '', 1)
                new_path = new_path.replace(id, '', 1)

            if skip:
                continue

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
                scan_id = None

            temp_fmap_dct = {"phase": fmap_phase}

            # TODO: fill these
            if site_id not in data_dct.keys():
                continue
            if sub_id not in data_dct[site_id].keys():
                continue
            if ses_id not in data_dct[site_id][sub_id].keys():
                continue

            if scan_id:
                # if the field map files are at scan level
                if 'fmap' not in data_dct[site_id][sub_id][ses_id][scan_id].keys():
                    data_dct[site_id][sub_id][ses_id][scan_id]['fmap'] = temp_fmap_dct
                else:
                    data_dct[site_id][sub_id][ses_id][scan_id]['fmap'].update(temp_fmap_dct)
            else:
                # if the field map fields are at session level
                if 'fmap' not in data_dct[site_id][sub_id][ses_id].keys():
                    data_dct[site_id][sub_id][ses_id]['fmap'] = temp_fmap_dct
                else:
                    data_dct[site_id][sub_id][ses_id]['fmap'].update(temp_fmap_dct)

        for fmap_mag in fmap_mag_pool:

            path_dct = {}

            site_parts = fmap_mag_template.split('{site}')

            partic_parts = []
            for part in site_parts:
                partic_parts = partic_parts + part.split('{participant}')
            ses_parts = []
            for part in partic_parts:
                ses_parts = ses_parts + part.split('{session}')
            scan_parts = []
            for part in ses_parts:
                scan_parts = scan_parts + part.split('{scan}')

            new_template = fmap_mag_template
            new_path = fmap_mag

            for idx in range(0, len(scan_parts)):
                part1 = scan_parts[idx]
                try:
                    part2 = scan_parts[idx + 1]
                except IndexError:
                    break

                label = new_template.split(part1, 1)[1]
                label = label.split(part2, 1)[0]

                id = new_path.split(part1, 1)[1]
                id = id.split(part2, 1)[0]

                if label not in path_dct.keys():
                    path_dct[label] = id
                    skip = False
                else:
                    if path_dct[label] != id:
                        warn = "\n\n[!] WARNING: While parsing your input data " \
                               "files, a file path was found with conflicting " \
                               "IDs for the same data level.\n\n" \
                               "File path: {0}\n" \
                               "Level: {1}\n" \
                               "Conflicting IDs: {2}, {3}\n\n" \
                               "This file has not been added to the data " \
                               "configuration.".format(fmap_mag, label,
                                                       path_dct[label], id)
                        print warn
                        skip = True
                        break

                new_template = new_template.replace(part1, '', 1)
                new_template = new_template.replace(label, '', 1)

                new_path = new_path.replace(part1, '', 1)
                new_path = new_path.replace(id, '', 1)

            if skip:
                continue

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
                scan_id = None

            temp_fmap_dct = {"magnitude": fmap_mag}

            # TODO: fill these
            if site_id not in data_dct.keys():
                continue
            if sub_id not in data_dct[site_id].keys():
                continue
            if ses_id not in data_dct[site_id][sub_id].keys():
                continue

            # TODO: reintroduce once scan-level nesting is implemented
            '''
            if scan_id:
                # if the field map files are at scan level
                if 'fmap' not in data_dct[site_id][sub_id][ses_id][scan_id].keys():
                    data_dct[site_id][sub_id][ses_id]['func'][scan_id]['fmap'] = temp_fmap_dct
                else:
                    data_dct[site_id][sub_id][ses_id]['func'][scan_id]['fmap'].update(temp_fmap_dct)

            else:
                # if the field map fields are at session level
                if 'fmap' not in data_dct[site_id][sub_id][ses_id].keys():
                    data_dct[site_id][sub_id][ses_id]['fmap'] = temp_fmap_dct
                else:
                    data_dct[site_id][sub_id][ses_id]['fmap'].update(temp_fmap_dct)
            '''

            # if the field map fields are at session level
            if 'fmap' not in data_dct[site_id][sub_id][ses_id].keys():
                data_dct[site_id][sub_id][ses_id]['fmap'] = temp_fmap_dct
            else:
                data_dct[site_id][sub_id][ses_id]['fmap'].update(
                    temp_fmap_dct)

    return data_dct


def run(data_settings_yml):

    import os
    import yaml

    print "\nGenerating data configuration file.."

    with open(data_settings_yml, "r") as f:
        settings_dct = yaml.load(f)

    incl_dct = format_incl_excl_dct(settings_dct['siteList'],
                                    settings_dct['subjectList'],
                                    settings_dct['sessionList'],
                                    settings_dct['scanList'])

    excl_dct = format_incl_excl_dct(settings_dct['exclusionSiteList'],
                                    settings_dct['exclusionSubjectList'],
                                    settings_dct['exclusionSessionList'],
                                    settings_dct['exclusionScanList'])

    if 'BIDS' in settings_dct['dataFormat'] or \
            'bids' in settings_dct['dataFormat']:
        # local (not on AWS S3 bucket), BIDS files

        #TODO: scan params

        if "s3://" not in settings_dct['bidsBaseDir']:
            if not os.path.isdir(settings_dct['bidsBaseDir']):
                err = "\n[!] The BIDS base directory you provided does not " \
                      "exist:\n{0}\n\n".format(settings_dct['bidsBaseDir'])
                raise Exception(err)

        params_dct = None

        data_dct = get_BIDS_data_dct(settings_dct['bidsBaseDir'],
                                     aws_creds_path=settings_dct['awsCredentialsFile'],
                                     inclusion_dct=incl_dct,
                                     exclusion_dct=excl_dct)

    elif 'Custom' in settings_dct['dataFormat'] or \
            'custom' in settings_dct['dataFormat']:
        # local (not on AWS S3 bucket), non-BIDS files

        params_dct = None
        if settings_dct['scanParametersCSV']:
            if '.csv' in settings_dct['scanParametersCSV']:
                params_dct = \
                    extract_scan_params_csv(settings_dct['scanParametersCSV'])

        data_dct = get_nonBIDS_data(settings_dct['anatomicalTemplate'],
                                    settings_dct['functionalTemplate'],
                                    params_dct,
                                    settings_dct['fieldMapPhase'],
                                    settings_dct['fieldMapMagnitude'],
                                    settings_dct['awsCredentialsFile'],
                                    incl_dct, excl_dct)

    #TODO: check data_dct for emptiness!

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
                         "data_config_{0}.yml"
                         "".format(settings_dct['subjectListName']))

        # put data_dct contents in an ordered list for the YAML dump
        data_list = []
        for site in sorted(data_dct.keys()):
            for sub in sorted(data_dct[site].keys()):
                for ses in sorted(data_dct[site][sub].keys()):
                    data_list.append(data_dct[site][sub][ses])

        with open(data_config_outfile, "wt") as f:
            # Make sure YAML doesn't dump aliases (so it's more human
            # read-able)
            f.write("# CPAC Data Configuration File\n# Version 1.0.3\n")
            f.write("#\n# http://fcp-indi.github.io for more info.\n#\n"
                    "# Tip: This file can be edited manually with "
                    "a text editor for quick modifications.\n\n")
            noalias_dumper = yaml.dumper.SafeDumper
            noalias_dumper.ignore_aliases = lambda self, data: True
            f.write(yaml.dump(data_list, default_flow_style=False,
                              Dumper=noalias_dumper))

        if os.path.exists(data_config_outfile):
            print "\nCPAC DATA SETTINGS file entered:" \
                  "\n{0}".format(data_settings_yml)
            print "\nCPAC DATA CONFIGURATION file created:" \
                  "\n{0}\n".format(data_config_outfile)
            print "Number of:"
            print "...sites: {0}".format(num_sites)
            print "...participants: {0}".format(num_subs)
            print "...participant-sessions: {0}".format(num_sess)
            print "...functional scans: {0}\n".format(num_scan)

    else:
        #TODO: fill this
        print "error nothing found"




