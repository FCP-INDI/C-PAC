

def gather_file_paths(base_directory, verbose=False):

    # this will go into core tools eventually

    # ideas: return number of paths, optionally instead
    #        or write out, optionally, a text file with all paths, for easy
    #        searching

    # test cases:
    #   - proper data directory
    #   - empty directory

    # TODO: this is not being used by the data config builder- do we need?

    import os

    path_list = []

    for root, dirs, files in os.walk(base_directory):
        for path in files:
            fullpath = os.path.join(root, path)
            path_list.append(fullpath)

    if verbose:
        print "Number of paths: {0}".format(len(path_list))

    return path_list


def download_single_s3_path(s3_path, download_dir=None, creds_path=None):
    """Download a single file from an AWS s3 bucket.

    :type s3_path: str
    :param s3_path: An "s3://" pre-pended path to a file stored on an
                    Amazon AWS s3 bucket.
    :type cfg_dict: dictionary
    :param cfg_dict: A dictionary containing the pipeline setup
                     parameters.
    :rtype: str
    :return: The local filepath of the downloaded s3 file.
    """

    import os
    from indi_aws import fetch_creds, aws_utils

    if not download_dir:
        download_dir = os.getcwd()

    if "s3://" in s3_path:
        s3_prefix = s3_path.replace("s3://", "")
    else:
        err = "[!] S3 file paths must be pre-pended with the 's3://' prefix."
        raise Exception(err)

    bucket_name = s3_prefix.split("/")[0]
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    data_dir = s3_path.split(bucket_name + "/")[1]
    local_dl = os.path.join(download_dir, data_dir)

    if os.path.isfile(local_dl):
        print "\nS3 bucket file already downloaded! Skipping download."
        print "S3 file: %s" % s3_path
        print "Local file already exists: %s\n" % local_dl
    else:
        aws_utils.s3_download(bucket, ([data_dir], [local_dl]))

    return local_dl


def pull_s3_sublist(data_folder, creds_path=None, keep_prefix=True):

    import os
    from indi_aws import fetch_creds

    if creds_path:
        creds_path = os.path.abspath(creds_path)

    s3_path = data_folder.split("s3://")[1]
    bucket_name = s3_path.split("/")[0]
    bucket_prefix = s3_path.split(bucket_name + "/")[1]

    print "Pulling from {0} ...".format(data_folder)

    s3_list = []
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    # ensure slash at end of bucket_prefix, so that if the final
    # directory name is a substring in other directory names, these
    # other directories will not be pulled into the file list
    if "/" not in bucket_prefix[-1]:
        bucket_prefix += "/"

    # Build S3-subjects to download
    for bk in bucket.objects.filter(Prefix=bucket_prefix):
        if keep_prefix:
            fullpath = os.path.join("s3://", bucket_name, str(bk.key))
            s3_list.append(fullpath)
        else:
            s3_list.append(str(bk.key).replace(bucket_prefix, ""))

    print "Finished pulling from S3. " \
          "{0} file paths found.".format(len(s3_list))

    if not s3_list:
        err = "\n\n[!] No input data found matching your data settings in " \
              "the AWS S3 bucket provided:\n{0}\n\n".format(data_folder)
        raise Exception(err)

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

    keys = {"TR (seconds)": "TR",
            "TE (seconds)": "TE",
            "Reference (slice no)": "reference",
            "Acquisition (pattern)": "acquisition",
            "FirstTR (start volume index)": "first_TR",
            "LastTR (final volume index)": "last_TR"}

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

            site_dict[site][sub][ses] = {keys[key]: val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and
                                         key != 'Participant' and
                                         key != 'Session' and key != 'Series'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            #site_dict[site][sub][ses]['tr'] = \
            #    site_dict[site][sub][ses].pop('tr (seconds)')

        elif sub != "All":
            # participant-specific scan parameters
            if site not in site_dict.keys():
                site_dict[site] = {}
            if sub not in site_dict[site].keys():
                site_dict[site][sub] = {}

            site_dict[site][sub][ses] = {keys[key]: val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and
                                         key != 'Participant' and
                                         key != 'Session' and key != 'Series'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            #site_dict[site][sub][ses]['tr'] =
            #    site_dict[site][sub][ses].pop('tr (seconds)')

        else:
            # site-specific scan parameters only
            if site not in site_dict.keys():
                site_dict[site] = {}
            if sub not in site_dict[site].keys():
                site_dict[site][sub] = {}

            site_dict[site][sub][ses] = {keys[key]: val
                                         for key, val in dict_row.items()
                                         if key != 'Site' and
                                         key != 'Participant' and
                                         key != 'Session' and key != 'Series'}

            # Assumes all other fields are formatted properly, but TR might
            # not be
            #site_dict[site][sub][ses]['tr'] = \
            #    site_dict[site][sub][ses].pop('tr (seconds)')

    return site_dict


def extract_scan_params_json(scan_params_json):

    # TODO: this is never used in data config builder- do we/will we need it?

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


def get_BIDS_data_dct(bids_base_dir, file_list=None, aws_creds_path=None,
                      inclusion_dct=None, exclusion_dct=None,
                      config_dir=None):

    import os
    import glob

    if not config_dir:
        config_dir = os.getcwd()

    anat_sess = os.path.join(bids_base_dir,
                             "sub-{participant}/ses-{session}/anat/sub-"
                             "{participant}_ses-{session}_T1w.nii.gz")
    anat = os.path.join(bids_base_dir,
                        "sub-{participant}/anat/sub-{participant}_T1w.nii.gz")

    func_sess = os.path.join(bids_base_dir,
                             "sub-{participant}"
                             "/ses-{session}/func/sub-"
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

    part_tsv_glob = os.path.join(bids_base_dir, "*participants.tsv")

    site_jsons_glob = os.path.join(bids_base_dir, "*bold.json")
    sub_jsons_glob = os.path.join(bids_base_dir, "*sub-*/*bold.json")
    ses_jsons_glob = os.path.join(bids_base_dir, "*sub-*/ses-*/*bold.json")

    site_dir_glob = os.path.join(bids_base_dir, "*", "sub-*/*.nii*")

    ses = False
    site_dir = False
    part_tsv = None

    if file_list:
        import fnmatch

        site_jsons = []
        sub_jsons = []
        ses_jsons = []

        for filepath in file_list:

            if fnmatch.fnmatch(filepath, site_dir_glob) and \
                    "derivatives" not in filepath:
                # check if there is a directory level encoding site ID, even
                # though that is not BIDS format
                site_dir = True
                sess_glob = os.path.join(bids_base_dir, "*", "sub-*/ses-*/*")

            if fnmatch.fnmatch(filepath, sess_glob):
                # check if there is a session level
                ses = True
            if fnmatch.fnmatch(filepath, part_tsv_glob):
                # check if there is a participants.tsv file
                part_tsv = filepath

            if fnmatch.fnmatch(filepath, ses_jsons_glob):
                ses_jsons.append(filepath)
            if fnmatch.fnmatch(filepath, sub_jsons_glob):
                sub_jsons.append(filepath)
            if fnmatch.fnmatch(filepath, site_jsons_glob):
                site_jsons.append(filepath)

    else:
        if len(glob.glob(site_dir_glob)) > 0:
            # check if there is a directory level encoding site ID, even
            # though that is not BIDS format
            site_dir = True
            sess_glob = os.path.join(bids_base_dir, "*", "sub-*/ses-*/*")

        ses = False
        if len(glob.glob(sess_glob)) > 0:
            # check if there is a session level
            ses = True

        # check if there is a participants.tsv file
        part_tsv = glob.glob(part_tsv_glob)
        site_jsons = glob.glob(site_jsons_glob)
        sub_jsons = glob.glob(sub_jsons_glob)
        ses_jsons = glob.glob(ses_jsons_glob)

    sites_dct = {}
    sites_subs_dct = {}

    if part_tsv:
        # handle participants.tsv file in BIDS dataset if one is present
        # this would contain site information if the dataset is multi-site
        import csv

        if file_list:
            print "\n\nFound a participants.tsv file in your BIDS data " \
                  "set on the S3 bucket. Downloading..\n"
            part_tsv = download_single_s3_path(part_tsv, config_dir,
                                               aws_creds_path)

        print "Checking participants.tsv file for site information:" \
              "\n{0}".format(part_tsv)

        with open(part_tsv, "r") as f:
            tsv = csv.DictReader(f)
            for row in tsv:
                if "site" in row.keys():
                    site = row["site"]
                    sub = row["participant_id"]
                    if site not in sites_dct.keys():
                        sites_dct[site] = []
                    sites_dct[site].append(sub)

        if sites_dct:
            # check for duplicates
            sites = sites_dct.keys()
            print "{0} sites found in the participant.tsv " \
                  "file.".format(len(sites))
            for site in sites:
                for other_site in sites:
                    if site == other_site:
                        continue
                    dups = set(sites_dct[site]) & set(sites_dct[other_site])
                    if dups:
                        err = "\n\n[!] There are duplicate participant IDs " \
                              "in different sites, as defined by your " \
                              "participants.tsv file! Consider pre-fixing " \
                              "the participant IDs with the site names.\n\n" \
                              "Duplicates:\n" \
                              "Sites: {0}, {1}\n" \
                              "Duplicate IDs: {2}" \
                              "\n\n".format(site, other_site, str(dups))
                        raise Exception(err)

                # now invert
                for sub in sites_dct[site]:
                    sites_subs_dct[sub] = site
        else:
            print "No site information found in the participants.tsv file."

    if not sites_subs_dct:
        # if there was no participants.tsv file, (or no site column in the
        # participants.tsv file), check for a directory level where multiple
        # sites might be encoded
        if site_dir:
            sub_dir = os.path.join(bids_base_dir, "sub-{participant}")
            new_dir = os.path.join(bids_base_dir, "{site}",
                                   "sub-{participant}")
            if ses:
                anat_sess = anat_sess.replace(sub_dir, new_dir)
                func_sess = func_sess.replace(sub_dir, new_dir)
                fmap_phase_sess = fmap_phase_sess.replace(sub_dir, new_dir)
                fmap_mag_sess = fmap_mag_sess.replace(sub_dir, new_dir)
            else:
                anat = anat.replace(sub_dir, new_dir)
                func = func.replace(sub_dir, new_dir)
                fmap_phase = fmap_phase.replace(sub_dir, new_dir)
                fmap_mag = fmap_mag.replace(sub_dir, new_dir)

    scan_params_dct = {}

    if len(ses_jsons) > 0:
        for json_file in ses_jsons:
            ids = os.path.basename(json_file).rstrip(".json").split("_")
            for id in ids:
                if "sub-" in id:
                    sub_id = id.replace("sub-", "")
                if "ses-" in id:
                    ses_id = id.replace("ses-", "")

            site_id = "All"
            if sites_subs_dct:
                if sub_id in sites_subs_dct.keys():
                    site_id = sites_subs_dct[sub_id]

            if site_id not in scan_params_dct.keys():
                scan_params_dct[site_id] = {}
            if sub_id not in scan_params_dct[site_id].keys():
                scan_params_dct[site_id][sub_id] = {}

            scan_params_dct[site_id][sub_id][ses_id] = json_file

    if len(sub_jsons) > 0:
        for json_file in sub_jsons:
            ids = os.path.basename(json_file).rstrip(".json").split("_")
            for id in ids:
                if "sub-" in id:
                    sub_id = id.replace("sub-", "")

            site_id = "All"
            if sites_subs_dct:
                if sub_id in sites_subs_dct.keys():
                    site_id = sites_subs_dct[sub_id]

            if site_id not in scan_params_dct.keys():
                scan_params_dct[site_id] = {}
            if sub_id not in scan_params_dct[site_id].keys():
                scan_params_dct[site_id][sub_id] = {}

            scan_params_dct[site_id][sub_id]["All"] = json_file

    if len(site_jsons) > 0:
        if len(site_jsons) > 1:
            if not site_dir:
                print "\nMore than one dataset-level JSON file. These may " \
                      "be scan-specific JSON files. Scan-specific scan " \
                      "parameters will be implemented soon. Files detected:"
                for json in site_jsons:
                    print "...{0}".format(json)
                print "\n"
            else:
                # if this kicks, then there is a site directory level, and the
                # site-specific JSON is sitting in it alongside the sub-*
                # participant ID folders
                for json in site_jsons:
                    json_levels = json.split("/")
                    site_name = json_levels[-2]
                    if site_name not in scan_params_dct.keys():
                        scan_params_dct[site_name] = {}
                    if "All" not in scan_params_dct[site_name].keys():
                        scan_params_dct[site_name]["All"] = {}
                    scan_params_dct[site_name]["All"]["All"] = json
        else:
            if "All" not in scan_params_dct.keys():
                scan_params_dct["All"] = {}
            if "All" not in scan_params_dct["All"].keys():
                scan_params_dct["All"]["All"] = {}

            for json_file in site_jsons:
                scan_params_dct["All"]["All"]["All"] = json_file

    if ses:
        # if there is a session level in the BIDS dataset
        data_dct = get_nonBIDS_data(anat_sess, func_sess, file_list=file_list,
                                    scan_params_dct=scan_params_dct,
                                    fmap_phase_template=fmap_phase_sess,
                                    fmap_mag_template=fmap_mag_sess,
                                    aws_creds_path=aws_creds_path,
                                    inclusion_dct=inclusion_dct,
                                    exclusion_dct=exclusion_dct,
                                    sites_dct=sites_subs_dct)
    else:
        # no session level
        data_dct = get_nonBIDS_data(anat, func, file_list=file_list,
                                    scan_params_dct=scan_params_dct,
                                    fmap_phase_template=fmap_phase,
                                    fmap_mag_template=fmap_mag,
                                    aws_creds_path=aws_creds_path,
                                    inclusion_dct=inclusion_dct,
                                    exclusion_dct=exclusion_dct,
                                    sites_dct=sites_subs_dct)

    return data_dct


def get_nonBIDS_data(anat_template, func_template, file_list=None,
                     scan_params_dct=None, fmap_phase_template=None,
                     fmap_mag_template=None, aws_creds_path=None,
                     inclusion_dct=None, exclusion_dct=None, sites_dct=None,
                     verbose=False):

    # go over the file paths, validate for nifti's?
    # work with the template

    # test cases:
    #   - no anats, no funcs
    #   - combination of .nii and .nii.gz?
    #   - throw error/warning if anat and func templates are identical
    #   - all permutations of scan parameters json/csv's at different levels

    import glob
    import fnmatch

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
    anat_pool = []
    func_pool = []

    if file_list:
        # mainly for AWS S3-stored data sets
        for filepath in file_list:
            if fnmatch.fnmatch(filepath, anat_glob):
                anat_pool.append(filepath)
            elif fnmatch.fnmatch(filepath, func_glob):
                func_pool.append(filepath)

    # run it anyway in case we're pulling anat from S3 and func from local or
    # vice versa - and if there is no file_list, this will run normally
    anat_local_pool = glob.glob(anat_glob)
    func_local_pool = glob.glob(func_glob)

    anat_pool = anat_pool + anat_local_pool
    func_pool = func_pool + func_local_pool

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
                warn = None
                if path_dct[label] != id:
                    if str(path_dct[label]) in id and "run-" in id:
                        # TODO: this is here only because we do not support
                        # TODO: multiple anat scans yet!!
                        if "run-1" in id:
                            warn = None
                            pass
                        else:
                            warn = "\n\n[!] WARNING: Multiple anatomical " \
                                   "scans found for a single participant. " \
                                   "Multiple anatomical scans are not yet " \
                                   "supported. Review the completed data " \
                                   "configuration file to ensure the " \
                                   "proper scans are included together.\n\n" \
                                   "Scan not included:\n{0}" \
                                   "\n\n".format(anat_path)
                    else:
                        warn = "\n\n[!] WARNING: While parsing your input data " \
                               "files, a file path was found with conflicting " \
                               "IDs for the same data level.\n\n" \
                               "File path: {0}\n" \
                               "Level: {1}\n" \
                               "Conflicting IDs: {2}, {3}\n\n" \
                               "This file has not been added to the data " \
                               "configuration.".format(anat_path, label,
                                                       path_dct[label], id)
                if warn:
                    print warn
                    skip = True
                    break
                else:
                    skip = False
                    pass

            new_template = new_template.replace(part1, '', 1)
            new_template = new_template.replace(label, '', 1)

            new_path = new_path.replace(part1, '', 1)
            new_path = new_path.replace(id, '', 1)

        if skip:
            continue

        sub_id = path_dct['{participant}']

        if '{site}' in path_dct.keys():
            site_id = path_dct['{site}']
        elif sites_dct:
            # mainly if we're pulling site info from a participants.tsv file
            # for a BIDS data set
            try:
                site_id = sites_dct[sub_id]
            except KeyError:
                site_id = "site-1"
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
            # doubt this ever happens, but just be safe
            warn = "\n\n[!] WARNING: Duplicate site-participant-session " \
                   "entry found in your input data directory!\n\nDuplicate " \
                   "sets:\n\n{0}\n\n{1}\n\nOnly adding the first one to " \
                   "the data configuration file." \
                   "\n\n".format(str(data_dct[site_id][sub_id][ses_id]),
                                 str(temp_sub_dct))
            print warn

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
        elif sites_dct:
            # mainly if we're pulling site info from a participants.tsv file
            # for a BIDS data set
            try:
                site_id = sites_dct[sub_id]
            except KeyError:
                site_id = "site-1"
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

            # TODO: implement scan-specific level once scan-level nesting is
            # TODO: available

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

        if site_id not in data_dct.keys():
            if verbose:
                print "No anatomical for functional for site:" \
                      "\n{0}\n{1}\n".format(func_path, site_id)
            continue
        if sub_id not in data_dct[site_id].keys():
            if verbose:
                print "No anatomical for functional for participant:" \
                      "\n{0}\n{1}\n".format(func_path, sub_id)
            continue
        if ses_id not in data_dct[site_id][sub_id].keys():
            if verbose:
                print "No anatomical for functional for session:" \
                      "\n{0}\n{1}\n".format(func_path, ses_id)
            continue

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
        if file_list:
            # mainly for AWS S3-stored data sets
            fmap_phase_pool = []
            fmap_mag_pool = []
            for filepath in file_list:
                if fnmatch.fnmatch(filepath, fmap_phase_glob):
                    fmap_phase_pool.append(filepath)
                elif fnmatch.fnmatch(filepath, fmap_mag_glob):
                    fmap_mag_pool.append(filepath)
        else:
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
                    warn = None
                    if path_dct[label] != id:
                        if str(path_dct[label]) in id and "run-" in id:
                            # TODO: this is here only because we do not
                            # TODO: support scan-level nesting yet!!!
                            warn = "\n\n[!] WARNING: Multiple field map " \
                                   "phase or magnitude files were found " \
                                   "for a single session. Multiple-run " \
                                   "field map files are not supported yet. " \
                                   "Review the completed data " \
                                   "configuration file to ensure the " \
                                   "proper scans are included together.\n\n" \
                                   "File not included:\n{0}" \
                                   "\n\n".format(fmap_phase)
                        else:
                            warn = "\n\n[!] WARNING: While parsing your input data " \
                                   "files, a file path was found with conflicting " \
                                   "IDs for the same data level.\n\n" \
                                   "File path: {0}\n" \
                                   "Level: {1}\n" \
                                   "Conflicting IDs: {2}, {3}\n\n" \
                                   "This file has not been added to the data " \
                                   "configuration.".format(fmap_phase, label,
                                                           path_dct[label], id)
                    if warn:
                        print warn
                        skip = True
                        break
                    else:
                        skip = False

                new_template = new_template.replace(part1, '', 1)
                new_template = new_template.replace(label, '', 1)

                new_path = new_path.replace(part1, '', 1)
                new_path = new_path.replace(id, '', 1)

            if skip:
                continue

            sub_id = path_dct['{participant}']

            if '{site}' in path_dct.keys():
                site_id = path_dct['{site}']
            elif sites_dct:
                # mainly if we're pulling site info from a participants.tsv
                # file for a BIDS data set
                try:
                    site_id = sites_dct[sub_id]
                except KeyError:
                    site_id = "site-1"
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

            if site_id not in data_dct.keys():
                if verbose:
                    print "Missing inputs (no anat/func) for field map for " \
                          "site:\n{0}\n{1}\n".format(fmap_phase, site_id)
                continue
            if sub_id not in data_dct[site_id].keys():
                if verbose:
                    print "Missing inputs (no anat/func) for field map for " \
                          "participant:" \
                          "\n{0}\n{1}\n".format(fmap_phase, sub_id)
                continue
            if ses_id not in data_dct[site_id][sub_id].keys():
                if verbose:
                    print "Missing inputs (no anat/func) for field map for " \
                          "session:\n{0}\n{1}\n".format(fmap_phase, ses_id)
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
                    warn = None
                    if path_dct[label] != id:
                        if str(path_dct[label]) in id and "run-" in id:
                            # TODO: this is here only because we do not
                            # TODO: support scan-level nesting yet!!!
                            warn = "\n\n[!] WARNING: Multiple field map " \
                                   "phase or magnitude files were found " \
                                   "for a single session. Multiple-run " \
                                   "field map files are not supported yet. " \
                                   "Review the completed data " \
                                   "configuration file to ensure the " \
                                   "proper scans are included together.\n\n" \
                                   "File not included:\n{0}" \
                                   "\n\n".format(fmap_phase)
                        else:
                            warn = "\n\n[!] WARNING: While parsing your input data " \
                                   "files, a file path was found with conflicting " \
                                   "IDs for the same data level.\n\n" \
                                   "File path: {0}\n" \
                                   "Level: {1}\n" \
                                   "Conflicting IDs: {2}, {3}\n\n" \
                                   "This file has not been added to the data " \
                                   "configuration.".format(fmap_mag, label,
                                                           path_dct[label], id)
                    if warn:
                        print warn
                        skip = True
                        break
                    else:
                        skip = False

                new_template = new_template.replace(part1, '', 1)
                new_template = new_template.replace(label, '', 1)

                new_path = new_path.replace(part1, '', 1)
                new_path = new_path.replace(id, '', 1)

            if skip:
                continue

            sub_id = path_dct['{participant}']

            if '{site}' in path_dct.keys():
                site_id = path_dct['{site}']
            elif sites_dct:
                # mainly if we're pulling site info from a participants.tsv
                # file for a BIDS data set
                try:
                    site_id = sites_dct[sub_id]
                except KeyError:
                    site_id = "site-1"
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

            if site_id not in data_dct.keys():
                if verbose:
                    print "Missing inputs (no anat/func) for field map for " \
                          "site:\n{0}\n{1}\n".format(fmap_mag, site_id)
                continue
            if sub_id not in data_dct[site_id].keys():
                if verbose:
                    print "Missing inputs (no anat/func) for field map for " \
                          "participant:" \
                          "\n{0}\n{1}\n".format(fmap_mag, sub_id)
                continue
            if ses_id not in data_dct[site_id][sub_id].keys():
                if verbose:
                    print "Missing inputs (no anat/func) for field map for " \
                          "session:\n{0}\n{1}\n".format(fmap_mag, ses_id)
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

    if "None" in settings_dct["awsCredentialsFile"]:
        settings_dct["awsCredentialsFile"] = None

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

        if "s3://" in settings_dct['bidsBaseDir']:
            # hosted on AWS S3 bucket
            file_list = pull_s3_sublist(settings_dct['bidsBaseDir'],
                                        settings_dct['awsCredentialsFile'])
        else:
            # local (not on AWS S3 bucket), BIDS files
            if not os.path.isdir(settings_dct['bidsBaseDir']):
                err = "\n[!] The BIDS base directory you provided does not " \
                      "exist:\n{0}\n\n".format(settings_dct['bidsBaseDir'])
                raise Exception(err)
            file_list = None

        data_dct = get_BIDS_data_dct(settings_dct['bidsBaseDir'],
                                     file_list=file_list,
                                     aws_creds_path=settings_dct['awsCredentialsFile'],
                                     inclusion_dct=incl_dct,
                                     exclusion_dct=excl_dct,
                                     config_dir=settings_dct["outputSubjectListLocation"])

    elif 'Custom' in settings_dct['dataFormat'] or \
            'custom' in settings_dct['dataFormat']:

        # keep as None if local data set (not on AWS S3 bucket)
        file_list = None
        base_dir = None

        if "s3://" in settings_dct["anatomicalTemplate"]:
            # hosted on AWS S3 bucket
            if '{site}' in settings_dct["anatomicalTemplate"]:
                base_dir = \
                    settings_dct["anatomicalTemplate"].split('{site}')[0]
            elif '{participant}' in settings_dct["anatomicalTemplate"]:
                base_dir = \
                    settings_dct["anatomicalTemplate"].split('{participant}')[0]

        elif "s3://" in settings_dct["functionalTemplate"]:
            # hosted on AWS S3 bucket
            if '{site}' in settings_dct["functionalTemplate"]:
                base_dir = \
                    settings_dct["functionalTemplate"].split('{site}')[0]
            elif '{participant}' in settings_dct["functionalTemplate"]:
                base_dir = \
                    settings_dct["functionalTemplate"].split('{participant}')[0]

        if base_dir:
            file_list = pull_s3_sublist(base_dir,
                                        settings_dct['awsCredentialsFile'])

        params_dct = None
        if settings_dct['scanParametersCSV']:
            if '.csv' in settings_dct['scanParametersCSV']:
                params_dct = \
                    extract_scan_params_csv(settings_dct['scanParametersCSV'])

        data_dct = get_nonBIDS_data(settings_dct['anatomicalTemplate'],
                                    settings_dct['functionalTemplate'],
                                    file_list=file_list,
                                    scan_params_dct=params_dct,
                                    fmap_phase_template=settings_dct['fieldMapPhase'],
                                    fmap_mag_template=settings_dct['fieldMapMagnitude'],
                                    aws_creds_path=settings_dct['awsCredentialsFile'],
                                    inclusion_dct=incl_dct,
                                    exclusion_dct=excl_dct)

    else:
        err = "\n\n[!] You must select a data format- either 'BIDS' or " \
              "'Custom', in the 'dataFormat' field in the data settings " \
              "YAML file.\n\n"
        raise Exception(err)

    if len(data_dct) > 0:

        # get some data
        num_sites = len(data_dct.keys())
        num_subs = num_sess = num_scan = 0
        for site in data_dct.keys():
            num_subs += len(data_dct[site])
            for sub in data_dct[site].keys():
                num_sess += len(data_dct[site][sub])
                for session in data_dct[site][sub].keys():
                    if 'func' in data_dct[site][sub][session].keys():
                        for scan in data_dct[site][sub][session]['func'].keys():
                            if scan != "scan_parameters":
                                num_scan += 1

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
        err = "\n\n[!] No anatomical input files were found given the data " \
              "settings provided.\n\n"
        raise Exception(err)




