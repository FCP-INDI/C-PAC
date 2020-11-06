import os
import yaml
import json


def bids_decode_fname(file_path, dbg=False, raise_error=True):
    import re

    f_dict = {}

    fname = os.path.basename(file_path)

    # first lets make sure that we know how to handle the file
    if 'nii' not in fname.lower() and 'json' not in fname.lower():
        raise IOError("File (%s) does not appear to be" % fname +
                      "a nifti or json file")

    if dbg:
        print("parsing %s" % file_path)

    # first figure out if there is a site directory level, this isn't
    # specified in BIDS currently, but hopefully will be in the future
    file_path_vals = os.path.dirname(file_path).split('/')
    sub = [s for s in file_path_vals if 'sub-' in s]
    if dbg:
        print("found subject %s in %s" % (sub, str(file_path_vals)))

    if len(sub) > 1:
        print("Odd that there is more than one subject directory" +
               "in (%s), does the filename conform to" % file_path +
               " BIDS format?")
    if sub:
        sub_ndx = file_path_vals.index(sub[0])
        if sub_ndx > 0 and file_path_vals[sub_ndx - 1]:
            if dbg:
                print("setting site to %s" % (file_path_vals[sub_ndx - 1]))
            f_dict["site"] = file_path_vals[sub_ndx - 1]
        else:
            f_dict["site"] = "none"
    elif file_path_vals[-1]:
        if dbg:
            print("looking for subject id didn't pan out settling for last"+
                   "subdir %s" % (str(file_path_vals[-1])))
        f_dict["site"] = file_path_vals[-1]
    else:
        f_dict["site"] = "none"

    f_dict["site"] = re.sub('[\s\-\_]+', '', f_dict["site"])

    fname = fname.split(".")[0]
    # convert the filename string into a dictionary to pull out the other
    # key value pairs

    for key_val_pair in fname.split("_"):
        # if the chunk has the shape key-val store key: val in f_dict
        if "-" in key_val_pair:
            chunks = key_val_pair.split("-")
            f_dict[chunks[0]] = "-".join(chunks[1:])
        else:
            f_dict["scantype"] = key_val_pair.split(".")[0]

    if "scantype" not in f_dict:
        msg = "Filename ({0}) does not appear to contain" \
              " scan type, does it conform to the BIDS format?".format(fname)
        if raise_error:
            raise ValueError(msg)
        else:
            print(msg)
    elif not f_dict["scantype"]:
        msg = "Filename ({0}) does not appear to contain" \
              " scan type, does it conform to the BIDS format?".format(fname)
        if raise_error:
            raise ValueError(msg)
        else:
            print(msg)
    else:
        if 'bold' in f_dict["scantype"] and not f_dict["task"]:
            msg = "Filename ({0}) is a BOLD file, but " \
                  "doesn't contain a task, does it conform to the" \
                  " BIDS format?".format(fname)
            if raise_error:
                raise ValueError(msg)
            else:
                print(msg)

    return f_dict


def bids_retrieve_params(bids_config_dict, f_dict, dbg=False):
    """

    Retrieve the BIDS parameters from bids_config_dict for BIDS file
    corresponding to f_dict. If an exact match for f_dict is not found
    the nearest match is returned, corresponding to the BIDS inheritance
    principle.

    :param bids_config_dict: BIDS configuration dictionary, this is a
      multi-level dictionary that maps the components of a bids filename
      (i.e. sub, ses, acq, run) to a dictionary that contains the BIDS
      parameters (RepetitionTime, EchoTime, etc). This information is
      extracted from sidecar json files using the principle of inheritance
      using the bids_parse_configs function
    :param f_dict: Dictionary built from the name of a file in the BIDS
      format. This is built using the bids_decode_fname by splitting on
      "-" and "_" delimeters
    :param dbg: boolean flag that indicates whether or not debug statements
      should be printed, defaults to "False"
    :return: returns a dictionary that contains the BIDS parameters
    """
    params = {}

    t_dict = bids_config_dict  # pointer to current dictionary
    # try to populate the configuration using information
    # already in the list
    for level in ['scantype', 'site', 'sub', 'ses', 'task', 'acq',
                  'rec', 'dir', 'run']:
        if level in f_dict:
            key = "-".join([level, f_dict[level]])
        else:
            key = "-".join([level, "none"])

        if dbg:
            print(key)
        # if the key doesn't exist in the config dictionary, check to see if
        # the generic key exists and return that
        if key in t_dict:
            t_dict = t_dict[key]
        else:
            if dbg:
                print("Couldn't find %s, so going with %s" % (key,
                        "-".join([level, "none"])))
            key = "-".join([level, "none"])
            if key in t_dict:
                t_dict = t_dict[key]

    # if we have an image parameter dictionary at this level, use it to
    # initialize our configuration we look for "RepetitionTime", because
    #  according to the spec it is a mandatory parameter for JSON
    # sidecare files

    if dbg:
        print(t_dict)

    for key in t_dict.keys():
        if 'RepetitionTime' in key:
            params = t_dict
            break

    for k, v in params.items():
        if isinstance(v, str):
            params[k] = v.encode('ascii', errors='ignore')

    return params


def bids_parse_sidecar(config_dict, dbg=False, raise_error=True):
    # type: (dict, bool) -> dict
    """
    Uses the BIDS principle of inheritance to build a data structure that
    maps parameters in side car .json files to components in the names of
    corresponding nifti files.

    :param config_dict: dictionary that maps paths of sidecar json files
       (the key) to a dictionary containing the contents of the files (the values)
    :param dbg: boolean flag that indicates whether or not debug statements
       should be printed
    :return: a dictionary that maps parameters to components from BIDS filenames
       such as sub, sess, run, acq, and scan type
    """

    # we are going to build a large-scale data structure, consisting of many
    # levels of dictionaries to hold the data.
    bids_config_dict = {}

    # initialize 'default' entries, this essentially is a pointer traversal
    # of the dictionary
    t_dict = bids_config_dict
    for level in ['scantype', 'site', 'sub', 'ses', 'task',
                  'acq', 'rec', 'dir', 'run']:
        key = '-'.join([level, 'none'])
        t_dict[key] = {}
        t_dict = t_dict[key]

    if dbg:
        print(bids_config_dict)

    # get the paths to the json yaml files in config_dict, the paths contain
    # the information needed to map the parameters from the jsons (the vals
    # of the config_dict) to corresponding nifti files. We sort the list
    # by the number of path components, so that we can iterate from the outer
    # most path to inner-most, which will help us address the BIDS inheritance
    # principle
    config_paths = sorted(
        list(config_dict.keys()),
        key=lambda p: len(p.split('/'))
    )

    if dbg:
        print(config_paths)

    for cp in config_paths:

        if dbg:
            print("processing %s" % (cp))

        # decode the filepath into its various components as defined by  BIDS
        f_dict = bids_decode_fname(cp, raise_error=raise_error)

        # handling inheritance is a complete pain, we will try to handle it by
        # build the key from the bottom up, starting with the most
        # parsimonious possible, incorporating configuration information that
        # exists at each level

        # first lets try to find any parameters that already apply at this
        # level using the information in the json's file path
        t_params = bids_retrieve_params(bids_config_dict, f_dict)

        # now populate the parameters
        bids_config = {}
        if t_params:
            bids_config.update(t_params)

        # add in the information from this config file
        t_config = config_dict[cp]
        if t_config is list:
            t_config = t_config[0]

        try:
            bids_config.update(t_config)
        except ValueError:
            err = "\n[!] Could not properly parse the AWS S3 path provided " \
                  "- please double-check the bucket and the path.\n\nNote: " \
                  "This could either be an issue with the path or the way " \
                  "the data is organized in the directory. You can also " \
                  "try providing a specific site sub-directory.\n\n"
            raise ValueError(err)

        # now put the configuration in the data structure, by first iterating
        # to the location of the key, and then inserting it. When a key isn't
        # defined we use the "none" value. A "none" indicates that the
        # corresponding parameters apply to all possible settings of that key
        # e.g. run-1, run-2, ... will all map to run-none if no jsons
        # explicitly define values for those runs
        t_dict = bids_config_dict  # pointer to current dictionary
        for level in ['scantype', 'site', 'sub', 'ses', 'task', 'acq',
                      'rec', 'dir', 'run']:
            if level in f_dict:
                key = "-".join([level, f_dict[level]])
            else:
                key = "-".join([level, "none"])

            if key not in t_dict:
                t_dict[key] = {}

            t_dict = t_dict[key]

        t_dict.update(bids_config)

    return(bids_config_dict)


def gen_bids_outputs_sublist(base_path, paths_list, key_list, creds_path):
    import copy

    func_keys = ["functional_to_anat_linear_xfm", "motion_params",
                 "movement_parameters", "motion_correct"]
    top_keys = list(set(key_list) - set(func_keys))
    bot_keys = list(set(key_list).intersection(func_keys))

    subjdict = {}

    if not base_path.endswith('/'):
        base_path = base_path + '/'

    # output directories are a bit different than standard BIDS, so
    # we handle things differently

    for p in paths_list:
        p = p.rstrip()

        # find the participant and session info which should be at
        # some level in the path
        path_base = p.replace(base_path, '')

        subj_info = path_base.split('/')[0]
        resource = path_base.split('/')[1]

        if resource not in key_list:
            continue

        if subj_info not in subjdict:
            subjdict[subj_info] = {"subj_info": subj_info}

        if creds_path:
            subjdict[subj_info]["creds_path"] = creds_path

        if resource in func_keys:
            run_info = path_base.split('/')[2]
            if "funcs" not in subjdict[subj_info]:
                subjdict[subj_info]["funcs"] = {}
            if run_info not in subjdict[subj_info]["funcs"]:
                subjdict[subj_info]["funcs"][run_info] = {'run_info': run_info}
            if resource in subjdict[subj_info]["funcs"][run_info]:
                print("warning resource %s already exists in subjdict ??" %
                      (resource))
            subjdict[subj_info]["funcs"][run_info][resource] = p
        else:
            subjdict[subj_info][resource] = p

    sublist = []
    for subj_info, subj_res in subjdict.items():
        missing = 0
        for tkey in top_keys:
            if tkey not in subj_res:
                print("%s not found for %s" % (tkey, subj_info))
                missing += 1
                break

        if missing == 0:
            for func_key, func_res in subj_res["funcs"].items():
                for bkey in bot_keys:
                    if bkey not in func_res:
                        print("%s not found for %s" % (bkey,
                                                       func_key))
                        missing += 1
                        break
                if missing == 0:
                    print("adding: %s, %s, %d" % (subj_info,
                                                  func_key,
                                                  len(sublist)))
                    tdict = copy.deepcopy(subj_res)
                    del tdict["funcs"]
                    tdict.update(func_res)
                    sublist.append(tdict)
    return sublist


def bids_gen_cpac_sublist(bids_dir, paths_list, config_dict, creds_path, dbg=False,
                          raise_error=True):
    """
    Generates a CPAC formatted subject list from information contained in a
    BIDS formatted set of data.

    :param bids_dir: base directory that contains all of the data, this could be
       a directory that contains data for a multiple BIDS datasets, in which
       case the intervening directories will be interpreted as site names
    :param paths_list: lists of all nifti files found in bids_dir, these paths
       are relative to bids_dir
    :param config_dict: dictionary that contains information from the json
       sidecars found in bids_dir, keys are relative paths and values are
       dictionaries containing all of the parameter information. if config_dict
       is None, the subject list will be built without the parameters
    :param creds_path: if using S3 bucket, this path credentials needed to
       access the bucket, if accessing anonymous bucket, this can be set
       to None
    :param dbg: boolean indicating whether or not the debug statements should
       be printed
    :return: a list of dictionaries suitable for use by CPAC to specify data
       to be processed
    """

    if dbg:
        print("gen_bids_sublist called with:")
        print("  bids_dir: {0}".format(bids_dir))
        print("  # paths: {0}".format(str(len(paths_list))))
        print("  config_dict: {0}".format(
            "missing" if not config_dict else "found")
        )
        print("  creds_path: {0}".format(creds_path))

    # if configuration information is not desired, config_dict will be empty,
    # otherwise parse the information in the sidecar json files into a dict
    # we can use to extract data for our nifti files
    if config_dict:
        bids_config_dict = bids_parse_sidecar(config_dict, raise_error=raise_error)

    subdict = {}

    for p in paths_list:
        p = p.rstrip()
        f = os.path.basename(p)

        if f.endswith(".nii") or f.endswith(".nii.gz"):

            f_dict = bids_decode_fname(p, raise_error=raise_error)

            if config_dict:
                t_params = bids_retrieve_params(bids_config_dict,
                                                f_dict)
                if not t_params:
                    print(f_dict)
                    print("Did not receive any parameters for %s," % (p) +
                          " is this a problem?")

                task_info = {"scan": os.path.join(bids_dir,p),
                             "scan_parameters": t_params.copy()}
            else:
                task_info = os.path.join(bids_dir ,p)

            if "ses" not in f_dict:
                f_dict["ses"] = "1"

            if "sub" not in f_dict:
                raise IOError("sub not found in %s," % (p) +
                              " perhaps it isn't in BIDS format?")

            if f_dict["sub"] not in subdict:
                subdict[f_dict["sub"]] = {}

            subjid = "-".join(["sub", f_dict["sub"]])

            if f_dict["ses"] not in subdict[f_dict["sub"]]:
                subdict[f_dict["sub"]][f_dict["ses"]] = \
                    {"creds_path": creds_path,
                     "site_id": "-".join(["site", f_dict["site"]]),
                     "subject_id": subjid,
                     "unique_id": "-".join(["ses", f_dict["ses"]])}

            if "T1w" in f_dict["scantype"]:
                if "lesion" in f_dict.keys() and "mask" in f_dict['lesion']:
                    if "lesion_mask" not in \
                            subdict[f_dict["sub"]][f_dict["ses"]]:
                        subdict[f_dict["sub"]][f_dict["ses"]]["lesion_mask"] = \
                            task_info["scan"]
                    else:
                        print("Lesion mask file (%s) already found" %
                              (subdict[f_dict["sub"]]
                               [f_dict["ses"]]
                               ["lesion_mask"]) + " for (%s:%s) discarding %s" %
                              (f_dict["sub"], f_dict["ses"], p))
                # TODO deal with scan parameters anatomical
                elif "anat" not in subdict[f_dict["sub"]][f_dict["ses"]]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["anat"] = \
                        task_info["scan"] if config_dict else task_info
                else:
                    print("Anatomical file (%s) already found" %
                          (subdict[f_dict["sub"]][f_dict["ses"]]["anat"]) +
                          " for (%s:%s) discarding %s" % (f_dict["sub"],
                                                          f_dict["ses"],
                                                          p))
            if "bold" in f_dict["scantype"]:
                task_key = f_dict["task"]
                if "run" in f_dict:
                    task_key = "_".join([task_key,
                                         "-".join(["run", f_dict["run"]])])
                if "acq" in f_dict:
                    task_key = "_".join([task_key,
                                         "-".join(["acq", f_dict["acq"]])])
                if "func" not in subdict[f_dict["sub"]][f_dict["ses"]]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["func"] = {}

                if task_key not in \
                    subdict[f_dict["sub"]][f_dict["ses"]]["func"]:

                    if not isinstance(task_info, dict):
                        task_info = {"scan": task_info}
                    subdict[f_dict["sub"]][f_dict["ses"]]["func"][task_key] = task_info

                else:
                    print("Func file (%s)" %
                        subdict[f_dict["sub"]][f_dict["ses"]]["func"][task_key] +
                        " already found for ( % s: %s: % s) discarding % s" % (
                               f_dict["sub"],
                               f_dict["ses"],
                               task_key,
                               p))

            if "epi" in f_dict["scantype"]:
                pe_dir = f_dict["dir"]
                if "acq" in f_dict:
                    if "fMRI" in f_dict["acq"]:
                        if "fmap" not in subdict[f_dict["sub"]][f_dict["ses"]]:
                            subdict[f_dict["sub"]][f_dict["ses"]]["fmap"] = {}
                        if "epi_{0}".format(
                            pe_dir
                        ) not in subdict[f_dict["sub"]][f_dict["ses"]]["fmap"]:
                            subdict[f_dict["sub"]][
                                f_dict["ses"]
                            ]["fmap"]["epi_{0}".format(pe_dir)] = task_info

    sublist = []
    for ksub, sub in subdict.items():
        for kses, ses in sub.items():
            if "anat" in ses and "func" in ses:
                sublist.append(ses)
            else:
                if "anat" not in ses:
                    print("%s %s %s is missing an anat" % (
                        ses["site_id"] if 'none' not in ses["site_id"] else '',
                        ses["subject_id"],
                        ses["unique_id"]
                    ))
                if "func" not in ses:
                    print("%s %s %s is missing an func" % (
                        ses["site_id"] if 'none' not in ses["site_id"] else '',
                        ses["subject_id"],
                        ses["unique_id"]
                    ))

    return sublist


def collect_bids_files_configs(bids_dir, aws_input_creds=''):
    """
    :param bids_dir:
    :param aws_input_creds:
    :return:
    """

    file_paths = []
    config_dict = {}

    suffixes = ['T1w', 'bold', '_epi', 'phasediff', 'magnitude',
                'magnitude1', 'magnitude2']

    if bids_dir.lower().startswith("s3://"):
        # s3 paths begin with s3://bucket/
        bucket_name = bids_dir.split('/')[2]
        s3_prefix = '/'.join(bids_dir.split('/')[:3])
        prefix = bids_dir.replace(s3_prefix, '').lstrip('/')

        if aws_input_creds:
            if not os.path.isfile(aws_input_creds):
                raise IOError("Could not find aws_input_creds (%s)" %
                              (aws_input_creds))

        from indi_aws import fetch_creds
        bucket = fetch_creds.return_bucket(aws_input_creds, bucket_name)

        print(f"gathering files from S3 bucket ({bucket}) for {prefix}")

        for s3_obj in bucket.objects.filter(Prefix=prefix):
            for suf in suffixes:
                if suf in str(s3_obj.key):
                    if str(s3_obj.key).endswith("json"):
                        try:
                            config_dict[s3_obj.key.replace(prefix, "").lstrip('/')] \
                                = json.loads(s3_obj.get()["Body"].read())
                        except Exception as e:
                            print("Error retrieving %s (%s)" %
                                  (s3_obj.key.replace(prefix, ""),
                                  e.message))
                            raise
                    elif 'nii' in str(s3_obj.key):
                        file_paths.append(str(s3_obj.key)
                                          .replace(prefix,'').lstrip('/'))

    else:
        for root, dirs, files in os.walk(bids_dir, topdown=False):
            if files:
                for f in files:
                    for suf in suffixes:
                        if 'nii' in f and suf in f:
                            file_paths += [os.path.join(root, f).replace(bids_dir,'')
                                   .lstrip('/')]
                        if f.endswith('json') and suf in f:
                            try:
                                config_dict.update(
                                    {os.path.join(root.replace(bids_dir, '').lstrip('/'), f):
                                         json.load(open(os.path.join(root, f), 'r'))})
                            except UnicodeDecodeError:
                                raise Exception("Could not decode {0}".format(os.path.join(root, f)))

    if not file_paths and not config_dict:
        raise IOError("Didn't find any files in {0}. Please verify that the "
                      "path is typed correctly, that you have read access to "
                      "the directory, and that it is not "
                      "empty.".format(bids_dir))

    return file_paths, config_dict


def test_gen_bids_sublist(bids_dir, test_yml, creds_path, dbg=False):

    (img_files, config) = collect_bids_files_configs(bids_dir, creds_path)
    print("Found %d config files for %d image files" % (len(config),
                                                        len(img_files)))

    sublist = bids_gen_cpac_sublist(bids_dir, img_files, config, creds_path, dbg)

    with open(test_yml, "w") as ofd:
        yaml.dump(sublist, ofd, encoding='utf-8')

    sublist = bids_gen_cpac_sublist(bids_dir, img_files, None, creds_path, dbg)

    test_yml = test_yml.replace(".yml","_no_param.yml")
    with open(test_yml, "w") as ofd:
        yaml.dump(sublist, ofd, encoding='utf-8')

    assert sublist

if __name__ == '__main__':

    test_gen_bids_sublist(
        "/Users/cameron.craddock/workspace/git_temp/CPAC"
        "/data/ADHD200/RawDataBIDS/",
        "/Users/cameron.craddock/workspace/git_temp/CPAC"
        "/test/rs_subject_list.yml",
        "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
        dbg=False)

    test_gen_bids_sublist(
        "/Users/cameron.craddock/workspace/git_temp/CPAC"
        "/data/ADHD200/RawDataBIDS/Peking_3",
        "/Users/cameron.craddock/workspace/git_temp/CPAC"
        "/test/rs_subject_list_pk3.yml",
        "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
        dbg=False)

    test_gen_bids_sublist(
        "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/",
        "/Users/cameron.craddock/workspace/git_temp/CPAC/test/"
           "rs_subject_list_s3.yml",
        "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
        dbg=False)

    test_gen_bids_sublist(
        "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/Peking_3",
        "/Users/cameron.craddock/workspace/git_temp/CPAC/test/"
           "rs_subject_list_pk3_s3.yml",
        "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
        dbg=False)
