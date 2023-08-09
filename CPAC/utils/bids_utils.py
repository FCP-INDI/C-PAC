# Copyright (C) 2016-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
import json
import os
import re
import sys
from warnings import warn

import yaml


def bids_decode_fname(file_path, dbg=False, raise_error=True):
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

    f_dict["site"] = re.sub(r'[\s\-\_]+', '', f_dict["site"])

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


def bids_entities_from_filename(filename):
    """Function to collect a list of BIDS entities from a given
    filename.

    Parameters
    ----------
    filename : str

    Returns
    -------
    entities : list

    Examples
    --------
    >>> bids_entities_from_filename(
    ...     's3://fake/data/sub-0001/ses-NFB3/func/'
    ...     'sub-0001_ses-NFB3_task-MSIT_bold.nii.gz')
    ['sub-0001', 'ses-NFB3', 'task-MSIT', 'bold']
    """
    return (
        filename.split('/')[-1] if '/' in filename else filename
    ).split('.')[0].split('_')


def bids_match_entities(file_list, entities, suffix):
    """Function to subset a list of filepaths by a passed BIDS entity.

    Parameters
    ----------
    file_list : list of str

    entities : str
        BIDS entities joined by underscores (e.g., 'ses-001_task-PEER1')

    suffix : str
        BIDS suffix (e.g., 'bold', 'T1w')

    Returns
    -------
    list of str

    Examples
    --------
    >>> bids_match_entities([
    ...     's3://fake/data/sub-001_ses-001_task-MSIT_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_task-PEER1_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_task-PEER2_bold.nii.gz'
    ... ], 'task-PEER1', 'bold')
    ['s3://fake/data/sub-001_ses-001_task-PEER1_bold.nii.gz']
    >>> bids_match_entities([
    ...     's3://fake/data/sub-001_ses-001_task-PEER1_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_task-PEER2_bold.nii.gz'
    ... ], 'PEER', 'bold')
    Traceback (most recent call last):
    LookupError: No match found for provided entity "PEER" in
    - s3://fake/data/sub-001_ses-001_task-PEER1_bold.nii.gz
    - s3://fake/data/sub-001_ses-001_task-PEER2_bold.nii.gz
    Perhaps you meant one of these?
    - task-PEER1
    - task-PEER2
    """
    matches = [
        file for file in file_list if (
            f'_{entities}_' in '_'.join(
                bids_entities_from_filename(file)
            ) and bids_entities_from_filename(file)[-1] == suffix
        ) or bids_entities_from_filename(file)[-1] != suffix
    ]
    if file_list and not matches:
        pp_file_list = '\n'.join([f'- {file}' for file in file_list])
        error_message = ' '.join([
            'No match found for provided',
            'entity' if len(entities.split('_')) == 1 else 'entities',
            f'"{entities}" in\n{pp_file_list}'
        ])
        partial_matches = [match.group() for match in [
            re.search(re.compile(f'[^_]*{entities}[^_]*'), file) for
            file in file_list
        ] if match is not None]
        if partial_matches:
            if len(partial_matches) == 1:
                error_message += f'\nPerhaps you meant "{partial_matches[0]}"?'
            else:
                error_message = '\n'.join([
                    error_message,
                    'Perhaps you meant one of these?',
                    *[f'- {match}' for match in partial_matches]
                ])
        raise LookupError(error_message)
    return matches


def bids_remove_entity(name, key):
    """Remove an entity from a BIDS string by key

    Parameters
    ----------
    name : str
        BIDS string to remove entity from
    key : str
        BIDS key of entity to remove

    Returns
    -------
    str
        BIDS name with entity removed

    Examples
    --------
    >>> bids_remove_entity('atlas-Yeo_space-MNI152NLin6_res-2x2x2', 'space')
    'atlas-Yeo_res-2x2x2'
    >>> bids_remove_entity('atlas-Yeo_space-MNI152NLin6_res-2x2x2', 'res')
    'atlas-Yeo_space-MNI152NLin6'
    """
    return '_'.join(entity for entity in bids_entities_from_filename(name)
                    if not entity.startswith(f'{key.rstrip("-")}-'))


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
    # sidecar files

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


def bids_shortest_entity(file_list):
    """Function to return the single file with the shortest chain of
    BIDS entities from a given list, returning the first if more than
    one have the same minimum length.

    Parameters
    ----------
    file_list : list of strings

    Returns
    -------
    str or None

    Examples
    --------
    >>> bids_shortest_entity([
    ...     's3://fake/data/sub-001_ses-001_task-MSIT_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_task-PEER1_bold.nii.gz',
    ...     's3://fake/data/sub-001_ses-001_task-PEER2_bold.nii.gz'
    ... ])
    's3://fake/data/sub-001_ses-001_bold.nii.gz'
    """
    entity_lists = [
        bids_entities_from_filename(filename) for filename in file_list
    ]

    if not entity_lists:
        return None

    shortest_len = min(len(entity_list) for entity_list in entity_lists)

    shortest_list = [
        file_list[i] for i in range(len(file_list)) if
        len(entity_lists[i]) == shortest_len
    ]

    return shortest_list[0] if len(shortest_list) == 1 else shortest_list


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


def bids_gen_cpac_sublist(bids_dir, paths_list, config_dict, creds_path,
                          dbg=False, raise_error=True, only_one_anat=True):
    """
    Generates a CPAC formatted subject list from information contained in a
    BIDS formatted set of data.

    Parameters
    ----------
    bids_dir : str
        base directory that contains all of the data, this could be a
        directory that contains data for a multiple BIDS datasets, in
        which case the intervening directories will be interpreted as
        site names

    paths_list : str
        lists of all nifti files found in bids_dir, these paths are
        relative to bids_dir

    config_dict : dict
        dictionary that contains information from the JSON sidecars
        found in bids_dir, keys are relative paths and values are
        dictionaries containing all of the parameter information. if
        config_dict is None, the subject list will be built without the
        parameters

    creds_path : str
        if using S3 bucket, this path credentials needed to access the
        bucket, if accessing anonymous bucket, this can be set to None

    dbg : bool
        indicating whether or not the debug statements should be
        printed

    raise_error : bool

    only_one_anat : bool
        The "anat" key for a subject expects a string value, but we can
        temporarily store a list instead by passing True here if we
        will be filtering that list down to a single string later

    Returns
    -------
    list
        a list of dictionaries suitable for use by CPAC to specify data
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
        bids_config_dict = bids_parse_sidecar(config_dict,
                                              raise_error=raise_error)

    subdict = {}
    for p in paths_list:
        if bids_dir in p:
            str_list = p.split(bids_dir)
            val = str_list[0]
            val = val.rsplit('/')
            val = val[0]
        else:
            str_list = p.split('/')
            val = str_list[0]

        if 'sub-' not in val:
            continue

        p = p.rstrip()
        f = os.path.basename(p)

        if f.endswith(".nii") or f.endswith(".nii.gz"):

            f_dict = bids_decode_fname(p, raise_error=raise_error)

            if config_dict:
                t_params = bids_retrieve_params(bids_config_dict,
                                                f_dict)
                if not t_params:
                    print("Did not receive any parameters for %s," % (p) +
                          " is this a problem?")

                task_info = {"scan": os.path.join(bids_dir, p),
                             "scan_parameters": t_params.copy()}
            else:
                task_info = os.path.join(bids_dir, p)

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

            if "T1w" in f_dict["scantype"] or "T2w" in f_dict["scantype"]:
                if "lesion" in f_dict.keys() and "mask" in f_dict['lesion']:
                    if "lesion_mask" not in \
                            subdict[f_dict["sub"]][f_dict["ses"]]:
                        subdict[f_dict["sub"]][f_dict["ses"]]["lesion_mask"] = \
                            task_info["scan"]
                    else:
                        print("Lesion mask file (%s) already found" %
                              (subdict[f_dict["sub"]]
                               [f_dict["ses"]]
                               ["lesion_mask"]) +
                              " for (%s:%s) discarding %s" %
                              (f_dict["sub"], f_dict["ses"], p))
                # TODO deal with scan parameters anatomical
                if "anat" not in subdict[f_dict["sub"]][f_dict["ses"]]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["anat"] = {}

                if f_dict["scantype"] not in subdict[f_dict["sub"]][
                    f_dict["ses"]
                ]["anat"]:
                    if only_one_anat:
                        subdict[f_dict["sub"]][f_dict["ses"]]["anat"][
                            f_dict["scantype"]
                        ] = task_info["scan"] if config_dict else task_info
                    else:
                        subdict[f_dict["sub"]][f_dict["ses"]]["anat"][
                            f_dict["scantype"]] = []
                if not only_one_anat:
                    subdict[f_dict["sub"]][f_dict["ses"]]["anat"][
                        f_dict["scantype"]].append(
                            task_info["scan"] if config_dict else task_info)

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

            if "phase" in f_dict["scantype"]:
                if "fmap" not in subdict[f_dict["sub"]][f_dict["ses"]]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["fmap"] = {}
                if f_dict["scantype"] not in subdict[f_dict["sub"]][f_dict["ses"]]["fmap"]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["fmap"][f_dict["scantype"]] = task_info

            if "magnitude" in f_dict["scantype"]:
                if "fmap" not in subdict[f_dict["sub"]][f_dict["ses"]]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["fmap"] = {}
                if f_dict["scantype"] not in subdict[f_dict["sub"]][f_dict["ses"]]["fmap"]:
                    subdict[f_dict["sub"]][f_dict["ses"]]["fmap"][f_dict["scantype"]] = task_info

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
            if "anat" in ses or "func" in ses:
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

    suffixes = ['T1w', 'T2w', 'bold', 'epi', 'phasediff', 'phase1',
                'phase2', 'magnitude', 'magnitude1', 'magnitude2']

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
                    if suf == 'epi' and 'acq-fMRI' not in s3_obj.key:
                        continue
                    if str(s3_obj.key).endswith("json"):
                        try:
                            config_dict[s3_obj.key.replace(prefix, "")
                                        .lstrip('/')] = json.loads(
                                            s3_obj.get()["Body"].read())
                        except Exception as e:
                            print("Error retrieving %s (%s)" %
                                  (s3_obj.key.replace(prefix, ""),
                                  e.message))
                            raise
                    elif 'nii' in str(s3_obj.key):
                        file_paths.append(str(s3_obj.key)
                                          .replace(prefix,'').lstrip('/'))

    else:
        for root, dirs, files in os.walk(bids_dir, topdown=False, followlinks=True):
            if files:
                for f in files:
                    for suf in suffixes:
                        if suf == 'epi' and 'acq-fMRI' not in f:
                            continue
                        if 'nii' in f and suf in f:
                            file_paths += [os.path.join(root, f)
                                           .replace(bids_dir, '').lstrip('/')]
                        if f.endswith('json') and suf in f:
                            try:
                                config_dict.update(
                                    {os.path.join(root.replace(bids_dir, '')
                                     .lstrip('/'), f):
                                         json.load(
                                             open(os.path.join(root, f), 'r')
                                         )})
                            except UnicodeDecodeError:
                                raise Exception("Could not decode {0}".format(
                                    os.path.join(root, f)))

    if not file_paths and not config_dict:
        raise IOError("Didn't find any files in {0}. Please verify that the "
                      "path is typed correctly, that you have read access to "
                      "the directory, and that it is not "
                      "empty.".format(bids_dir))

    return file_paths, config_dict


def camelCase(string: str) -> str:  # pylint: disable=invalid-name
    """Convert a hyphenated string to camelCase

    Parameters
    ----------
    string : str
        string to convert to camelCase

    Returns
    -------
    str

    Examples
    --------
    >>> camelCase('PearsonNilearn-aCompCor')
    'PearsonNilearnACompCor'
    >>> camelCase('mean-Pearson-Nilearn-aCompCor')
    'meanPearsonNilearnACompCor'
    """
    pieces = string.split('-')
    for i in range(1, len(pieces)):  # don't change case of first piece
        if pieces[i]:  # don't do anything to falsy pieces
            pieces[i] = f'{pieces[i][0].upper()}{pieces[i][1:]}'
    return ''.join(pieces)


def combine_multiple_entity_instances(bids_str: str) -> str:
    """Combines mutliple instances of a key in a BIDS string to a single
    instance by camelCasing and concatenating the values

    Parameters
    ----------
    bids_str : str

    Returns
    -------
    str

    Examples
    --------
    >>> combine_multiple_entity_instances(
    ...     'sub-1_ses-HBN_site-RU_task-rest_atlas-AAL_'
    ...     'desc-Nilearn_desc-36-param_suffix.ext')
    'sub-1_ses-HBN_site-RU_task-rest_atlas-AAL_desc-Nilearn36Param_suffix.ext'
    >>> combine_multiple_entity_instances(
    ...     'sub-1_ses-HBN_site-RU_task-rest_'
    ...     'run-1_framewise-displacement-power.1D')
    'sub-1_ses-HBN_site-RU_task-rest_run-1_framewiseDisplacementPower.1D'
    """
    _entity_list = bids_str.split('_')
    entity_list = _entity_list[:-1]
    suffixes = [camelCase(_entity_list[-1])]
    entities = {}
    for entity in entity_list:
        if '-' in entity:
            key, value = entity.split('-', maxsplit=1)
            if key not in entities:
                entities[key] = []
            entities[key].append(value)
    for key, value in entities.items():
        entities[key] = camelCase('-'.join(value))
    if 'desc' in entities:  # make 'desc' final entity
        suffixes.insert(0, f'desc-{entities.pop("desc")}')
    return '_'.join([f'{key}-{value}' for key, value in entities.items()
                     ] + suffixes)


def insert_entity(resource, key, value):
    """Insert a `f'{key}-{value}'` BIDS entity before `desc-` if
    present or before the suffix otherwise

    Parameters
    ----------
    resource, key, value : str

    Returns
    -------
    str

    Examples
    --------
    >>> insert_entity('run-1_desc-preproc_bold', 'reg', 'default')
    'run-1_reg-default_desc-preproc_bold'
    >>> insert_entity('run-1_bold', 'reg', 'default')
    'run-1_reg-default_bold'
    >>> insert_entity('run-1_desc-preproc_bold', 'filt', 'notch4c0p31bw0p12')
    'run-1_filt-notch4c0p31bw0p12_desc-preproc_bold'
    >>> insert_entity('run-1_reg-default_bold', 'filt', 'notch4c0p31bw0p12')
    'run-1_reg-default_filt-notch4c0p31bw0p12_bold'
    """
    entities = resource.split('_')[:-1]
    suff = resource.split('_')[-1]
    new_entities = [[], []]
    for entity in entities:
        if entity.startswith('desc-'):
            new_entities[1].append(entity)
        else:
            new_entities[0].append(entity)
    return '_'.join([*new_entities[0], f'{key}-{value}', *new_entities[1],
                     suff])


def load_yaml_config(config_filename, aws_input_creds):

    if config_filename.lower().startswith('data:'):
        try:
            header, encoded = config_filename.split(",", 1)
            config_content = b64decode(encoded)
            config_data = yaml.safe_load(config_content)
            return config_data
        except:
            print("Error! Could not find load config from data URI")
            raise

    if config_filename.lower().startswith("s3://"):
        # s3 paths begin with s3://bucket/
        bucket_name = config_filename.split('/')[2]
        s3_prefix = '/'.join(config_filename.split('/')[:3])
        prefix = config_filename.replace(s3_prefix, '').lstrip('/')

        if aws_input_creds:
            if not os.path.isfile(aws_input_creds):
                raise IOError("Could not find aws_input_creds (%s)" %
                              (aws_input_creds))

        from indi_aws import fetch_creds
        bucket = fetch_creds.return_bucket(aws_input_creds, bucket_name)
        downloaded_config = '/tmp/' + os.path.basename(config_filename)
        bucket.download_file(prefix, downloaded_config)
        config_filename = downloaded_config

    config_filename = os.path.realpath(config_filename)

    try:
        config_data = yaml.safe_load(open(config_filename, 'r'))
        return config_data
    except IOError:
        print("Error! Could not find config file {0}".format(config_filename))
        raise


def cl_strip_brackets(arg_list):
    """Removes '[' from before first and ']' from after final
    arguments in a list of commandline arguments

    Parameters
    ----------
    arg_list : list

    Returns
    -------
    list

    Examples
    --------
    >>> cl_strip_brackets('[a b c]'.split(' '))
    ['a', 'b', 'c']
    >>> cl_strip_brackets('a b c'.split(' '))
    ['a', 'b', 'c']
    >>> cl_strip_brackets('[ a b c ]'.split(' '))
    ['a', 'b', 'c']
    """
    arg_list[0] = arg_list[0].lstrip('[')
    arg_list[-1] = arg_list[-1].rstrip(']')
    return [arg for arg in arg_list if arg]


def create_cpac_data_config(bids_dir, participant_labels=None,
                            aws_input_creds=None, skip_bids_validator=False,
                            only_one_anat=True):
    """
    Create a C-PAC data config YAML file from a BIDS directory.

    Parameters
    ----------
    bids_dir : str

    participant_labels : list or None

    aws_input_creds

    skip_bids_validator : bool

    only_one_anat : bool
        The "anat" key for a subject expects a string value, but we
        can temporarily store a list instead by passing True here if
        we will be filtering that list down to a single string later

    Returns
    -------
    list
    """
    print("Parsing {0}..".format(bids_dir))

    (file_paths, config) = collect_bids_files_configs(bids_dir,
                                                      aws_input_creds)

    if participant_labels and file_paths:
        file_paths = [
            file_path for file_path in file_paths if any(
                participant_label in file_path
                for participant_label in participant_labels
            )
        ]

    if not file_paths:
        print("Did not find data for {0}".format(
            ", ".join(participant_labels)
        ))
        sys.exit(1)

    raise_error = not skip_bids_validator

    sub_list = bids_gen_cpac_sublist(
        bids_dir,
        file_paths,
        config,
        aws_input_creds,
        raise_error=raise_error,
        only_one_anat=only_one_anat
    )

    if not sub_list:
        print("Did not find data in {0}".format(bids_dir))
        sys.exit(1)

    return sub_list


def load_cpac_data_config(data_config_file, participant_labels,
                          aws_input_creds):
    """
    Loads the file as a check to make sure it is available and readable

    Parameters
    ----------
    data_config_file : str
        path to data config

    participants_labels : list or None

    aws_input_creds

    Returns
    -------
    list
    """
    sub_list = load_yaml_config(data_config_file, aws_input_creds)

    if participant_labels:

        sub_list = [
            d
            for d in sub_list
            if (
                d["subject_id"]
                if d["subject_id"].startswith('sub-')
                else 'sub-' + d["subject_id"]
            ) in participant_labels
        ]

        if not sub_list:
            print("Did not find data for {0} in {1}".format(
                ", ".join(participant_labels),
                (
                    data_config_file
                    if not data_config_file.startswith("data:")
                    else "data URI"
                )
            ))
            sys.exit(1)

    return sub_list


def res_in_filename(cfg, label):
    """Specify resolution in filename

    Parameters
    ----------
    cfg : CPAC.utils.configuration.Configuration

    label : str

    Returns
    -------
    label : str

    Examples
    --------
    >>> from CPAC.utils.configuration import Configuration
    >>> res_in_filename(Configuration({
    ...     'registration_workflows': {
    ...         'anatomical_registration': {'resolution_for_anat': '2x2x2'}}}),
    ...     'sub-1_res-anat_bold')
    'sub-1_res-2x2x2_bold'
    >>> res_in_filename(Configuration({
    ...     'registration_workflows': {
    ...         'anatomical_registration': {'resolution_for_anat': '2x2x2'}}}),
    ...     'sub-1_res-3mm_bold')
    'sub-1_res-3mm_bold'
    """
    if '_res-' in label:
        # replace resolution text with actual resolution
        resolution = label.split('_res-', 1)[1].split('_', 1)[0]
        resolution = {
            'anat': cfg['registration_workflows', 'anatomical_registration',
                        'resolution_for_anat'],
            'bold': cfg['registration_workflows', 'functional_registration',
                        'func_registration_to_template', 'output_resolution',
                        'func_preproc_outputs'],
            'derivative': cfg['registration_workflows',
                              'functional_registration',
                              'func_registration_to_template',
                              'output_resolution', 'func_derivative_outputs']
        }.get(resolution, resolution)
        label = re.sub('_res-[A-Za-z0-9]*_', f'_res-{resolution}_', label)
    return label


def sub_list_filter_by_labels(sub_list, labels):
    """Function to filter a sub_list by provided BIDS labels for
    specified suffixes

    Parameters
    ----------
    sub_list : list

    labels : dict

    labels['T1w'] : str or None
        C-PAC currently only uses a single T1w image

    labels['bold'] : str, list, or None

    Returns
    -------
    list
    """
    if labels.get('T1w'):
        sub_list = _sub_list_filter_by_label(sub_list, 'T1w', labels['T1w'])
    if labels.get('bold'):
        labels['bold'] = cl_strip_brackets(labels['bold'])
        sub_list = _sub_list_filter_by_label(sub_list, 'bold', labels['bold'])
    return sub_list


def with_key(entity: str, key: str) -> str:
    """Return a keyed BIDS entity

    Parameters
    ----------
    entity, key : str

    Returns
    -------
    str

    Examples
    --------
    >>> with_key('sub-1', 'sub')
    'sub-1'
    >>> with_key('1', 'sub')
    'sub-1'
    """
    if not isinstance(entity, str):
        entity = str(entity)
    if not entity.startswith(f'{key}-'):
        entity = '-'.join((key, entity))
    return entity


def without_key(entity: str, key: str) -> str:
    """Return a BIDS entity value

    Parameters
    ----------
    entity, key : str

    Returns
    -------
    str

    Examples
    --------
    >>> without_key('sub-1', 'sub')
    '1'
    >>> without_key('1', 'sub')
    '1'
    """
    if not isinstance(entity, str):
        entity = str(entity)
    if entity.startswith(f'{key}-'):
        entity = entity.replace(f'{key}-', '')
    return entity


def _t1w_filter(anat, shortest_entity, label):
    """Helper function to filter T1w paths

    Parameters
    ----------
    anat: list or str

    shortest_entity: bool

    label: str

    Returns
    -------
    anat: list
    """
    if not isinstance(anat, list):
        anat = [anat]
    if shortest_entity:
        anat = bids_shortest_entity(anat)
    else:
        anat = bids_match_entities(anat, label, 'T1w')
        # pylint: disable=invalid-name
        try:
            anat_T2 = bids_match_entities(anat, label, 'T2w')
        except LookupError:
            anat_T2 = None
        if anat_T2 is not None:
            anat = anat_T2
    return anat


def _sub_anat_filter(anat, shortest_entity, label):
    """Helper function to filter anat paths in sub_list

    Parameters
    ----------
    anat : list or dict

    shortest_entity : bool

    label : str

    Returns
    -------
    list or dict
        same type as 'anat' parameter
    """
    if isinstance(anat, dict):
        if 'T1w' in anat:
            anat['T1w'] = _t1w_filter(anat['T1w'],
                                      shortest_entity,
                                      label)
        return anat
    return _t1w_filter(anat, shortest_entity, label)


def _sub_list_filter_by_label(sub_list, label_type, label):
    """Function to filter a sub_list by a CLI-provided label.

    Parameters
    ----------
    sub_list : list

    label_type : str
        'T1w' or 'bold'

    label : str or list

    Returns
    -------
    list

    Examples
    --------
    >>> from CPAC.pipeline.test.sample_data import sub_list
    >>> _sub_list_filter_by_label(sub_list, 'bold', 'task-PEER1')[
    ...     0]['func'].keys()
    dict_keys(['PEER1'])
    """
    label_list = [label] if isinstance(label, str) else list(label)
    new_sub_list = []
    if label_type in label_list:
        shortest_entity = True
        label_list.remove(label_type)
    else:
        shortest_entity = False
    if label_type == 'T1w':
        for sub in [sub for sub in sub_list if 'anat' in sub]:
            try:
                sub['anat'] = _sub_anat_filter(sub['anat'],
                                               shortest_entity,
                                               label_list[0] if not
                                               shortest_entity else None)
                if sub['anat']:
                    new_sub_list.append(sub)
            except LookupError as lookup_error:
                warn(str(lookup_error))

    elif label_type == 'bold':
        for sub in [sub for sub in sub_list if 'func' in sub]:
            try:
                all_scans = [sub['func'][scan].get('scan') for
                             scan in sub['func']]
                new_func = {}
                for entities in label_list:
                    matched_scans = bids_match_entities(all_scans, entities,
                                                        label_type)
                    for scan in matched_scans:
                        new_func = {
                            **new_func,
                            **_match_functional_scan(sub['func'], scan)
                        }
                if shortest_entity:
                    new_func = {
                        **new_func,
                        **_match_functional_scan(
                            sub['func'], bids_shortest_entity(all_scans)
                        )
                    }
                sub['func'] = new_func
                new_sub_list.append(sub)
            except LookupError as lookup_error:
                warn(str(lookup_error))
    return new_sub_list


def _match_functional_scan(sub_list_func_dict, scan_file_to_match):
    """Function to subset a scan from a sub_list_func_dict by a scan filename

    Parameters
    ---------
    sub_list_func_dict : dict
        sub_list[sub]['func']

    scan_file_to_match : str

    Returns
    -------
    dict

    Examples
    --------
    >>> from CPAC.pipeline.test.sample_data import sub_list
    >>> matched = _match_functional_scan(
    ...     sub_list[0]['func'],
    ...     '/fake/data/sub-0001/ses-NFB3/func/'
    ...     'sub-0001_ses-NFB3_task-PEER1_bold.nii.gz')
    >>> matched.keys()
    dict_keys(['PEER1'])
    >>> all([key in matched['PEER1'] for key in [
    ...     'fmap_mag', 'fmap_phase', 'scan', 'scan_parameters'
    ... ]])
    True
    """
    return {
        entity: sub_list_func_dict[entity] for entity in
        sub_list_func_dict if
        sub_list_func_dict[entity].get('scan') == scan_file_to_match
    }
