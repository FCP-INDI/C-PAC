# Copyright (C) 2022-2024  C-PAC Developers

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
"""Run a group-level analysis with C-PAC."""

import fnmatch
import os

from CPAC.utils.monitoring import WFLOGGER


def load_config_yml(config_file: str) -> dict:
    """Load a configuration YAML file.

    Parameters
    ----------
    config_file : str
        full filepath to YAML (.yml) file

    Returns
    -------
    config : dict
    """
    import os

    import yaml

    try:
        config_path = os.path.realpath(config_file)
        config_dict = yaml.safe_load(open(config_path, "r"))
        config = config_dict

    except (OSError, TypeError, yaml.YAMLError) as e:
        err = (
            "\n\n[!] CPAC says: Could not load or read the configuration "
            f"YAML file:\n{config_file}\n"
        )
        raise yaml.YAMLError(err) from e

    return config


def load_text_file(filepath: str, label: str = "file") -> list[str]:
    """Load a text file and returns the lines in a list.

    Parameters
    ----------
    filepath : str
        full filepath to the text file

    label : str

    Returns
    -------
    lines_list : list of str
        list of lines from text file
    """
    if not filepath.endswith(".txt"):
        err = (
            "\n\n[!] CPAC says: The %s should be a text file (.txt).\n"
            "Path provided: %s\n\n" % (label, filepath)
        )
        raise Exception(err)

    try:
        with open(filepath, "r") as f:
            lines_list = f.readlines()
    except Exception as e:
        err = (
            "\n\n[!] CPAC says: Could not load or read the %s:\n%s\n"
            "Details: %s\n\n" % (label, filepath, e)
        )
        raise Exception(err)

    # get rid of those \n's that love to show up everywhere
    return [i.rstrip("\n") for i in lines_list]


def grab_pipeline_dir_subs(pipeline_dir, ses=False):
    """Grab the subject IDs from the pipeline output directory."""
    import os

    inclusion_list = []
    if ses:
        pipeline_list = [
            x
            for x in os.listdir(pipeline_dir)
            if os.path.isdir(os.path.join(pipeline_dir, x))
        ]
    else:
        pipeline_list = [
            x.split("_")[0]
            for x in os.listdir(pipeline_dir)
            if os.path.isdir(os.path.join(pipeline_dir, x))
        ]
    for sub_id in pipeline_list:
        if sub_id not in inclusion_list:
            inclusion_list.append(sub_id)
    return sorted(inclusion_list)


def read_pheno_csv_into_df(pheno_csv, id_label=None):
    """Read the phenotypic file CSV or TSV into a Pandas DataFrame."""
    import pandas as pd

    with open(pheno_csv, "r") as f:
        if id_label:
            if ".tsv" in pheno_csv or ".TSV" in pheno_csv:
                pheno_df = pd.read_table(f, dtype={id_label: object})
            else:
                pheno_df = pd.read_csv(f, dtype={id_label: object})
        elif ".tsv" in pheno_csv or ".TSV" in pheno_csv:
            pheno_df = pd.read_table(f)
        else:
            pheno_df = pd.read_csv(f)

    return pheno_df


def gather_nifti_globs(pipeline_output_folder, resource_list, pull_func=False):
    """Gather the NIFTI file globs for the derivatives selected.

    The number of directory levels under each participant's output folder
    can vary depending on what preprocessing strategies were chosen, and
    there may be several output filepaths with varying numbers of directory
    levels.

    This parses them quickly while also catching each preprocessing strategy.
    """
    import glob
    import os

    import pandas as pd
    import pkg_resources as p

    exts = ".nii"
    nifti_globs = []

    keys_tsv = p.resource_filename("CPAC", "resources/cpac_outputs.tsv")
    try:
        keys = pd.read_csv(keys_tsv, delimiter="\t")
    except Exception as e:
        err = (
            "\n[!] Could not access or read the cpac_outputs.tsv "
            f"resource file:\n{keys_tsv}\n\nError details {e}\n"
        )
        raise Exception(err)

    derivative_list = list(keys[keys["Sub-Directory"] == "func"]["Resource"])
    derivative_list = derivative_list + list(
        keys[keys["Sub-Directory"] == "anat"]["Resource"]
    )

    if pull_func:
        derivative_list = derivative_list + list(
            keys[keys["Space"] == "functional"]["Resource"]
        )

    if len(resource_list) == 0:
        err = "\n\n[!] No derivatives selected!\n\n"
        raise Exception(err)

    # remove any extra /'s
    pipeline_output_folder = pipeline_output_folder.rstrip("/")

    WFLOGGER.info(
        "\n\nGathering the output file paths from %s...", pipeline_output_folder
    )

    # this is just to keep the fsl feat config file derivative_list entries
    # nice and lean
    dirs_to_grab = []
    for derivative_name in derivative_list:
        for resource_name in resource_list:
            if resource_name in derivative_name:
                dirs_to_grab.append(derivative_name)

    # grab MeanFD_Jenkinson just in case
    dirs_to_grab.append("framewise-displacement-jenkinson")

    for resource_name in dirs_to_grab:
        glob_string = os.path.join(
            pipeline_output_folder, "*", "*", f"*{resource_name}*"
        )
        # get all glob strings that result in a list of paths where every path
        # ends with a NIFTI file
        prog_string = ".."
        while len(glob.glob(glob_string)) != 0:
            if any(exts in x for x in glob.glob(glob_string)):
                nifti_globs.append(glob_string)

            glob_string = os.path.join(glob_string, "*")
            prog_string = prog_string + "."

    if len(nifti_globs) == 0:
        err = (
            "\n\n[!] No output filepaths found in the pipeline output "
            "directory provided for the derivatives selected!\n\nPipeline "
            "output directory provided: %s\nDerivatives selected:%s\n\n"
            % (pipeline_output_folder, resource_list)
        )
        raise Exception(err)

    return nifti_globs


def grab_raw_score_filepath(filepath, resource_id):
    """Grab the filepath for the raw score of the resource."""
    # this lives in the output path collector

    import glob
    import os

    if "vmhc" in resource_id:
        raw_score_path = filepath.replace(resource_id, "vmhc_raw_score")
        raw_score_path = raw_score_path.replace(raw_score_path.split("/")[-1], "")
        raw_score_path = glob.glob(os.path.join(raw_score_path, "*"))[0]
    else:
        raw_score_path = filepath.replace("_zstd", "")
        raw_score_path = raw_score_path.replace("_fisher", "")
        raw_score_path = raw_score_path.replace("_zstat", "")

        if "sca_roi_files_to_standard" in resource_id:
            sub_folder = raw_score_path.split("/")[-2] + "/"
            if "z_score" in sub_folder:
                raw_score_path = raw_score_path.replace(sub_folder, "")
        elif "sca_tempreg_maps_zstat" in resource_id:
            sca_filename = raw_score_path.split("/")[-1]
            globpath = raw_score_path.replace(sca_filename, "*")
            globpath = os.path.join(globpath, sca_filename)
            raw_score_path = glob.glob(globpath)[0]
        elif "dr_tempreg_maps" in resource_id:
            raw_score_path = raw_score_path.replace("map_z_", "map_")
            raw_filename = raw_score_path.split("/")[-1]
            raw_score_path = raw_score_path.replace(raw_filename, "")
            raw_score_path = glob.glob(os.path.join(raw_score_path, "*", raw_filename))[
                0
            ]
        else:
            # in case filenames are different between z-standardized and raw
            raw_score_path = raw_score_path.replace(raw_score_path.split("/")[-1], "")
            try:
                raw_score_path = glob.glob(os.path.join(raw_score_path, "*"))[0]
            except (FileNotFoundError, IndexError, TypeError):
                raw_score_path = os.path.join(raw_score_path, "*")

    if (raw_score_path is None) or (not os.path.exists(raw_score_path)):
        err = (
            "\n\n[!] The filepath for the raw score of "
            "%s can not be found.\nFilepath: %s\n\nThis "
            "is needed for the Measure Mean calculation."
            "\n\n" % (resource_id, raw_score_path)
        )
        raise Exception(err)

    return raw_score_path


def find_power_params_file(filepath, resource_id, series_id):
    """Find the power parameters file for the participant and series."""
    import os

    try:
        power_path = filepath.replace(resource_id, "power_params", 1)
        series_id_string = "_scan_%s" % series_id
        power_first_half = power_path.split(series_id_string)[0]
        power_first_half = os.path.join(power_first_half, series_id_string)
        participant_id = power_first_half.split("/")[-3]
    except Exception as e:
        err = (
            "\n\n[!] Something went wrong with finding the power "
            "parameters file for at least one of the participants.\n\n"
            "Error details: %s\n\n" % e
        )
        raise Exception(err)

    power_params_file = None
    for root, dirs, files in os.walk(power_first_half):
        for filename in files:
            filepath = os.path.join(root, filename)
            if "pow_params.txt" in filepath:
                power_params_file = filepath

    if not power_params_file:
        err = (
            "\n\n[!] Could not find the power parameters file for the "
            "following participant and series..\nParticipant: %s\n"
            "Series: %s\n\nIt should be available here: %s\n\n"
            % (participant_id, series_id, power_first_half)
        )
        raise Exception(err)

    return power_params_file


def extract_power_params(power_params_lines, power_params_filepath):
    """Extract the power parameters from the power parameters file."""
    # check formatting
    if len(power_params_lines) != 2:  # noqa: PLR2004
        err = (
            "\n\n[!] There is something wrong with the formatting of the "
            f"power parameters file.\nFilepath: {power_params_filepath}\n\n"
        )
        raise ValueError(err)

    names_list = power_params_lines[0].split(",")
    values_list = power_params_lines[1].split(",")

    # let's make extra sure
    if (values_list[0].replace(" ", "") not in power_params_filepath) or (
        values_list[1].replace(" ", "") not in power_params_filepath
    ):
        err = (
            "\n\n[!] There is a mismatch between the contents of the "
            "power parameters file and where it is located!\n"
            "Filepath: %s\n\n" % power_params_filepath
        )
        raise Exception(err)

    if (
        (names_list[2].replace(" ", "") != "MeanFD_Power")
        or (names_list[3].replace(" ", "") != "MeanFD_Jenkinson")
        or (names_list[-1].replace(" ", "") != "MeanDVARS")
    ):
        err = (
            "\n\n[!] There is a mismatch between the power parameters "
            "format and what is expected!!\nFilepath: %s\n\n" % power_params_filepath
        )
        raise Exception(err)

    meanfd_power = values_list[2]
    meanfd_jenk = values_list[3]
    meandvars = values_list[-1]

    return meanfd_power, meanfd_jenk, meandvars


def create_output_dict_list(
    nifti_globs,
    pipeline_output_folder,
    resource_list,
    get_motion=False,
    get_raw_score=False,
    pull_func=False,
    derivatives=None,
    exts=["nii", "nii.gz"],
):
    """Create a dictionary of output filepaths and their associated information."""
    import os

    import pandas as pd
    import pkg_resources as p

    if len(resource_list) == 0:
        err = "\n\n[!] No derivatives selected!\n\n"
        raise Exception(err)

    if derivatives is None:
        keys_tsv = p.resource_filename("CPAC", "resources/cpac_outputs.tsv")
        try:
            keys = pd.read_csv(keys_tsv, delimiter="\t")
        except Exception as e:
            err = (
                "\n[!] Could not access or read the cpac_outputs.csv "
                f"resource file:\n{keys_tsv}\n\nError details {e}\n"
            )
            raise Exception(err)

        derivatives = list(keys[keys["Sub-Directory"] == "func"]["Resource"])
        derivatives = derivatives + list(
            keys[keys["Sub-Directory"] == "anat"]["Resource"]
        )

        if pull_func:
            derivatives = derivatives + list(
                keys[keys["Space"] == "functional"]["Resource"]
            )

    # remove any extra /'s
    pipeline_output_folder = pipeline_output_folder.rstrip("/")

    WFLOGGER.info(
        "\n\nGathering the output file paths from %s...", pipeline_output_folder
    )

    # this is just to keep the fsl feat config file derivatives entries
    # nice and lean
    search_dirs = []
    for derivative_name in derivatives:
        for resource_name in resource_list:
            if resource_name in derivative_name:
                search_dirs.append(derivative_name)

    """
    search_dirs = [
        resource_name
        for resource_name in resource_list
        if any([resource_name in derivative_name
                for derivative_name in derivatives])
    ]
    """

    # grab MeanFD_Jenkinson just in case
    search_dirs += ["framewise-displacement-jenkinson"]

    exts = ["." + ext.lstrip(".") for ext in exts]

    # parse each result of each "valid" glob string
    output_dict_list = {}

    for root, _, files in os.walk(pipeline_output_folder):
        for filename in files:
            filepath = os.path.join(root, filename)
            if not any(fnmatch.fnmatch(filepath, pattern) for pattern in nifti_globs):
                continue

            if not any(filepath.endswith(ext) for ext in exts):
                continue
            relative_filepath = filepath.split(pipeline_output_folder)[1]
            filepath_pieces = [_f for _f in relative_filepath.split("/") if _f]

            resource_id = "_".join(filepath_pieces[2].split(".")[0].split("_")[3:])

            if resource_id not in search_dirs:
                continue

            series_id_string = filepath_pieces[2].split("_")[1]
            strat_info = "_".join(filepath_pieces[2].split("_")[2:3])

            unique_resource_id = (resource_id, strat_info)

            if unique_resource_id not in output_dict_list.keys():
                output_dict_list[unique_resource_id] = []

            unique_id = filepath_pieces[0]

            series_id = series_id_string.replace("_scan_", "")
            series_id = series_id.replace("_rest", "")

            new_row_dict = {}
            new_row_dict["participant_session_id"] = unique_id
            new_row_dict["participant_id"], new_row_dict["Sessions"] = unique_id.split(
                "_"
            )

            new_row_dict["Series"] = series_id
            new_row_dict["Filepath"] = filepath

            WFLOGGER.info(
                "%s - %s - %s", unique_id.split("_")[0], series_id, resource_id
            )

            if get_motion:
                # if we're including motion measures
                power_params_file = find_power_params_file(
                    filepath, resource_id, series_id
                )
                power_params_lines = load_text_file(
                    power_params_file, "power parameters file"
                )
                meanfd_p, meanfd_j, meandvars = extract_power_params(
                    power_params_lines, power_params_file
                )
                new_row_dict["MeanFD_Power"] = meanfd_p
                new_row_dict["MeanFD_Jenkinson"] = meanfd_j
                new_row_dict["MeanDVARS"] = meandvars

            if get_raw_score:
                # grab raw score for measure mean just in case
                raw_score_path = grab_raw_score_filepath(filepath, resource_id)
                new_row_dict["Raw_Filepath"] = raw_score_path

            # unique_resource_id is tuple (resource_id,strat_info)
            output_dict_list[unique_resource_id].append(new_row_dict)

    return output_dict_list


def create_output_df_dict(output_dict_list, inclusion_list=None):
    """Create a dictionary of dataframes from the output_dict_list."""
    import pandas as pd

    output_df_dict = {}

    # unique_resource_id is tuple (resource_id,strat_info)
    for unique_resource_id in output_dict_list.keys():
        # NOTE: this dataframe reflects what was found in the C-PAC output
        #       directory for individual-level analysis outputs,
        #       NOT what is in the pheno file
        new_df = pd.DataFrame(output_dict_list[unique_resource_id])

        # drop whatever is not in the inclusion lists
        if inclusion_list:
            new_df = new_df[new_df.participant_id.isin(inclusion_list)]

        if new_df.empty:
            WFLOGGER.warning(
                "No outputs found for %s for the participants listed in the the group"
                " analysis participant list you used. Skipping generating a model for"
                " this output.",
                unique_resource_id,
            )
            continue

        # unique_resource_id is tuple (resource_id,strat_info)
        if unique_resource_id not in output_df_dict.keys():
            output_df_dict[unique_resource_id] = new_df

    return output_df_dict


def gather_outputs(
    pipeline_folder,
    resource_list,
    inclusion_list,
    get_motion,
    get_raw_score,
    get_func=False,
    derivatives=None,
):
    """Gather the output filepaths and their associated information."""
    nifti_globs = gather_nifti_globs(pipeline_folder, resource_list, get_func)

    output_dict_list = create_output_dict_list(
        nifti_globs,
        pipeline_folder,
        resource_list,
        get_motion,
        get_raw_score,
        get_func,
        derivatives,
    )

    return create_output_df_dict(output_dict_list, inclusion_list)


def pheno_sessions_to_repeated_measures(pheno_df, sessions_list):
    """Match participant-session IDs for FEAT repeated measures analysis.

    Take in the selected session names, and match them to the unique
    participant-session IDs appropriately for an FSL FEAT repeated measures
    analysis.

    More info:
      https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FEAT/
          UserGuide#Paired_Two-Group_Difference_.28Two-Sample_Paired_T-Test.29
    Sample input:
        pheno_df
          sub01
          sub02
        sessions_list
          [ses01, ses02]
    Expected output:
        pheno_df      Sessions  participant_sub01  participant_sub02
          sub01          ses01                  1                  0
          sub02          ses01                  0                  1
          sub01          ses02                  1                  0
          sub02          ses02                  0                  1
    """
    import pandas as pd

    # first, check to see if this design matrix setup has already been done
    # in the pheno CSV file
    #     NOTE: this is mainly for PRESET GROUP ANALYSIS MODELS!!!
    num_partic_cols = 0
    for col_names in pheno_df.columns:
        if "participant_" in col_names:
            num_partic_cols += 1

    if num_partic_cols > 1 and (
        "Sessions" in pheno_df.columns or "Sessions_column_one" in pheno_df.columns
    ):
        for part_id in pheno_df["participant_id"]:
            if f"participant_{part_id}" in pheno_df.columns:
                continue
            break
        else:
            # if it's already set up properly, then just send the pheno_df
            # back and bypass all the machinery below
            return pheno_df
    else:
        # if not an FSL model preset, continue as normal
        new_rows = []
        another_new_row = []

        # grab the ordered sublist before we double the rows
        sublist = pheno_df["participant_id"]

        for session in sessions_list:
            sub_pheno_df = pheno_df.copy()
            sub_pheno_df["Sessions"] = session
            sub_pheno_df["participant_session_id"] = (
                pheno_df.participant_id + "_ses-%s" % session
            )
            new_rows.append(sub_pheno_df)
            another_new_row.append(sub_pheno_df)
            pheno_df = pd.concat(new_rows)
            pheno_df = pd.concat(another_new_row)

    sessions_col = []
    part_ids_col = []

    # participant IDs new columns
    participant_id_cols = {}
    i = 0

    for participant_unique_id in pheno_df["participant_session_id"]:
        part_col = [0] * len(pheno_df["participant_session_id"])

        for session in sessions_list:
            if session in participant_unique_id.split("_")[1]:
                # print(participant_unique_id)# generate/update sessions categorical column
                part_id = participant_unique_id.split("_")[0]

                part_ids_col.append(str(part_id))
                sessions_col.append(str(session))

                header_title = "participant_%s" % part_id

                # generate/update participant ID column (1's or 0's)
                if header_title not in participant_id_cols.keys():
                    part_col[i] = 1
                    participant_id_cols[header_title] = part_col
                else:
                    participant_id_cols[header_title][i] = 1
        i += 1

    pheno_df["Sessions"] = sessions_col
    pheno_df["participant"] = part_ids_col

    # add new participant ID columns
    for sub_id in sublist:
        new_col = f"participant_{sub_id}"
        pheno_df[new_col] = participant_id_cols[new_col]

    return pheno_df.astype("object")


def pheno_series_to_repeated_measures(pheno_df, series_list, repeated_sessions=False):
    """Take in the selected series/scans, and create all of the permutations...

    ...of unique participant IDs (participant_site_session) and series/scans
    and populate the pheno

    This is so the user does not have to have a specially-formatted version of the
    phenotype CSV for repeated measures; they can just enter the regular one.
    """
    import pandas as pd

    # first, check to see if this design matrix setup has already been done
    # in the pheno CSV file
    num_partic_cols = 0
    for col_names in pheno_df.columns:
        if "participant_" in col_names:
            num_partic_cols += 1
    if num_partic_cols > 1 and "Series" in pheno_df.columns:
        for part_id in pheno_df["participant_id"]:
            if f"participant_{part_id}" in pheno_df.columns:
                continue
            break
        else:
            # if it's already set up properly, then just send the pheno_df
            # back and bypass all the machinery below
            return pheno_df

    new_rows = []
    for series in series_list:
        sub_pheno_df = pheno_df.copy()
        sub_pheno_df["Series"] = series
        new_rows.append(sub_pheno_df)
    pheno_df = pd.concat(new_rows)

    if not repeated_sessions:
        # participant IDs new columns
        participant_id_cols = {}
        i = 0

        for participant_unique_id in pheno_df["participant_id"]:
            part_col = [0] * len(pheno_df["participant_id"])
            header_title = "participant_%s" % participant_unique_id

            if header_title not in participant_id_cols.keys():
                part_col[i] = 1
                participant_id_cols[header_title] = part_col
            else:
                participant_id_cols[header_title][i] = 1

            i += 1

        for new_col in participant_id_cols.keys():
            pheno_df[new_col] = participant_id_cols[new_col]

    return pheno_df.astype("object")


def balance_repeated_measures(pheno_df, sessions_list, series_list=None):
    """Balance the repeated measures design matrix.

    This is for repeated measures only.
    If the user selects a participant list like this:
       sub001_session_1
       sub001_session_2
       sub002_session_1
       sub002_session_2
       sub003_session_1
    then have this drop "sub003_session_1", because repeated measures
    requires a uniform balance of repeats
    """
    from collections import Counter

    part_ID_count = Counter(pheno_df["participant_id"])

    sessions_x_series = len(sessions_list)
    if series_list is not None:
        sessions_x_series *= len(series_list)

    dropped_parts = []

    for part_ID in part_ID_count.keys():
        if part_ID_count[part_ID] != sessions_x_series:
            pheno_df = pheno_df[pheno_df.participant_id != part_ID]
            try:
                del pheno_df[f"participant_{part_ID}"]
            except (KeyError, TypeError):
                pass
            dropped_parts.append(part_ID)

    return pheno_df, dropped_parts


def prep_feat_inputs(group_config_file: str) -> dict:
    """Prep group analysis run.

    Parameters
    ----------
    config_file : str
        filepath to the C-PAC group-level config file
    """
    import os

    import pandas as pd
    import pkg_resources as p

    keys_tsv = p.resource_filename("CPAC", "resources/cpac_outputs.tsv")
    try:
        keys = pd.read_csv(keys_tsv, delimiter="\t")
    except Exception as e:
        err = (
            "\n[!] Could not access or read the cpac_outputs.tsv "
            f"resource file:\n{keys_tsv}\n\nError details {e}\n"
        )
        raise Exception(err)

    derivatives = list(
        keys[keys["Derivative"] == "yes"][keys["Space"] == "template"][
            keys["Values"] == "z-score"
        ]["Resource"]
    )
    derivatives = derivatives + list(
        keys[keys["Derivative"] == "yes"][keys["Space"] == "template"][
            keys["Values"] == "z-stat"
        ]["Resource"]
    )

    group_model = load_config_yml(group_config_file)
    pipeline_dir = group_model.pipeline_dir

    #   - create participant list
    #   - get output measure list
    #   - see if any of the models will require the raw scores

    get_motion = False
    get_raw_score = False

    if not group_model.participant_list:
        inclusion_list = grab_pipeline_dir_subs(pipeline_dir)
    elif "." in group_model.participant_list:
        if not os.path.isfile(group_model.participant_list):
            msg = (
                "\n[!] C-PAC says: Your participant "
                "inclusion list is not a valid file!\n\n"
                f"File path: {group_model.participant_list}"
                "\n"
            )
            if os.path.exists(group_model.participant_list):
                raise ValueError(msg)
            raise FileNotFoundError(msg)
        inclusion_list = load_text_file(
            group_model.participant_list, "group-level analysis participant list"
        )
    else:
        inclusion_list = grab_pipeline_dir_subs(pipeline_dir)

    output_measure_list = group_model.derivative_list

    # if any of the models will require motion parameters
    if ("MeanFD" in group_model.design_formula) or (
        "MeanDVARS" in group_model.design_formula
    ):
        get_motion = True

    # make sure "None" gets processed properly here...
    if group_model.custom_roi_mask in ("None", "none"):
        custom_roi_mask = None
    else:
        custom_roi_mask = group_model.custom_roi_mask

    if ("Measure_Mean" in group_model.design_formula) or (custom_roi_mask is not None):
        get_raw_score = True

    # sammin sammin mmmm samin in gray v

    # create encompassing output dataframe dictionary
    #     note, it is still limited to the lowest common denominator of all
    #     group model choices- it does not pull in the entire output directory
    # - there will be a dataframe for each combination of output measure
    #   type and preprocessing strategy
    # - each dataframe will contain output filepaths and their associated
    #   information, and each dataframe will include ALL SERIES/SCANS
    # - the dataframes will be pruned for each model LATER
    output_df_dict = gather_outputs(
        pipeline_dir, output_measure_list, inclusion_list, get_motion, get_raw_score
    )

    # alright, group model processing time
    #   going to merge the phenotype DFs with the output file DF
    analysis_dict = {}

    model_name = group_model.model_name

    if len(group_model.derivative_list) == 0:
        err = (
            "\n\n[!] There are no derivatives listed in the "
            "derivative_list field of your group analysis "
            "configuration file.\n\nConfiguration file: "
            f"{group_config_file}\n"
        )
        raise Exception(err)

    # load original phenotype CSV into a dataframe
    pheno_df = read_pheno_csv_into_df(
        group_model.pheno_file, group_model.participant_id_label
    )

    # enforce the sub ID label to "Participant"
    pheno_df.rename(
        columns={group_model.participant_id_label: "participant_id"}, inplace=True
    )
    pheno_df["participant_id"] = pheno_df["participant_id"].astype(str)

    # unique_resource = (output_measure_type, preprocessing strategy)
    # output_df_dict[unique_resource] = dataframe
    for unique_resource in output_df_dict.keys():
        resource_id = unique_resource[0]

        # do this backwards, because the group_model.derivative_list is a
        # list of substrings that would be in a derivative name
        # for example:
        #     group_model.derivative_list = ['centrality']
        #     this would include both 'centrality_zstd' and
        #     'centrality_smooth_zstd', both of which could be the current
        #     value of 'resource_id'
        # also, 'derivatives' is a list of group-analysis eligible
        # derivatives (standard space, z-score standardized)
        for derivative in group_model.derivative_list:
            if derivative in resource_id and resource_id in derivatives:
                break
        else:
            continue

        strat_info = unique_resource[1]

        # output_df has the information for ALL of the output files for
        # this unique_resource_id- all series, and if applicable, motion
        # params numbers, and paths to raw outputs (for measure mean or
        # custom ROI means)
        #   then cut it down and merge with the phenotype DF as needed
        #   depending on the analysis
        output_df = output_df_dict[unique_resource]

        # prune the output_df for this specific group model and output +
        # preprocessing strategy
        if not group_model.participant_list:
            inclusion_list = grab_pipeline_dir_subs(pipeline_dir)
            output_df = output_df[output_df["participant_id"].isin(inclusion_list)]
        elif os.path.isfile(group_model.participant_list):
            inclusion_list = load_text_file(
                group_model.participant_list, "group-level analysis participant list"
            )
            output_df = output_df[output_df["participant_id"].isin(inclusion_list)]
        else:
            msg = "\nCannot read group-level analysis participant list.\n"
            raise Exception(msg)

        new_pheno_df = pheno_df.copy()

        # check for inconsistency with leading zeroes
        # (sometimes, the sub_ids from individual will be something like
        #  '0002601' and the phenotype will have '2601')
        sublist_subs = output_df["participant_id"]
        pheno_subs = list(new_pheno_df["participant_id"])

        for sub in sublist_subs:
            if sub in pheno_subs:
                # okay, there's at least one match
                break
        else:
            new_sublist_subs = [str(x).lstrip("0") for x in sublist_subs]
            for sub in new_sublist_subs:
                if sub in pheno_subs:
                    # that's better
                    output_df["participant_id"] = new_sublist_subs
                    break
            else:
                msg = (
                    "the participant IDs in your group "
                    "analysis participant list and the "
                    "participant IDs in your phenotype file "
                    "do not match"
                )
                raise Exception(msg)

        repeated_measures = False
        repeated_sessions = False
        repeated_series = False

        if len(group_model.sessions_list) > 0:
            repeated_sessions = True

        if len(group_model.series_list) > 0:
            repeated_series = True

        if repeated_sessions or repeated_series:
            repeated_measures = True

        if repeated_measures:
            if repeated_sessions:
                # IF USING FSL PRESETS: new_pheno_df will get passed
                #                       through unchanged
                new_pheno_df = pheno_sessions_to_repeated_measures(
                    new_pheno_df, group_model.sessions_list
                )

            # create new rows for all of the series, if applicable
            #   ex. if 10 subjects and two sessions, 10 rows -> 20 rows
            if repeated_series:
                # IF USING FSL PRESETS: new_pheno_df will get passed
                #                       through unchanged
                new_pheno_df = pheno_series_to_repeated_measures(
                    new_pheno_df, group_model.series_list, repeated_sessions
                )

            # drop the pheno rows - if there are participants missing in
            # the output files (ex. if ReHo did not complete for 2 of the
            # participants, etc.), then drop these rows from the phenotype
            #   we are dropping all instances of this participant, all
            #   sessions and all series, because in repeated measures/
            #   within-subject, if one goes, they all have to go
            new_pheno_df = new_pheno_df[
                pheno_df["participant_id"].isin(output_df["participant_id"])
            ]

            if len(new_pheno_df) == 0:
                err = (
                    "\n\n[!] There is a mis-match between the "
                    "participant IDs in the output directory/particip"
                    "ant list and the phenotype file.\n\n"
                )
                raise Exception(err)

            join_columns = ["participant_id"]

            if "Series" in new_pheno_df:
                # if Series is one of the categorically-encoded covariates
                # make sure we only are including the series the user has
                # selected to include in the repeated measures analysis

                # check in case the pheno has series IDs that doesn't
                # exist in the output directory, first
                new_pheno_df = new_pheno_df[
                    new_pheno_df["Series"].isin(output_df["Series"])
                ]

                # okay, now check against the user-specified series list
                new_pheno_df = new_pheno_df[
                    new_pheno_df["Series"].isin(group_model.series_list)
                ]
                join_columns.append("Series")

                # pull together the pheno DF and the output files DF!
                new_pheno_df = pd.merge(
                    new_pheno_df, output_df, how="inner", on=join_columns
                )

                if repeated_sessions:
                    # this can be removed/modified once sessions are no
                    # longer integrated in the full unique participant IDs
                    new_pheno_df, dropped_parts = balance_repeated_measures(
                        new_pheno_df, group_model.sessions_list, group_model.series_list
                    )

                    session = "repeated_measures_multiple_sessions"
                    series = "repeated_measures_multiple_series"
                else:
                    series = "repeated_measures_multiple_series"
                    if "session" in output_df:
                        for ses_df_tuple in new_pheno_df.groupby("Sessions"):
                            session = f"ses-{ses_df_tuple[0]}"
                            ses_df = ses_df_tuple[1]

                            # send it in
                            analysis_dict[
                                (
                                    model_name,
                                    group_config_file,
                                    resource_id,
                                    strat_info,
                                    session,
                                    series,
                                )
                            ] = new_pheno_df
                    else:
                        # default a session
                        session = "ses-1"
                        # send it in
                        analysis_dict[
                            (
                                model_name,
                                group_config_file,
                                resource_id,
                                strat_info,
                                session,
                                series,
                            )
                        ] = new_pheno_df

            else:
                # this runs if there are repeated sessions but not
                # repeated series
                #   split up the series here
                #   iterate over the Series/Scans
                for series_df_tuple in output_df.groupby("Series"):
                    session = "repeated_measures_multiple_sessions"
                    series = series_df_tuple[0]

                    # series_df is output_df but with only one of the
                    # Series
                    series_df = series_df_tuple[1]
                    series_df = series_df.astype("str")

                    # trim down the pheno DF to match the output DF, remove any
                    # extra participant_<ID> identity columns, and merge
                    newer_pheno_df = new_pheno_df[
                        pheno_df["participant_id"].isin(series_df["participant_id"])
                    ]
                    for col in newer_pheno_df.columns:
                        if (
                            "participant_" in col
                            and "id" not in col
                            and "session" not in col
                        ):
                            if col.replace("participant_", "") not in list(
                                series_df["participant_id"]
                            ):
                                newer_pheno_df = newer_pheno_df.drop(labels=col, axis=1)

                    newer_pheno_df = newer_pheno_df.astype("str")
                    newer_pheno_df = pd.merge(
                        newer_pheno_df,
                        series_df,
                        how="inner",
                        on=["participant_id", "Sessions"],
                    )

                    # this can be removed/modified once sessions are no
                    # longer integrated in the full unique participant IDs
                    if "Sessions" in newer_pheno_df.columns:
                        newer_pheno_df, dropped_parts = balance_repeated_measures(
                            newer_pheno_df, group_model.sessions_list, None
                        )

                    # unique_resource =
                    #        (output_measure_type, preprocessing strategy)
                    analysis_dict[
                        (
                            model_name,
                            group_config_file,
                            resource_id,
                            strat_info,
                            session,
                            series,
                        )
                    ] = newer_pheno_df

        else:
            # no repeated measures

            # split up the output files list DataFrame by series, then
            # merge with the pheno DataFrame and send it off for analysis

            # essentially, make sure each series combination goes into its
            # own model (and dataframe) for this unique_resource_id

            # iterate over the Series/Scans
            for series_df_tuple in output_df.groupby("Series"):
                series = series_df_tuple[0]

                # series_df - this is output_df but with only one of the Series
                series_df = series_df_tuple[1]

                # trim down the pheno DF to match the output DF and merge
                newer_pheno_df = new_pheno_df[
                    pheno_df["participant_id"].isin(series_df["participant_id"])
                ]

                # multiple sessions?
                if "Sessions" in series_df:
                    for ses_df_tuple in series_df.groupby("Sessions"):
                        session = f"ses-{ses_df_tuple[0]}"
                        ses_df = ses_df_tuple[1]
                        newer_ses_pheno_df = pd.merge(
                            newer_pheno_df, ses_df, how="inner", on=["participant_id"]
                        )
                        # send it in
                        analysis_dict[
                            (
                                model_name,
                                group_config_file,
                                resource_id,
                                strat_info,
                                session,
                                series,
                            )
                        ] = newer_ses_pheno_df
                else:
                    # default a session
                    session = "ses-1"
                    newer_pheno_df = pd.merge(
                        newer_pheno_df, series_df, how="inner", on=["participant_id"]
                    )
                    # send it in
                    analysis_dict[
                        (
                            model_name,
                            group_config_file,
                            resource_id,
                            strat_info,
                            session,
                            series,
                        )
                    ] = newer_pheno_df

    if len(analysis_dict) == 0:
        err = (
            "\n\n[!] C-PAC says: Could not find a match between the "
            "participants in your pipeline output directory that were "
            "included in your analysis, and the participants in the "
            "phenotype file provided.\n\n"
        )
        raise Exception(err)

    return analysis_dict


def build_feat_models(group_config_file):
    """Build FSL's FEAT models for group analysis."""
    import os

    from CPAC.pipeline.cpac_ga_model_generator import build_feat_model

    analysis_dict = prep_feat_inputs(group_config_file)

    for unique_resource_id in analysis_dict.keys():
        # unique_resource_id is a 6-long tuple:
        #    ( model name, group model config file, output measure name,
        #          preprocessing strategy string, session_id,
        #          series_id or "repeated_measures" )

        model_name = unique_resource_id[0]
        group_config_file = unique_resource_id[1]
        resource_id = unique_resource_id[2]
        preproc_strat = unique_resource_id[3]
        session = unique_resource_id[4]
        series = unique_resource_id[5]
        model_df = analysis_dict[unique_resource_id]

        dmat_csv_path, new_sub_file, empty_csv = build_feat_model(
            model_df,
            model_name,
            group_config_file,
            resource_id,
            preproc_strat,
            session,
            series,
        )

    if os.path.isfile(empty_csv):
        return 0
    return -1


def run_feat(group_config_file, feat=True):
    """Run FSL's FEAT tool for group analysis."""
    from multiprocessing import Process
    import os

    import numpy as np
    import pandas as pd

    from CPAC.pipeline.cpac_ga_model_generator import create_dir
    from CPAC.utils.create_flame_model_files import create_flame_model_files

    # let's get the show on the road
    procss = []

    # get group pipeline config loaded
    c = load_config_yml(group_config_file)

    pipeline_dir = c["pipeline_setup"]["output_directory"]["source_outputs_path"]
    model_name = c["fsl_feat"]["model_name"]
    out_dir = c["pipeline_setup"]["output_directory"]["output_path"]

    pipeline_name = pipeline_dir.rstrip("/").split("/")[-1]

    model_dir = os.path.join(
        out_dir,
        "cpac_group_analysis",
        "FSL_FEAT",
        f"{pipeline_name}",
        f"group_model_{model_name}",
    )

    custom_contrasts_csv = os.path.join(model_dir, "contrasts.csv")

    contrasts_df = pd.read_csv(custom_contrasts_csv)
    if contrasts_df.shape[0] == 1 and np.count_nonzero(contrasts_df.values[0][1:]) == 0:
        err = (
            "\n\n[!] C-PAC says: It appears you haven't defined any "
            f"contrasts in your contrasts CSV file.\n\nContrasts file:\n{custom_contrasts_csv}"
            "\n\nDefine your contrasts in this file and run again."
            "\n\n"
        )
        raise Exception(err)

    models = {}
    for root, dirs, files in os.walk(model_dir):
        for filename in files:
            filepath = os.path.join(root, filename)
            second_half = filepath.split(model_dir)[1].split("/")
            second_half.remove("")

            try:
                id_tuple = (
                    second_half[0],
                    second_half[1],
                    second_half[2],
                    second_half[3],
                )
            except IndexError:
                # not a file we are interested in
                continue

            if id_tuple not in models.keys():
                models[id_tuple] = {}

            if "group_sublist" in filepath:
                models[id_tuple]["group_sublist"] = filepath
            elif "design_matrix" in filepath:
                models[id_tuple]["design_matrix"] = filepath
                models[id_tuple]["dir_path"] = filepath.replace(
                    "model_files/design_matrix.csv", ""
                )
            elif "groups" in filepath:
                models[id_tuple]["group_vector"] = filepath
            elif "merged_mask" in filepath:
                models[id_tuple]["merged_mask"] = filepath
            elif "merged" in filepath:
                models[id_tuple]["merged"] = filepath

    if len(models) == 0:
        err = (
            "\n\n[!] C-PAC says: Cannot find the FSL-FEAT/Randomise model "
            f"files.\n\nI am looking here:\n{model_dir}\n\nIf that doesn't sound "
            "right, double-check your group configuration file.\n\nDid you "
            "build the model beforehand?\n\n"
        )
        raise Exception(err)

    for id_tuple in models.keys():
        # generate working/log directory for this sub-model's group analysis run
        work_dir = os.path.join(
            c.work_dir, models[id_tuple]["dir_path"].replace(out_dir, "").lstrip("/")
        )
        work_dir = work_dir.replace(
            "cpac_group_analysis", "cpac_group_analysis_workdir"
        )
        work_dir = work_dir.replace("model_files/", "")
        log_dir = os.path.join(
            c["pipeline_setup"]["log_directory"]["path"],
            models[id_tuple]["dir_path"].replace(out_dir, "").lstrip("/"),
        )
        log_dir = log_dir.replace("cpac_group_analysis", "cpac_group_analysis_logdir")
        log_dir = log_dir.replace("model_files/", "")

        model_out_dir = os.path.join(models[id_tuple]["dir_path"], "fsl-feat_results")
        model_out_dir = model_out_dir.replace("model_files/", "")
        input_files_dir = os.path.join(
            models[id_tuple]["dir_path"], "flame_input_files"
        )

        create_dir(work_dir, "group analysis working")
        create_dir(log_dir, "group analysis logfile")
        create_dir(input_files_dir, "FSL-FEAT FLAME tool input files")
        create_dir(model_out_dir, "FSL-FEAT output files")

        design_matrix = pd.read_csv(models[id_tuple]["design_matrix"])
        design_matrix = design_matrix.drop(labels="participant_id", axis=1)

        grp_vector = load_text_file(
            models[id_tuple]["group_vector"], "group vector file"
        )

        mat, grp, con, fts = create_flame_model_files(
            design_matrix,
            design_matrix.columns,
            None,
            None,
            custom_contrasts_csv,
            None,
            c["fsl_feat"]["group_sep"],
            grp_vector,
            c["fsl_feat"]["coding_scheme"],
            model_name,
            id_tuple[0],
            input_files_dir,
        )
        if fts:
            f_test = True
        else:
            f_test = False

        if not con:
            WFLOGGER.warning(
                "\n\n################## MODEL NOT BEING INCLUDED ##################"
                "\n\n[!] C-PAC says: There is a mismatch between the design matrix and"
                " contrasts matrix for this model:\n\nDerivative: %s\nSession: %s"
                "\nScan: %s\nPreprocessing strategy:\n    %s\n\nThe model is not"
                " proceeding into the FSL-FEAT FLAME run."
                "\n\n#########################################################\n",
                id_tuple[0],
                id_tuple[1],
                id_tuple[2],
                id_tuple[3],
            )
            continue

        if feat:
            from CPAC.group_analysis.group_analysis import run_feat_pipeline

            procss.append(
                Process(
                    target=run_feat_pipeline,
                    args=(
                        c,
                        models[id_tuple]["merged"],
                        models[id_tuple]["merged_mask"],
                        f_test,
                        mat,
                        con,
                        grp,
                        model_out_dir,
                        work_dir,
                        log_dir,
                        model_name,
                        fts,
                    ),
                )
            )
        else:
            from CPAC.randomise.randomise import prep_randomise_workflow

            model_out_dir = model_out_dir.replace("feat_results", "randomise_results")
            procss.append(
                Process(
                    target=prep_randomise_workflow,
                    args=(
                        c,
                        models[id_tuple]["merged"],
                        models[id_tuple]["merged_mask"],
                        f_test,
                        mat,
                        con,
                        grp,
                        model_out_dir,
                        work_dir,
                        log_dir,
                        model_name,
                        fts,
                    ),
                )
            )

    manage_processes(procss, out_dir, c["fsl_feat"]["num_models_at_once"])


def run_cwas_group(
    pipeline_dir,
    out_dir,
    working_dir,
    crash_dir,
    roi_file,
    regressor_file,
    participant_column,
    columns,
    permutations,
    parallel_nodes,
    plugin_args,
    z_score,
    inclusion=None,
):
    """Run a group-level CWAS analysis."""
    import os

    from CPAC.cwas.pipeline import create_cwas

    pipeline_dir = os.path.abspath(pipeline_dir)

    out_dir = os.path.join(
        out_dir, "cpac_group_analysis", "MDMR", os.path.basename(pipeline_dir)
    )

    working_dir = os.path.join(
        working_dir, "cpac_group_analysis", "MDMR", os.path.basename(pipeline_dir)
    )

    crash_dir = os.path.join(
        crash_dir, "cpac_group_analysis", "MDMR", os.path.basename(pipeline_dir)
    )

    inclusion_list = None

    if inclusion:
        inclusion_list = load_text_file(inclusion, "MDMR participant inclusion list")

    output_df_dct = gather_outputs(
        pipeline_dir,
        ["space-template_desc-preproc_bold"],
        inclusion_list,
        False,
        False,
        get_func=True,
    )

    for preproc_strat in output_df_dct.keys():
        # go over each preprocessing strategy

        df_dct = {}
        strat_df = output_df_dct[preproc_strat]

        if len(set(strat_df["Series"])) > 1:
            # more than one scan/series ID
            for strat_scan in list(set(strat_df["Series"])):
                # make a list of sub-dataframes, each one with only file paths
                # from one scan ID each
                df_dct[strat_scan] = strat_df[strat_df["Series"] == strat_scan]
        else:
            df_dct[next(iter(set(strat_df["Series"])))] = strat_df
        for df_scan in df_dct.keys():
            func_paths = {
                p.split("_")[0]: f
                for p, f in zip(
                    df_dct[df_scan].participant_id, df_dct[df_scan].Filepath
                )
            }

            if plugin_args["n_procs"] == 1:
                plugin = "Linear"
            else:
                plugin = "MultiProc"

            cwas_wf = create_cwas(
                name=f"MDMR_{df_scan}",
                working_dir=working_dir,
                crash_dir=crash_dir,
            )
            cwas_wf.inputs.inputspec.subjects = func_paths
            cwas_wf.inputs.inputspec.roi = roi_file
            cwas_wf.inputs.inputspec.regressor = regressor_file
            cwas_wf.inputs.inputspec.participant_column = participant_column
            cwas_wf.inputs.inputspec.columns = columns
            cwas_wf.inputs.inputspec.permutations = permutations
            cwas_wf.inputs.inputspec.parallel_nodes = parallel_nodes
            cwas_wf.inputs.inputspec.z_score = z_score
            cwas_wf.run(plugin=plugin, plugin_args=plugin_args)


def run_cwas(pipeline_config):
    """Run CWAS."""
    import os

    import yaml

    pipeline_config = os.path.abspath(pipeline_config)

    pipeconfig_dct = yaml.safe_load(open(pipeline_config, "r"))

    num_cpus = pipeconfig_dct["pipeline_setup"]["system_config"]["num_cpus"]
    mem_gb = pipeconfig_dct["pipeline_setup"]["system_config"]["num_memory"]

    plugin_args = {"n_procs": num_cpus, "memory_gb": mem_gb}

    pipeline = pipeconfig_dct["pipeline_setup"]["output_directory"][
        "source_outputs_path"
    ]
    output_dir = pipeconfig_dct["pipeline_setup"]["output_directory"]["output_path"]
    working_dir = pipeconfig_dct["pipeline_setup"]["working_directory"]["path"]
    crash_dir = pipeconfig_dct["pipeline_setup"]["log_directory"]["path"]

    roi_file = pipeconfig_dct["mdmr"]["roi_file"]
    regressor_file = pipeconfig_dct["mdmr"]["regressor_file"]
    participant_column = pipeconfig_dct["mdmr"]["regressor_participant_column"]
    columns = pipeconfig_dct["mdmr"]["regressor_columns"]
    permutations = pipeconfig_dct["mdmr"]["permutations"]
    parallel_nodes = pipeconfig_dct["mdmr"]["parallel_nodes"]
    inclusion = pipeconfig_dct["mdmr"]["inclusion_list"]
    z_score = pipeconfig_dct["mdmr"]["zscore"]

    if not inclusion or "None" in inclusion or "none" in inclusion:
        inclusion = None

    run_cwas_group(
        pipeline,
        output_dir,
        working_dir,
        crash_dir,
        roi_file,
        regressor_file,
        participant_column,
        columns,
        permutations,
        parallel_nodes,
        plugin_args,
        z_score,
        inclusion=inclusion,
    )


def find_other_res_template(template_path: str, new_resolution: int) -> str:
    """Find the same template/standard file in another resolution, if it exists.

    Parameters
    ----------
    template_path : str
        file path to the template NIfTI file.

    new_resolution : int
        the resolution of the template file you need

    Notes
    -----
    Makes an assumption regarding the filename format of the files.
    """
    # TODO: this is assuming there is a mm resolution in the file path - not
    # TODO: robust to varying templates - look into alternatives

    ref_file = None

    if "mm" in template_path:
        parts = {}
        try:
            parts["left"], parts["right"] = template_path.rsplit("mm", 1)
        except ValueError:
            # TODO: better message
            msg = "no resolution in the file path!"
            raise Exception(msg)

        parts["left"] = str(new_resolution).join(
            parts["left"].rsplit(parts["left"][-1], 1)
        )
        ref_file = f"{parts['left']}{parts['right']}"

    elif "${resolution_for_func_preproc}" in template_path:
        ref_file = template_path.replace(
            "${resolution_for_func_preproc}", f"{new_resolution}mm"
        )

    if ref_file:
        WFLOGGER.info(
            "\nAttempting to find %smm version of the template:\n%s\n\n",
            new_resolution,
            ref_file,
        )

    return ref_file


def check_cpac_output_image(image_path, reference_path, out_dir=None, roi_file=False):
    """Check if the image needs to be resampled. If so, return the command."""
    import os

    import nibabel as nib

    if not out_dir:
        out_dir = os.getcwd()

    # we want to preserve the original directory structure of the input image,
    # but place that sub-tree into the BASC working directory (in this case,
    # 'out_dir')
    try:
        orig_dir = "pipeline_{0}".format(image_path.split("pipeline_")[1])
    except IndexError:
        if roi_file:
            orig_dir = os.path.join("ROI_files", os.path.basename(image_path))
        else:
            raise IndexError(image_path)
    out_path = os.path.join(out_dir, "resampled_input_images", orig_dir)

    # if this was already done
    if os.path.isfile(out_path):
        image_path = out_path

    resample = False

    image_nb = nib.load(image_path)
    ref_nb = nib.load(reference_path)

    # check: do we even need to resample?
    if int(image_nb.header.get_zooms()[0]) != int(ref_nb.header.get_zooms()[0]):
        WFLOGGER.info(
            "Input image resolution is %smm\nTemplate image resolution is %smm\n",
            image_nb.header.get_zooms()[0],
            ref_nb.header.get_zooms()[0],
        )
        resample = True
    if image_nb.shape != ref_nb.shape:
        WFLOGGER.info(
            "Input image shape is %s\nTemplate image shape is %s\n",
            image_nb.shape,
            ref_nb.shape,
        )
        resample = True

    if resample:
        if not os.path.isdir(out_path.replace(os.path.basename(out_path), "")):
            try:
                os.makedirs(out_path.replace(os.path.basename(out_path), ""))
            except (OSError, TypeError) as os_error:
                # TODO: better message
                msg = "couldn't make the dirs!"
                raise OSError(msg) from os_error

        WFLOGGER.info(
            "Resampling input image:\n%s\n\n..to this reference:\n%s\n\n..and writing"
            " this file here:\n%s\n",
            image_path,
            reference_path,
            out_path,
        )
        cmd = ["flirt", "-in", image_path, "-ref", reference_path, "-out", out_path]
        if roi_file:
            cmd.append("-interp")
            cmd.append("nearestneighbour")
        return cmd
    return resample


def resample_cpac_output_image(cmd_args):
    """Run resampling command and return the output file path."""
    import subprocess

    WFLOGGER.info("Running:\n%s\n\n", " ".join(cmd_args))

    flag = "resampled_input_images"

    for arg in cmd_args:
        if flag in arg:
            out_file = arg
            subprocess.check_output(cmd_args)
            return out_file
    msg = f"Missing required argument '{flag}'"
    raise ValueError(msg)


def launch_PyBASC(pybasc_config):
    """Run PyBASC."""
    import subprocess

    WFLOGGER.info("Running PyBASC with configuration file:\n%s", pybasc_config)

    cmd_args = ["PyBASC", pybasc_config]
    return subprocess.check_output(cmd_args)


def run_basc(pipeline_config):
    """
    Run the PyBASC module.

    PyBASC is a separate Python package built and maintained by Aki Nikolaidis
    which implements the BASC analysis via Python.
    PyBASC is based off of the following work:
        - Garcia-Garcia, M., Nikolaidis, A., Bellec, P., Craddock, R. C., Cheung, B., Castellanos, F. X., & Milham, M. P. (2017).
              Detecting stable individual differences in the functional organization of the human basal ganglia. NeuroImage.
        - Bellec, P., Rosa-Neto, P., Lyttelton, O. C., Benali, H., & Evans, A. C. (2010).
              Multi-level bootstrap analysis of stable clusters in resting-state fMRI. Neuroimage, 51(3), 1126-1139.
        - Bellec, P., Marrelec, G., & Benali, H. (2008).
              A bootstrap test to investigate changes in brain connectivity for functional MRI. Statistica Sinica, 1253-1268.
    PyBASC GitHub repository:
        https://github.com/AkiNikolaidis/PyBASC
    PyBASC author:
        https://www.researchgate.net/profile/Aki_Nikolaidis
    Inputs
        pipeline_config: path to C-PAC pipeline configuration YAML file
    Steps (of the C-PAC interface for PyBASC, not PyBASC itself)
        1. Read in the PyBASC-relevant pipeline config items and create a new
           PyBASC config dictionary.
        2. Ensure the functional template is in the selected resolution.
        3. If ROI mask files are on AWS S3, download these into the group
           level analysis working directory.
        4. Resample the ROI mask files if they are not in the selected
           output resolution.
        5. Create sub-directories for each C-PAC pipeline the user has
           selected to run PyBASC for (preprocessed and template-space
           functional time series are pulled from each pipeline output
           directory, for input into PyBASC).
        6. Gather space-template_bold outputs from each pipeline.
        7. Create further sub-directories for each nuisance regression
           strategy and functional scan within each C-PAC pipeline, and
           separate the functional outputs by strategy and scan as well.
        8. Resample functional time series if they are not in the selected
           output resolution.
        9. Finish populating the PyBASC config dictionary, and write it out
           into a config YAML file for each pipeline-strategy-scan we are
           running.
        10. Launch PyBASC for each configuration generated.
    """
    from multiprocessing import pool
    import os

    import yaml

    from CPAC.utils.datasource import check_for_s3

    pipeline_config = os.path.abspath(pipeline_config)

    pipeconfig_dct = yaml.safe_load(open(pipeline_config, "r"))

    output_dir = os.path.abspath(
        pipeconfig_dct["pipeline_setup"]["output_directory"]["output_path"]
    )
    working_dir = os.path.abspath(
        pipeconfig_dct["pipeline_setup"]["working_directory"]["path"]
    )
    if pipeconfig_dct["pipeline_setup"]["Amazon-AWS"]["aws_output_bucket_credentials"]:
        creds_path = os.path.abspath(
            pipeconfig_dct["pipeline_setup"]["Amazon-AWS"][
                "aws_output_bucket_credentials"
            ]
        )

    func_template = pipeconfig_dct["basc"]["template_brain_only_for_func"]
    if "$FSLDIR" in func_template:
        if os.environ.get("FSLDIR"):
            func_template = func_template.replace("$FSLDIR", os.environ["FSLDIR"])

    basc_inclusion = pipeconfig_dct["pipeline_setup"]["output_directory"][
        "participant_list"
    ]
    basc_scan_inclusion = pipeconfig_dct["basc"]["scan_inclusion"]
    basc_resolution = pipeconfig_dct["basc"]["resolution"]

    basc_config_dct = {"run": True, "reruns": 1}

    for key in pipeconfig_dct.keys():
        if "basc" in key:
            basc_config_dct = pipeconfig_dct[key]

    iterables = [
        "dataset_bootstrap_list",
        "timeseries_bootstrap_list",
        "blocklength_list",
        "n_clusters_list",
        "output_sizes",
    ]
    for iterable in iterables:
        basc_config_dct[iterable] = [
            int(x) for x in str(basc_config_dct[iterable]).split(",")
        ]

    basc_config_dct["proc_mem"] = [basc_config_dct["proc"], basc_config_dct["memory"]]
    del basc_config_dct["proc"]
    del basc_config_dct["memory"]

    if "None" in basc_inclusion or "none" in basc_inclusion:
        basc_inclusion = None

    # we have the functional template only for potential resampling - to make
    # sure everything is the same resolution and shape (as what the user has
    # selected)
    if "mm" in basc_resolution:
        basc_resolution = basc_resolution.replace("mm", "")

    # get the functional template, but in the specified resolution for BASC
    ref_file = find_other_res_template(func_template, basc_resolution)

    # did that actually work?
    if not os.path.isfile(ref_file):
        # TODO: better message
        msg = f"\n[!] The reference file could not be found.\nPath: {ref_file}\n"
        raise FileNotFoundError(msg)

    working_dir = os.path.join(
        working_dir,
        "cpac_group_analysis",
        "PyBASC",
        f"{basc_resolution}mm_resolution",
        "working_dir",
    )

    # check for S3 of ROI files here!
    for key in ["roi_mask_file", "cross_cluster_mask_file"]:
        if "s3://" in basc_config_dct[key]:
            basc_config_dct[key] = check_for_s3(
                basc_config_dct[key], creds_path=creds_path, dl_dir=working_dir
            )

    # resample ROI files if necessary
    roi_cmd_args = check_cpac_output_image(
        basc_config_dct["roi_mask_file"], ref_file, out_dir=working_dir, roi_file=True
    )
    roi_two_cmd_args = check_cpac_output_image(
        basc_config_dct["cross_cluster_mask_file"],
        ref_file,
        out_dir=working_dir,
        roi_file=True,
    )

    if roi_cmd_args:
        roi_file = resample_cpac_output_image(roi_cmd_args)
        basc_config_dct["roi_mask_file"] = roi_file
    if roi_two_cmd_args:
        roi_file_two = resample_cpac_output_image(roi_two_cmd_args)
        basc_config_dct["cross_cluster_mask_file"] = roi_file_two

    pipeline_dir = os.path.abspath(
        pipeconfig_dct["pipeline_setup"]["output_directory"]["source_outputs_path"]
    )

    out_dir = os.path.join(
        output_dir,
        "cpac_group_analysis",
        "PyBASC",
        f"{basc_resolution}mm_resolution",
        os.path.basename(pipeline_dir),
    )
    working_dir = os.path.join(working_dir, os.path.basename(pipeline_dir))

    inclusion_list = None
    scan_inclusion = None

    if basc_inclusion:
        inclusion_list = load_text_file(
            basc_inclusion, "BASC participant inclusion list"
        )

    if "none" in basc_scan_inclusion.lower():
        basc_scan_inclusion = None
    if basc_scan_inclusion:
        scan_inclusion = basc_scan_inclusion.split(",")

    # create encompassing output dataframe dictionary
    #     note, it is still limited to the lowest common denominator of all
    #     group model choices- it does not pull in the entire output directory
    # - there will be a dataframe for each combination of output measure
    #   type and preprocessing strategy
    # - each dataframe will contain output filepaths and their associated
    #   information, and each dataframe will include ALL SERIES/SCANS
    output_df_dct = gather_outputs(
        pipeline_dir,
        ["space-template_bold"],
        inclusion_list,
        False,
        False,
        get_func=True,
    )

    for preproc_strat in output_df_dct.keys():
        # go over each preprocessing strategy

        df_dct = {}
        strat_df = output_df_dct[preproc_strat]

        nuisance_string = preproc_strat[1].replace(
            os.path.basename(preproc_strat[1]), ""
        )

        if len(set(strat_df["Series"])) > 1:
            # more than one scan/series ID
            for strat_scan in list(set(strat_df["Series"])):
                # make a list of sub-dataframes, each one with only file paths
                # from one scan ID each
                df_dct[strat_scan] = strat_df[strat_df["Series"] == strat_scan]
        else:
            df_dct[next(iter(set(strat_df["Series"])))] = strat_df

        # TODO: need a catch for if none of the df_scans below are in scan_inclusion

        for df_scan in df_dct.keys():
            # do only the selected scans
            if scan_inclusion:
                if df_scan not in scan_inclusion:
                    continue

            basc_config_dct["analysis_ID"] = (
                f"{os.path.basename(pipeline_dir)}_{df_scan}"
            )

            # add scan label and nuisance regression strategy label to the
            # output directory path
            scan_out_dir = os.path.join(out_dir, df_scan, nuisance_string.lstrip("/"))
            scan_working_dir = os.path.join(
                working_dir, df_scan, nuisance_string.lstrip("/")
            )

            basc_config_dct["home"] = scan_out_dir
            basc_config_dct["cluster_methods"] = ["ward"]

            func_paths = list(df_dct[df_scan]["Filepath"])

            # affinity threshold is an iterable, and must match the number of
            # functional file paths for the MapNodes
            pipeconfig_dct["basc"]["affinity_thresh"] * len(func_paths)

            # resampling if necessary
            #     each run should take the file, resample it and write it
            #     into the BASC sub-dir of the working directory
            #         should end up with a new "func_paths" list with all of
            #         these file paths in it
            ref_file_iterable = [ref_file] * len(func_paths)
            working_dir_iterable = [scan_working_dir] * len(func_paths)
            func_cmd_args_list = map(
                check_cpac_output_image,
                func_paths,
                ref_file_iterable,
                working_dir_iterable,
            )

            # resample them now
            if func_cmd_args_list[0]:
                p = pool.Pool(int(basc_config_dct["proc_mem"][0]))
                func_paths = p.map(resample_cpac_output_image, func_cmd_args_list)

            # TODO: add list into basc_config here
            basc_config_dct["subject_file_list"] = func_paths

            basc_config_outfile = os.path.join(scan_working_dir, "PyBASC_config.yml")
            WFLOGGER.info(
                "\nWriting PyBASC configuration file for %s scan in\n%s",
                df_scan,
                basc_config_outfile,
            )
            with open(basc_config_outfile, "wt") as f:
                noalias_dumper = yaml.dumper.SafeDumper
                noalias_dumper.ignore_aliases = lambda self, data: True
                f.write(
                    yaml.dump(
                        basc_config_dct, default_flow_style=False, Dumper=noalias_dumper
                    )
                )

            # go!
            launch_PyBASC(basc_config_outfile)


def run_isc_group(
    pipeline_dir,
    out_dir,
    working_dir,
    crash_dir,
    isc,
    isfc,
    levels=[],
    permutations=1000,
    std_filter=None,
    scan_inclusion=None,
    roi_inclusion=None,
    num_cpus=1,
):
    """Run the ISC pipeline for group-level analysis."""
    import os

    from CPAC.isc.pipeline import create_isc, create_isfc

    pipeline_dir = os.path.abspath(pipeline_dir)

    out_dir = os.path.join(
        out_dir, "cpac_group_analysis", "ISC", os.path.basename(pipeline_dir)
    )

    working_dir = os.path.join(
        working_dir, "cpac_group_analysis", "ISC", os.path.basename(pipeline_dir)
    )

    crash_dir = os.path.join(
        crash_dir, "cpac_group_analysis", "ISC", os.path.basename(pipeline_dir)
    )

    output_df_dct = gather_outputs(
        pipeline_dir,
        ["space-template_bold", "space-template_desc-Mean_timeseries"],
        inclusion_list=None,
        get_motion=False,
        get_raw_score=False,
        get_func=True,
        derivatives=["space-template_bold", "space-template_desc-Mean_timeseries"],
        # exts=['nii', 'nii.gz', 'csv']
    )

    for preproc_strat in output_df_dct.keys():
        # go over each preprocessing strategy

        derivative, _ = preproc_strat

        if "voxel" not in levels and derivative == "space-template_bold":
            continue

        if "roi" not in levels and derivative == "space-template_desc-Mean_timeseries":
            continue

        if derivative == "space-template_desc-Mean_timeseries":
            if roi_inclusion:
                # backwards because ROI labels embedded as substrings
                for roi_label in roi_inclusion:
                    if roi_label in _:
                        break
                else:
                    WFLOGGER.warning(
                        "ROI label '%s' not found in\n%s/%s\n", roi_label, derivative, _
                    )
                    continue

        df_dct = {}
        strat_df = output_df_dct[preproc_strat]

        if len(set(strat_df["Series"])) > 1:
            # more than one scan/series ID
            for strat_scan in list(set(strat_df["Series"])):
                # make a list of sub-dataframes, each one with only file paths
                # from one scan ID each
                df_dct[strat_scan] = strat_df[strat_df["Series"] == strat_scan]
        else:
            df_dct[next(iter(set(strat_df["Series"])))] = strat_df

        if isc:
            for df_scan in df_dct.keys():
                if scan_inclusion:
                    if df_scan not in scan_inclusion:
                        continue
                func_paths = {
                    p.split("_")[0]: f
                    for p, f in zip(
                        df_dct[df_scan].participant_id, df_dct[df_scan].Filepath
                    )
                }

                unique_out_dir = os.path.join(out_dir, "ISC", derivative, _, df_scan)

                it_id = "ISC_{0}_{1}_{2}".format(
                    df_scan, derivative, _.replace(".", "").replace("+", "")
                )

                isc_wf = create_isc(
                    name=it_id,
                    output_dir=unique_out_dir,
                    working_dir=working_dir,
                    crash_dir=crash_dir,
                )
                isc_wf.inputs.inputspec.subjects = func_paths
                isc_wf.inputs.inputspec.permutations = permutations
                isc_wf.inputs.inputspec.std = std_filter
                isc_wf.inputs.inputspec.collapse_subj = False
                isc_wf.run(plugin="MultiProc", plugin_args={"n_procs": num_cpus})

        if isfc:
            for df_scan in df_dct.keys():
                if scan_inclusion:
                    if df_scan not in scan_inclusion:
                        continue
                func_paths = {
                    p.split("_")[0]: f
                    for p, f in zip(
                        df_dct[df_scan].participant_id, df_dct[df_scan].Filepath
                    )
                }

                unique_out_dir = os.path.join(out_dir, "ISFC", derivative, _, df_scan)

                it_id = "ISFC_{0}_{1}_{2}".format(
                    df_scan, derivative, _.replace(".", "").replace("+", "")
                )

                isfc_wf = create_isfc(
                    name=it_id,
                    output_dir=unique_out_dir,
                    working_dir=working_dir,
                    crash_dir=crash_dir,
                )
                isfc_wf.inputs.inputspec.subjects = func_paths
                isfc_wf.inputs.inputspec.permutations = permutations
                isfc_wf.inputs.inputspec.std = std_filter
                isfc_wf.inputs.inputspec.collapse_subj = False
                isfc_wf.run(plugin="MultiProc", plugin_args={"n_procs": num_cpus})


def run_isc(pipeline_config):
    """Run the ISC pipeline."""
    import os

    import yaml

    pipeline_config = os.path.abspath(pipeline_config)

    pipeconfig_dct = yaml.safe_load(open(pipeline_config, "r"))

    pipeline_dir = pipeconfig_dct["pipeline_setup"]["output_directory"][
        "source_outputs_path"
    ]

    output_dir = pipeconfig_dct["pipeline_setup"]["output_directory"]["output_path"]
    working_dir = pipeconfig_dct["pipeline_setup"]["working_directory"]["path"]
    crash_dir = pipeconfig_dct["pipeline_setup"]["log_directory"]["path"]

    scan_inclusion = None
    if "scan_inclusion" in pipeconfig_dct.keys():
        scan_inclusion = pipeconfig_dct["pipeline_setup"]["system_config"][
            "scan_inclusion"
        ]

    roi_inclusion = None
    if "isc_roi_inclusion" in pipeconfig_dct.keys():
        roi_inclusion = pipeconfig_dct["isc_isfc"]["roi_inclusion"]

    num_cpus = 1
    if "num_cpus" in pipeconfig_dct.keys():
        num_cpus = pipeconfig_dct["pipeline_setup"]["system_config"]["num_cpus"]

    isc = 1 in pipeconfig_dct.get("runISC", [])
    isfc = 1 in pipeconfig_dct.get("runISFC", [])
    permutations = pipeconfig_dct.get("isc_permutations", 1000)
    std_filter = pipeconfig_dct.get("isc_level_voxel_std_filter", None)

    if std_filter == 0.0:
        std_filter = None

    levels = []
    if 1 in pipeconfig_dct.get("isc_level_voxel", []):
        levels += ["voxel"]

    if 1 in pipeconfig_dct.get("isc_level_roi", []):
        levels += ["roi"]

    if len(levels) == 0:
        return

    if not isc and not isfc:
        WFLOGGER.info(
            "\nISC and ISFC are not enabled to run in the group-level analysis"
            " configuration YAML file, and will not run.\n"
        )
        return

    pipeline_dirs = []
    if "pipeline_" in pipeline_dir:
        pipeline_dirs.append(pipeline_dir)
    for dirname in os.listdir(pipeline_dir):
        if "pipeline_" in dirname:
            pipeline_dirs.append(os.path.join(pipeline_dir, dirname))

    if not pipeline_dirs:
        WFLOGGER.error(
            "\nNo pipeline output directories found- make sure your 'pipeline_dir'"
            " field in the group configuration YAML file is pointing to a C-PAC"
            " pipeline output directory populated with a folder or folders that begin"
            " with the 'pipeline_' prefix.\n\nPipeline directory provided:\n%s\n",
            pipeline_dir,
        )

    for pipeline in pipeline_dirs:
        run_isc_group(
            pipeline,
            output_dir,
            working_dir,
            crash_dir,
            isc=isc,
            isfc=isfc,
            levels=levels,
            permutations=permutations,
            std_filter=std_filter,
            scan_inclusion=scan_inclusion,
            roi_inclusion=roi_inclusion,
            num_cpus=num_cpus,
        )


def run_qpp(group_config_file):
    """Run the QPP pipeline."""
    from CPAC.qpp.pipeline import create_qpp

    c = load_config_yml(group_config_file)

    pipeline_dir = os.path.abspath(
        c["pipeline_setup"]["output_directory"]["source_outputs_path"]
    )
    out_dir = os.path.join(
        c["pipeline_setup"]["output_directory"]["output_path"],
        "cpac_group_analysis",
        "QPP",
        os.path.basename(pipeline_dir),
    )
    working_dir = os.path.join(
        c["pipeline_setup"]["working_directory"]["path"],
        "cpac_group_analysis",
        "QPP",
        os.path.basename(pipeline_dir),
    )
    crash_dir = os.path.join(
        c["pipeline_setup"]["crash_log_directory"]["path"],
        "cpac_group_analysis",
        "QPP",
        os.path.basename(pipeline_dir),
    )

    try:
        os.makedirs(out_dir)
        os.makedirs(working_dir)
        os.makedirs(crash_dir)
    except OSError:
        pass

    outputs = gather_outputs(
        pipeline_dir,
        ["space-template_bold"],
        inclusion_list=c["pipeline_setup"]["output_directory"]["participant_list"],
        get_motion=False,
        get_raw_score=False,
        get_func=True,
        derivatives=["space-template_bold"],
        # exts=['nii', 'nii.gz']
    )

    if c["qpp"]["stratification"] == "Scan":
        qpp_stratification = ["Series"]
    elif c["qpp"]["stratification"] == "Session":
        qpp_stratification = ["Sessions"]
    elif c["qpp"]["stratification"] in ["Session and Scan", "Scan and Session"]:
        qpp_stratification = ["Sessions", "Series"]
    else:
        qpp_stratification = []

    for _output_df in outputs.values():
        output_df = _output_df
        if c["qpp"]["session_inclusion"]:
            output_df = output_df[
                output_df["Sessions"].isin(c["qpp"]["session_inclusion"])
            ]
        if c["qpp"]["scan_inclusion"]:
            output_df = output_df[output_df["Series"].isin(c["qpp"]["scan_inclusion"])]

        if qpp_stratification:
            output_df_groups = output_df.groupby(by=qpp_stratification)
        else:
            output_df_groups = [([], output_df)]

        for _group_id, _output_df_group in output_df_groups:
            group_id, output_df_group = _group_id, _output_df_group
            group = list(zip(qpp_stratification, group_id))

            group_id = "_".join(
                [
                    "%s-%s"
                    % (
                        {
                            "Sessions": "ses",
                            "Series": "scan",
                        }[k],
                        v,
                    )
                    for k, v in group
                ]
            )

            group_working_dir = os.path.join(working_dir, group_id)
            group_crash_dir = os.path.join(crash_dir, group_id)

            output_df_group, _ = balance_repeated_measures(
                output_df_group,
                output_df_group.Sessions.unique(),
                output_df_group.Series.unique(),
            )

            output_df_group = output_df_group.sort_values(by="participant_session_id")

            wf = create_qpp(
                name="QPP", working_dir=group_working_dir, crash_dir=group_crash_dir
            )

            wf.inputs.inputspec.window_length = c["qpp"]["window"]
            wf.inputs.inputspec.permutations = c["qpp"]["permutations"]
            wf.inputs.inputspec.lower_correlation_threshold = c["qpp"][
                "initial_threshold"
            ]
            wf.inputs.inputspec.higher_correlation_threshold = c["qpp"][
                "final_threshold"
            ]
            wf.inputs.inputspec.iterations = c["qpp"]["iterations"]
            wf.inputs.inputspec.correlation_threshold_iteration = c["qpp"][
                "initial_threshold_iterations"
            ]
            wf.inputs.inputspec.convergence_iterations = 1

            wf.inputs.inputspec.datasets = output_df_group.Filepath.tolist()

            wf.run()


def manage_processes(procss, output_dir, num_parallel=1):
    """Manage multiple processes in parallel."""
    import os

    # start kicking it off
    pid = open(os.path.join(output_dir, "pid_group.txt"), "w")

    jobQueue = []
    if len(procss) <= num_parallel:
        """
        Stream all the subjects as sublist is
        less than or equal to the number of
        subjects that need to run
        """
        for p in procss:
            p.start()
            print(p.pid, file=pid)

    else:
        """
        Stream the subject workflows for preprocessing.
        At Any time in the pipeline c.numSubjectsAtOnce
        will run, unless the number remaining is less than
        the value of the parameter stated above
        """
        idx = 0
        while idx < len(procss):
            if len(jobQueue) == 0 and idx == 0:
                idc = idx
                for p in procss[idc : idc + num_parallel]:
                    p.start()
                    print(p.pid, file=pid)
                    jobQueue.append(p)
                    idx += 1
            else:
                for job in jobQueue:
                    if not job.is_alive():
                        WFLOGGER.warning("found dead job %s", job)
                        loc = jobQueue.index(job)
                        del jobQueue[loc]
                        procss[idx].start()
                        jobQueue.append(procss[idx])
                        idx += 1

    pid.close()


def run(config_file):
    """Run all group analyses.

    This function only really exists for
    the "Run Group-Level Analysis" command on the GUI
    """
    # get MAIN pipeline config loaded
    c = load_config_yml(config_file)

    # Run MDMR, if selected
    if 1 in c["mdmr"]["run"]:
        run_cwas(config_file)

    # Run ISC, if selected
    if 1 in c["isc_isfc"]["runISC"] or 1 in c["isc_isfc"]["runISFC"]:
        run_isc(config_file)

    # Run PyBASC, if selected
    if 1 in c["basc"]["run"]:
        run_basc(config_file)

    # Run FSL FEAT group analysis, if selected
    if 1 in c["fsl_feat"]["run"]:
        run_feat(config_file)

    # Run randomise, if selected
    if 1 in c["fsl_randomise"]["run"]:
        run_feat(config_file, feat=False)

    # Run QPP, if selected
    if 1 in c["qpp"]["run"]:
        run_qpp(config_file)
