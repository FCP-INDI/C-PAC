# Copyright (C) 2012-2024  C-PAC Developers

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
import glob
import logging
import os
from pathlib import Path
import string
import sys
from typing import BinaryIO, Optional

import yaml

logger = logging.getLogger("extract_data_logs")
if logger.handlers:
    for handler in logger.handlers:
        logger.removeHandler(handler)
logging.basicConfig(
    filename=os.path.join(os.getcwd(), "extract_data_logs.log"),
    filemode="w",
    level=logging.DEBUG,
    format="%(levelname)s %(asctime)s %(lineno)d %(message)s",
)


def extract_data(c, param_map):
    """
    Generate a CPAC input subject list Python file.

    The method extracts anatomical and functional data for each site (if multiple site)
    and/or scan parameters for each site and put it into a data structure read by Python.

    Examples
    --------
    subjects_list =[
       {
        'subject_id' : '0050386',
        'unique_id' : 'session_1',
        'anat': '/Users/home/data/NYU/0050386/session_1/anat_1/anat.nii.gz',
        'rest':{
            'rest_1_rest' : '/Users/home/data/NYU/0050386/session_1/rest_1/rest.nii.gz',
            'rest_2_rest' : '/Users/home/data/NYU/0050386/session_1/rest_2/rest.nii.gz',
            }
        'scan_parameters':{
            'tr': '2',
            'acquisition': 'alt+z2',
            'reference': '17',
            'first_tr': '',
            'last_tr': '',
            }
        },
    ]

    or

    subjects_list =[
       {
        'subject_id' : '0050386',
        'unique_id' : 'session_1',
        'anat': '/Users/home/data/NYU/0050386/session_1/anat_1/anat.nii.gz',
        'rest':{
            'rest_1_rest' : '/Users/home/data/NYU/0050386/session_1/rest_1/rest.nii.gz',
            'rest_2_rest' : '/Users/home/data/NYU/0050386/session_1/rest_2/rest.nii.gz',
            }
          },
    ]

    """

    def get_list(arg) -> list:
        """Read each line of the file into list."""
        if isinstance(arg, list):
            ret_list = arg
        else:
            ret_list = [fline.rstrip("\r\n") for fline in open(arg, "r").readlines()]

        return ret_list

    exclusion_list = []
    if c.exclusionSubjectList is not None:
        exclusion_list = get_list(c.exclusionSubjectList)

    subject_list = []
    if c.subjectList is not None:
        subject_list = get_list(c.subjectList)

    def checkTemplate(template) -> None:
        """Check if `template` is correct."""
        if template.count("%s") != 2:
            msg = (
                "Please provide '%s' in the template"
                "where your site and subjects are present"
                "Please see examples"
            )
            logger.exception(msg)
            raise Exception(msg)

        filename, ext = os.path.splitext(os.path.basename(template))
        ext = os.path.splitext(filename)[1] + ext

        if ext not in [".nii", ".nii.gz"]:
            msg = "Invalid file name", os.path.basename(template)
            logger.exception(msg)
            raise Exception(msg)

    def get_site_list(path):
        base, relative = path.split("%s")
        return os.listdir(base)

    def check_length(scan_name, file_name):
        if len(file_name) > 30:
            msg = (
                "filename- %s is too long."
                "It should not be more than 30 characters." % (file_name)
            )
            logger.exception(msg)
            raise Exception(msg)

        if (
            len(scan_name) - len(os.path.splitext(os.path.splitext(file_name)[0])[0])
            >= 40
        ):
            msg = (
                "scan name %s is too long."
                "It should not be more than 20 characters"
                % (
                    scan_name.replace(
                        "_" + os.path.splitext(os.path.splitext(file_name)[0])[0], ""
                    )
                )
            )
            logger.exception(msg)
            raise Exception(msg)

    def create_site_subject_mapping(base, relative):
        """Create mapping between site and subject."""
        site_subject_map = {}
        base_path_list = []

        if c.siteList is not None:
            site_list = get_list(c.siteList)
        else:
            site_list = get_site_list(base)

        for site in site_list:
            paths = glob.glob(string.replace(base, "%s", site))
            base_path_list.extend(paths)
            for path in paths:
                for sub in os.listdir(path):
                    # check if subject is present in subject_list
                    if subject_list:
                        if sub in subject_list and sub not in exclusion_list:
                            site_subject_map[sub] = site
                    elif sub not in exclusion_list:
                        if sub not in ".DS_Store":
                            site_subject_map[sub] = site

        return base_path_list, site_subject_map

    def getPath(template):
        """Split the input template path...

        ...into base, path before subject directory and relative, path after subject directory.
        """
        checkTemplate(template)
        base, relative = template.rsplit("%s", 1)
        base, subject_map = create_site_subject_mapping(base, relative)
        base.sort()
        relative = relative.lstrip("/")
        return base, relative, subject_map

    # get anatomical base path and anatomical relative path
    anat_base, anat_relative = getPath(c.anatomicalTemplate)[:2]

    # get functional base path, functional relative path and site-subject map
    func_base, func_relative, subject_map = getPath(c.functionalTemplate)

    if not anat_base:
        msg = (
            "Anatomical Data template incorrect. No such file or directory %s",
            anat_base,
        )
        logger.exception(msg)
        raise Exception(msg)

    if not func_base:
        msg = "Functional Data template incorrect. No such file or directory %s, func_base"
        logger.exception(msg)
        raise Exception(msg)

    if len(anat_base) != len(func_base):
        msg1 = (
            "Some sites are missing, Please check your template",
            anat_base,
            "!=",
            func_base,
        )
        logger.exception(msg1)

        msg2 = (
            " Base length Unequal. Some sites are missing."
            "extract_data doesn't script support this.Please"
            "Provide your own subjects_list file"
        )
        logger.exception(msg2)
        raise Exception(msg2)

    # calculate the length of relative paths(path after subject directory)
    func_relative_len = len(func_relative.split("/"))
    anat_relative_len = len(anat_relative.split("/"))

    def check_for_sessions(relative_path, path_length):
        """Check if there are sessions present."""
        # default
        session_present = False
        session_path = "session_1"

        # session present if path_length is equal to 3
        if path_length == 3:  # noqa: PLR2004
            relative_path_list = relative_path.split("/")
            session_path = relative_path_list[0]
            relative_path = string.join(relative_path_list[1:], "/")
            session_present = True
        elif path_length > 3:  # noqa: PLR2004
            msg = (
                "extract_data script currently doesn't support this directory structure."
                "Please provide the subjects_list file to run CPAC."
                "For more information refer to manual"
            )
            logger.exception(msg)
            raise Exception(msg)
        return session_present, session_path, relative_path

    func_session_present, func_session_path, func_relative = check_for_sessions(
        func_relative, func_relative_len
    )

    anat_session_present, anat_session_path, anat_relative = check_for_sessions(
        anat_relative, anat_relative_len
    )

    f = open(
        os.path.join(
            c.outputSubjectListLocation, "CPAC_subject_list_%s.yml" % c.subjectListName
        ),
        "wb",
    )

    def fetch_path(i, anat_sub, func_sub, session_id):
        """
        Extract anatomical and functional path for a session and print to file.

        Parameters
        ----------
        i : int
            index of site
        anat_sub : string
            string containing subject/ concatenated
            subject-session path for anatomical file
        func_sub : string
            string containing subject/ concatenated
            subject-session path for functional file
        session_id : string
            session

        Raises
        ------
        Exception
        """
        try:

            def print_begin_of_file(sub, session_id):
                print("-", file=f)
                print("    subject_id: '" + sub + "'", file=f)
                print("    unique_id: '" + session_id + "'", file=f)

            def print_end_of_file(sub):
                if param_map is not None:
                    try:
                        logger.debug("site for sub %s -> %s", sub, subject_map.get(sub))
                        logger.debug(
                            "scan parameters for the above site %s",
                            param_map.get(subject_map.get(sub)),
                        )
                        print("    scan_parameters:", file=f)
                        print(
                            "        tr: '"
                            + param_map.get(subject_map.get(sub))[4]
                            + "'",
                            file=f,
                        )
                        print(
                            "        acquisition: '"
                            + param_map.get(subject_map.get(sub))[0]
                            + "'",
                            file=f,
                        )
                        print(
                            "        reference: '"
                            + param_map.get(subject_map.get(sub))[3]
                            + "'",
                            file=f,
                        )
                        print(
                            "        first_tr: '"
                            + param_map.get(subject_map.get(sub))[1]
                            + "'",
                            file=f,
                        )
                        print(
                            "        last_tr: '"
                            + param_map.get(subject_map.get(sub))[2]
                            + "'",
                            file=f,
                        )
                    except:
                        msg = (
                            " No Parameter values for the %s site is defined in the scan"
                            " parameters csv file" % subject_map.get(sub)
                        )
                        raise ValueError(msg)

            # get anatomical file
            anat_base_path = os.path.join(anat_base[i], anat_sub)
            func_base_path = os.path.join(func_base[i], func_sub)

            anat = None
            func = None

            anat = glob.glob(os.path.join(anat_base_path, anat_relative))
            func = glob.glob(os.path.join(func_base_path, func_relative))

            if anat and func:
                print_begin_of_file(anat_sub.split("/")[0], session_id)
                print("    anat: '" + os.path.realpath(anat[0]) + "'", file=f)
                print("    rest: ", file=f)

                # iterate for each rest session
                for _iter in func:
                    # get scan_id
                    iterable = os.path.splitext(
                        os.path.splitext(_iter.replace(func_base_path, "").lstrip("/"))[
                            0
                        ]
                    )[0]
                    iterable = iterable.replace("/", "_")
                    check_length(iterable, os.path.basename(os.path.realpath(_iter)))
                    print(
                        "      " + iterable + ": '" + os.path.realpath(_iter) + "'",
                        file=f,
                    )

                print_end_of_file(anat_sub.split("/")[0])

            else:
                logger.debug("skipping subject %s", anat_sub.split("/")[0])

        except ValueError:
            logger.exception(ValueError.message)
            raise

        except Exception as e:
            err_msg = (
                "Exception while felching anatomical and functional "
                "paths: \n" + str(e)
            )

            logger.exception(err_msg)
            raise Exception(err_msg)

    def walk(index, sub):
        """
        Walk across each subject path in the data site path.

        Parameters
        ----------
        index : int
            index of site
        sub : string
            subject_id

        Raises
        ------
        Exception
        """
        try:
            if func_session_present:
                # if there are sessions
                if "*" in func_session_path:
                    session_list = glob.glob(
                        os.path.join(
                            func_base[index], os.path.join(sub, func_session_path)
                        )
                    )
                else:
                    session_list = [func_session_path]

                if session_list:
                    for session in session_list:
                        session_id = os.path.basename(session)
                        if anat_session_present:
                            if func_session_path == anat_session_path:
                                fetch_path(
                                    index,
                                    os.path.join(sub, session_id),
                                    os.path.join(sub, session_id),
                                    session_id,
                                )
                            else:
                                fetch_path(
                                    index,
                                    os.path.join(sub, anat_session_path),
                                    os.path.join(sub, session_id),
                                    session_id,
                                )
                        else:
                            fetch_path(
                                index, sub, os.path.join(sub, session_id), session_id
                            )
                else:
                    logger.debug("Skipping subject %s", sub)

            else:
                logger.debug("No sessions")
                session_id = ""
                fetch_path(index, sub, sub, session_id)

        except Exception:
            logger.exception(Exception.message)
            raise

        except:
            err_msg = "Please make sessions are consistent across all subjects.\n\n"

            logger.exception(err_msg)
            raise Exception(err_msg)

    try:
        for i in range(len(anat_base)):
            for sub in os.listdir(anat_base[i]):
                # check if subject is present in subject_list
                if subject_list:
                    if sub in subject_list and sub not in exclusion_list:
                        logger.debug("extracting data for subject: %s", sub)
                        walk(i, sub)
                # check that subject is not in exclusion list
                elif sub not in exclusion_list and sub not in ".DS_Store":
                    logger.debug("extracting data for subject: %s", sub)
                    walk(i, sub)

        _name = os.path.join(c.outputSubjectListLocation, "CPAC_subject_list.yml")
        logger.info(
            "Extraction Successfully Completed...Input Subjects_list for CPAC - %s",
            _name,
        )

    except Exception:
        logger.exception(Exception.message)
        raise

    finally:
        f.close()


def generate_supplementary_files(data_config_outdir, data_config_name):
    """Generate phenotypic template file and subject list for group analysis."""
    import csv
    import os

    data_config_path = os.path.join(data_config_outdir, data_config_name)

    try:
        subjects_list = yaml.safe_load(open(data_config_path, "r"))
    except:
        f"\n\n[!] Data configuration file couldn't be read!\nFile path: {data_config_path}\n"

    subject_scan_set = set()
    subID_set = set()
    session_set = set()
    subject_set = set()
    scan_set = set()
    data_list = []

    try:
        for sub in subjects_list:
            if sub["unique_id"]:
                subject_id = sub["subject_id"] + "_" + sub["unique_id"]
            else:
                subject_id = sub["subject_id"]

            try:
                for scan in sub["func"]:
                    subject_scan_set.add((subject_id, scan))
                    subID_set.add(sub["subject_id"])
                    session_set.add(sub["unique_id"])
                    subject_set.add(subject_id)
                    scan_set.add(scan)
            except KeyError:
                try:
                    for scan in sub["rest"]:
                        subject_scan_set.add((subject_id, scan))
                        subID_set.add(sub["subject_id"])
                        session_set.add(sub["unique_id"])
                        subject_set.add(subject_id)
                        scan_set.add(scan)
                except KeyError:
                    # one of the participants in the subject list has no
                    # functional scans
                    subID_set.add(sub["subject_id"])
                    session_set.add(sub["unique_id"])
                    subject_set.add(subject_id)

    except TypeError:
        err_str = (
            "Subject list could not be populated!\nThis is most likely due to a"
            " mis-formatting in your inclusion and/or exclusion subjects txt file or"
            " your anatomical and/or functional path templates.\nCheck formatting of"
            " your anatomical/functional path templates and inclusion/exclusion"
            " subjects text files"
        )
        raise TypeError(err_str)

    for item in subject_scan_set:
        list1 = []
        list1.append(item[0] + "/" + item[1])
        for val in subject_set:
            if val in item:
                list1.append(1)
            else:
                list1.append(0)

        for val in scan_set:
            if val in item:
                list1.append(1)
            else:
                list1.append(0)

        data_list.append(list1)

    # generate the phenotypic file templates for group analysis
    file_name = os.path.join(
        data_config_outdir, "phenotypic_template_%s.csv" % data_config_name
    )

    f = _sassy_try_open_wb(file_name)

    writer = csv.writer(f)

    writer.writerow(["participant", "EV1", ".."])
    for sub in sorted(subID_set):
        writer.writerow([sub, ""])

    f.close()

    logger.info("Template Phenotypic file for group analysis - %s", file_name)

    """
    # generate the phenotypic file templates for repeated measures
    if (len(session_set) > 1) and (len(scan_set) > 1):

        file_name = os.path.join(data_config_outdir, 'phenotypic_template_repeated' \
                '_measures_mult_sessions_and_scans_%s.csv' \
                % data_config_name)

        f = _sassy_try_open_wb(file_name)

        writer = csv.writer(f)
        writer.writerow(['participant', 'session', 'series', 'EV1', '..'])

        for session in sorted(session_set):
            for scan in sorted(scan_set):
                for sub in sorted(subID_set):
                    writer.writerow([sub, session, scan, ''])

        f.close()

        logger.info(
            "Template Phenotypic file for group analysis with repeated "
            "measures (multiple sessions and scans) - %s", file_name
        )

    if (len(session_set) > 1):

        file_name = os.path.join(data_config_outdir, 'phenotypic_template_repeated' \
                '_measures_multiple_sessions_%s.csv' % data_config_name)

        f = _sassy_try_open_wb(file_name)

        writer = csv.writer(f)

        writer.writerow(['participant', 'session', 'EV1', '..'])

        for session in sorted(session_set):
            for sub in sorted(subID_set):
                writer.writerow([sub, session, ''])

        f.close()

        logger.info(
            "Template Phenotypic file for group analysis with repeated "
            "measures (multiple sessions) - %s", file_name
        )

    if (len(scan_set) > 1):

        file_name = os.path.join(data_config_outdir, 'phenotypic_template_repeated' \
                '_measures_multiple_scans_%s.csv' % data_config_name)

        f = _sassy_try_open_wb(file_name)

        writer = csv.writer(f)

        writer.writerow(['participant', 'series', 'EV1', '..'])

        for scan in sorted(scan_set):
            for sub in sorted(subID_set):
                writer.writerow([sub, scan, ''])

        f.close()

    logger.info("Template Phenotypic file for group analysis with repeated "
        "measures (multiple scans) - %s", file_name
    )
    """

    # generate the group analysis subject lists
    file_name = os.path.join(
        data_config_outdir, "participant_list_group_analysis_%s.txt" % data_config_name
    )

    try:
        with open(file_name, "w") as f:
            for sub in sorted(subID_set):
                print(sub, file=f)
    except:
        _sassy_oserror(file_name)

    logger.info(
        "Participant list required later for group analysis - %s\n\n", file_name
    )


def read_csv(csv_input):
    """Read CSV file.

    'Acquisition'
    'Reference'
    'Site'
    'TR (seconds)'
    """
    from collections import defaultdict
    import csv

    try:
        reader = csv.DictReader(open(csv_input, "U"))

        dict_labels = defaultdict(list)
        for line in reader:
            csv_dict = {k.lower(): v for k, v in line.items()}
            dict_labels[csv_dict.get("site")] = [
                csv_dict[key]
                for key in sorted(csv_dict.keys())
                if key not in ("site", "scan")
            ]

        if len(dict_labels) < 1:
            msg = "Scan Parameters File is either empty or missing header"
            logger.exception(msg)
            raise Exception(msg)

        return dict_labels

    except IOError:
        msg = "Error reading the csv file %s", csv_input
        logger.exception(msg)
        raise Exception(msg)
    except:
        msg = "Error reading scan parameters csv. Make sure you are using the correct template"
        logger.exception(msg)
        raise Exception(msg)


def _sassy_oserror(file_name: str) -> None:
    """Raise a sassy OSError."""
    msg = (
        f"\n\nCPAC says: I couldn't save this file to your drive:\n {file_name}"
        "\n\nMake sure you have write access? Then come back. Don't worry.. I'll"
        " wait.\n\n"
    )
    raise OSError(msg)


def _sassy_try_open_wb(file_name: str) -> Optional[BinaryIO]:
    """Open a file in 'wb' mode or raise a sassy OSError if a file can't be saved."""
    f = None
    try:
        f = open(file_name, "wb")
    except (OSError, TypeError):
        _sassy_oserror(file_name)
    return f


class Configuration(object):
    """Set dictionary keys as map attributes."""

    def __init__(self, config_map):
        for key in config_map:
            if config_map[key] == "None":
                config_map[key] = None
            setattr(self, key, config_map[key])


def run(data_config: Path | str) -> None:
    """Run a data config.

    Parameters
    ----------
    data_config : ~pathlib.Path or str
        path to data_config file
    """
    logger.info(
        "For any errors or messages check the log file - %s",
        os.path.join(os.getcwd(), "extract_data_logs.log"),
    )

    c = Configuration(yaml.safe_load(open(os.path.realpath(data_config), "r")))

    if c.scanParametersCSV is not None:
        read_csv(c.scanParametersCSV)
    else:
        logger.debug(
            "no scan parameters csv included\n"
            "make sure you turn off slice timing correction option\n"
            "in CPAC configuration\n"
        )

    generate_supplementary_files(c.outputSubjectListLocation, c.subjectListName)


if __name__ == "__main__":
    if len(sys.argv) != 2:  # noqa: PLR2004
        print("Usage: python extract_data.py data_config.yml")  # noqa T201
        sys.exit()
    else:
        run(sys.argv[1])
