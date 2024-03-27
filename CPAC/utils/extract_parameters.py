# Copyright (C) 2013-2024  C-PAC Developers

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
def merge(output_dir, scan_name, threshold, motion_f, power_f, flag):
    """Merge power parameters and motion parameters file."""
    import os
    import re

    if threshold is None:
        filename = scan_name + "_all_params.csv"
        filename = filename.lstrip("_")
        outfile = os.path.join(output_dir, filename)
    else:
        filename = scan_name + threshold + "_all_params.csv"
        filename = filename.lstrip("_")
        outfile = os.path.join(output_dir, filename)
        float(re.sub(r"[a-zA-Z_]", "", threshold))

    # Read in the motion and power parameters files
    try:
        motion = open(motion_f, "r").readlines()
    except Exception as e:
        err_string = (
            "\n\n[!] CPAC says: Could not read the motion "
            "parameters file.\n\nFilepath: %s\n\nError details: %s"
            "\n\n" % (motion_f, e)
        )
        raise Exception(err_string)

    try:
        power = open(power_f, "r").readlines()
    except Exception as e:
        err_string = (
            "\n\n[!] CPAC says: Could not read the power "
            "parameters file.\n\nFilepath: %s\n\nError details: %s"
            "\n\n" % (power_f, e)
        )
        raise Exception(err_string)

    # Write the combined motion and power parameters CSV file
    try:
        if flag:
            f = open(outfile, "w")

            m = motion[0].strip("\n")
            p = ",".join(power[0].split(",")[1:])

            f.write(m + p)
        else:
            f = open(outfile, "a")

        m = motion[1]
        p = ",".join(power[1].split(",")[2:])

        f.write(m + p + "\n")
        f.close()
    except Exception as e:
        err_string = (
            "\n\n[!] CPAC says: Could not create or open the motion "
            "and power parameters CSV file. Ensure you have write "
            "permissions for the directory it is writing to.\n\n"
            "Attempted write path: %s\n\nError details: %s\n\n" % (outfile, e)
        )
        raise Exception(err_string)


def grab(output_dir, scrubbing):
    """
    Method to grab all the motion parameters
    and power parameters file from each subject
    for each pipeline and merge them.

    Parameters
    ----------
    output_dir : string
        Path to the datasink output directory of CPAC
    """
    import glob
    import os
    import re

    pipelines = glob.glob(os.path.join(output_dir, "pipeline*"))

    for p in pipelines:
        scan_list = []
        threshold_list = []

        pattern1 = re.compile(r"(\w)*scan(\w)*(\d)*(\w)*[/]")
        pattern2 = re.compile(r"(\w)*threshold_[-+]?([0-9]*\.[0-9]+|[0-9]+)")

        scans = glob.glob(os.path.join(p, "*/power_params/*/*"))

        # get the unique scans and threshold value
        for s in scans:
            val = re.search(pattern1, s)
            if val:
                scan_list.append(val.group(0).rstrip("/"))

            val = re.search(pattern2, s)
            if val:
                threshold_list.append(val.group(0))

        scan_list = set(scan_list)
        threshold_list = set(threshold_list)

        for scan in scan_list:
            for threshold in threshold_list:
                Flag = 1
                # merge files for each subject
                for sub in os.listdir(p):
                    sub = os.path.join(p, sub)
                    motion_file = os.path.join(
                        sub, "motion_params", scan, "motion_parameters.txt"
                    )
                    power_file = os.path.join(
                        sub, "power_params", scan, threshold, "pow_params.txt"
                    )
                    if os.path.exists(motion_file) and os.path.exists(power_file):
                        merge(p, scan, threshold, motion_file, power_file, Flag)
                        Flag = 0

            if 0 in scrubbing:
                for sub in os.listdir(p):
                    sub = os.path.join(p, sub)
                    motion_file = os.path.join(
                        sub, "motion_params", scan, "motion_parameters.txt"
                    )
                    power_file = os.path.join(
                        sub, "power_params", scan, "pow_params.txt"
                    )

                    if os.path.exists(motion_file) and os.path.exists(power_file):
                        threshold = None
                        merge(p, scan, threshold, motion_file, power_file, Flag)
                        Flag = 0

    return threshold


def run(output_path, scrubbing):
    return grab(output_path, scrubbing)


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 2:
        grab(sys.argv[1], [0])
    else:
        print("Usage: python extract_parameters.py /path/to/output/dir")  # noqa: T201
