#!/usr/bin/env python
# Copyright (C) 2018-2024  C-PAC Developers

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
"""Run C-PAC in a container."""

import argparse
import datetime
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Optional
from warnings import simplefilter

import yaml

from CPAC import __version__, license_notice
from CPAC.pipeline import AVAILABLE_PIPELINE_CONFIGS
from CPAC.pipeline.random_state import set_up_random_state
from CPAC.pipeline.schema import str_to_bool1_1
from CPAC.utils.bids_utils import (
    cl_strip_brackets,
    create_cpac_data_config,
    load_cpac_data_config,
    load_yaml_config,
    sub_list_filter_by_labels,
)
from CPAC.utils.configuration import Configuration, preconfig_yaml, set_subject
from CPAC.utils.configuration.yaml_template import (
    create_yaml_from_template,
    hash_data_config,
    upgrade_pipeline_to_1_8,
)
from CPAC.utils.docs import DOCS_URL_PREFIX
from CPAC.utils.monitoring import failed_to_start, FMLOGGER, log_nodes_cb, WFLOGGER
from CPAC.utils.utils import update_nested_dict

from bids2table import bids2table

simplefilter(action="ignore", category=FutureWarning)
DEFAULT_TMP_DIR = "/tmp"


def run(command: str, env: Optional[dict] = None) -> None:
    """Run a command in the shell."""
    if env is None:
        env = {}
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=env
    )
    while True:
        line = process.stdout.readline()
        line = line.decode()[:-1]
        if line == "" and process.poll() is not None:
            break


def parse_yaml(value: str) -> dict:
    """Parse a string as a YAML dictionary."""
    try:
        config = yaml.safe_load(value)
        if not isinstance(config, dict):
            msg = "config must be a dictionary"
            raise TypeError(msg)
        return config
    except Exception:
        # pylint: disable=raise-missing-from
        msg = f"Invalid configuration: '{value}'"
        raise argparse.ArgumentTypeError(msg)


def resolve_aws_credential(source: Path | str) -> str:
    """Set AWS credentials from a file or environment variable."""
    if source == "env":
        from urllib.request import urlopen

        aws_creds_address = "169.254.170.2{}".format(
            os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"]
        )
        aws_creds = urlopen(aws_creds_address).read()

        aws_input_creds = "/tmp/aws_input_creds_%d.csv" % int(round(time.time() * 1000))
        with open(aws_input_creds) as ofd:
            for key, vname in [
                ("AccessKeyId", "AWSAcessKeyId"),
                ("SecretAccessKey", "AWSSecretKey"),
            ]:
                ofd.write(f"{vname}={aws_creds[key]}")

        return aws_input_creds

    if os.path.isfile(source):
        return source
    msg = f"Could not find aws credentials {source}"
    raise IOError(msg)


def run_main():
    """Run this function if not importing as a script."""
    parser = argparse.ArgumentParser(
        description="C-PAC Pipeline Runner. " + license_notice
    )
    parser.add_argument(
        "bids_dir",
        help="The directory with the input dataset "
        "formatted according to the BIDS standard. "
        "Use the format s3://bucket/path/to/bidsdir to "
        "read data directly from an S3 bucket. This may "
        "require AWS S3 credentials specified via the "
        "--aws_input_creds option.",
    )
    parser.add_argument(
        "output_dir",
        help="The directory where the output files should be "
        "stored. If you are running group level analysis "
        "this folder should be prepopulated with the "
        "results of the participant level analysis. Use "
        "the format s3://bucket/path/to/bidsdir to "
        "write data directly to an S3 bucket. This may "
        "require AWS S3 credentials specified via the "
        "--aws_output_creds option.",
    )
    parser.add_argument(
        "analysis_level",
        help="Level of the analysis that will be performed. "
        "Multiple participant level analyses can be run "
        "independently (in parallel) using the same "
        "output_dir. test_config will run through the "
        "entire configuration process but will not "
        "execute the pipeline.",
        choices=["participant", "group", "test_config", "cli"],
        type=lambda choice: choice.replace("-", "_").lower(),
    )

    parser.add_argument(
        "--pipeline-file",
        "--pipeline_file",
        help="Path for the pipeline configuration file to "
        "use. Use the format s3://bucket/path/to/"
        "pipeline_file to read data directly from an "
        "S3 bucket. This may require AWS S3 credentials "
        "specified via the --aws_input_creds option.",
        default=preconfig_yaml("default"),
    )
    parser.add_argument(
        "--group-file",
        "--group_file",
        help="Path for the group analysis configuration file "
        "to use. Use the format s3://bucket/path/to/"
        "pipeline_file to read data directly from an S3 "
        "bucket. This may require AWS S3 credentials "
        "specified via the --aws_input_creds option. "
        "The output directory needs to refer to the "
        "output of a preprocessing individual pipeline.",
        default=None,
    )
    parser.add_argument(
        "--data-config-file",
        "--data_config_file",
        help="Yaml file containing the location of the data "
        "that is to be processed. This file is not "
        "necessary if the data in bids_dir is organized "
        "according to the BIDS format. This enables "
        "support for legacy data organization and cloud "
        "based storage. A bids_dir must still be "
        "specified when using this option, but its "
        "value will be ignored. Use the format s3://"
        "bucket/path/to/data_config_file to read data "
        "directly from an S3 bucket. This may require "
        "AWS S3 credentials specified via the "
        "--aws_input_creds option.",
        default=None,
    )

    parser.add_argument(
        "--preconfig",
        help="Name of the preconfigured pipeline to run. "
        "Available preconfigured pipelines: "
        + str(AVAILABLE_PIPELINE_CONFIGS)
        + ". See "
        f"{DOCS_URL_PREFIX}/user/pipelines/preconfig "
        "for more information about the preconfigured "
        "pipelines.",
        default=None,
    )
    if [
        _ for _ in ["--pipeline-override", "--pipeline_override"] if _ in sys.argv
    ]:  # secret option
        parser.add_argument(
            "--pipeline-override",
            "--pipeline_override",
            type=parse_yaml,
            action="append",
            help="Override specific options from the "
            "pipeline configuration. E.g.: "
            "\"{'pipeline_setup': {'system_config': "
            "{'maximum_memory_per_participant': 1}}}\"",
        )

    parser.add_argument(
        "--aws-input-creds",
        "--aws_input_creds",
        help="Credentials for reading from S3. If not "
        "provided and s3 paths are specified in the "
        "data config we will try to access the bucket "
        'anonymously use the string "env" to indicate '
        "that input credentials should read from the "
        "environment. (E.g. when using AWS iam roles).",
        default=None,
    )
    parser.add_argument(
        "--aws-output-creds",
        "--aws_output_creds",
        help="Credentials for writing to S3. If not provided "
        "and s3 paths are specified in the output "
        "directory we will try to access the bucket "
        'anonymously use the string "env" to indicate '
        "that output credentials should read from the "
        "environment. (E.g. when using AWS iam roles).",
        default=None,
    )
    # TODO: restore <default=3> for <--n_cpus> once we remove
    #       <max_cores_per_participant> from config file
    #       <https://github.com/FCP-INDI/C-PAC/pull/1264#issuecomment-631643708>
    parser.add_argument(
        "--n-cpus",
        "--n_cpus",
        type=int,
        default=0,
        help="Number of execution resources per participant "
        "available for the pipeline. This flag takes "
        "precidence over max_cores_per_participant in "
        "the pipeline configuration file.",
    )
    parser.add_argument(
        "--mem-mb",
        "--mem_mb",
        type=float,
        help="Amount of RAM available per participant in "
        "megabytes. Included for compatibility with "
        "BIDS-Apps standard, but mem_gb is preferred. "
        "This flag takes precedence over "
        "maximum_memory_per_participant in the pipeline "
        "configuration file.",
    )
    parser.add_argument(
        "--mem-gb",
        "--mem_gb",
        type=float,
        help="Amount of RAM available per participant in "
        "gigabytes. If this is specified along with "
        "mem_mb, this flag will take precedence. This "
        "flag also takes precedence over "
        "maximum_memory_per_participant in the pipeline "
        "configuration file.",
    )
    parser.add_argument(
        "--runtime-usage",
        "--runtime_usage",
        type=str,
        help="Path to a callback.log from a prior run of the "
        "same pipeline configuration (including any "
        "resource-management parameters that will be "
        "applied in this run, like 'n_cpus' and "
        "'num_ants_threads'). This log will be used to "
        "override per-node memory estimates with "
        "observed values plus a buffer.",
    )
    parser.add_argument(
        "--runtime-buffer",
        "--runtime_buffer",
        type=float,
        help="Buffer to add to per-node memory estimates if "
        "--runtime_usage is specified. This number is a "
        "percentage of the observed memory usage.",
    )
    parser.add_argument(
        "--num-ants-threads",
        "--num_ants_threads",
        type=int,
        default=0,
        help="The number of cores to allocate to ANTS-"
        "based anatomical registration per "
        "participant. Multiple cores can greatly "
        "speed up this preprocessing step. This "
        "number cannot be greater than the number of "
        "cores per participant.",
    )
    parser.add_argument(
        "--random-seed",
        "--random_seed",
        type=str,
        help="Random seed used to fix the state of execution. "
        "If unset, each process uses its own default. If "
        "set, a `random.log` file will be generated "
        "logging the random state used by each process. "
        "If set to a positive integer (up to 2147483647"
        "), that integer will be used to seed each "
        "process. If set to 'random', a random seed "
        "will be generated and recorded for each "
        "process.",
    )
    parser.add_argument(
        "--save-working-dir",
        "--save_working_dir",
        nargs="?",
        help="Save the contents of the working directory.",
        default=False,
    )
    parser.add_argument(
        "--fail-fast",
        "--fail_fast",
        type=str.title,
        help="Stop worklow execution on first crash?",
    )
    parser.add_argument(
        "--participant-label",
        "--participant_label",
        help="The label of the participant that should be "
        "analyzed. The label corresponds to "
        "sub-<participant_label> from the BIDS spec "
        '(so it does not include "sub-"). If this '
        "parameter is not provided all participants "
        "should be analyzed. Multiple participants "
        "can be specified with a space separated "
        "list.",
        nargs="+",
    )
    parser.add_argument(
        "--participant-ndx",
        "--participant_ndx",
        help="The index of the participant that should be "
        "analyzed. This corresponds to the index of "
        "the participant in the data config file. "
        "This was added to make it easier to "
        "accommodate SGE array jobs. Only a single "
        "participant will be analyzed. Can be used "
        "with participant label, in which case it is "
        "the index into the list that follows the "
        'participant_label flag. Use the value "-1" '
        "to indicate that the participant index "
        "should be read from the "
        "AWS_BATCH_JOB_ARRAY_INDEX environment "
        "variable.",
        default=None,
        type=int,
    )

    parser.add_argument(
        "--T1w-label",
        "--T1w_label",
        help="C-PAC only runs one T1w per participant-"
        "session at a time, at this time. Use this "
        'flag to specify any BIDS entity (e.g., "acq-'
        'VNavNorm") or sequence of BIDS entities ('
        'e.g., "acq-VNavNorm_run-1") to specify '
        "which of multiple T1w files to use. Specify "
        '"--T1w_label T1w" to choose the T1w file '
        "with the fewest BIDS entities (i.e., the "
        "final option of [*_acq-VNavNorm_T1w.nii.gz, "
        '*_acq-HCP_T1w.nii.gz, *_T1w.nii.gz"]). '
        "C-PAC will choose the first T1w it finds if "
        "the user does not provide this flag, or "
        "if multiple T1w files match the --T1w_label "
        "provided.\nIf multiple T2w files are present "
        "and a comparable filter is possible, T2w "
        "files will be filtered as well. If no T2w files "
        "match this --T1w_label, T2w files will be "
        "processed as if no --T1w_label were provided.",
    )
    parser.add_argument(
        "--bold-label",
        "--bold_label",
        help="To include a specified subset of available "
        "BOLD files, use this flag to specify any "
        'BIDS entity (e.g., "task-rest") or sequence '
        'of BIDS entities (e.g. "task-rest_run-1"). '
        "To specify the bold file with the fewest "
        "BIDS entities in the file name, specify "
        '"--bold_label bold". Multiple `--bold_'
        "label`s can be specified with a space-"
        "separated list. If multiple `--bold_label`s "
        'are provided (e.g., "--bold_label task-rest_'
        'run-1 task-rest_run-2", each scan that '
        "includes all BIDS entities specified in any "
        "of the provided `--bold_label`s will be "
        "analyzed. If this parameter is not provided "
        "all BOLD scans should be analyzed.",
        nargs="+",
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"C-PAC BIDS-App version {__version__}",
    )
    parser.add_argument(
        "--bids-validator-config",
        "--bids_validator_config",
        help="JSON file specifying configuration of "
        "bids-validator: See https://github.com/bids-"
        "standard/bids-validator for more info.",
    )
    parser.add_argument(
        "--skip-bids-validator",
        "--skip_bids_validator",
        help="Skips bids validation.",
        action="store_true",
    )

    parser.add_argument(
        "--anat-only",
        "--anat_only",
        help="run only the anatomical preprocessing",
        action="store_true",
    )

    parser.add_argument(
        "--user_defined",
        type=str,
        help="Arbitrary user defined string that will be "
        "included in every output sidecar file.",
    )

    parser.add_argument(
        "--tracking-opt-out",
        "--tracking_opt-out",
        action="store_true",
        help="Disable usage tracking. Only the number of "
        "participants on the analysis is tracked.",
        default=False,
    )

    parser.add_argument(
        "--monitoring",
        help="Enable monitoring server on port 8080. You "
        "need to bind the port using the Docker "
        'flag "-p".',
        action="store_true",
    )

    # get the command line arguments
    args = parser.parse_args(
        sys.argv[1 : (sys.argv.index("--") if "--" in sys.argv else len(sys.argv))]
    )

    bids_dir_is_s3 = args.bids_dir.lower().startswith("s3://")
    bids_dir = args.bids_dir if bids_dir_is_s3 else os.path.realpath(args.bids_dir)
    output_dir_is_s3 = args.output_dir.lower().startswith("s3://")
    output_dir = (
        args.output_dir if output_dir_is_s3 else os.path.realpath(args.output_dir)
    )
    exitcode = 0
    if args.analysis_level == "cli":
        from CPAC.__main__ import main

        main.main(args=sys.argv[sys.argv.index("--") + 1 :])
        sys.exit(0)

    elif args.analysis_level == "group":
        if not args.group_file or not os.path.exists(args.group_file):
            import pkg_resources as p

            WFLOGGER.warning("\nNo group analysis configuration file was supplied.\n")

            args.group_file = p.resource_filename(
                "CPAC",
                os.path.join("resources", "configs", "group_config_template.yml"),
            )

            output_group = os.path.join(output_dir, "group_config.yml")

            try:
                if output_dir.lower().startswith("s3://"):
                    raise Exception

                if not os.path.exists(output_group):
                    shutil.copyfile(args.group_file, output_group)
            except (Exception, IOError):
                FMLOGGER.warning(
                    "Could not create group analysis configuration file.\nPlease refer to the C-PAC documentation for group analysis setup."
                )
            else:
                WFLOGGER.warning(
                    "Please refer to the output directory for a template of  the file and, after customizing to your analysis, add the flag\n\n    --group_file %s\n\nto your `docker run` command\n",
                    output_group,
                )

            sys.exit(1)

        else:
            import CPAC.pipeline.cpac_group_runner as cgr

            WFLOGGER.info(
                "Starting group level analysis of data in %s using %s",
                bids_dir,
                args.group_file,
            )
            cgr.run(args.group_file)

            sys.exit(0)

    elif args.analysis_level in ["test_config", "participant"]:
        # check to make sure that the input directory exists
        if (
            not args.data_config_file
            and not bids_dir_is_s3
            and not os.path.exists(bids_dir)
        ):
            msg = f"Error! Could not find {bids_dir}"
            raise FileNotFoundError(msg)

        # check to make sure that the output directory exists
        if not output_dir_is_s3 and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except Exception as e:
                msg = f"Error! Could not find/create output dir {output_dir}"
                raise FileNotFoundError(msg) from e

        # validate input dir (if skip_bids_validator is not set)
        if not args.data_config_file:
            if args.bids_validator_config:
                WFLOGGER.info("Running BIDS validator...")
                run(f"bids-validator --config {args.bids_validator_config} {bids_dir}")
            elif args.skip_bids_validator:
                WFLOGGER.info("Skipping BIDS validator...")
            elif bids_dir_is_s3:
                WFLOGGER.info("Skipping BIDS validator for S3 datasets...")
            else:
                WFLOGGER.info("Running BIDS validator...")
                run(f"bids-validator {bids_dir}")

        if args.preconfig:
            args.pipeline_file = preconfig_yaml(args.preconfig)

        # otherwise, if we are running group, participant, or dry run we
        # begin by conforming the configuration
        if isinstance(args.pipeline_file, dict):
            c = args.pipeline_file
        else:
            c = load_yaml_config(args.pipeline_file, args.aws_input_creds)

        if "pipeline_setup" not in c:
            _url = f"{DOCS_URL_PREFIX}/user/pipelines/1.7-1.8-nesting-mappings"

            WFLOGGER.warning(
                "\nC-PAC changed its pipeline configuration "
                "format in v1.8.0.\nSee %s for details.\n",
                _url,
            )

            updated_config = os.path.join(
                output_dir, "updated_config", os.path.basename(args.pipeline_file)
            )
            os.makedirs(os.path.join(output_dir, "updated_config"), exist_ok=True)

            open(updated_config, "w").write(yaml.dump(c))

            upgrade_pipeline_to_1_8(updated_config)
            c = load_yaml_config(updated_config, args.aws_input_creds)

        overrides = {}
        if hasattr(args, "pipeline_override") and args.pipeline_override:
            overrides = {k: v for d in args.pipeline_override for k, v in d.items()}
            c = update_nested_dict(c, overrides)

        if args.anat_only:
            c = update_nested_dict(c, {"FROM": "anat-only"})

        if args.user_defined:
            c["pipeline_setup"]["output_directory"]["user_defined"] = args.user_defined

        c = Configuration(c)

        # get the aws_input_credentials, if any are specified
        if args.aws_input_creds:
            c["awsCredentialsFile"] = resolve_aws_credential(args.aws_input_creds)

        if args.aws_output_creds:
            c["pipeline_setup"]["Amazon-AWS"]["aws_output_bucket_credentials"] = (
                resolve_aws_credential(args.aws_output_creds)
            )

        c["pipeline_setup"]["output_directory"]["path"] = os.path.join(
            output_dir, "output"
        )

        if not output_dir_is_s3:
            c["pipeline_setup"]["log_directory"]["path"] = os.path.join(
                output_dir, "log"
            )
        else:
            c["pipeline_setup"]["log_directory"]["path"] = os.path.join(
                DEFAULT_TMP_DIR, "log"
            )

        if args.mem_gb:
            c["pipeline_setup"]["system_config"]["maximum_memory_per_participant"] = (
                float(args.mem_gb)
            )
        elif args.mem_mb:
            c["pipeline_setup"]["system_config"]["maximum_memory_per_participant"] = (
                float(args.mem_mb) / 1024.0
            )
        else:
            try:
                c[
                    "pipeline_setup", "system_config", "maximum_memory_per_participant"
                ] = float(
                    c[
                        "pipeline_setup",
                        "system_config",
                        "maximum_memory_per_participant",
                    ]
                )
            except KeyError:
                c[
                    "pipeline_setup", "system_config", "maximum_memory_per_participant"
                ] = 6.0

        # Preference: n_cpus if given, override if present, else from config if
        # present, else n_cpus=3
        if int(args.n_cpus) == 0:
            try:
                args.n_cpus = c[
                    "pipeline_setup", "system_config", "max_cores_per_participant"
                ]
            except KeyError:
                args.n_cpus = 3
        c["pipeline_setup", "system_config", "max_cores_per_participant"] = int(
            args.n_cpus
        )

        c["pipeline_setup"]["system_config"]["num_participants_at_once"] = int(
            c["pipeline_setup"]["system_config"].get("num_participants_at_once", 1)
        )
        # Reduce cores per participant if cores times participants is more than
        # available CPUS. n_cpus is a hard upper limit.
        if (
            c["pipeline_setup"]["system_config"]["max_cores_per_participant"]
            * c["pipeline_setup"]["system_config"]["num_participants_at_once"]
        ) > int(args.n_cpus):
            c["pipeline_setup"]["system_config"]["max_cores_per_participant"] = (
                int(args.n_cpus)
                // c["pipeline_setup"]["system_config"]["num_participants_at_once"]
            )
            if c["pipeline_setup"]["system_config"]["max_cores_per_participant"] == 0:
                c["pipeline_setup"]["system_config"]["max_cores_per_participant"] = (
                    args.n_cpus
                )
                c["pipeline_setup"]["system_config"]["num_participants_at_once"] = 1

        if int(args.num_ants_threads) == 0:
            try:
                args.num_ants_threads = c[
                    "pipeline_setup", "system_config", "num_ants_threads"
                ]
            except KeyError:
                args.num_ants_threads = 3
        c["pipeline_setup", "system_config", "num_ants_threads"] = int(
            args.num_ants_threads
        )

        c["pipeline_setup"]["system_config"]["num_ants_threads"] = min(
            c["pipeline_setup"]["system_config"]["max_cores_per_participant"],
            int(c["pipeline_setup"]["system_config"]["num_ants_threads"]),
        )

        if args.random_seed:
            c["pipeline_setup"]["system_config"]["random_seed"] = args.random_seed

        if c["pipeline_setup"]["system_config"]["random_seed"] is not None:
            c["pipeline_setup"]["system_config"]["random_seed"] = set_up_random_state(
                c["pipeline_setup"]["system_config"]["random_seed"]
            )

        if args.runtime_usage is not None:
            c["pipeline_setup"]["system_config"]["observed_usage"]["callback_log"] = (
                args.runtime_usage
            )
        if args.runtime_buffer is not None:
            c["pipeline_setup"]["system_config"]["observed_usage"]["buffer"] = (
                args.runtime_buffer
            )

        if args.save_working_dir is not False:
            c["pipeline_setup"]["working_directory"]["remove_working_dir"] = False
        if isinstance(args.save_working_dir, str):
            c["pipeline_setup"]["working_directory"]["path"] = os.path.abspath(
                args.save_working_dir
            )
        elif not output_dir_is_s3:
            c["pipeline_setup"]["working_directory"]["path"] = os.path.join(
                output_dir, "working"
            )
        else:
            FMLOGGER.warning(
                "Cannot write working directory to S3 bucket. "
                "Either change the output directory to something "
                "local or turn off the --save_working_dir flag"
            )

        if args.fail_fast is not None:
            c["pipeline_setup", "system_config", "fail_fast"] = str_to_bool1_1(
                args.fail_fast
            )

        if c["pipeline_setup"]["output_directory"]["quality_control"][
            "generate_xcpqc_files"
        ]:
            c["functional_preproc"]["motion_estimates_and_correction"][
                "motion_estimates"
            ]["calculate_motion_first"] = True
            c["functional_preproc"]["motion_estimates_and_correction"][
                "motion_estimates"
            ]["calculate_motion_after"] = True

        if args.participant_label:
            WFLOGGER.info(
                "#### Running C-PAC for %s", ", ".join(args.participant_label)
            )
        else:
            WFLOGGER.info("#### Running C-PAC")

        WFLOGGER.info(
            "Number of participants to run in parallel: %s",
            c["pipeline_setup", "system_config", "num_participants_at_once"],
        )

        if not args.data_config_file:
            WFLOGGER.info("Input directory: %s", bids_dir)

        WFLOGGER.info(
            "Output directory: %s\nWorking directory: %s\nLog directory: %s\n"
            "Remove working directory: %s\nAvailable memory: %s (GB)\n"
            "Available threads: %s\nNumber of threads for ANTs: %s",
            c["pipeline_setup", "output_directory", "path"],
            c["pipeline_setup", "working_directory", "path"],
            c["pipeline_setup", "log_directory", "path"],
            c["pipeline_setup", "working_directory", "remove_working_dir"],
            c["pipeline_setup", "system_config", "maximum_memory_per_participant"],
            c["pipeline_setup", "system_config", "max_cores_per_participant"],
            c["pipeline_setup", "system_config", "num_ants_threads"],
        )

        # create a timestamp for writing config files
        # pylint: disable=invalid-name
        st = datetime.datetime.now().strftime("%Y-%m-%dT%H-%M-%SZ")

        if args.participant_label:
            args.participant_label = cl_strip_brackets(args.participant_label)
            args.participant_label = [
                "sub-" + pt if not pt.startswith("sub-") else pt
                for pt in args.participant_label
            ]

        # otherwise we move on to conforming the data configuration
        if not args.data_config_file:
            sub_list = create_cpac_data_config(
                bids_dir,
                args.participant_label,
                args.aws_input_creds,
                args.skip_bids_validator,
                only_one_anat=False,
            )
            # Initializing the bidstable on the bids_directory
            bids_table = bids2table(bids_dir, workers=10)

            try:
                # fillna
                bids_table['ent__ses'] = bids_table['ent__ses'].fillna('None')
                grouped_tab = bids_table.groupby(["ent__sub", "ent__ses"])
            except Exception as e:
                WFLOGGER.warning("Could not create bids table: %s", e)
                print("Could not create bids table: %s", e)
                sys.exit(1)
        # else:
        #     sub_list = load_cpac_data_config(
        #         args.data_config_file, args.participant_label, args.aws_input_creds
        #     )
        list(sub_list)
        sub_list = sub_list_filter_by_labels(
            sub_list, {"T1w": args.T1w_label, "bold": args.bold_label}
        )

        # C-PAC only handles single anatomical images (for now)
        # so we take just the first as a string if we have a list
        for i, sub in enumerate(sub_list):
            if isinstance(sub.get("anat"), dict):
                for anat_key in sub["anat"]:
                    if isinstance(sub["anat"][anat_key], list) and len(
                        sub["anat"][anat_key]
                    ):
                        sub_list[i]["anat"][anat_key] = sub["anat"][anat_key][0]
            if isinstance(sub.get("anat"), list) and len(sub["anat"]):
                sub_list[i]["anat"] = sub["anat"][0]

        if args.participant_ndx is not None:
            participant_ndx = int(args.participant_ndx)
            if participant_ndx == -1:
                args.participant_ndx = os.environ["AWS_BATCH_JOB_ARRAY_INDEX"]

            if 0 <= participant_ndx < len(sub_list):
                WFLOGGER.info(
                    "Processing data for participant %s (%s)",
                    args.participant_ndx,
                    sub_list[participant_ndx]["subject_id"],
                )
                sub_list = [sub_list[participant_ndx]]
                data_hash = hash_data_config(sub_list)
                data_config_file = (
                    f"cpac_data_config_{data_hash}_idx-"
                    f"{args.participant_ndx}_{st}.yml"
                )
            else:
                msg = f"Participant ndx {participant_ndx} is out of bounds [0, {len(sub_list)})"
                raise IndexError(msg)
        else:
            data_hash = hash_data_config(sub_list)
            data_config_file = f"cpac_data_config_{data_hash}_{st}.yml"

        sublogdirs = [set_subject(sub, c)[2] for sub in grouped_tab]
        # write out the data configuration file
        data_config_file = os.path.join(sublogdirs[0], data_config_file)
        with open(data_config_file, "w", encoding="utf-8") as _f:
            noalias_dumper = yaml.dumper.SafeDumper
            noalias_dumper.ignore_aliases = lambda self, data: True
            yaml.dump(sub_list, _f, default_flow_style=False, Dumper=noalias_dumper)

        # update and write out pipeline config file
        pipeline_config_file = os.path.join(
            sublogdirs[0], f"cpac_pipeline_config_{data_hash}_{st}.yml"
        )
        with open(pipeline_config_file, "w", encoding="utf-8") as _f:
            _f.write(create_yaml_from_template(c))
        minimized_config = f"{pipeline_config_file[:-4]}_min.yml"
        with open(minimized_config, "w", encoding="utf-8") as _f:
            _f.write(create_yaml_from_template(c, import_from="blank"))
        for config_file in (data_config_file, pipeline_config_file, minimized_config):
            os.chmod(config_file, 0o444)  # Make config files readonly

        if len(sublogdirs) > 1:
            # If more than one run is included in the given data config
            # file, an identical copy of the data and pipeline config
            # will be included in the log directory for each run
            for sublogdir in sublogdirs[1:]:
                for config_file in (
                    data_config_file,
                    pipeline_config_file,
                    minimized_config,
                ):
                    try:
                        os.link(
                            config_file, config_file.replace(sublogdirs[0], sublogdir)
                        )
                    except FileExistsError:
                        pass

        if args.analysis_level in ["participant", "test_config"]:
            # build pipeline easy way
            import CPAC.pipeline.cpac_runner
            from CPAC.utils.monitoring import monitor_server

            monitoring = None
            if args.monitoring:
                from json import JSONDecodeError

                try:
                    monitoring = monitor_server(
                        c["pipeline_setup"]["pipeline_name"],
                        c["pipeline_setup"]["log_directory"]["path"],
                    )
                except (
                    AttributeError,
                    FileNotFoundError,
                    JSONDecodeError,
                    KeyError,
                    OSError,
                    PermissionError,
                    TypeError,
                    ValueError,
                ) as e:
                    WFLOGGER.warning(
                        "The run will continue without monitoring. Monitoring was configured to be enabled, but the monitoring server failed to start, so : %s\n",
                        e,
                    )

            plugin_args = {
                "n_procs": int(
                    c["pipeline_setup"]["system_config"]["max_cores_per_participant"]
                ),
                "memory_gb": int(
                    c["pipeline_setup"]["system_config"][
                        "maximum_memory_per_participant"
                    ]
                ),
                "raise_insufficient": c["pipeline_setup"]["system_config"][
                    "raise_insufficient"
                ],
                "status_callback": log_nodes_cb,
            }
            if (
                c["pipeline_setup"]["system_config"]["observed_usage"]["callback_log"]
                is not None
            ):
                plugin_args["runtime"] = {
                    "usage": c["pipeline_setup"]["system_config"]["observed_usage"][
                        "callback_log"
                    ],
                    "buffer": c["pipeline_setup"]["system_config"]["observed_usage"][
                        "buffer"
                    ],
                }

            WFLOGGER.info("Starting participant level processing")
            exitcode = CPAC.pipeline.cpac_runner.run(
                grouped_tab,
                pipeline_config_file,
                plugin="MultiProc" if plugin_args["n_procs"] > 1 else "Linear",
                plugin_args=plugin_args,
                tracking=not args.tracking_opt_out,
                test_config=args.analysis_level == "test_config",
            )

            if monitoring:
                monitoring.join(10)

            if args.analysis_level == "test_config":
                if exitcode == 0:
                    WFLOGGER.info(
                        "\nPipeline and data configuration files should"
                        " have been written to %s and %s respectively.\n",
                        pipeline_config_file,
                        data_config_file,
                    )

            # wait to import `LOGTAIL` here so it has any runtime updates
            from CPAC.utils.monitoring import LOGTAIL

            for warning in LOGTAIL["warnings"]:
                WFLOGGER.warning("%s\n", warning.rstrip())

    sys.exit(exitcode)


if __name__ == "__main__":
    try:
        run_main()
    except Exception as exception:
        # if we hit an exception before the pipeline starts to build but
        # we're still able to create a logfile, log the error in the file
        failed_to_start(
            sys.argv[2]
            if len(sys.argv) > 2  # noqa: PLR2004
            else os.getcwd(),
            exception,
        )
        raise exception
