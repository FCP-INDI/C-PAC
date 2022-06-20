#!/usr/bin/env python
"""Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>."""
import argparse
import datetime
import os
import subprocess
import sys
import time
import shutil
import yaml

from CPAC import license_notice, __version__
from CPAC.pipeline import AVAILABLE_PIPELINE_CONFIGS
from CPAC.pipeline.random_state import set_up_random_state
from CPAC.utils.bids_utils import create_cpac_data_config, \
                                  load_cpac_data_config, \
                                  load_yaml_config, \
                                  sub_list_filter_by_labels
from CPAC.utils.configuration import Configuration
from CPAC.utils.docs import DOCS_URL_PREFIX
from CPAC.utils.monitoring import log_nodes_cb
from CPAC.utils.yaml_template import create_yaml_from_template, \
                                     upgrade_pipeline_to_1_8
from CPAC.utils.utils import cl_strip_brackets, load_preconfig, \
                             update_nested_dict

import yamlordereddictloader
from warnings import simplefilter, warn
simplefilter(action='ignore', category=FutureWarning)

DEFAULT_TMP_DIR = "/tmp"
DEFAULT_PIPELINE = "/cpac_resources/default_pipeline.yml"
if not os.path.exists(DEFAULT_PIPELINE):
    DEFAULT_PIPELINE = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "default_pipeline.yml"
    )


def run(command, env={}):
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               shell=True, env=env)
    while True:
        line = process.stdout.readline()
        line = line.decode()[:-1]
        if line == '' and process.poll() is not None:
            break


def parse_yaml(value):
    try:
        config = yaml.safe_load(value)
        if type(config) != dict:
            raise
        return config
    except:
        raise argparse.ArgumentTypeError("Invalid configuration: '%s'" % value)


def resolve_aws_credential(source):

    if source == "env":
        from urllib.request import urlopen
        aws_creds_address = "169.254.170.2{}".format(
            os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"]
        )
        aws_creds = urlopen(aws_creds_address).read()

        aws_input_creds = "/tmp/aws_input_creds_%d.csv" % int(
            round(time.time() * 1000)
        )
        with open(aws_input_creds) as ofd:
            for key, vname in [
                ("AccessKeyId", "AWSAcessKeyId"),
                ("SecretAccessKey", "AWSSecretKey")
            ]:
                ofd.write("{0}={1}".format(vname, aws_creds[key]))

        return aws_input_creds

    if os.path.isfile(source):
        return source
    else:
        raise IOError(
            "Could not find aws credentials {0}"
            .format(source)
        )


def run_main():
    """Run this function if not importing as a script"""
    parser = argparse.ArgumentParser(description='C-PAC Pipeline Runner. ' +
                                     license_notice)
    parser.add_argument('bids_dir',
                        help='The directory with the input dataset '
                             'formatted according to the BIDS standard. '
                             'Use the format s3://bucket/path/to/bidsdir to '
                             'read data directly from an S3 bucket. This may '
                             'require AWS S3 credentials specified via the '
                             '--aws_input_creds option.')
    parser.add_argument('output_dir',
                        help='The directory where the output files should be '
                             'stored. If you are running group level analysis '
                             'this folder should be prepopulated with the '
                             'results of the participant level analysis. Use '
                             'the format s3://bucket/path/to/bidsdir to '
                             'write data directly to an S3 bucket. This may '
                             'require AWS S3 credentials specified via the '
                             '--aws_output_creds option.')
    parser.add_argument('analysis_level',
                        help='Level of the analysis that will be performed. '
                             'Multiple participant level analyses can be run '
                             'independently (in parallel) using the same '
                             'output_dir. test_config will run through the '
                             'entire configuration process but will not '
                             'execute the pipeline.',
                        choices=['participant', 'group', 'test_config', 'cli'],
                        type=str.lower)

    parser.add_argument('--pipeline_file',
                        help='Path for the pipeline configuration file to '
                             'use. Use the format s3://bucket/path/to/'
                             'pipeline_file to read data directly from an '
                             'S3 bucket. This may require AWS S3 credentials '
                             'specified via the --aws_input_creds option.',
                        default=DEFAULT_PIPELINE)
    parser.add_argument('--group_file',
                        help='Path for the group analysis configuration file '
                             'to use. Use the format s3://bucket/path/to/'
                             'pipeline_file to read data directly from an S3 '
                             'bucket. This may require AWS S3 credentials '
                             'specified via the --aws_input_creds option. '
                             'The output directory needs to refer to the '
                             'output of a preprocessing individual pipeline.',
                        default=None)
    parser.add_argument('--data_config_file',
                        help='Yaml file containing the location of the data '
                             'that is to be processed. This file is not '
                             'necessary if the data in bids_dir is organized '
                             'according to the BIDS format. This enables '
                             'support for legacy data organization and cloud '
                             'based storage. A bids_dir must still be '
                             'specified when using this option, but its '
                             'value will be ignored. Use the format s3://'
                             'bucket/path/to/data_config_file to read data '
                             'directly from an S3 bucket. This may require '
                             'AWS S3 credentials specified via the '
                             '--aws_input_creds option.',
                        default=None)

    parser.add_argument('--preconfig',
                        help='Name of the preconfigured pipeline to run. '
                             'Available preconfigured pipelines: ' +
                             str(AVAILABLE_PIPELINE_CONFIGS) + '. See '
                             f'{DOCS_URL_PREFIX}/user/pipelines/preconfig '
                             'for more information about the preconfigured '
                             'pipelines.',
                        default=None)

    if '--pipeline_override' in sys.argv:  # secret option
        parser.add_argument('--pipeline_override', type=parse_yaml,
                            action='append',
                            help='Override specific options from the '
                                 'pipeline configuration. E.g.: '
                                 '"{\'pipeline_setup\': {\'system_config\': '
                                 '{\'maximum_memory_per_participant\': 1}}}"')

    parser.add_argument('--aws_input_creds',
                        help='Credentials for reading from S3. If not '
                             'provided and s3 paths are specified in the '
                             'data config we will try to access the bucket '
                             'anonymously use the string "env" to indicate '
                             'that input credentials should read from the '
                             'environment. (E.g. when using AWS iam roles).',
                        default=None)
    parser.add_argument('--aws_output_creds',
                        help='Credentials for writing to S3. If not provided '
                             'and s3 paths are specified in the output '
                             'directory we will try to access the bucket '
                             'anonymously use the string "env" to indicate '
                             'that output credentials should read from the '
                             'environment. (E.g. when using AWS iam roles).',
                        default=None)
    # TODO: restore <default=3> for <--n_cpus> once we remove
    #       <max_cores_per_participant> from config file
    #       <https://github.com/FCP-INDI/C-PAC/pull/1264#issuecomment-631643708>
    parser.add_argument('--n_cpus', type=int, default=0,
                        help='Number of execution resources per participant '
                             'available for the pipeline. This flag takes '
                             'precidence over max_cores_per_participant in '
                             'the pipeline configuration file.')
    parser.add_argument('--mem_mb', type=float,
                        help='Amount of RAM available per participant in '
                             'megabytes. Included for compatibility with '
                             'BIDS-Apps standard, but mem_gb is preferred. '
                             'This flag takes precedence over '
                             'maximum_memory_per_participant in the pipeline '
                             'configuration file.')
    parser.add_argument('--mem_gb', type=float,
                        help='Amount of RAM available per participant in '
                             'gigabytes. If this is specified along with '
                             'mem_mb, this flag will take precedence. This '
                             'flag also takes precedence over '
                             'maximum_memory_per_participant in the pipeline '
                             'configuration file.')
    parser.add_argument('--runtime_usage', type=str,
                        help='Path to a callback.log from a prior run of the '
                             'same pipeline configuration (including any '
                             'resource-management parameters that will be '
                             "applied in this run, like 'n_cpus' and "
                             "'num_ants_threads'). This log will be used to "
                             'override per-node memory estimates with '
                             'observed values plus a buffer.')
    parser.add_argument('--runtime_buffer', type=float,
                        help='Buffer to add to per-node memory estimates if '
                             '--runtime_usage is specified. This number is a '
                             'percentage of the observed memory usage.')
    parser.add_argument('--num_ants_threads', type=int, default=0,
                        help='The number of cores to allocate to ANTS-'
                             'based anatomical registration per '
                             'participant. Multiple cores can greatly '
                             'speed up this preprocessing step. This '
                             'number cannot be greater than the number of '
                             'cores per participant.')
    parser.add_argument('--random_seed', type=str,
                        help='Random seed used to fix the state of execution. '
                             'If unset, each process uses its own default. If '
                             'set, a `random.log` file will be generated '
                             'logging the random state used by each process. '
                             'If set to a positive integer (up to 2147483647'
                             '), that integer will be used to seed each '
                             'process. If set to \'random\', a random seed '
                             'will be generated and recorded for each '
                             'process.')
    parser.add_argument('--save_working_dir', nargs='?',
                        help='Save the contents of the working directory.',
                        default=False)
    parser.add_argument('--disable_file_logging', action='store_true',
                        help='Disable file logging, this is useful for '
                             'clusters that have disabled file locking.',
                        default=False)

    parser.add_argument('--participant_label',
                        help='The label of the participant that should be '
                             'analyzed. The label corresponds to '
                             'sub-<participant_label> from the BIDS spec '
                             '(so it does not include "sub-"). If this '
                             'parameter is not provided all participants '
                             'should be analyzed. Multiple participants '
                             'can be specified with a space separated '
                             'list.',
                        nargs="+")
    parser.add_argument('--participant_ndx',
                        help='The index of the participant that should be '
                             'analyzed. This corresponds to the index of '
                             'the participant in the data config file. '
                             'This was added to make it easier to '
                             'accommodate SGE array jobs. Only a single '
                             'participant will be analyzed. Can be used '
                             'with participant label, in which case it is '
                             'the index into the list that follows the '
                             'participant_label flag. Use the value "-1" '
                             'to indicate that the participant index '
                             'should be read from the '
                             'AWS_BATCH_JOB_ARRAY_INDEX environment '
                             'variable.',
                        default=None, type=int)

    parser.add_argument('--T1w_label',
                        help='C-PAC only runs one T1w per participant-'
                             'session at a time, at this time. Use this '
                             'flag to specify any BIDS entity (e.g., "acq-'
                             'VNavNorm") or sequence of BIDS entities ('
                             'e.g., "acq-VNavNorm_run-1") to specify '
                             'which of multiple T1w files to use. Specify '
                             '"--T1w_label T1w" to choose the T1w file '
                             'with the fewest BIDS entities (i.e., the '
                             'final option of [*_acq-VNavNorm_T1w.nii.gz, '
                             '*_acq-HCP_T1w.nii.gz, *_T1w.nii.gz"]). '
                             'C-PAC will choose the first T1w it finds if '
                             'the user does not provide this flag, or '
                             'if multiple T1w files match the --T1w_label '
                             'provided.\nIf multiple T2w files are present '
                             'and a comparable filter is possible, T2w '
                             'files will be filtered as well. If no T2w files '
                             'match this --T1w_label, T2w files will be '
                             'processed as if no --T1w_label were provided.')
    parser.add_argument('--bold_label',
                        help='To include a specified subset of available '
                             'BOLD files, use this flag to specify any '
                             'BIDS entity (e.g., "task-rest") or sequence '
                             'of BIDS entities (e.g. "task-rest_run-1"). '
                             'To specify the bold file with the fewest '
                             'BIDS entities in the file name, specify '
                             '"--bold_label bold". Multiple `--bold_'
                             'label`s can be specified with a space-'
                             'separated list. If multiple `--bold_label`s '
                             'are provided (e.g., "--bold_label task-rest_'
                             'run-1 task-rest_run-2", each scan that '
                             'includes all BIDS entities specified in any '
                             'of the provided `--bold_label`s will be '
                             'analyzed. If this parameter is not provided '
                             'all BOLD scans should be analyzed.',
                        nargs="+")

    parser.add_argument('-v', '--version', action='version',
                        version=f'C-PAC BIDS-App version {__version__}')
    parser.add_argument('--bids_validator_config',
                        help='JSON file specifying configuration of '
                             'bids-validator: See https://github.com/bids-'
                             'standard/bids-validator for more info.')
    parser.add_argument('--skip_bids_validator',
                        help='Skips bids validation.',
                        action='store_true')

    parser.add_argument('--anat_only',
                        help='run only the anatomical preprocessing',
                        action='store_true')

    parser.add_argument('--tracking_opt-out', action='store_true',
                        help='Disable usage tracking. Only the number of '
                             'participants on the analysis is tracked.',
                        default=False)

    parser.add_argument('--monitoring',
                        help='Enable monitoring server on port 8080. You '
                             'need to bind the port using the Docker '
                             'flag "-p".',
                        action='store_true')

    # get the command line arguments
    args = parser.parse_args(
        sys.argv[
            1:(
                sys.argv.index('--')
                if '--' in sys.argv
                else len(sys.argv)
            )
        ]
    )

    bids_dir_is_s3 = args.bids_dir.lower().startswith("s3://")
    bids_dir = args.bids_dir if bids_dir_is_s3 else os.path.realpath(
        args.bids_dir)
    output_dir_is_s3 = args.output_dir.lower().startswith("s3://")
    output_dir = args.output_dir if output_dir_is_s3 else os.path.realpath(
        args.output_dir)

    if args.analysis_level == "cli":
        from CPAC.__main__ import main
        main.main(args=sys.argv[sys.argv.index('--') + 1:])
        sys.exit(0)

    elif args.analysis_level == "group":
        if not args.group_file or not os.path.exists(args.group_file):

            print()
            print("No group analysis configuration file was supplied.")
            print()

            import pkg_resources as p
            args.group_file = \
                p.resource_filename(
                    "CPAC",
                    os.path.join(
                        "resources",
                        "configs",
                        "group_config_template.yml"
                    )
                )

            output_group = os.path.join(output_dir, "group_config.yml")

            try:
                if output_dir.lower().startswith("s3://"):
                    raise Exception

                if not os.path.exists(output_group):
                    shutil.copyfile(args.group_file, output_group)
            except (Exception, IOError):
                print("Could not create group analysis configuration file.")
                print("Please refer to the C-PAC documentation for group "
                      "analysis setup.")
                print()
            else:
                print(
                    "Please refer to the output directory for a template of "
                    "the file and, after customizing to your analysis, add "
                    "the flag"
                    "\n\n"
                    "    --group_file {0}"
                    "\n\n"
                    "to your `docker run` command"
                    "\n"
                    .format(output_group)
                )

            sys.exit(1)

        else:
            import CPAC.pipeline.cpac_group_runner as cgr
            print("Starting group level analysis of data in {0} using "
                  "{1}".format(bids_dir, args.group_file))
            cgr.run(args.group_file)

            sys.exit(0)

    elif args.analysis_level in ["test_config", "participant"]:

        # check to make sure that the input directory exists
        if (
            not args.data_config_file and
            not bids_dir_is_s3 and
            not os.path.exists(bids_dir)
        ):

            print(f"Error! Could not find {bids_dir}")
            sys.exit(1)

        # check to make sure that the output directory exists
        if not output_dir_is_s3 and not os.path.exists(output_dir):

            try:
                os.makedirs(output_dir)
            except Exception:
                print(f"Error! Could not find/create output dir {output_dir}")
                sys.exit(1)

        # validate input dir (if skip_bids_validator is not set)
        if not args.data_config_file:
            print()
            if args.bids_validator_config:
                print("Running BIDS validator")
                run("bids-validator --config {config} {bids_dir}".format(
                    config=args.bids_validator_config,
                    bids_dir=bids_dir
                ))
            elif args.skip_bids_validator:
                print('Skipping bids-validator...')
            elif bids_dir_is_s3:
                print('Skipping bids-validator for S3 datasets...')
            else:
                print("Running BIDS validator")
                run(f"bids-validator {bids_dir}")

        if args.preconfig:
            args.pipeline_file = load_preconfig(args.preconfig)

        # otherwise, if we are running group, participant, or dry run we
        # begin by conforming the configuration
        c = load_yaml_config(args.pipeline_file, args.aws_input_creds)

        if 'pipeline_setup' not in c:
            _url = (f'{DOCS_URL_PREFIX}/user/pipelines/'
                    '1.7-1.8-nesting-mappings')

            warn('\nC-PAC changed its pipeline configuration format in '
                 f'v1.8.0.\nSee {_url} for details.\n',
                 category=DeprecationWarning)

            updated_config = os.path.join(
                output_dir,
                'updated_config',
                os.path.basename(args.pipeline_file)
            )
            os.makedirs(
                os.path.join(output_dir, 'updated_config'), exist_ok=True)

            open(updated_config, 'w').write(yaml.dump(c))

            upgrade_pipeline_to_1_8(updated_config)
            c = load_yaml_config(updated_config, args.aws_input_creds)

        overrides = {}
        if hasattr(args, 'pipeline_override') and args.pipeline_override:
            overrides = {
                k: v for d in args.pipeline_override for k, v in d.items()}
            c = update_nested_dict(c, overrides)

        if args.anat_only:
            c = update_nested_dict(c, {'FROM': 'anat-only'})

        c = Configuration(c)

        # get the aws_input_credentials, if any are specified
        if args.aws_input_creds:
            c['awsCredentialsFile'] = resolve_aws_credential(
                args.aws_input_creds)

        if args.aws_output_creds:
            c['pipeline_setup']['Amazon-AWS'][
                'aws_output_bucket_credentials'
            ] = resolve_aws_credential(
                args.aws_output_creds
            )

        c['pipeline_setup']['output_directory']['path'] = os.path.join(
            output_dir, "output")

        if not output_dir_is_s3:
            c['pipeline_setup']['log_directory']['path'] = os.path.join(
                output_dir, "log")
        else:
            c['pipeline_setup']['log_directory']['path'] = os.path.join(
                DEFAULT_TMP_DIR, "log")

        if args.mem_gb:
            c['pipeline_setup']['system_config'][
                'maximum_memory_per_participant'] = float(args.mem_gb)
        elif args.mem_mb:
            c['pipeline_setup']['system_config'][
                'maximum_memory_per_participant'] = float(args.mem_mb) / 1024.0
        else:
            try:
                c['pipeline_setup', 'system_config',
                  'maximum_memory_per_participant'] = float(
                    c['pipeline_setup', 'system_config',
                      'maximum_memory_per_participant'])
            except KeyError:
                c['pipeline_setup', 'system_config',
                  'maximum_memory_per_participant'] = 6.0

        # Preference: n_cpus if given, override if present, else from config if
        # present, else n_cpus=3
        if int(args.n_cpus) == 0:
            try:
                args.n_cpus = c['pipeline_setup', 'system_config',
                                'max_cores_per_participant']
            except KeyError:
                args.n_cpus = 3
        c['pipeline_setup', 'system_config',
          'max_cores_per_participant'] = int(args.n_cpus)

        c['pipeline_setup']['system_config']['num_participants_at_once'] = int(
            c['pipeline_setup']['system_config'].get(
                'num_participants_at_once', 1))
        # Reduce cores per participant if cores times participants is more than
        # available CPUS. n_cpus is a hard upper limit.
        if (
            c['pipeline_setup']['system_config']['max_cores_per_participant'] *
            c['pipeline_setup']['system_config']['num_participants_at_once']
        ) > int(args.n_cpus):
            c['pipeline_setup']['system_config'][
                'max_cores_per_participant'
            ] = int(args.n_cpus) // c['pipeline_setup']['system_config'][
                'num_participants_at_once'
            ]
            if c['pipeline_setup']['system_config'][
                'max_cores_per_participant'
            ] == 0:
                c['pipeline_setup']['system_config'][
                    'max_cores_per_participant'] = args.n_cpus
                c['pipeline_setup']['system_config'][
                    'num_participants_at_once'] = 1

        if int(args.num_ants_threads) == 0:
            try:
                args.num_ants_threads = c['pipeline_setup', 'system_config',
                                          'num_ants_threads']
            except KeyError:
                args.num_ants_threads = 3
        c['pipeline_setup', 'system_config', 'num_ants_threads'] = int(
            args.num_ants_threads)

        c['pipeline_setup']['system_config']['num_ants_threads'] = min(
            c['pipeline_setup']['system_config']['max_cores_per_participant'],
            int(c['pipeline_setup']['system_config']['num_ants_threads'])
        )

        if args.random_seed:
            c['pipeline_setup']['system_config']['random_seed'] = \
                args.random_seed

        if c['pipeline_setup']['system_config']['random_seed'] is not None:
            c['pipeline_setup']['system_config']['random_seed'] =  \
                set_up_random_state(c['pipeline_setup']['system_config'][
                    'random_seed'])

        if args.runtime_usage is not None:
            c['pipeline_setup']['system_config']['observed_usage'][
                'callback_log'] = args.runtime_usage
        if args.runtime_buffer is not None:
            c['pipeline_setup']['system_config']['observed_usage'][
                'buffer'] = args.runtime_buffer

        c['disable_log'] = args.disable_file_logging

        if args.save_working_dir is not False:
            c['pipeline_setup']['working_directory'][
                'remove_working_dir'] = False
        if isinstance(args.save_working_dir, str):
            c['pipeline_setup']['working_directory']['path'] = \
                os.path.abspath(args.save_working_dir)
        elif not output_dir_is_s3:
            c['pipeline_setup']['working_directory']['path'] = \
                os.path.join(output_dir, "working")
        else:
            warn('Cannot write working directory to S3 bucket. '
                 'Either change the output directory to something '
                 'local or turn off the --save_working_dir flag',
                 category=UserWarning)

        if c['pipeline_setup']['output_directory']['quality_control'][
                'generate_xcpqc_files']:
            c['functional_preproc']['motion_estimates_and_correction'][
                'motion_estimates']['calculate_motion_first'] = True
            c['functional_preproc']['motion_estimates_and_correction'][
                'motion_estimates']['calculate_motion_after'] = True

        if args.participant_label:
            print(
                "#### Running C-PAC for {0}"
                .format(", ".join(args.participant_label))
            )
        else:
            print("#### Running C-PAC")

        print("Number of participants to run in parallel: {0}"
              .format(c['pipeline_setup']['system_config'][
                  'num_participants_at_once']))

        if not args.data_config_file:
            print("Input directory: {0}".format(bids_dir))

        print("Output directory: {0}".format(
            c['pipeline_setup']['output_directory']['path']))
        print("Working directory: {0}".format(
            c['pipeline_setup']['working_directory']['path']))
        print("Log directory: {0}".format(
            c['pipeline_setup']['log_directory']['path']))
        print("Remove working directory: {0}".format(
            c['pipeline_setup']['working_directory']['remove_working_dir']))
        print("Available memory: {0} (GB)".format(
            c['pipeline_setup']['system_config'][
                'maximum_memory_per_participant']))
        print("Available threads: {0}".format(
            c['pipeline_setup']['system_config']['max_cores_per_participant']))
        print("Number of threads for ANTs: {0}".format(
            c['pipeline_setup']['system_config']['num_ants_threads']))

        # create a timestamp for writing config files
        st = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%SZ')

        # update config file
        if not output_dir_is_s3:
            pipeline_config_file = os.path.join(
                output_dir, "cpac_pipeline_config_{0}.yml".format(st)
            )
        else:
            pipeline_config_file = os.path.join(
                DEFAULT_TMP_DIR, "cpac_pipeline_config_{0}.yml".format(st)
            )

        open(pipeline_config_file, 'w').write(
            create_yaml_from_template(c, DEFAULT_PIPELINE, True))
        open(f'{pipeline_config_file[:-4]}_min.yml', 'w').write(
            create_yaml_from_template(c, DEFAULT_PIPELINE, False))

        if args.participant_label:
            args.participant_label = cl_strip_brackets(args.participant_label)
            args.participant_label = [
                'sub-' + pt if not pt.startswith('sub-') else pt
                for pt in args.participant_label
            ]

        # otherwise we move on to conforming the data configuration
        if not args.data_config_file:
            sub_list = create_cpac_data_config(bids_dir,
                                               args.participant_label,
                                               args.aws_input_creds,
                                               args.skip_bids_validator,
                                               only_one_anat=False)
        else:
            sub_list = load_cpac_data_config(args.data_config_file,
                                             args.participant_label,
                                             args.aws_input_creds)
        prefilter = list(sub_list)
        sub_list = sub_list_filter_by_labels(sub_list,
                                             {'T1w': args.T1w_label,
                                              'bold': args.bold_label})

        # C-PAC only handles single anatomical images (for now)
        # so we take just the first as a string if we have a list
        for i, sub in enumerate(sub_list):
            if isinstance(sub.get('anat'), dict):
                for anat_key in sub['anat']:
                    if(
                        isinstance(sub['anat'][anat_key], list) and
                        len(sub['anat'][anat_key])
                    ):
                        sub_list[i]['anat'][
                            anat_key] = sub['anat'][anat_key][0]
            if isinstance(sub.get('anat'), list) and len(sub['anat']):
                sub_list[i]['anat'] = sub['anat'][0]

        if args.participant_ndx is not None:

            participant_ndx = int(args.participant_ndx)
            if participant_ndx == -1:
                args.participant_ndx = os.environ['AWS_BATCH_JOB_ARRAY_INDEX']

            if 0 <= participant_ndx < len(sub_list):
                print('Processing data for participant {0} ({1})'.format(
                    args.participant_ndx,
                    sub_list[participant_ndx]["subject_id"]
                ))
                sub_list = [sub_list[participant_ndx]]
                data_config_file = "cpac_data_config_idx-%s_%s.yml" % (
                    args.participant_ndx, st)
            else:
                print("Participant ndx {0} is out of bounds [0, {1})".format(
                    participant_ndx,
                    str(len(sub_list))
                ))
                sys.exit(1)
        else:
            # write out the data configuration file
            data_config_file = "cpac_data_config_{0}.yml".format(st)

        if not output_dir_is_s3:
            data_config_file = os.path.join(output_dir, data_config_file)
        else:
            data_config_file = os.path.join(DEFAULT_TMP_DIR, data_config_file)

        with open(data_config_file, 'w') as f:
            noalias_dumper = yaml.dumper.SafeDumper
            noalias_dumper.ignore_aliases = lambda self, data: True
            yaml.dump(sub_list, f, default_flow_style=False,
                      Dumper=noalias_dumper)

        if args.analysis_level in ["participant", "test_config"]:
            # build pipeline easy way
            from CPAC.utils.monitoring import monitor_server
            import CPAC.pipeline.cpac_runner

            monitoring = None
            if args.monitoring:
                try:
                    monitoring = monitor_server(
                        c['pipeline_setup']['pipeline_name'],
                        c['pipeline_setup']['log_directory']['path']
                    )
                except:
                    pass

            plugin_args = {
                'n_procs': int(c['pipeline_setup']['system_config'][
                    'max_cores_per_participant']),
                'memory_gb': int(c['pipeline_setup']['system_config'][
                    'maximum_memory_per_participant']),
                'raise_insufficient': c['pipeline_setup']['system_config'][
                                        'raise_insufficient'],
                'status_callback': log_nodes_cb
            }
            if c['pipeline_setup']['system_config']['observed_usage'][
                    'callback_log'] is not None:
                plugin_args['runtime'] = {
                    'usage': c['pipeline_setup']['system_config'][
                        'observed_usage']['callback_log'],
                    'buffer': c['pipeline_setup']['system_config'][
                        'observed_usage']['buffer']}

            print("Starting participant level processing")
            CPAC.pipeline.cpac_runner.run(
                data_config_file,
                pipeline_config_file,
                plugin='MultiProc' if plugin_args[
                    'n_procs'
                ] > 1 else 'Linear',
                plugin_args=plugin_args,
                tracking=not args.tracking_opt_out,
                test_config=(1 if args.analysis_level == "test_config" else 0)
            )

            if monitoring:
                monitoring.join(10)

            if args.analysis_level == "test_config":
                print(
                    '\nPipeline and data configuration files should'
                    ' have been written to {0} and {1} respectively.'.format(
                        pipeline_config_file,
                        data_config_file
                    )
                )

    sys.exit(0)


if __name__ == '__main__':
    run_main()
