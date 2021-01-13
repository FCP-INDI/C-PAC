#!/usr/bin/env python

import argparse
import datetime
import os
import subprocess
import sys
import time
import shutil
import yaml
from base64 import b64decode
from urllib import request
from urllib.error import HTTPError

from CPAC import __version__
from CPAC.utils.configuration import Configuration
from CPAC.utils.yaml_template import create_yaml_from_template, \
                                     upgrade_pipeline_to_1_8
from CPAC.utils.utils import load_preconfig

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

parser = argparse.ArgumentParser(description='C-PAC Pipeline Runner')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                                     'formatted according to the BIDS standard. '
                                     'Use the format'
                                     ' s3://bucket/path/to/bidsdir to read data directly from an S3 bucket.'
                                     ' This may require AWS S3 credentials specified via the'
                                     ' --aws_input_creds option.')
parser.add_argument('output_dir', help='The directory where the output files '
                                       'should be stored. If you are running group level analysis '
                                       'this folder should be prepopulated with the results of the '
                                       'participant level analysis. Use the format '
                                       ' s3://bucket/path/to/bidsdir to write data directly to an S3 bucket.'
                                       ' This may require AWS S3 credentials specified via the'
                                       ' --aws_output_creds option.')
parser.add_argument('analysis_level', help='Level of the analysis that will '
                                           ' be performed. Multiple participant level analyses can be run '
                                           ' independently (in parallel) using the same output_dir. '
                                           ' GUI will open the CPAC gui (currently only works with singularity) and'
                                           ' test_config will run through the entire configuration process but will'
                                           ' not execute the pipeline.',
                    choices=['participant', 'group', 'test_config', 'gui', 'cli'], type=str.lower)

parser.add_argument('--pipeline_file', help='Path for the pipeline '
                                            ' configuration file to use. '
                                            'Use the format'
                                            ' s3://bucket/path/to/pipeline_file to read data directly from an S3 bucket.'
                                            ' This may require AWS S3 credentials specified via the'
                                            ' --aws_input_creds option.',
                    default=DEFAULT_PIPELINE)
parser.add_argument('--group_file', help='Path for the group analysis configuration file to use. '
                                         'Use the format s3://bucket/path/to/pipeline_file to read data directly from an S3 bucket. '
                                         'This may require AWS S3 credentials specified via the --aws_input_creds option. '
                                         'The output directory needs to refer to the output of a preprocessing individual pipeline.',
                    default=None)
parser.add_argument('--data_config_file', help='Yaml file containing the location'
                                               ' of the data that is to be processed. Can be generated from the CPAC'
                                               ' gui. This file is not necessary if the data in bids_dir is organized'
                                               ' according to'
                                               ' the BIDS format. This enables support for legacy data organization'
                                               ' and cloud based storage. A bids_dir must still be specified when'
                                               ' using this option, but its value will be ignored.'
                                               ' Use the format'
                                               ' s3://bucket/path/to/data_config_file to read data directly from an S3 bucket.'
                                               ' This may require AWS S3 credentials specified via the'
                                               ' --aws_input_creds option.',
                    default=None)


parser.add_argument('--preconfig', help='Name of the pre-configured pipeline to run.',
                    default=None)

parser.add_argument('--pipeline_override', type=parse_yaml, action='append',
                    help='Override specific options from the pipeline configuration. E.g.: "maximumMemoryPerParticipant: 10"')

parser.add_argument('--aws_input_creds', help='Credentials for reading from S3.'
                                              ' If not provided and s3 paths are specified in the data config'
                                              ' we will try to access the bucket anonymously'
                                              ' use the string "env" to indicate that input credentials should'
                                              ' read from the environment. (E.g. when using AWS iam roles).',
                    default=None)
parser.add_argument('--aws_output_creds', help='Credentials for writing to S3.'
                                               ' If not provided and s3 paths are specified in the output directory'
                                               ' we will try to access the bucket anonymously'
                    ' use the string "env" to indicate that output credentials should'
                    ' read from the environment. (E.g. when using AWS iam roles).',
                    default=None)
# TODO: restore <default=3> for <--n_cpus> once we remove
#       <maxCoresPerParticipant> from config file
#       <https://github.com/FCP-INDI/C-PAC/pull/1264#issuecomment-631643708>
parser.add_argument('--n_cpus', type=int, default=0,
                    help='Number of execution resources per participant '
                         ' available for the pipeline.')
parser.add_argument('--mem_mb', type=float,
                    help='Amount of RAM available to the pipeline in megabytes.'
                         ' Included for compatibility with BIDS-Apps standard, but mem_gb is preferred')
parser.add_argument('--mem_gb', type=float,
                    help='Amount of RAM available to the pipeline in gigabytes.'
                         ' if this is specified along with mem_mb, this flag will take precedence.')

parser.add_argument('--save_working_dir', nargs='?',
                    help='Save the contents of the working directory.', default=False)
parser.add_argument('--disable_file_logging', action='store_true',
                    help='Disable file logging, this is useful for clusters that have disabled file locking.',
                    default=False)
parser.add_argument('--participant_label', help='The label of the participant'
                                                ' that should be analyzed. The label '
                                                'corresponds to sub-<participant_label> from the BIDS spec '
                                                '(so it does not include "sub-"). If this parameter is not '
                                                'provided all participants should be analyzed. Multiple '
                                                'participants can be specified with a space separated list. To work'
                                                ' correctly this should come at the end of the command line.',
                    nargs="+")
parser.add_argument('--participant_ndx', help='The index of the participant'
                                              ' that should be analyzed. This corresponds to the index of the'
                                              ' participant in the data config file. This was added to make it easier'
                                              ' to accommodate SGE array jobs. Only a single participant will be'
                                              ' analyzed. Can be used with participant label, in which case it is the'
                                              ' index into the list that follows the participant_label flag.'
                                              ' Use the value "-1" to indicate that the participant index should'
                                              ' be read from the AWS_BATCH_JOB_ARRAY_INDEX environment variable.',
                    default=None, type=int)

parser.add_argument('-v', '--version', action='version',
                    version='C-PAC BIDS-App version {}'.format(__version__))
parser.add_argument('--bids_validator_config', help='JSON file specifying configuration of '
                    'bids-validator: See https://github.com/bids-standard/bids-validator for more info.')
parser.add_argument('--skip_bids_validator',
                    help='Skips bids validation.',
                    action='store_true')

parser.add_argument('--anat_only', help='run only the anatomical preprocessing',
                    action='store_true')

parser.add_argument('--tracking_opt-out', action='store_true',
                    help='Disable usage tracking. Only the number of participants on the analysis is tracked.',
                    default=False)

parser.add_argument('--monitoring',
                    help='Enable monitoring server on port 8080. You need to bind the port using the Docker flag "-p".',
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

# if we are running the GUI, then get to it
if args.analysis_level == "gui":
    print("Starting CPAC GUI")
    import CPAC.GUI

    CPAC.GUI.run()
    sys.exit(0)

elif args.analysis_level == "cli":
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

        output_group = os.path.join(args.output_dir, "group_config.yml")

        try:
            if args.output_dir.lower().startswith("s3://"):
                raise Exception

            if not os.path.exists(output_group):
                shutil.copyfile(args.group_file, output_group)
        except (Exception, IOError) as e:
            print("Could not create group analysis configuration file.")
            print("Please refer to the C-PAC documentation for group analysis "
                  "setup.")
            print()
        else:
            print(
                "Please refer to the output directory for a template of the "
                "file and, after customizing to your analysis, add the flag"
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
        print("Starting group level analysis of data in {0} using {1}".format(
            args.bids_dir, args.group_file
        ))
        cgr.run(args.group_file)

        sys.exit(0)

elif args.analysis_level in ["test_config", "participant"]:
    # check to make sure that the input directory exists
    if not args.data_config_file and \
        not args.bids_dir.lower().startswith("s3://") and \
        not os.path.exists(args.bids_dir):

        print("Error! Could not find {0}".format(args.bids_dir))
        sys.exit(1)

    # check to make sure that the output directory exists
    if not args.output_dir.lower().startswith("s3://") and \
        not os.path.exists(args.output_dir):

        try:
            os.makedirs(args.output_dir)
        except:
            print("Error! Could not find/create output dir {0}".format(
                args.output_dir
            ))
            sys.exit(1)

    # validate input dir (if skip_bids_validator is not set)
    if not args.data_config_file:
        print()
        if args.bids_validator_config:
            print("Running BIDS validator")
            run("bids-validator --config {config} {bids_dir}".format(
                config=args.bids_validator_config,
                bids_dir=args.bids_dir
            ))
        elif args.skip_bids_validator:
            print('Skipping bids-validator...')
        elif args.bids_dir.lower().startswith("s3://"):
            print('Skipping bids-validator for S3 datasets...')
        else:
            print("Running BIDS validator")
            run("bids-validator {bids_dir}".format(bids_dir=args.bids_dir))

    if args.preconfig:
        args.pipeline_file = load_preconfig(args.preconfig)

    # otherwise, if we are running group, participant, or dry run we
    # begin by conforming the configuration
    c = load_yaml_config(args.pipeline_file, args.aws_input_creds)

    if 'pipeline_setup' not in c:
        url_version = f'v{__version__}'
        _url = (f'https://fcp-indi.github.io/docs/{url_version}/'
            'user/pipelines/1.7-1.8-nesting-mappings')
        try:
            request.urlopen(_url)

        except HTTPError:
            if 'dev' in url_version:
                url_version = 'nightly'
            else:
                url_version = 'latest'

        _url = (f'https://fcp-indi.github.io/docs/{url_version}/'
            'user/pipelines/1.7-1.8-nesting-mappings')

        warn('\nC-PAC changed its pipeline configuration format in v1.8.0.\n'
             f'See {_url} for details.\n',
             category=DeprecationWarning)

        updated_config = os.path.join(
            args.output_dir,
            'updated_config',
            os.path.basename(args.pipeline_file)
        )
        os.makedirs(
            os.path.join(args.output_dir, 'updated_config'), exist_ok=True)

        open(updated_config, 'w').write(yaml.dump(c))

        upgrade_pipeline_to_1_8(updated_config)
        c = load_yaml_config(updated_config, args.aws_input_creds)

    c = Configuration(c).dict()

    overrides = {}
    if args.pipeline_override:
        overrides = {k: v for d in args.pipeline_override for k, v in d.items()}
        c.update(overrides)

    if args.anat_only:
        # TODO
        # c.update({'functional_preproc': 'run': Off})?
        c.update({"runFunctional": [0]})

    # get the aws_input_credentials, if any are specified
    if args.aws_input_creds:
        c['awsCredentialsFile'] = resolve_aws_credential(args.aws_input_creds)

    if args.aws_output_creds:
        c['pipeline_setup']['Amazon-AWS']['aws_output_bucket_credentials'] = resolve_aws_credential(
            args.aws_output_creds
        )

    c['pipeline_setup']['output_directory']['path'] = os.path.join(args.output_dir, "output")

    if "s3://" not in args.output_dir.lower():
        c['pipeline_setup']['crash_log_directory']['path'] = os.path.join(args.output_dir, "crash")
        c['pipeline_setup']['log_directory']['path'] = os.path.join(args.output_dir, "log")
    else:
        c['pipeline_setup']['crash_log_directory']['path'] = os.path.join(DEFAULT_TMP_DIR, "crash")
        c['pipeline_setup']['log_directory']['path'] = os.path.join(DEFAULT_TMP_DIR, "log")

    if args.mem_gb:
        c['pipeline_setup']['system_config']['maximum_memory_per_participant'] = float(args.mem_gb)
    elif args.mem_mb:
        c['pipeline_setup']['system_config']['maximum_memory_per_participant'] = float(args.mem_mb) / 1024.0
    else:
        c['pipeline_setup']['system_config']['maximum_memory_per_participant'] = 6.0

    # Preference: n_cpus if given, override if present, else from config if
    # present, else n_cpus=3
    if args.n_cpus == 0:
        c['pipeline_setup']['system_config']['max_cores_per_participant'] = int(c['pipeline_setup']['system_config'].get('max_cores_per_participant', 3))
        args.n_cpus = 3
    else:
        c['pipeline_setup']['system_config']['max_cores_per_participant'] = args.n_cpus
    c['pipeline_setup']['system_config']['num_participants_at_once'] = int(c['pipeline_setup']['system_config'].get('num_participants_at_once', 1))
    # Reduce cores per participant if cores times particiapants is more than
    # available CPUS. n_cpus is a hard upper limit.
    if (c['pipeline_setup']['system_config']['max_cores_per_participant']  * c['pipeline_setup']['system_config']['num_participants_at_once']) > int(
        args.n_cpus
    ):
        c['pipeline_setup']['system_config']['max_cores_per_participant']  = int(
            args.n_cpus
        ) // c['pipeline_setup']['system_config']['num_participants_at_once']
    c['pipeline_setup']['system_config']['num_ants_threads'] = min(
        c['pipeline_setup']['system_config']['max_cores_per_participant'] , int(c['pipeline_setup']['system_config']['num_ants_threads'])
    )

    c['disable_log'] = args.disable_file_logging

    if args.save_working_dir is not False:
        c['pipeline_setup']['working_directory']['remove_working_dir'] = False
        if args.save_working_dir is not None:
            c['pipeline_setup']['working_directory']['path'] = \
                os.path.abspath(args.save_working_dir)
        elif "s3://" not in args.output_dir.lower():
            c['pipeline_setup']['working_directory']['path'] = \
                os.path.join(args.output_dir, "working")
        else:
            print('Cannot write working directory to S3 bucket.'
                ' Either change the output directory to something'
                ' local or turn off the --save_working_dir flag')


    if args.participant_label:
        print(
            "#### Running C-PAC for {0}"
                .format(", ".join(args.participant_label))
        )
    else:
        print("#### Running C-PAC")

    print("Number of participants to run in parallel: {0}"
          .format(c['pipeline_setup']['system_config']['num_participants_at_once']))

    if not args.data_config_file:
        print("Input directory: {0}".format(args.bids_dir))

    print("Output directory: {0}".format(c['pipeline_setup']['output_directory']['path']))
    print("Working directory: {0}".format(c['pipeline_setup']['working_directory']['path']))
    print("Crash directory: {0}".format(c['pipeline_setup']['crash_log_directory']['path']))
    print("Log directory: {0}".format(c['pipeline_setup']['log_directory']['path']))
    print("Remove working directory: {0}".format(c['pipeline_setup']['working_directory']['remove_working_dir']))
    print("Available memory: {0} (GB)".format(c['pipeline_setup']['system_config']['maximum_memory_per_participant']))
    print("Available threads: {0}".format(c['pipeline_setup']['system_config']['max_cores_per_participant']))
    print("Number of threads for ANTs: {0}".format(c['pipeline_setup']['system_config']['num_ants_threads']))

    # create a timestamp for writing config files
    st = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%SZ')

    # update config file
    if "s3://" not in args.output_dir.lower():
        pipeline_config_file = os.path.join(
            args.output_dir, "cpac_pipeline_config_{0}.yml".format(st)
        )
    else:
        pipeline_config_file = os.path.join(
            DEFAULT_TMP_DIR, "cpac_pipeline_config_{0}.yml".format(st)
        )

    with open(pipeline_config_file, 'w') as f:
        f.write(create_yaml_from_template(c, DEFAULT_PIPELINE))

    participant_labels = []
    if args.participant_label:
        participant_labels = [
            'sub-' + pt if not pt.startswith('sub-') else pt
            for pt in args.participant_label
        ]

    # otherwise we move on to conforming the data configuration
    if not args.data_config_file:

        from bids_utils import collect_bids_files_configs, bids_gen_cpac_sublist

        print("Parsing {0}..".format(args.bids_dir))

        (file_paths, config) = collect_bids_files_configs(
            args.bids_dir, args.aws_input_creds)

        if args.participant_label:
            file_paths = [
                file_path
                for file_path in file_paths
                if any(
                    participant_label in file_path
                    for participant_label in participant_labels
                )
            ]

        if not file_paths:
            print("Did not find data for {0}".format(
                ", ".join(args.participant_label)
            ))
            sys.exit(1)

        raise_error = not args.skip_bids_validator

        sub_list = bids_gen_cpac_sublist(
            args.bids_dir,
            file_paths,
            config,
            args.aws_input_creds,
            raise_error=raise_error
        )

        if not sub_list:
            print("Did not find data in {0}".format(args.bids_dir))
            sys.exit(1)

    else:
        # load the file as a check to make sure it is available and readable
        sub_list = load_yaml_config(args.data_config_file, args.aws_input_creds)

        if args.participant_label:

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
                    ", ".join(args.participant_label),
                    (
                        args.data_config_file
                        if not args.data_config_file.startswith("data:")
                        else "data URI"
                    )
                ))
                sys.exit(1)

    if args.participant_ndx:

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

    if "s3://" not in args.output_dir.lower():
        data_config_file = os.path.join(args.output_dir, data_config_file)
    else:
        data_config_file = os.path.join(DEFAULT_TMP_DIR, data_config_file)

    with open(data_config_file, 'w') as f:
        noalias_dumper = yaml.dumper.SafeDumper
        noalias_dumper.ignore_aliases = lambda self, data: True
        yaml.dump(sub_list, f, default_flow_style=False, Dumper=noalias_dumper)

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
            'n_procs': int(c['pipeline_setup']['system_config']['max_cores_per_participant']),
            'memory_gb': int(c['pipeline_setup']['system_config']['maximum_memory_per_participant']),
        }

        print ("Starting participant level processing")
        CPAC.pipeline.cpac_runner.run(
            data_config_file,
            pipeline_config_file,
            plugin='MultiProc' if plugin_args['n_procs'] > 1 else 'Linear',
            plugin_args=plugin_args,
            tracking=not args.tracking_opt_out,
            test_config = 1 if args.analysis_level == "test_config" else 0
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
