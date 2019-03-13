#!/usr/bin/env python
from __future__ import print_function

import argparse
import datetime
import os
import subprocess
import sys
import time
from base64 import b64decode
import shutil
import yaml

from CPAC.utils.yaml_template import create_yaml_from_template

__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

DEFAULT_PIPELINE = "/cpac_resources/default_pipeline.yaml"
if not os.path.exists(DEFAULT_PIPELINE):
    DEFAULT_PIPELINE = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "default_pipeline.yaml"
    )


def load_yaml_config(config_filename, aws_input_creds):

    if config_filename.lower().startswith('data:'):
        try:
            header, encoded = config_filename.split(",", 1)
            config_content = b64decode(encoded)
            config_data = yaml.load(config_content)
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
        with open(config_filename, 'r') as f:
            config_data = yaml.load(f)
            return config_data
    except IOError:
        print("Error! Could not find config file {0}".format(config_filename))
        raise


def run(command, env={}):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               shell=True, env=env)
    while True:
        line = process.stdout.readline()
        line = str(line)[:-1]
        print(line)
        if line == '' and process.poll() is not None:
            break


def parse_yaml(value):
    try:
        config = yaml.load(value)
        if type(config) != dict:
            raise
        return config
    except:
         raise argparse.ArgumentTypeError("Invalid configuration: '%s'" % value)


parser = argparse.ArgumentParser(description='C-PAC Pipeline Runner')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                                     'formatted according to the BIDS standard. '
                                     'Use the format'
                                     ' s3://bucket/path/to/bidsdir to read data directly from an S3 bucket.'
                                     ' This may require AWS S3 credentials specificied via the'
                                     ' --aws_input_creds option.')
parser.add_argument('output_dir', help='The directory where the output files '
                                       'should be stored. If you are running group level analysis '
                                       'this folder should be prepopulated with the results of the '
                                       'participant level analysis. Us the format '
                                       ' s3://bucket/path/to/bidsdir to write data directly to an S3 bucket.'
                                       ' This may require AWS S3 credentials specificied via the'
                                       ' --aws_output_creds option.')
parser.add_argument('analysis_level', help='Level of the analysis that will '
                                           ' be performed. Multiple participant level analyses can be run '
                                           ' independently (in parallel) using the same output_dir. '
                                           ' GUI will open the CPAC gui (currently only works with singularity) and'
                                           ' test_config will run through the entire configuration process but will'
                                           ' not execute the pipeline.',
                    choices=['participant', 'group', 'test_config', 'gui'], type=str.lower)
parser.add_argument('--pipeline_file', help='Path for the pipeline '
                                            ' configuration file to use. '
                                            'Use the format'
                                            ' s3://bucket/path/to/pipeline_file to read data directly from an S3 bucket.'
                                            ' This may require AWS S3 credentials specificied via the'
                                            ' --aws_input_creds option.',
                    default=DEFAULT_PIPELINE)
parser.add_argument('--group_file', help='Path for the group analysis configuration file to use. '
                                         'Use the format s3://bucket/path/to/pipeline_file to read data directly from an S3 bucket. '
                                         'This may require AWS S3 credentials specificied via the --aws_input_creds option. '
                                         'The output directory needs to refer to the output of a preprocessing individual pipeline.',
                    default=None)

parser.add_argument('--pipeline_override', type=parse_yaml, action='append',
                    help='Override specific options from the pipeline configuration. E.g.: "maximumMemoryPerParticipant: 10"')

parser.add_argument('--data_config_file', help='Yaml file containing the location'
                                               ' of the data that is to be processed. Can be generated from the CPAC'
                                               ' gui. This file is not necessary if the data in bids_dir is organized'
                                               ' according to'
                                               ' the BIDS format. This enables support for legacy data organization'
                                               ' and cloud based storage. A bids_dir must still be specified when'
                                               ' using this option, but its value will be ignored.'
                                               ' Use the format'
                                               ' s3://bucket/path/to/data_config_file to read data directly from an S3 bucket.'
                                               ' This may require AWS S3 credentials specificied via the'
                                               ' --aws_input_creds option.',
                    default=None)
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
parser.add_argument('--n_cpus', help='Number of execution '
                                     ' resources available for the pipeline', default="1")
parser.add_argument('--mem_mb', help='Amount of RAM available to the pipeline in megabytes.'
                                     ' Included for compatibility with BIDS-Apps standard, but mem_gb is preferred')
parser.add_argument('--mem_gb', help='Amount of RAM available to the pipeline in gigabytes.'
                                     ' if this is specified along with mem_mb, this flag will take precedence.')
parser.add_argument('--save_working_dir', action='store_true',
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
                                                ' correctly this should come at the end of the command line',
                    nargs="+")
parser.add_argument('--participant_ndx', help='The index of the participant'
                                              ' that should be analyzed. This corresponds to the index of the'
                                              ' participant in the data config file. This was added to make it easier'
                                              ' to accomodate SGE array jobs. Only a single participant will be'
                                              ' analyzed. Can be used with participant label, in which case it is the'
                                              ' index into the list that follows the particpant_label flag.'
                                              ' Use the value "-1" to indicate that the participant index should'
                                              ' be read from the AWS_BATCH_JOB_ARRAY_INDEX environment variable.',
                    default=None)
parser.add_argument('-v', '--version', action='version',
                    version='C-PAC BIDS-App version {}'.format(__version__))
parser.add_argument('--bids_validator_config', help='JSON file specifying configuration of '
                    'bids-validator: See https://github.com/INCF/bids-validator for more info')
parser.add_argument('--skip_bids_validator',
                    help='skips bids validation',
                    action='store_true')
parser.add_argument('--ndmg_mode', help='produce ndmg connectome graphs and '
                    'write out in the ndmg output format',
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
args = parser.parse_args()

# if we are running the GUI, then get to it
if args.analysis_level == "gui":
    print("Starting CPAC GUI")
    import CPAC

    CPAC.GUI.run()
    sys.exit(1)


if args.analysis_level == "cli":
    from CPAC.__main__ import main
    main()
    sys.exit(0)

# check to make sure that the input directory exists
if not args.bids_dir.lower().startswith("s3://") and not os.path.exists(args.bids_dir):
    print("Error! Could not find {0}".format(args.bids_dir))
    sys.exit(0)

# check to make sure that the output directory exists
if not args.output_dir.lower().startswith("s3://") and not os.path.exists(args.output_dir):
    print("Error! Could not find {0}".format(args.output_dir))
    sys.exit(0)

# validate input dir (if skip_bids_validator is not set)
if args.bids_validator_config:
    print("\nRunning BIDS validator")
    run("bids-validator --config {config} {bids_dir}".format(
        config=args.bids_validator_config,
        bids_dir=args.bids_dir))
elif args.skip_bids_validator:
    print('\nSkipping bids-validator...')
elif args.bids_dir.lower().startswith("s3://"):
    print('\nSkipping bids-validator for S3 datasets...')
else:
    print("\nRunning BIDS validator")
    run("bids-validator {bids_dir}".format(bids_dir=args.bids_dir))

if args.ndmg_mode:
    print('\nRunning ndmg mode')
    import os
    import pkg_resources as p
    args.pipeline_file = \
        p.resource_filename("CPAC",
                            os.path.join("resources",
                                         "configs",
                                         "pipeline_config_ndmg.yml"))

# otherwise, if we are running group, participant, or dry run we
# begin by conforming the configuration
c = load_yaml_config(args.pipeline_file, args.aws_input_creds)
if args.pipeline_override:
    overrides = {k: v for d in args.pipeline_override for k, v in d.items()}
    c.update(overrides)

if args.anat_only:
    c.update({ "runFunctional": [0] })

# get the aws_input_credentials, if any are specified
if args.aws_input_creds:
    if args.aws_input_creds is "env":
        import urllib2
        aws_creds_address = "169.254.170.2" + os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"]
        aws_creds = urllib2.urlopen(aws_creds_address).read()

        args.aws_input_creds = "/tmp/aws_input_creds.csv"

        with open(args.aws_input_creds) as ofd:
            for key, vname in [("AccessKeyId","AWSAcessKeyId"), ("SecretAccessKey","AWSSecretKey")]:
                ofd.write("{0}={1}".format(vname,aws_creds[key])) 

    if os.path.isfile(args.aws_input_creds):
        c['awsCredentialsFile'] = args.aws_input_creds
    else:
        raise IOError("Could not find aws credentials {0}".format(args.aws_input_creds))

# set the parameters using the command line arguements
# TODO: we will need to check that the directories exist, and
# make them if they do not
c['outputDirectory'] = os.path.join(args.output_dir, "output")

if "s3://" not in args.output_dir.lower():
    c['crashLogDirectory'] = os.path.join(args.output_dir, "crash")
    c['logDirectory'] = os.path.join(args.output_dir, "log")
else:
    c['crashLogDirectory'] = os.path.join("/scratch", "crash")
    c['logDirectory'] = os.path.join("/scratch", "log")

if args.mem_gb:
    c['maximumMemoryPerParticipant'] = float(args.mem_gb)
elif args.mem_mb:
    c['maximumMemoryPerParticipant'] = float(args.mem_mb) / 1024.0
else:
    c['maximumMemoryPerParticipant'] = 6.0

c['maxCoresPerParticipant'] = int(args.n_cpus)
c['numParticipantsAtOnce'] = 1
c['num_ants_threads'] = min(int(args.n_cpus), int(c['num_ants_threads']))

if args.aws_output_creds:
    if args.aws_output_creds == "env":
        import urllib2
        aws_creds_address = "169.254.170.2" + os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"]
        aws_creds = urllib2.urlopen(aws_creds_address).read()

        args.aws_output_creds = "/tmp/aws_output_creds.csv"

        with open(args.aws_output_creds) as ofd:
            for key, vname in [("AccessKeyId","AWSAcessKeyId"), ("SecretAccessKey","AWSSecretKey")]:
                ofd.write("{0}={1}".format(vname, aws_creds[key])) 

    if os.path.isfile(args.aws_output_creds):
        c['awsOutputBucketCredentials'] = args.aws_output_creds
    else:
        raise IOError("Could not find aws credentials {0}".format(args.aws_output_creds))

if args.disable_file_logging is True:
    c['disable_log'] = True
else:
    c['disable_log'] = False

if args.save_working_dir is True:
    if "s3://" not in args.output_dir.lower():
        c['removeWorkingDir'] = False
        c['workingDirectory'] = os.path.join(args.output_dir, "working")
    else:
        print('Cannot write working directory to S3 bucket.'
               ' Either change the output directory to something'
               ' local or turn off the --save_working_dir flag')
else:
    c['removeWorkingDir'] = True
    c['workingDirectory'] = os.path.join('/scratch', "working")

if args.participant_label:
    print("#### Running C-PAC on {0}".format(args.participant_label))
else:
    print("#### Running C-PAC")

print("Number of participants to run in parallel: {0}".format(c['numParticipantsAtOnce']))
print("Input directory: {0}".format(args.bids_dir))
print("Output directory: {0}".format(c['outputDirectory']))
print("Working directory: {0}".format(c['workingDirectory']))
print("Crash directory: {0}".format(c['crashLogDirectory']))
print("Log directory: {0}".format(c['logDirectory']))
print("Remove working directory: {0}".format(c['removeWorkingDir']))
print("Available memory: {0} (GB)".format(c['maximumMemoryPerParticipant']))
print("Available threads: {0}".format(c['maxCoresPerParticipant']))
print("Number of threads for ANTs: {0}".format(c['num_ants_threads']))

# create a timestamp for writing config files
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d%H%M%S')

# update config file
if "s3://" not in args.output_dir.lower():
    config_file = os.path.join(args.output_dir, "cpac_pipeline_config_{0}.yml".format(st))
else:
    config_file = os.path.join("/scratch", "cpac_pipeline_config_{0}.yml".format(st))


with open(config_file, 'w') as f:
    f.write(create_yaml_from_template(c, DEFAULT_PIPELINE))


# we have all we need if we are doing a group level analysis
if args.analysis_level == "group":

    if not args.group_file or not os.path.exists(args.group_file):

        print()
        print()
        print("No group analysis configuration file was supplied.")
        print()

        import pkg_resources as p
        args.group_file = \
            p.resource_filename("CPAC",
                                os.path.join("resources",
                                            "configs",
                                            "group_config_template.yml"))

        output_group = os.path.join(args.output_dir, "group_config.yml")

        try:
            if args.output_dir.lower().startswith("s3://"):
                raise Exception

            if not os.path.exists(output_group):
                shutil.copyfile(args.group_file, output_group)
        except Exception, IOError:
            print("Could not create group analysis configuration file.")
            print("Please refer to the C-PAC documentation for group analysis set up.")
            print()
        else:
            print(
                "Please refer to the output directory for a template of the file "
                "and, after customizing to your analysis, add the flag\n\n"
                "    --group_file %s"
                "\n\nto your `docker run` command\n"
                % output_group
            )

        sys.exit(1)

    else:
        import CPAC.pipeline.cpac_group_runner as cgr
        print("Starting group level analysis of data in {0} using {1}".format(args.bids_dir, args.group_file))
        cgr.run(args.group_file)
        
        sys.exit(0)

# otherwise we move on to conforming the data configuration
if not args.data_config_file:

    from bids_utils import collect_bids_files_configs, bids_gen_cpac_sublist

    (file_paths, config) = collect_bids_files_configs(args.bids_dir, args.aws_input_creds)

    if args.participant_label:

        pt_file_paths = []
        for pt in args.participant_label:

            if 'sub-' not in pt:
                pt = 'sub-' + pt

            pt_file_paths += [fp for fp in file_paths if pt in fp]

        file_paths = pt_file_paths

    if not file_paths:
        print("Did not find any files to process")
        sys.exit(1)

    sub_list = bids_gen_cpac_sublist(args.bids_dir, file_paths, config,
                                     args.aws_input_creds)

    if not sub_list:
        print("Did not find data in {0}".format(args.bids_dir))
        sys.exit(1)

else:
    # load the file as a check to make sure it is available and readable
    sub_list = load_yaml_config(args.data_config_file, args.aws_input_creds)

    if args.participant_label:
        t_sub_list = []
        for sub_dict in sub_list:
            if sub_dict["subject_id"] in args.participant_label or \
                            sub_dict["subject_id"].replace("sub-", "") in args.participant_label:
                t_sub_list.append(sub_dict)

        sub_list = t_sub_list

        if not sub_list:
            print("Did not find data for {0} in {1}".format(", ".join(args.participant_label),
                                                             args.data_config_file
                                                             if not args.data_config_file.startswith("data:")
                                                             else "data URI"))
            sys.exit(1)


if args.participant_ndx:

    if int(args.participant_ndx) == -1:
        args.participant_ndx = os.environ['AWS_BATCH_JOB_ARRAY_INDEX']

    if 0 <= int(args.participant_ndx) < len(sub_list):
        # make sure to keep it a list
        print('Processing data for participant {0} ({1})'.format(args.participant_ndx, sub_list[int(args.participant_ndx)]["subject_id"]))
        sub_list = [sub_list[int(args.participant_ndx)]]
        data_config_file = "cpac_data_config_pt%s_%s.yml" % (args.participant_ndx, st)
    else:
        print("Participant ndx {0} is out of bounds [0,{1})".format(int(args.participant_ndx),
                                                                     len(sub_list)))
        sys.exit(1)
else:
    # write out the data configuration file
    data_config_file = "cpac_data_config_{0}.yml".format(st)


if "s3://" not in args.output_dir.lower():
    data_config_file = os.path.join(args.output_dir, data_config_file)
else:
    data_config_file = os.path.join("/scratch", data_config_file)

with open(data_config_file, 'w') as f:

    # Avoid dict/list references
    noalias_dumper = yaml.dumper.SafeDumper
    noalias_dumper.ignore_aliases = lambda self, data: True
    yaml.dump(sub_list, f, default_flow_style=False, Dumper=noalias_dumper)


if args.analysis_level == "participant":
    # build pipeline easy way
    import CPAC
    from CPAC.utils.monitoring import log_nodes_cb, monitor_server

    monitoring = None
    if args.monitoring:
        try:
            monitoring = monitor_server(c['pipelineName'], c['logDirectory'])
        except:
            pass

    plugin_args = {'n_procs': int(c['maxCoresPerParticipant']),
                   'memory_gb': int(c['maximumMemoryPerParticipant']),
                   'status_callback': log_nodes_cb}

    print ("Starting participant level processing")
    CPAC.pipeline.cpac_runner.run(config_file, data_config_file,
                                  plugin='MultiProc', plugin_args=plugin_args,
                                  tracking=not args.tracking_opt_out)

    if monitoring:
        monitoring.join(10)
else:
    print ('This has been a test run, the pipeline and data configuration files should'
           ' have been written to {0} and {1} respectively.'
           ' CPAC will not be run.'.format(config_file, data_config_file))

sys.exit(0)
