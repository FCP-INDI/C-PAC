#!/usr/bin/env python
import argparse
import os
import subprocess
import yaml
import sys

import datetime
import time

__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

def load_yaml_config(config_filename, aws_input_creds):

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

        bucket.download_file(prefix, '/tmp/'+os.path.basename(config_filename))

        config_filename = '/tmp/'+os.path.basename(config_filename)

    config_filename = os.path.realpath(config_filename)
    if os.path.isfile(config_filename):
        with open(config_filename,'r') as infd:
            config_data = yaml.load(infd)

    return(config_data)



def run(command, env={}):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               shell=True, env=env)
    while True:
        line = process.stdout.readline()
        line = str(line)[:-1]
        print(line)
        if line == '' and process.poll() is not None:
            break


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
                    choices=['participant', 'group', 'test_config', 'GUI'])
parser.add_argument('--pipeline_file', help='Name for the pipeline '
                                            ' configuration file to use. '
                                            'Use the format'
                                            ' s3://bucket/path/to/pipeline_file to read data directly from an S3 bucket.'
                                            ' This may require AWS S3 credentials specificied via the'
                                            ' --aws_input_creds option.',
                    default="/cpac_resources/default_pipeline.yaml")
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

# get the command line arguments
args = parser.parse_args()

print(args)

# if we are running the GUI, then get to it
if args.analysis_level == "GUI":
    print "Starting CPAC GUI"
    import CPAC

    CPAC.GUI.run()
    sys.exit(1)

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
    print('skipping bids-validator...')
else:
    print("\nRunning BIDS validator")
    run("bids-validator {bids_dir}".format(bids_dir=args.bids_dir))

# get the aws_input_credentials, if any are specified
if args.aws_input_creds:
    if args.aws_input_creds is "env":
        import urllib2
        aws_creds_address = "169.254.170.2"+os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"]
        aws_creds = urllib2.urlopen(aws_creds_address).read()

        args.aws_input_creds = "/tmp/aws_input_creds.csv"

        with open(args.aws_input_creds) as ofd:
            for key, vname in [("AccessKeyId","AWSAcessKeyId"), ("SecretAccessKey","AWSSecretKey")]:
                ofd.write("{0}={1}".format(vname,aws_creds[key])) 

    if os.path.isfile(args.aws_input_creds):
        c['awsCredentialsFile'] = args.aws_input_creds
    else:
        raise IOError("Could not find aws credentials {0}".format(args.aws_input_creds))

# otherwise, if we are running group, participant, or dry run we
# begin by conforming the configuration
c = load_yaml_config(args.pipeline_file, args.aws_input_creds)

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

    if args.aws_output_creds is "env":
        import urllib2
        aws_creds_address = "169.254.170.2"+os.environ["AWS_CONTAINER_CREDENTIALS_RELATIVE_URI"]
        aws_creds = urllib2.urlopen(aws_creds_address).read()

        args.aws_output_creds = "/tmp/aws_output_creds.csv"

        with open(args.aws_output_creds) as ofd:
            for key, vname in [("AccessKeyId","AWSAcessKeyId"), ("SecretAccessKey","AWSSecretKey")]:
                ofd.write("{0}={1}".format(vname,aws_creds[key])) 

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
        print ('Cannot write working directory to S3 bucket.'
               ' Either change the output directory to something'
               ' local or turn off the --removeWorkingDir flag')
else:
    c['removeWorkingDir'] = True
    c['workingDirectory'] = os.path.join('/scratch', "working")

if args.participant_label:
    print ("#### Running C-PAC on {0}".format(args.participant_label))
else:
    print ("#### Running C-PAC")

print ("Number of participants to run in parallel: {0}".format(c['numParticipantsAtOnce']))
print ("Input directory: {0}".format(args.bids_dir))
print ("Output directory: {0}".format(c['outputDirectory']))
print ("Working directory: {0}".format(c['workingDirectory']))
print ("Crash directory: {0}".format(c['crashLogDirectory']))
print ("Log directory: {0}".format(c['logDirectory']))
print ("Remove working directory: {0}".format(c['removeWorkingDir']))
print ("Available memory: {0} (GB)".format(c['maximumMemoryPerParticipant']))
print ("Available threads: {0}".format(c['maxCoresPerParticipant']))
print ("Number of threads for ANTs: {0}".format(c['num_ants_threads']))

# create a timestamp for writing config files
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d%H%M%S')

# update config file
if "s3://" not in args.output_dir.lower():
    config_file = os.path.join(args.output_dir, "cpac_pipeline_config_{0}.yml".format(st))
else:
    config_file = os.path.join("/scratch", "cpac_pipeline_config_{0}.yml".format(st))

with open(config_file, 'w') as f:
    yaml.dump(c, f)

# we have all we need if we are doing a group level analysis
if args.analysis_level == "group":
    # print ("Starting group level analysis of data in %s using %s"%(args.bids_dir, config_file))
    # import CPAC
    # CPAC.pipeline.cpac_group_runner.run(config_file, args.bids_dir)
    # sys.exit(1)
    print ("Starting group level analysis of data in {0} using {1}".format(args.bids_dir, config_file))
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
        print ("Did not find any files to process")
        sys.exit(1)

    # TODO: once CPAC is updated to use per-scan parameters from subject list,
    # change the 3rd arguement to the config dict returned from
    # collect_bids_files_configs
    sub_list = bids_gen_cpac_sublist(args.bids_dir, file_paths, [], args.aws_input_creds)

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
            print ("Did not find data for {0} in {1}".format(", ".join(args.participant_label),
                                                             args.data_config_file))
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
        print ("Participant ndx {0} is out of bounds [0,{1})".format(int(args.participant_ndx),
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
    yaml.dump(sub_list, f)

if args.analysis_level == "participant":
    # build pipeline easy way
    import CPAC
    from nipype.pipeline.plugins.callback_log import log_nodes_cb

    plugin_args = {'n_procs': int(c['maxCoresPerParticipant']),
                   'memory_gb': int(c['maximumMemoryPerParticipant']),
                   'callback_log': log_nodes_cb}

    print ("Starting participant level processing")
    CPAC.pipeline.cpac_runner.run(config_file, data_config_file,
                                  plugin='MultiProc', plugin_args=plugin_args)
else:
    print ('This has been a test run, the pipeline and data configuration files should'
           ' have been written to {0} and {1} respectively.'
           ' CPAC will not be run.'.format(config_file, data_config_file))

sys.exit(0)
