#!/usr/bin/python

import sys
from os import path

def check_inputs(*pathstrs):
    for pathstr in pathstrs:
        if not path.exists(pathstr):
            print "ERROR: input '%s' doesn't exist" % pathstr
            raise SystemExit(2)
        _,ext = path.splitext(pathstr)
        if ext != ".yml":
            print "ERROR: input '%s' is not a YAML file (*.yml)" % pathstr
            raise SystemExit(2)
    return

if len(sys.argv) != 3 and len(sys.argv) != 5:
    print "Usage: %s /path/to/pipeline_config.yml /path/to/CPAC_subject_list.yml" % path.basename(sys.argv[0])
    print "Alternate usage: %s /path/to/pipeline_config.yml /path/to/CPAC_subject_list.yml nipype=/path/to/custom/Nipype cpac=/path/to/custom/cpac" % path.basename(sys.argv[0])
    print "Will run C-PAC"
    raise SystemExit(1)


def check_custom_path(custom_path):

    if custom_path != None:

        if "nipype=" in custom_path:
            nipype_path = custom_path.replace("nipype=","")
            sys.path.insert(0,nipype_path)
        elif "cpac=" in custom_path:
            cpac_path = custom_path.replace("cpac=","")
            sys.path.insert(0,cpac_path)
        else:
            print "Improper format for custom paths.\n"
            print "Usage: %s /path/to/pipeline_config.yml /path/to/CPAC_subject_list.yml" % path.basename(sys.argv[0])
            print "Alternate usage: %s /path/to/pipeline_config.yml /path/to/CPAC_subject_list.yml nipype=/path/to/custom/Nipype cpac=/path/to/custom/cpac" % path.basename(sys.argv[0])
            print "Will run C-PAC"
            raise SystemExit(1)


config_file         = sys.argv[1]
subject_list_file   = sys.argv[2]

if len(sys.argv)==5:
    check_custom_path(sys.argv[3])
    check_custom_path(sys.argv[4])

check_inputs(config_file, subject_list_file)

import CPAC
CPAC.pipeline.cpac_runner.run(config_file, subject_list_file)
