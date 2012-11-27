#!/usr/bin/env python

import sys
from os import path

def check_inputs(*pathstrs):
    for pathstr in pathstrs:
        if not path.exists(pathstr):
            print "ERROR: input '%s' doesn't exist" % pathstr
            raise SystemExit(2)
        _,ext = path.splitext(pathstr):
        if ext != ".py":
            print "ERROR: input '%s' is not a python file (*.py)" % pathstr
            raise SystemExit(2)
    return

if len(sys.argv) != 3:
    print "Usage: %s /path/to/config.py /path/to/CPAC_subject_list.py" % path.basename(sys.argv[0])
    print "Will run C-PAC"
    raise SystemExit(1)

config_file         = sys.argv[1]
subject_list_file   = sys.argv[2]

check_inputs(config_file, subject_list_file)

import CPAC
CPAC.pipeline.cpac_runner.run(config_file, subject_list_file)
