#!/usr/bin/env python

# Goal of this script is to update the scan.js file of this subject
# with the supplied python log-file containing workflow info

import os, sys
from os import path as op
from lockfile import FileLock
import yaml

if len(sys.argv) != 3:
    print "usage: %s log-file log-dir" % sys.argv[0]
    sys.exit(127)

logfile = sys.argv[1]
logdir  = sys.argv[2]

# Get dictionary containing log values
wf_info = yaml.load(open(logfile))

# Fix Nones
for k,v in wf_info.iteritems():
    if v is None:
        wf_info[k] = ''

# Add ID
wf_id = wf_info['workflow_name'] + "_" + wf_info['strategy'].replace(".", "_")
wf_info['wf_id'] = wf_id

# Output path
jsfile = op.join(logdir, "reports", "%s.js" % wf_info["scan_id"])

# Lock file
lock = FileLock(jsfile)
with lock:
    # Append workflow info for this pipeline
    js = open(jsfile, 'a')
    js.write("wf_info[%(pipeline_index)s].push({\n" % wf_info)
    js.write("\tsubject_id: '%(subject_id)s', \n" % wf_info) 
    js.write("\tscan_id: '%(scan_id)s', \n" % wf_info) 
    js.write("\twf_id: '%(wf_id)s', \n" % wf_info) 
    js.write("\tstrategy: '%(strategy)s', \n" % wf_info) 
    js.write("\twf_name: '%(workflow_name)s', \n" % wf_info)
    js.write("\twf_status: '%(wf_status)s', \n" % wf_info)
    js.write("\ttimestamp: '%(timestamp)s'\n" % wf_info)
    js.write("});\n\n")
    js.close()
