#!/usr/bin/env python

import argparse
import os
from os import path

import sys
sys.path.insert(0, '/home2/data/Projects/CPAC_Regression_Test/nipype-installs/fcp-indi-nipype/running-install/lib/python2.7/site-packages')
sys.path.insert(1, '/home2/data/Projects/CPAC_Regression_Test/zarrar/centrality_tests/lib/C-PAC')
sys.path.append("/home/data/PublicProgram/epd-7.2-2-rh5-x86_64/lib/python2.7/site-packages")

from CPAC.network_centrality import calc_centrality


###
# Load command-line argument
###

parser = argparse.ArgumentParser(description='Compute centrality for a given timeseries.')

# Inputs
parser.add_argument('-i', '--input', help='Input timeseries data file', 
                    required=True)
parser.add_argument('-m', '--mask', help='Brain mask (by default the program will create a mask where voxels have non-zero variance regardeless of the user specified mask).')

# Option: Method
parser.add_argument('--degree', help='Calculate degree centrality', 
                    action="store_true")
parser.add_argument('--eigen', help='Calculate eigen centrality', 
                    action="store_true")

# Option: Outputs
parser.add_argument('--binarize', action="store_true", 
                    help='For a given voxel, save the number of connections that pass a threshold')
parser.add_argument('--weighted', action="store_true", 
                    help='For a given voxel, save the sum of all connection weights that pass a threshold.')

# Option: Threshold
# TODO: explicitly check that only one is specified.
parser.add_argument('--no-threshold', action='store_true', 
                    help="The raw correlations will be used and the approach will be very fast and low memory usage.")
parser.add_argument('--sparsity', type=float, 
                    help="Sparsity based threshold. (Only one threshold option can be specified.)")
parser.add_argument('--pvalue', type=float, 
                    help='P-value threshold for each connection. (Only one threshold option can be specified.)')
parser.add_argument('--rho', type=float, 
                    help="Regular correlation threshold. (Only one threshold option can be specified.)")

# Option: Zcore
parser.add_argument('--zscore', action='store_true', 
                    help="Z-score each of the output centrality images (will create additional outputs)")

# Option: Memory
parser.add_argument('--memlimit', type=float, 
                    help="Memory limit that should be set.")

# Option: Threads
parser.add_argument('--nthreads', type=int, default=1, 
                    help="Number of threads to parallel processes for Intel MKL")

# Output
parser.add_argument('-o', '--outdir', default=os.getcwd(), help="Output directory")



###
# Parse and Read User Args
###

args = parser.parse_args()

try:
    import mkl
    mkl.set_num_threads(args.nthreads)
except ImportError:
    pass

if not args.degree and not args.eigen:
    raise SystemExit("--degree and/or --eigen must be specified")
method_options = [args.degree, args.eigen]

if not args.binarize and not args.weighted:
    raise SystemExit("--binarize and/or --weighted must be specified")
weight_options = [args.binarize, args.weighted]

if args.pvalue is not None:
    option = 0
    threshold = args.pvalue
elif args.sparsity is not None:
    option = 1
    threshold = args.sparsity
elif args.rho is not None:
    option = 2
    threshold = args.rho
elif args.no_threshold:
    option = 3
    threshold = None
else:
    raise SystemExit("You must specify one threshold option: --pvalue, --sparsity, --rho, or --no-threshold.")


###
# Call on the Big Guy/Gal (CPAC)
###

curdir = os.getcwd()
os.chdir(args.outdir)

out_list = calc_centrality(args.input, args.mask, method_options, weight_options, option, 
                           threshold, args.memlimit)

os.chdir(curdir)


###
# Z-Score
###

def zscore_image(infile, mask, outfile):
    import os, commands
    
    cmd     = "fslstats %s -k %s -m" % (infile, mask)
    print cmd
    mean    = commands.getstatusoutput(cmd)
    mean    = float(mean[1])
    
    cmd     = "fslstats %s -k %s -s" % (infile, mask)
    print cmd
    sd      = commands.getstatusoutput(cmd)
    sd      = float(sd[1])
    
    cmd = "3dcalc -a %s -b %s -expr '((a-%f)/%f)*step(b)' -prefix %s" % (infile, mask, mean, sd, outfile)
    print cmd
    os.system(cmd)
    
    return

print ''
for centfile in out_list:
    print "Z-score: %s" % os.path.basename(centfile)
    
    # New output filename
    prefix,ext  = os.path.splitext(centfile)
    if ext == '.gz':
        prefix,ext  = os.path.splitext(prefix)
        ext = ext + '.gz'
    zfile       = "%s_zscore%s" % (prefix, ext)
    
    # Create zscore image
    zscore_image(centfile, args.mask, zfile)
