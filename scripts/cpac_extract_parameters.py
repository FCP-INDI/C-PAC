#! /usr/bin/env python
import sys
from CPAC.utils.extract_parameters import grab

if __name__ == '__main__':
    if (len(sys.argv) == 2):
        grab(sys.argv[1], [0])
    else:
        print 'Usage: cpac_extract_parameters.py /path/to/output/dir'
