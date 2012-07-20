import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def set_gauss(fwhm):

    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string


def pToFile(time_series):

    import os
    import re
    import commands
    import sys

    dir = os.path.dirname(time_series)


    dir1 = re.sub(r'(.)+_seeds_', '', dir)

    dir1 = dir1.split('.nii.gz')
    dir1 = dir1[0]

    ts_oneD = os.path.join(dir, dir1 + '.1D')
    cmd = "cp %s %s" % (time_series, ts_oneD)
    print cmd

    sys.stderr.write('\n' + commands.getoutput(cmd))
    return os.path.abspath(ts_oneD)
