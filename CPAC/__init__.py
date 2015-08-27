"""
Configurable Pipeline for the Analysis of Connectomes
=====================================================

CPAC is a configurable, open-source, Nipype-based, automated processing pipeline for resting state functional MRI (R-fMRI) data, for use by both novice and expert users.
"""

import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs': [np.get_include()]})

from numpy.testing import nosetester
test = nosetester.NoseTester.test

class _NoseTester(nosetester.NoseTester):
    def test(self, label='fast', verbose=1, extra_argv=['--exe'], doctests = False, coverage=False):
        return super(_NoseTester, self).test(label=label, verbose=verbose, extra_argv=extra_argv, doctests=doctests, coverage=coverage)

test = _NoseTester().test

import anat_preproc, \
       func_preproc, \
       reho, \
       seg_preproc, \
       registration, \
       sca, \
       basc, \
       nuisance, \
       generate_motion_statistics, \
       alff, \
       qc, \
       seg_preproc, \
       vmhc, \
       median_angle, \
       timeseries, \
       network_centrality, \
       scrubbing, \
       group_analysis, \
       easy_thresh,\
       utils, \
       pipeline, \
       cwas, \
       GUI

__all__ = ['GUI', 'pipeline', 'anat_preproc', 'func_preproc', 'registration', 'seg_preproc', 'reho', 'sca', 'basc', 'nuisance', 'alff', 'vmhc', 'median_angle', 'generate_motion_statistics', 'timeseries', 'network_centrality', 'scrubbing', 'utils', 'group_analysis', 'easy_thresh', 'qc', 'cwas']

from subprocess import Popen, PIPE
import re
#__version__ = '0.1-git'

try:
    version = '0.3.9.2'


#    gitproc = Popen(['git', 'log', '--oneline'], stdout = PIPE)
#    (stdout, stderr) = gitproc.communicate()
#    rows = stdout.split("\n")
#    v_num = re.search( r'(?<=(version_|release_))(.)*', rows[0])
#    if v_num:
#        version = v_num.group(0).strip("'")
#    else:
#        version = 'unknown_version'
except OSError:
    version = 'unknown_version'

__version__ =  str(version)
