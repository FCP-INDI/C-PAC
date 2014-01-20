"""
This tests the functions in network_centrality/core.py
"""

import os, sys
import numpy as np
from numpy.testing import *

from nose.tools import ok_, eq_, raises, with_setup
from nose.plugins.attrib import attr    # http://nose.readthedocs.org/en/latest/plugins/attrib.html

import sys
sys.path.insert(0, '/home2/data/Projects/CPAC_Regression_Test/nipype-installs/fcp-indi-nipype/running-install/lib/python2.7/site-packages')
sys.path.insert(1, "/home2/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/C-PAC")
sys.path.append("/home/data/PublicProgram/epd-7.2-2-rh5-x86_64/lib/python2.7/site-packages")

from CPAC.network_centrality import degree_centrality, fast_degree_centrality
from CPAC.network_centrality import eigenvector_centrality, fast_eigenvector_centrality


def test_degree_centrality_binarize():
    print "testing centrality binarize"
    
    method      = "binarize"
    nblock      = 20
    nvoxs       = 100
    r_value     = 0.2
    corr_matrix = np.random.random((nblock, nvoxs))
    
    ref  = np.sum(corr_matrix>r_value, axis=1)
    comp = degree_centrality(corr_matrix, r_value, method)
    
    assert_equal(ref, comp)

def test_degree_centrality_weighted():
    print "testing centrality weighted"
    
    method      = "weighted"
    nblock      = 20
    nvoxs       = 100
    r_value     = 0.2
    corr_matrix = np.random.random((nblock, nvoxs))
    
    ref  = np.sum(corr_matrix*(corr_matrix>r_value), axis=1)
    comp = degree_centrality(corr_matrix, r_value, method)
    
    assert_equal(ref, comp)


def test_fast_eigenvector_centrality(ntpts=100, nvoxs=1000):
    print "testing fast_eigenvector_centrality"
    
    # Simulate Data
    import numpy as np
    from CPAC.cwas.subdist import norm_cols
    # Normalize Random Time-Series Data
    m = np.random.random((ntpts,nvoxs))
    m = norm_cols(m)
    # Correlation Data with Range 0-1
    mm = m.T.dot(m) # note that need to generate connectivity matrix here
    
    # Execute
    #from CPAC.network_centrality.core import fast_eigenvector_centrality,slow_eigenvector_centrality
    
    ref  = eigenvector_centrality(mm, verbose=False)  # we need to transform mm to be a distance
    comp = fast_eigenvector_centrality(m, verbose=False)
    
    diff = np.abs(ref-comp).mean()  # mean diff
    print(diff)
    
    ok_(diff < np.spacing(1e2)) # allow minimal difference


def test_fast_on_real_data():
    from pandas import read_table
    from os import path as op
    
    k       = 200
    subdir  = "/home2/data/Projects/CWAS/share/nki/subinfo/40_Set1_N104"
    ffile   = op.join(subdir, "short_compcor_rois_random_k%04i.txt" % k)
    fpaths  = read_table(ffile, header=None)
    fpath   = fpaths.ix[0,0]

    import nibabel as nib
    img = nib.load(fpath)
    dat = img.get_data()

    import numpy as np
    from CPAC.cwas.subdist import norm_cols, ncor
    norm_dat = norm_cols(dat)
    corr_dat = norm_dat.T.dot(norm_dat)
    
    ref    = eigenvector_centrality(corr_dat)
    comp   = fast_eigenvector_centrality(norm_dat)
    
    diff = np.abs(ref-comp).mean()  # mean diff
    print(diff)
    
    ok_(diff < np.spacing(1e10)) # allow some differences

