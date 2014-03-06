"""
This tests the functions in network_centrality/thresh_and_sum.pyx
"""

import os, sys
import numpy as np
from numpy.testing import *

from nose.tools import ok_, eq_, raises, with_setup
from nose.plugins.attrib import attr    # http://nose.readthedocs.org/en/latest/plugins/attrib.html

import sys
sys.path.insert(0, '/home2/data/Projects/CPAC_Regression_Test/zarrar/centrality_tests/lib/nipype')
sys.path.insert(1, "/home2/data/Projects/CPAC_Regression_Test/zarrar/centrality_tests/lib/C-PAC")

# For eigen centrality
from CPAC.network_centrality.thresh_and_sum import \
        thresh_binarize_float, thresh_binarize_double, \
        thresh_weighted_float, thresh_weighted_double, \
        thresh_transform_weighted_float, thresh_transform_weighted_double

# For degree centrality
from CPAC.network_centrality.thresh_and_sum import \
        centrality_binarize_float, centrality_binarize_double, \
        centrality_weighted_float, centrality_weighted_double, \
        centrality_both_float, centrality_both_double    # these aren't currently used


###
# TEST thresholding of matrices for eigenvector centrality
###

class TestThresholding:
    @attr('threshold', 'binarize')
    def test_thresh_binarize():
        print "testing threshold binarize"
    
        nvoxs       = 1000
        r_value     = 0.2
        corr_matrix = np.random.random((nvoxs, nvoxs)).astype('float32')

        ref  = 1*(corr_matrix>r_value)

        comp = corr_matrix.copy()
        thresh_binarize_float(comp, r_value)

        assert_equal(ref, comp)

    @attr('threshold', 'weighted')
    def test_thresh_weighted():
        print "testing threshold weighted"
    
        nvoxs       = 1000
        r_value     = 0.2
        corr_matrix = np.random.random((nvoxs, nvoxs)).astype('float32')
    
        ref  = corr_matrix*(corr_matrix>r_value)

        comp = corr_matrix.copy()
        thresh_weighted_float(comp, r_value)
    
        assert_equal(ref, comp)
    
    @attr('threshold', 'transform', 'weighted')
    def test_thresh_transform_weighted():
        print "testing threshold weighted"
    
        nvoxs       = 1000
        r_value     = 0.2
        corr_matrix = np.random.random((nvoxs, nvoxs)).astype('float32')
    
        ref  = ((1.0+corr_matrix)/2.0)*(corr_matrix>r_value)

        comp = corr_matrix.copy()
        thresh_transform_weighted_float(comp, r_value)
    
        assert_equal(ref, comp)


###
# TEST centrality functions
###
