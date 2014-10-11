"""
This tests the functions in network_centrality/core.py
"""

import os, sys
import numpy as np
from numpy.testing import *

from nose.tools import ok_, eq_, raises, with_setup
from nose.plugins.attrib import attr    # http://nose.readthedocs.org/en/latest/plugins/attrib.html

import sys
sys.path.insert(0, '/home2/data/Projects/CPAC_Regression_Test/zarrar/centrality_tests/lib/nipype')
sys.path.insert(1, "/home2/data/Projects/CPAC_Regression_Test/zarrar/centrality_tests/lib/C-PAC")

from CPAC.network_centrality import degree_centrality, fast_degree_centrality
from CPAC.network_centrality import eigenvector_centrality, fast_eigenvector_centrality

class TestDegreeCentrality:
    @attr('degree', 'centrality', 'binarize')
    def test_degree_centrality_binarize():
        from CPAC.cwas.subdist import norm_cols, ncor
        
        # Settings
        nvoxblocks  = 50
        nvoxs       = 200
        r_value     = 0.2
        method      = "binarize"
        
        print "testing centrality binarize - float"
        corr_matrix = np.random.random((nvoxblocks, nvoxs)).astype('float32')
        ref  = np.sum(corr_matrix>r_value, axis=1)
        comp = degree_centrality(corr_matrix, r_value, method)
        assert_equal(ref, comp)
        
        print "testing centrality binarize - double"
        corr_matrix = np.random.random((nvoxblocks, nvoxs)).astype('float64')
        ref  = np.sum(corr_matrix>r_value, axis=1)
        comp = degree_centrality(corr_matrix, r_value, method)
        assert_equal(ref, comp)
        
        return
    
    @attr('degree', 'centrality', 'weighted')
    def test_degree_centrality_weighted():
        from CPAC.cwas.subdist import norm_cols, ncor
        
        # Settings
        nvoxblocks  = 50
        nvoxs       = 200
        r_value     = 0.2
        method      = "weighted"
        
        print "testing centrality binarize - float"
        corr_matrix = np.random.random((nvoxblocks, nvoxs)).astype('float32')
        ref  = np.sum(corr_matrix>r_value, axis=1)
        comp = degree_centrality(corr_matrix, r_value, method)
        assert_equal(ref, comp)
        
        print "testing centrality binarize - double"
        corr_matrix = np.random.random((nvoxblocks, nvoxs)).astype('float64')
        ref  = np.sum(corr_matrix>r_value, axis=1)
        comp = degree_centrality(corr_matrix, r_value, method)
        assert_equal(ref, comp)
        
        return
    
    def test_degree_on_real_data():
        # TODO: Replace the mask and func with a standard testing one
        mpath = "/home2/data/Projects/CPAC_Regression_Test/centrality_template/mask-thr50-3mm.nii.gz"
        mask_inds = nib.load(mpath).get_data().nonzero()
    
        fpath   = "/home/data/Projects/CPAC_Regression_Test/2014-02-24_v-0-3-4/run/w/resting_preproc_0010042_session_1/_scan_rest_1_rest/_scan_rest_1_rest/_csf_threshold_0.98/_gm_threshold_0.7/_wm_threshold_0.98/_compcor_ncomponents_5_selector_pc10.linear1.wm0.global0.motion1.quadratic0.gm0.compcor1.csf0/_bandpass_freqs_0.009.0.1/_mask_mask-thr50-3mm/resample_functional_to_template_0/bandpassed_demeaned_filtered_wtsimt_flirt.nii.gz"
    
        import nibabel as nib
        img = nib.load(fpath)
        dat = img.get_data()
        dat = dat[mask_inds]

        import numpy as np
        from CPAC.cwas.subdist import norm_cols, ncor
        norm_dat = norm_cols(dat.T)
        corr_dat = norm_dat.T.dot(norm_dat)
    
        r_value = 0.2
    
        ref     = np.sum(corr_dat[:5,:5]>r_value, axis=1)
        comp    = degree_centrality(corr_dat[:5,:5], r_value, "binarize")
        assert_equal(ref, comp)
    
        ref     = np.sum(corr_dat*(corr_dat>r_value), axis=1)
        comp    = degree_centrality(corr_dat, r_value, "weighted")
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

