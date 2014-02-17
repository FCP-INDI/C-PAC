#!/usr/bin/env python

"""
This are tests to compare CPAC's degree centrality output 
to a reference in this case AFNI's 3dTcorr...

Also needs to be cleaned up into specific functions.
"""

import os, sys
import numpy as np
from numpy.testing import *

from nose.tools import ok_, eq_, raises, with_setup
from nose.plugins.attrib import attr    # http://nose.readthedocs.org/en/latest/plugins/attrib.html

import nibabel as nib

import sys
sys.path.append("/home2/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/C-PAC")

#def ref_centrality():
import os
templ = {
    "input": "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/_mask_mask_abide_90percent_gm_4mm/_scan_rest_1_rest/resample_functional_to_template_0/residual_wtsimt_flirt.nii.gz", 
    "mask": "/home2/data/Projects/ABIDE_Initiative/CPAC/abide/templates/masks/mask_abide_90percent_gm_4mm.nii.gz", 
    "thresh": 0.266037302736, 
    "out_wt": "tmp_weighted.nii.gz", 
    "out_uwt": "tmp_binarize.nii.gz"
}
cmd = "3dTcorrMap -input %(input)s -mask %(mask)s -Sexpr 'step(r-%(thresh)s)' %(out_uwt)s"
cmd = cmd % templ
ret = os.system(cmd)
    
ref_file = "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/_mask_mask_abide_90percent_gm_4mm/_scan_rest_1_rest/resample_functional_to_template_0/tmp.nii.gz"
comp_file = "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/network_centrality_0/_mask_mask_abide_90percent_gm_4mm/_scan_rest_1_rest/calculate_centrality/degree_centrality_weighted.nii.gz"

ref_img = nib.load(ref_file)
comp_img = nib.load(comp_file)

ref_data = ref_img.get_data()
comp_data = comp_img.get_data()

mask = nib.load(templ["mask"]).get_data().nonzero()
nvoxs = float(len(mask[0]))

assert_equal(ref_data[mask], comp_data[mask])
assert_allclose(ref_data[mask], comp_data[mask])

np.mean(ref_data[mask] - comp_data[mask])




# Binary!

ref_file = "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/_mask_mask_abide_90percent_gm_4mm/_scan_rest_1_rest/resample_functional_to_template_0/tmp_binarize.nii.gz"
comp_file = "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/network_centrality_0/_mask_mask_abide_90percent_gm_4mm/_scan_rest_1_rest/calculate_centrality/degree_centrality_binarize.nii.gz"

ref_img = nib.load(ref_file)
comp_img = nib.load(comp_file)

ref_data = ref_img.get_data()
comp_data = comp_img.get_data()

mask = nib.load(templ["mask"]).get_data().nonzero()
nvoxs = float(len(mask[0]))

diff_data = np.abs(ref_data[mask] - comp_data[mask])
w = diff_data.nonzero()

# We allow for some differences but the majority of voxels are assumed to be the same
ok_(w[0].shape[0] < 25)
ok_(diff_data[w].mean() <= 1)




from CPAC.cwas.subdist import norm_cols, ncor

ts_img = nib.load(templ['input'])
ts_data = ts_img.get_data()
ts = ts_data[mask]
norm_ts = norm_cols(ts.T)
corr_matrix = np.nan_to_num(ts.T.dot(ts))

binarize = np.sum(corr_matrix > templ['thresh'], axis=1)
binarize[binarize!=0] = binarize[binarize!=0] - 1
