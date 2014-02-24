#!/usr/bin/env python

"""
This is a pit of a hodge podge of tests that needs to be cleaned up!
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

from CPAC.network_centrality import calc_centrality
# from .. import calc_centrality

@attr('subject', 'centrality')
def test_calc_centrality():
    
    infile = "/home2/data/Projects/ABIDE_Initiative/CPAC/test_qp/All_Output/sym_links/pipeline_MerrittIsland/_compcor_ncomponents_5_linear1.global1.motion1.quadratic1.compcor1.CSF_0.96_GM_0.7_WM_0.96/0051466_session_1/scan_rest_1_rest/func/functional_mni.nii.gz"
    maskfile = "/home2/data/Projects/ABIDE_Initiative/CPAC/abide/templates/masks/mask_abide_90percent_gm.nii.gz"
    
    tmp1 = calc_centrality(infile, maskfile, [1,0], [1,1], 0, 0.01, 12)
    
    # datafile : string (nifti file)
    #     path to subject data file
    # template : string (nifti file)
    #     path to mask/parcellation unit
    # method_options : list (boolean)
    #     list of two booleans for binarize and weighted options respectively
    # weight_options : list (boolean)
    #     list of two booleans for binarize and weighted options respectively
    # option : an integer
    #     0 for probability p_value, 1 for sparsity threshold, 
    #     any other for threshold value
    # threshold : a float
    #     pvalue/sparsity_threshold/threshold value
    # allocated_memory : string
    #     amount of memory allocated to degree centrality
    
    
    
    
    
    
    config_file     = "files/config_fsl.yml"
    subject_infos   = common.gen_file_map(CPAC_OUTPUT)
    
    curdir = os.getcwd()
    os.chdir(__file__)
    
    conf             = common.load_configuration(config_file)
    old_subject_file = common.load_subject_list(conf.subjectListFile)
    new_subject_file = setup_group_subject_list(config_file, subject_infos)
    
    os.chdir(curdir)
    
    # All should be right here
    assert_equal(np.array(old_subject_file), np.array(new_subject_file))
    
    # Case where only one path is correct
    
    
    """
    I want to test a couple of cases:
    - If nothing changes
    - If the subject list actually has subjects with missing data, it detects that
        (here I can actually just modify the subject_infos or actually recall gen_file_map with a list of subject's desired)
    """
    



def compare_correlation_approaches():
    import time
    
    nvoxs = 100
    ntpts = 25
    timeseries = np.random.random((nvoxs,ntpts))
    
    print "Norming"
    start = time.clock()
    norm_timeseries = norm_cols(timeseries.T)
    print (time.clock() - start)
    
    i=50; j = 1
    
    print "\nRegular"
    start = time.clock()
    ref = calc_corrcoef(timeseries[j:i].T, timeseries.T)
    print (time.clock() - start)
    
    print "\nNew"
    start = time.clock()
    comp = norm_timeseries[:,j:i].T.dot(norm_timeseries)
    print (time.clock() - start)
    
    return np.allclose(ref, comp)
    
    

# First test is to see if the regular output from CPAC's processing matches with AFNI's approach
#def compare_more():
import os
templ = {
    "input": "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/_mask_mask_abide_90percent_gm_4mm/_scan_rest_1_rest/resample_functional_to_template_0/residual_wtsimt_flirt.nii.gz", 
    "mask": "/home2/data/Projects/ABIDE_Initiative/CPAC/abide/templates/masks/mask_abide_90percent_gm_4mm.nii.gz", 
    "thresh": 0.266037302736, 
    "output": "tmp.nii.gz"
}
cmd = "3dTcorrMap -input %(input)s -mask %(mask)s -Cexpr 'step(r-%(thresh)s)*r' %(output)s"
cmd = cmd % templ
os.system(cmd)



def call_calc_centrality():
import time
    
datafile = "/data/Projects/ABIDE_Initiative/CPAC/test_qp/Centrality_Working/resting_preproc_0051466_session_1/_mask_mask_abide_90percent_gm_3mm/_scan_rest_1_rest/resample_functional_to_template_0/residual_wtsimt_flirt.nii.gz"
template = "/home2/data/Projects/ABIDE_Initiative/CPAC/abide/templates/masks/mask_abide_90percent_gm_3mm.nii.gz"
method_options = [True, False]
weight_options = [True, True]
option = 0
threshold = 0.266037302736
allocated_memory = 12

start = time.clock()
calc_centrality(datafile, template, method_options, weight_options, option, threshold, allocated_memory)
end = time.clock()
print (end-start)
    