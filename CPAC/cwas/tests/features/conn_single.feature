Feature: Connectivity via Pearson Correlations for a single subject

@connectivity, @simulated
Scenario: Correlation for simulated data
    Given simulated time-series data 
    When we norm the data
    and we compute the connectivity on normed data
    and we compute the connectivity with numpy
    Then the correlation values for cpac should be like numpy

@connectivity, @adhd200
Scenario: Correlation for nki data
    Given subject data from "/home2/data/Projects/CPAC_Regression_Test/sink_beta_June3/sym_links/pipeline_HackettCity/_compcor_ncomponents_5_linear1.motion1.compcor1.CSF_0.98_GM_0.7_WM_0.98/2014113_session_1/scan_rest_1_rest/func/bandpass_freqs_0.009.0.1/functional_mni_4mm.nii.gz"
    and mask data from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/rois/adhd04_mask_gray_4mm.nii.gz"
    and the subject data is masked
    When we norm the data
    and we compute the connectivity on normed data
    and we compute the connectivity with numpy
    Then the correlation values for cpac should be like numpy
