Feature: Connectivity via Pearson Correlations for a single subject

@connectivity
Scenario: Correlation for simulated data
    Given simulated time-series data 
    When we norm the data
    and we compute the connectivity on normed data
    Then the correlation values should be like the standard correlation function

@connectivity, @nki
Scenario: Correlation for nki data
    Given subject data from "/home2/data/Projects/NKI_ROCKLAND_CPAC_test/Sink/sym_links/pipeline_OakhurstCity/_compcor_ncomponents_5_linear1.motion1.compcor1.CSF_0.98_GM_0.7_WM_0.98/M10902157_session_1/scan_RfMRI_mx_645_rest_RPI/func/bandpass_freqs_0.009.0.1/functional_mni_4mm_smoothed.nii.gz"
    and mask data from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/rois/nki_mask_gray_4mm.nii.gz"
    and the subject data is masked
    When we norm the data
    and we compute the connectivity on normed data
    Then the correlation values should be like the standard correlation function
