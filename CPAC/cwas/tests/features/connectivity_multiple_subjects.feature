Feature: Connectivity via Pearson Correlations for multiple subjects

@connectivity
Scenario: Correlation for simulated multi-subject data
    Given simulated subjects time-series data
    When we norm the subjects data
    and we compute the connectivity on normed subjects data
    Then the correlation values for multiple subjects should be like the standard correlation function

@connectivity, @nki
Scenario: Correlation for nki data
    Given the subject list "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_funcpaths_n4.txt"
    and mask data from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/rois/nki_mask_gray_4mm.nii.gz"
    and the subjects data are masked
    When we norm the subjects data
    and we compute the connectivity on normed subjects data
    Then the correlation values for multiple subjects should be like the standard correlation function
