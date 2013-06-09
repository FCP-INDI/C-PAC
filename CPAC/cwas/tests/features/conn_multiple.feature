Feature: Connectivity via Pearson Correlations for multiple subjects

@connectivity
Scenario: Correlation for simulated multi-subject data
    Given simulated subjects time-series data
    When we norm the subjects data
    and we compute the connectivity for multiple subjects on normed data
    and we compute the connectivity for multiple subjects with numpy
    Then the correlation values for cpac should be like numpy

@connectivity, @adhd200
Scenario: Correlation for nki data
    Given the subject list "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd04_funcpaths_4mm.txt"
    and mask data from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/rois/adhd04_mask_gray_4mm.nii.gz"
    and the subjects data are masked
    When we norm the subjects data
    and we compute the connectivity for multiple subjects on normed data
    Then the correlation values for cpac should be like numpy
