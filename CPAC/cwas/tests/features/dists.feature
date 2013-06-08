Feature: Distances Between Subjects From Connectivity Data

@distances, @simulation, @single_voxel
Scenario: Compute distances for multi-subject data with a single voxel
    Given simulated subjects connectivity data for 1 seed voxels
    When we compute the distances with cpac
    and we compute the distances with numpy
    Then the cpac distances should be like numpy

@distances, @simulation, @multiple_voxels
Scenario: Compute distances for multi-subject data with multiple voxels
    Given simulated subjects connectivity data for 12 seed voxels
    When we compute the distances with cpac
    and we compute the distances with numpy
    Then the cpac distances should be like numpy

#@distances, @nki
#Scenario: Correlation for nki data
#    Given the subject list "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_funcpaths_n4.txt"
#    and mask data from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/rois/nki_mask_gray_4mm.nii.gz"
#    and the subjects data are masked
#    When we norm the subjects data
#    and we compute the connectivity on normed subjects data
#    Then the correlation values for multiple subjects should be like the standard correlation function
#