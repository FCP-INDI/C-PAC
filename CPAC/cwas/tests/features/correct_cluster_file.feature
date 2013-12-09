Feature: Get clusters for a given file

@cluster, @wip
Scenario: Cluster imaging data from a file
    Given image file '/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r/diagnosis.mdmr/zstats_diagnosis.nii.gz'
    and mask file '/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r/mask_r2py.nii.gz'
    and voxel threshold of 1.65
    When we calculate the cluster sizes using CPAC
    and we calculate the cluster sizes using FSL
    Then the cluster sizes derived from CPAC should be like FSL
