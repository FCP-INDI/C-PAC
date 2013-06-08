Feature: Computing many distances between connectivity maps

@cwas, @nki, @single_voxel
Scenario: Calculate distances for a single voxel
    Given a connectir-based distance folder "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/results_nki.r"
    and the subject list "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_funcpaths.txt"
    and the subjects data are masked
    When we calculate distances for 1 seeds with cpac
    and we calculate distances for 1 seeds with numpy
    Then the cpac distances should be the same as numpy
    and the cpac distances should be like connectir

@cwas, @nki, @multiple_voxels
Scenario: Calculate distances for multiple voxels
    Given a connectir-based distance folder "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/results_nki.r"
    and the subject list "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_funcpaths.txt"
    and the subjects data are masked
    When we calculate distances for 12 seeds with cpac
    and we calculate distances for 12 seeds with numpy
    Then the cpac distances should be the same as numpy
    and the cpac distances should be like connectir


