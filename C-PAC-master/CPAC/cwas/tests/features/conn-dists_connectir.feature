Feature: Computing many distances between connectivity maps

@cwas, @adhd200, @single_voxel
Scenario: Calculate distances for a single voxel
    Given a connectir-based distance folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r"
    and the subject list "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd04_funcpaths_4mm.txt"
    and the subjects data are masked
    When we calculate distances for 1 seeds with cpac
    and we calculate distances for 1 seeds with numpy
    Then the cpac distances should be the same as numpy
    and the cpac distances should be like connectir

@cwas, @adhd200, @multiple_voxels
Scenario: Calculate distances for multiple voxels
    Given a connectir-based distance folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r"
    and the subject list "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd04_funcpaths_4mm.txt"
    and the subjects data are masked
    When we calculate distances for 12 seeds with cpac
    and we calculate distances for 12 seeds with numpy
    Then the cpac distances should be the same as numpy
    and the cpac distances should be like connectir


