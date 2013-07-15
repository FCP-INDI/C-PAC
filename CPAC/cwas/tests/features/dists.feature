Feature: Distances Between Subjects From Connectivity Data

@distances, @simulated, @single_voxel
Scenario: Compute distances for multi-subject data with a single voxel
    Given simulated subjects connectivity data for 1 seed voxels
    When we compute the distances with cpac
    and we compute the distances with numpy
    Then the cpac distances should be like numpy

@distances, @simulated, @multiple_voxels
Scenario: Compute distances for multi-subject data with multiple voxels
    Given simulated subjects connectivity data for 12 seed voxels
    When we compute the distances with cpac
    and we compute the distances with numpy
    Then the cpac distances should be like numpy
