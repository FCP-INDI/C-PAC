Feature: Distances Between Subjects From Functional Time-Series Data

@distances, @simulated
Scenario: Compute distances for multi-subject data with a single voxel blocks
    Given simulated subjects time-series data
    When we compute the distances from functional data with cpac for 1 voxel block
    and we compute the distances from functional data with numpy
    Then the cpac distances should be like numpy

@distances, @simulated
Scenario: Compute distances for multi-subject data with a multi-voxel blocks
    Given simulated subjects time-series data
    When we compute the distances from functional data with cpac for 10 voxel block
    and we compute the distances from functional data with numpy
    Then the cpac distances should be like numpy

@distances, @simulated
Scenario: Compute distances for multi-subject data using a voxel range
    Given simulated subjects time-series data
    When we compute the distances from functional data with cpac using a voxel range
    and we compute the distances from functional data with numpy
    Then the cpac distances should be like numpy
