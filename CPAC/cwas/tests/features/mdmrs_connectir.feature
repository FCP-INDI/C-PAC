Feature: Many MDMRs

@cwas, @mdmr, @adhd200
Scenario: Single MDMR for ADHD200 N of 4
    Given a connectir-based MDMR folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r/diagnosis.mdmr" looking at the factor "diagnosis"
    and regressors from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd04_regressors.txt" with columns of interest [1]
    When we compute a single MDMR in python
    Then the many pseudo-F values should be like R
    and the many F permutations should be like R
    and the many p-values should be like R

@cwas, @mdmr, @adhd200
Scenario: Many MDMRs for ADHD200 N of 4
    Given a connectir-based MDMR folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r/diagnosis.mdmr" looking at the factor "diagnosis"
    and regressors from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd04_regressors.txt" with columns of interest [1]
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be like R

@cwas, @mdmr, @adhd200
Scenario: Vary the permutations of Many MDMRs for ADHD200 N of 4
    Given a connectir-based MDMR folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd04.r/diagnosis.mdmr" looking at the factor "diagnosis"
    and regressors from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd04_regressors.txt" with columns of interest [1]
    and 1000 permutations
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be similar to R

@cwas, @mdmr, @adhd200
Scenario: Single MDMR for ADHD200 N of 40
    Given a connectir-based MDMR folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd40.r/diagnosis.mdmr" looking at the factor "diagnosis"
    and regressors from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd40_regressors.txt" with columns of interest [1]
    When we compute a single MDMR in python
    Then the many pseudo-F values should be like R
    and the many F permutations should be like R
    and the many p-values should be like R

@cwas, @mdmr, @adhd200
Scenario: Many MDMRs for ADHD200 N of 40
    Given a connectir-based MDMR folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd40.r/diagnosis.mdmr" looking at the factor "diagnosis"
    and regressors from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd40_regressors.txt" with columns of interest [1]
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be like R

@cwas, @mdmr, @adhd200
Scenario: Vary the permutations of Many MDMRs for ADHD200 N of 40
    Given a connectir-based MDMR folder "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/results_adhd40.r/diagnosis.mdmr" looking at the factor "diagnosis"
    and regressors from "/home/data/Projects/CPAC_Regression_Test/2013-05-30_cwas/configs/adhd40_regressors.txt" with columns of interest [1]
    and 1000 permutations
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be similar to R
