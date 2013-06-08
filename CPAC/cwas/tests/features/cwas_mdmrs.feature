Feature: Many MDMRs

@cwas, @nki
Scenario: Single MDMR for NKI RS
    Given a connectir-based MDMR folder "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/results_nki.r/iq_meanFD+age+sex.mdmr" looking at the factor "FSIQ"
    and regressors from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_regressors.txt" with columns of interest [1]
    When we compute a single MDMR in python
    Then the many pseudo-F values should be like R
    and the many F permutations should be like R
    and the many p-values should be like R

@cwas, @nki
Scenario: Many MDMRs for NKI RS
    Given a connectir-based MDMR folder "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/results_nki.r/iq_meanFD+age+sex.mdmr" looking at the factor "FSIQ"
    and regressors from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_regressors.txt" with columns of interest [1]
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be like R

@cwas, @nki
Scenario: Vary the permutations of Many MDMRs for NKI RS
    Given a connectir-based MDMR folder "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/results_nki.r/iq_meanFD+age+sex.mdmr" looking at the factor "FSIQ"
    and regressors from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_regressors.txt" with columns of interest [1]
    and 1000 permutations
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be similar to R
    