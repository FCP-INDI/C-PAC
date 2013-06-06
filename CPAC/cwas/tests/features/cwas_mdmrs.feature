Feature: Many MDMRs or MDMR Wrapper

Scenario: 
    Given a connectir-based MDMR folder "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/results_nki.r/iq_meanFD+age+sex.mdmr" looking at the factor "FSIQ"
    and regressors from "/home/data/Projects/Z_CPAC_Regression_Test/2013-05-30_cwas/configs/nki_regressors.txt" with columns of interest [1]
    When we compute many mdmrs in python
    Then the many pseudo-F values should be like R
    and the many p-values should be similar to R

@when('we compute many mdmrs in python')
def step(context):
    context.nperms = 1000
    context.Fs, context.ps = calc_mdmrs(context.dmats, context.regressors, context.cols, context.nperms)

@then('the many pseudo-F values should be like R')
def step(context):
    comp = np.allclose(context.Fs, context.r_Fs)
    assert_that(comp, "pseudo-F stats")

@then('the many p-values should be similar to R')
def step(context):
    comp = np.corrcoef(context.ps, context.r_ps)
    assert_that(comp, is greater than 0.98)
