Feature: Deriving a pseudo-F statistic from a p-value based on distributions of pseudo-F statistics and p-values

@pseudoF, @simulated
Scenario: Compute pseudo-F from a specified p-value
    Given simulated fstats
    and simulated pvals
    and a pval of 0.05
    When we infer an fstat given a pval using CPAC
    and we infer an fstat given a pval using R
    Then the fstat derived from CPAC should be like R
