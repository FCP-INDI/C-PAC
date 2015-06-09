#%% tsdata_to_autocov
#%
#% Calculate sample autocovariance sequence from time series data
#%
#% <matlab:open('tsdata_to_autocov.m') code>
#%
#%% Syntax
#%
#%     G = tsdata_to_autocov(X,q)
#%
#%% Arguments
#%
#% See also <mvgchelp.html#4 Common variable names and data structures>.
#%
#% _input_
#%
#%     X          multi-trial time series data
#%     q          number of lags
#%
#% _output_
#%
#%     G          sample autocovariance sequence
#%
#%% Description
#%
#% Returns |q|-lag sample autocovariance sequence |G| defined as
#% [[ii_acseq.png]] for the (presumed stationary) multivariate process |X|.
#% |X| may contain single- or multi-trial time series data.
#%
#% _*Note 1:*_ This routine is discouraged for VAR numerical modelling, and is
#% only included for completeness; sample autocovariances are notoriously noisy
#% and biased (but see the experimental <tsdata_to_autocov_debias.html
#% |tsdata_to_autocov_debias|>). The recommended practice is to estimate a VAR
#% model via <tsdata_to_var.html |tsdata_to_var|> and then calculate
#% autocovariance via <var_to_autocov.html |var_to_autocov|>.
#%
#% _*Note 2:*_ For multi-trial data we don't calculate autocovariance on a
#% per-trial basis, since this doesn't really make sense... trials in multi-trial
#% data must be assumed to be from the same distribution. If you feel you
#% absolutely have to calculate per-trial autocovariance (not recommended), call
#% this function for each trial series |X(:,:,r)| and average the results over
#% trials. Alternatively, if you feel you have to at least _demean_ per-trial
#% (not recommended), call <demean.html |demean|> for each trial series
#% |X(:,:,r)| _before_ calling this routine.
#%
#%% References
#%
#% [1] L. Barnett and A. K. Seth,
#% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
#%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
#% Inference>, _J. Neurosci. Methods_ 223, 2014
#% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
#%
#%% See also
#%
#% <demean.html |demean|> |
#% <tsdata_to_var.html |tsdata_to_var|> |
#% <var_to_autocov.html |var_to_autocov|> |
#% <tsdata_to_autocov_debias.html |tsdata_to_autocov_debias|>
#%
#% (C) Lionel Barnett and Anil K. Seth, 2012. 




def tsdata_to_autocov(X, q):
    
    import numpy as np
    from matplotlib import pylab

    # Local Variables: G, M, k, m, N, q, X, n
    # Function calls: zeros, tsdata_to_autocov, demean, reshape, size
    if len(X.shape) == 2:
        X = np.expand_dims(X, axis=2)
        [n, m, N] = np.shape(X)
    else:
        [n, m, N] = np.shape(X)
    #%assert(q < m-1,'too many lags');
    X = pylab.demean(X, axis=1) ## This is correct
    G = np.zeros((n, n, (q+1)))
    
    for k in range(0,q+1):
        M = N * m-k
        G[:,:,k] = np.dot(np.reshape(X[:,k:m,:], (n, M)), np.reshape(X[:,0:m-k,:], (n, M)).conj().T) / M-1
        #c = numpy.linalg.lstsq(b.T, a.T)[0].T
    return [G]