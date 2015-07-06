#%% tsdata_to_autocov
#
# Calculate sample autocovariance sequence from time series data
#
# <matlab:open('tsdata_to_autocov.m') code>
#
#% Syntax
#
#     G = tsdata_to_autocov(X,q)
#
#%% Arguments
#
# See also <mvgchelp.html#4 Common variable names and data structures>.
#
# _input_
#
#     X          multi-trial time series data
#     q          number of lags
#
# _output_
#
#     G          sample autocovariance sequence
#
#%% Description
#
# Returns |q|-lag sample autocovariance sequence |G| defined as
# [[ii_acseq.png]] for the (presumed stationary) multivariate process |X|.
# |X| may contain single- or multi-trial time series data.
#
#%% _*Note 1:*_ This routine is discouraged for VAR numerical modelling, and is
# only included for completeness; sample autocovariances are notoriously noisy
# and biased (but see the experimental <tsdata_to_autocov_debias.html
# |tsdata_to_autocov_debias|>). The recommended practice is to estimate a VAR
# model via <tsdata_to_var.html |tsdata_to_var|> and then calculate
# autocovariance via <var_to_autocov.html |var_to_autocov|>.
#
# _*Note 2:*_ For multi-trial data we don't calculate autocovariance on a
# per-trial basis, since this doesn't really make sense... trials in multi-trial
# data must be assumed to be from the same distribution. If you feel you
# absolutely have to calculate per-trial autocovariance (not recommended), call
# this function for each trial series |X(:,:,r)| and average the results over
# trials. Alternatively, if you feel you have to at least _demean_ per-trial
# (not recommended), call <demean.html |demean|> for each trial series
# |X(:,:,r)| _before_ calling this routine.
#
#%% References
#
# [1] L. Barnett and A. K. Seth,
# <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
#     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
# Inference>, _J. Neurosci. Methods_ 223, 2014
# [ <matlab:open('mvgc_preprint.pdf') preprint> ].
#
#%% See also
#
# <demean.html |demean|> |
# <tsdata_to_var.html |tsdata_to_var|> |
# <var_to_autocov.html |var_to_autocov|> |
# <tsdata_to_autocov_debias.html |tsdata_to_autocov_debias|>
#
# (C) Lionel Barnett and Anil K. Seth, 2012. 




def tsdata_to_autocov(X, q):
    
    import numpy as np
    from matplotlib import pylab

    if len(X.shape) == 2:
        X = np.expand_dims(X, axis=2)
        [n, m, N] = np.shape(X)
    else:
        [n, m, N] = np.shape(X)

    X = pylab.demean(X, axis=1) ## This is correct
    G = np.zeros((n, n, (q+1)))
    
    for k in range(q+1):
        M = N * (m-k)
        G[:,:,k] = np.dot(np.reshape(X[:,k:m,:], (n, M)), np.reshape(X[:,0:m-k,:], (n, M)).conj().T) / (M-1)
    return G
    
def autocov_to_pwcgc(G):
    
    import numpy as np

    n = G.shape[0]
    F = np.zeros([n,n])
    
    # full regression   
    [AF,SIG] = autocov_to_var(G)    
    LSIG = np.log(abs(np.diag(SIG)))
    
    for j_ in range(n):
    
        # reduced regression
        
        jo = np.arange(n) # omit j
        jo = np.delete(jo,j_)

        ixgrid1 = np.ix_(jo,jo)
        [AF,SIGj] = autocov_to_var(G[ixgrid1])
          
        LSIGj = np.log(abs(np.diag(SIGj)))
    
        for ii_ in range(n-1):
            i_ = jo[ii_]
            F[i_,j_] = LSIGj[ii_]-LSIG[i_]
        
    return F
    

#%% autocov_to_mvgc
#%
#% Calculate conditional time-domain MVGC (multivariate Granger causality)
#%
#% <matlab:open('autocov_to_mvgc.m') code>
#%
#%% Syntax
#%
#%     F = autocov_to_mvgc(G,x,y)
#%
#%% Arguments
#%
#% See also <mvgchelp.html#4 Common variable names and data structures>.
#%
#% _input_
#%
#%     G          autocovariance sequence
#%     x          vector of indices of target (causee) multi-variable
#%     y          vector of indices of source (causal) multi-variable
#%
#% _output_
#%
#%     F          Granger causality
#%
#%% Description
#%
#% Returns the time-domain MVGC
#%
#% <<eq_mvgc.png>>
#%
#% from the variable |Y| (specified by the vector of indices |y|) to the
#% variable |X| (specified by the vector of indices |x|), conditional on all
#% other variables |Z| represented in |G|, for a stationary VAR process with
#% autocovariance sequence |G|. See ref. [1] for details.
#%
#% The caller should take note of any warnings issued by this function and test
#% results with a call <isbad.html |isbad|>|(F,false)|.
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
#% <autocov_to_var.html |autocov_to_var|> |
#% <isbad.html |isbad|>
#%
#% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
#% installation directory for licensing terms.
#%
#%%
def autocov_to_mvgc(G, x, y):
    
    import numpy as np
    from series_mod import autocov_to_var

    # Local Variables: G, F, xzy, n, xz, SIGR, SIG, y, x, z
    # Function calls: autocov_to_var, log, det, NaN, length, autocov_to_mvgc, size
    n = G.shape[0]
    
    
     #WARNING PART 1!!
#    x = x.flatten(0).conj() #be carefull for multi variate case
#    #% vectorise
#    y = y.flatten(0).conj()
#    #% vectorise
    
    
    z = np.arange(n)
    z = np.delete(z,[np.array(np.hstack((x, y)))])
    #% indices of other variables (to condition out)
    xz = np.array(np.hstack((x, z)))
    xzy = np.array(np.hstack((xz, y)))
    F = 0
    #% full regression
    #%owstate = warn_supp;
    ixgrid1 = np.ix_(xzy,xzy)
    [AF,SIG] = autocov_to_var(G[ixgrid1]) #G[ixgrid1,:])
    #%warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
    #%if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!
    #% reduced regression
    #%owstate = warn_supp;
    ixgrid2 = np.ix_(xz,xz)
    [AF,SIGR] = autocov_to_var(G[ixgrid2]) #G[ixgrid2,:])
    #% reduced regression
    #%warn_test(owstate,     'in reduced regression - bad autocovariance matrix? Check output of ''var_info''');
    #%if warn_if(isbad(SIGR),'in reduced regression - regression failed'), return; end % show-stopper!
    
    #WARNING PART 2!!
    #x = np.arange(np.size(x,axis=1)+1)
    
    ###########
    ixgrid3 = np.ix_(x,x)
    F = np.log(abs(np.linalg.det(SIGR[ixgrid3])))-np.log(abs(np.linalg.det(SIG[ixgrid3]))) 
    #####   not probed
    return F
    
    
    
    
    
    
    
    
    
    
    
# autocov_to_var
#
# Calculate VAR parameters from autocovariance sequence
#
# <matlab:open('autocov_to_var.m') code>
#
#%% Syntax
#
#     [A,SIG] = autocov_to_var(G)
#
#%% Arguments
#
# See also <mvgchelp.html#4 Common variable names and data structures>.
#
# _input_
#
#     G          autocovariance sequence
#
# _output_
#
#     A          VAR coefficients matrix
#     SIG        residuals covariance matrix
#
#%% Description
#
# Calculates regression coefficients |A| and residuals covariance matrix
# |SIG| from the autocovariance sequence |G| defined as [[ii_acseq.png]]
# by solving the Yule-Walker equations
#
# <<eq_yweqs.png>>
#
# (where  [[ii_Sigma.png]] = |SIG|). For a |q|-lag autocovariance sequence,
# this routine corresponds to an autoregression of |q| lags. It also
# effects an efficient spectral factorisation if called with the
# autocovariance sequence derived from the cross-power spectral density
# (_e.g._ as calculated by <cpsd_to_autocov.html |cpsd_to_autocov|>).
#
# This routine implements Whittle's recursive LWR algorithm [2] which, for |n|
# variables, performs |2q| separate |n x n| matrix inversions as compared with a
# single |nq x nq| matrix inversion for the conventional "OLS" solution of the
# Yule-Walker equations (see [1]). The LWR algorithm also (unlike OLS)
# guarantees that if the "true" regression model is stable, then the estimated
# model is also stable, even if not of the correct order.
#
# *_Note_*: If the regressions are rank-deficient or ill-conditioned then A may
# be "bad" (i.e. will contain a |NaN| or |Inf|; see <isbad.html |isbad|>) and/or
# warnings will may be issued. The caller should test for both these
# possibilities, by calls to <isbad.html |isbad|> and <warn_supp.html
# |warn_supp|> ... <warn_test.html |warn_test|> respectively. The likely cause
# is that something went wrong in <var_to_autocov.html |var_to_autocov|>, which
# is typically called prior to this function; check the results of the latter
# with a <var_info.html |var_info|> call.
#
#%% References
#
# [1] L. Barnett and A. K. Seth,
# <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
#     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
# Inference>, _J. Neurosci. Methods_ 223, 2014
# [ <matlab:open('mvgc_preprint.pdf') preprint> ].
#
# [2] P. Whittle, "On the fitting of multivariate autoregressions, and the
# approximate canonical factorization of a spectral density matrix",
# _Biometrika_, 50, 1963.
#
#%% See also
#
# <var_to_autocov.html |var_to_autocov|> |
# <cpsd_to_autocov.html |cpsd_to_autocov|> |
# <warn_supp.html |warn_supp|> |
# <warn_test.html |warn_test|> |
# <isbad.html |isbad|> |
# <var_info.html |var_info|>
#
# (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
# installation directory for licensing terms.
#


def autocov_to_var(G):
    
    import numpy as np
    
    # Local Variables: kb, AB, kf, G, G0, AFPREV, ABPREV, k, AF, AAB, q, GF, r, SIG, GB, AAF, qn
    # Function calls: q1, autocov_to_var, reshape, nargout, n, zeros, flipdim, permute
    [n,m,q1] = G.shape; 
    q = q1 - 1
    qn = q * n
    G0 = G[:,:,0]
    # covariance
    GF = np.reshape(G[:,:,1:], (n, qn)).conj().T
  

    # TO DO: Solve FLIPDIM behaviour in python. In Matlab:
    # GB = reshape(permute(flipdim(G(:,:,2:end),3),[1 3 2]),qn,n);

    # forward  autocov sequence
    #GB = np.reshape(np.transpose(flipdim(G[:,:,1:], 2), (0, 2, 1)), (qn, n))

    # for lag=1 (our case in fMRI), flipdim has no effect
    GB = np.reshape(np.transpose(G[:,:,1:], (0, 2, 1)), (qn, n))
   
   
    # backward autocov sequence np.transpose(x, (1, 0, 2))
    AF = np.zeros([n, qn])
    # forward  coefficients
    AB = np.zeros([n, qn])
    # backward coefficients (reversed compared with Whittle's treatment)
    # initialise recursion
    k = 1 # model order
    
    r = q-k
    kf = np.arange(k*n)
    # forward  indices
    kb = np.arange(r*n, qn)
    # backward indices
    AF[:,kf] = np.linalg.lstsq(G0.T,GB.T)[0].T
    AB[:,kb] = np.linalg.lstsq(G0.T,GF.T)[0].T
    
#    a= np.linalg.lstsq(G0.T,GB.T)[0].T
#    b= np.dot(GB,np.linalg.pinv(G0))
        
    SIG = G0-np.dot(AF, GF)
    AF = np.reshape(AF, (n, n, q))
    return AF, SIG
    
    
    
    
    
# TO DO: Test if this loop computes the same as the equivalent in Matlab
# for lag>1 cases
# and loop for lag>1
#   for k in np.arange(2, q+1):
#        
#        
#        DF = GB[(r-1)*n+1:r*n,:] - np.dot(AF[:,kf],GB[kb,:])
#        VB = G0 - np.dot(AB[:,kb],GB[kb,:])
#        
#        AAF = np.dot(DF,np.linalg.inv(VB)); # DF/VB
#        
#        DB = GF[(k-1)*n+1:k*n,:] - np.dot(AB[:,kb],GF[kf,:])
#        VF = np.dot(G0-AF[:,kf],GF[kf,:])
#        
#        AAB = np.dot(DB,np.linalg.inv(VF)); # DB/VF
#        
#        AFPREV = AF[:,kf-1]
#        ABPREV = AB[:,kb-1]
#        r = q-k
#        kf = np.arange(1, (np.dot(k, n))+1)
#        kb = np.arange(np.dot(r, n)+1, (qn)+1)
#        AF[:,kf-1] = np.array(np.hstack((AFPREV-np.dot(AAF, ABPREV), AAF)))
#        AB[:,kb-1] = np.array(np.hstack((AAB, ABPREV-np.dot(AAB, AFPREV))))



