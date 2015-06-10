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
    AF[:,kf] = np.dot(GB[kb,:],np.linalg.inv(G0))
    AB[:,kb] = np.dot(GF[kf,:],np.linalg.inv(G0)) 
    
    # TO DO: Test if this loop computes the same as the equivalent in Matlab
    # for lag>1 cases
    # and loop for lag>1
#    for k in np.arange(2, q+1):
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
        
    
    SIG = G0-np.dot(AF, GF)
    AF = np.reshape(AF, (n, n, q))
    return [AF, SIG]