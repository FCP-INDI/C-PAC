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
    import scipy

    # Local Variables: G, F, xzy, n, xz, SIGR, SIG, y, x, z
    # Function calls: autocov_to_var, log, det, NaN, length, autocov_to_mvgc, size
    n = G.shape[0]
    
    x = x.flatten(0).conj() #be carefull for multi variate case
    #% vectorise
    y = y.flatten(0).conj()
    #% vectorise
    
    
    z = np.arange(n)
    z = np.delete(z,[np.array(np.hstack((x, y)))])
    #% indices of other variables (to condition out)
    xz = np.array(np.hstack((x, z)))
    xzy = np.array(np.hstack((xz, y)))
    F = NaN
    #% full regression
    #%owstate = warn_supp;
    [~,SIG] = autocov_to_var(G[xzy,xzy,:])
    #%warn_test(owstate,    'in full regression - bad autocovariance matrix? Check output of ''var_info''');
    #%if warn_if(isbad(SIG),'in full regression - regression failed'), return; end % show-stopper!
    #% reduced regression
    #%owstate = warn_supp;
    [~,SIGR] = autocov_to_var(G[xz,xz,:])
    #% reduced regression
    #%warn_test(owstate,     'in reduced regression - bad autocovariance matrix? Check output of ''var_info''');
    #%if warn_if(isbad(SIGR),'in reduced regression - regression failed'), return; end % show-stopper!
    x = np.arange(1, (length(x))+1)
    
    ###########
    F = np.log(plt.det(SIGR[x-1,x-1]))-np.log(plt.det(SIG[x-1,x-1])) 
    #####   not probed
    
    
    return [F]
