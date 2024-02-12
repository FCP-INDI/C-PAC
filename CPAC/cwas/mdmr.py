import numpy as np

def check_rank(X):
    k    = X.shape[1]
    rank = np.linalg.matrix_rank(X)
    if rank < k:
        raise Exception("matrix is rank deficient (rank %i vs cols %i)" % (rank, k))

def hat(X):
    Q1, _ = np.linalg.qr(X)
    return Q1.dot(Q1.T)

def gower(D):
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    I = np.eye(n, n)
    uno = np.ones((n, 1))
    C = I - (1.0 / n) * uno.dot(uno.T)
    G = C.dot(A).dot(C)
    return G

def gen_h2(x, cols, indexperm):
    H = gen_h(x, cols, indexperm)
    other_cols = [i for i in range(x.shape[1]) if i not in cols]
    Xj = x[:,other_cols]
    H2 = H - hat(Xj)
    return H2

def permute_design(x, cols, indexperm):
    Xj = x.copy()
    Xj[:, cols] = Xj[indexperm][:, cols]
    return Xj

def gen_h(x, cols, indexperm):
    x = permute_design(x, cols, indexperm)
    H = hat(x)
    return H

def gen_h2_perms(x, cols, perms):
    nperms, nobs = perms.shape
    H2perms = np.zeros((nobs**2, nperms))
    for i in range(nperms):
        H2 = gen_h2(x, cols, perms[i,:])
        H2perms[:,i] = H2.flatten()

    return H2perms

def gen_ih_perms(x, cols, perms):
    nperms, nobs = perms.shape
    I = np.eye(nobs, nobs)

    IHperms = np.zeros((nobs ** 2, nperms))
    for i in range(nperms):
        IH = I - gen_h(x, cols, perms[i, :])
        IHperms[:, i] = IH.flatten()

    return IHperms

def calc_ssq_fast(Hs, Gs, transpose=True):
    if transpose:
        ssq = Hs.T.dot(Gs)
    else:
        ssq = Hs.dot(Gs)
    return ssq

def ftest_fast(Hs, IHs, Gs, df_among, df_resid, **ssq_kwrds):
    SS_among = calc_ssq_fast(Hs, Gs, **ssq_kwrds)
    SS_resid = calc_ssq_fast(IHs, Gs, **ssq_kwrds)
    F = (SS_among / df_among) / (SS_resid / df_resid)
    return F

def mdmr(D, X, columns, permutations):

    check_rank(X)
    
    subjects = X.shape[0]
    if subjects != D.shape[1]:
        raise Exception("# of subjects incompatible between X and D")
    
    voxels = D.shape[0]
    Gs = np.zeros((subjects ** 2, voxels))
    for di in range(voxels):
        Gs[:, di] = gower(D[di]).flatten()
    
    X1 = np.hstack((np.ones((subjects, 1)), X))
    columns = columns.copy() #removed a +1

    regressors = X1.shape[1]

    permutation_indexes = np.zeros((permutations, subjects), dtype=int)
    permutation_indexes[0, :] = range(subjects)
    for i in range(1, permutations):
        permutation_indexes[i,:] = np.random.permutation(subjects)
    
    H2perms = gen_h2_perms(X1, columns, permutation_indexes)
    IHperms = gen_ih_perms(X1, columns, permutation_indexes)

    df_among = len(columns)
    df_resid = subjects - regressors

    F_perms = ftest_fast(H2perms, IHperms, Gs, df_among, df_resid)

    p_vals = (F_perms[1:, :] >= F_perms[0, :]) \
                .sum(axis=0) \
                .astype('float')
    p_vals /= permutations

    return F_perms[0, :], p_vals

