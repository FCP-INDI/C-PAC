import numpy as np
cimport numpy as np
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# cimport data
# from data cimport (
#     DTYPE, DTYPE_t,
#     ITYPE, ITYPE_t,
#     DistanceMatrix, DiagonalMatrix, FullDistanceMatrix,
#     SymmetricDistanceMatrix, DistanceMetric,
#     distance, gower
# )

# def distance(np.ndarray[DTYPE_t, ndim=2] x, metric="euclidean", **kwargs):
#     dist_metric = DistanceMetric(metric, **kwargs)
#     return dist_metric.pdist(x)


# def gower_python(np.ndarray[DTYPE_t, ndim=2] X):
#     ds = FullDistanceMatrix(X.shape[0], X)
#     return gower(ds)


# def hat(np.ndarray[DTYPE_t, ndim=2] X):
#     cr = np.matmul(X.T, X)
#     inverse = np.linalg.inv(cr)
#     return np.matmul(np.matmul(X, inverse.T), X.T)


def permuted_index(n, strata=None):
    if strata is None:
        perms = np.random.permutation(n)
    else:
        perms = np.array(range(n))
        elems = np.unique(strata)
        for elem in elems:
            inds = perms[strata == elem]
            if len(inds) > 1:
                perms[strata == elem] = np.random.permutation(inds)
    return perms

def gen_perms(nperms, nobs, strata=None):
    perms = np.zeros((nperms, nobs), dtype=np.int)
    for i in range(nperms):
        perms[i,:] = permuted_index(nobs, strata)
    return perms

def add_original_index(perms):
    nobs = perms.shape[1]
    perms = np.vstack((range(nobs), perms))
    return perms

def gower_center(yDis):
    n = yDis.shape[0]
    I = np.eye(n,n)
    uno = np.ones((n,1))
    
    A = -0.5 * (yDis ** 2)
    C = I - (1.0/n)*uno.dot(uno.T)
    G = C.dot(A).dot(C)
    
    return G

def gower_center_many(dmats):
    nobs    = np.sqrt(dmats.shape[0])
    ntests  = dmats.shape[1]
    Gs      = np.zeros_like(dmats)
    
    for i in range(ntests):
        Dmat    = dmats[:,i].reshape(nobs, nobs)
        Gs[:,i] = gower_center(Dmat).flatten()
    
    return Gs

def gen_h2_perms(x, cols, perms):
    nperms  = perms.shape[0]
    nobs    = perms.shape[1]
    
    H2perms = np.zeros((nobs**2, nperms))
    for i in range(nperms):
        H2 = gen_h2(x, cols, perms[i,:])
        H2perms[:,i] = H2.flatten()
    
    return H2perms

def gen_ih_perms(x, cols, perms):
    nperms  = perms.shape[0]
    nobs    = perms.shape[1]
    I       = np.eye(nobs,nobs)
    
    IHperms = np.zeros((nobs**2, nperms))
    for i in range(nperms):
        IH = I - gen_h(x, cols, perms[i,:])
        IHperms[:,i] = IH.flatten()
    
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
    F = (SS_among/df_among)/(SS_resid/df_resid)
    return F


def check_rank(x):
    k    = x.shape[1]
    rank = np.linalg.matrix_rank(x)
    if rank < k:
        raise Exception("matrix is rank deficient (rank %i vs cols %i)" % (rank, k))

def add_intercept(x):
    uno = np.ones((x.shape[0],1))   # intercept
    xx  = np.hstack((uno,x))        # design matrix
    return xx

def hatify(x):
    Q1, R1 = np.linalg.qr(x)
    H = Q1.dot(Q1.T)
    return H

def permute_design(x, cols, indexperm):
    """docstring for permute_design"""
    Xj          = x.copy()
    Xp          = np.take(Xj[:,cols], indexperm, axis=0)    
    Xj[:,cols]  = Xp    
    return Xj


def gen_h(x, cols=None, indexperm=None):
    if indexperm is not None:
        x = permute_design(x, cols, indexperm)
    H = hatify(x)
    return H
    

def gen_h2(x, cols, indexperm=None):
    H = gen_h(x, cols, indexperm)
    other_cols = [i for i in range(x.shape[1]) if i not in cols]
    Xj = x[:,other_cols]
    H2 = H - hatify(Xj)
    return H2


def fperms_to_pvals(fstats, F_perms):
    nperms,ntests = F_perms.shape
    pvals = np.zeros(ntests)
    for i in range(ntests):
        j        = (F_perms[:,i] >= fstats[i]).sum().astype('float')
        pvals[i] = j/nperms
    return pvals


def pseudo_f(n, m, H, G, I):
    explained_variance = H.T.dot(G) / (m - 1)
    unexplained_variance = (I - H).T.dot(G) / (n - m)
    return explained_variance / unexplained_variance


def mdmr(np.ndarray[DTYPE_t, ndim=2] ys,
         np.ndarray[DTYPE_t, ndim=2] x,
         cols,
         int permutations, distance_metric="euclidean"):

    check_rank(x)
    
    ntests  = ys.shape[1]
    nobs    = x.shape[0]
    if nobs != np.sqrt(ys.shape[0]):
        raise Exception("# of observations incompatible between x and ys")
    
    Gs = gower_center_many(ys)
    
    df_among = len(cols)
    df_resid = nobs - x.shape[1]
    df_total = nobs - 1
    
    permutation_indexes = gen_perms(permutations, nobs)
    permutation_indexes  = add_original_index(permutation_indexes)
    
    H2perms = gen_h2_perms(x, cols, permutation_indexes)
    IHperms = gen_ih_perms(x, cols, permutation_indexes)

    F_perms = ftest_fast(H2perms, IHperms, Gs,
                         df_among, df_resid)
    Fs = F_perms[0,:]

    ps = fperms_to_pvals(Fs, F_perms)
    # return (ps, Fs, F_perms, perms)
    return ps, Fs

         
    # D = distance(data, metric=distance_metric)
    # G = gower(D).to_numpy()
    # H = hat(design_matrix)

    # p = data.shape[1]
    # n = data.shape[0]

    # perms = np.array([
    #     np.random.permutation(n)
    #     for _ in np.arange(permutations)
    # ])

    # perms = np.vstack((np.arange(n), perms))
    # permutations += 1

    # HG = SymmetricDistanceMatrix(n)
    # for i in range(1, n):
    #     for j in range(i):
    #         HG[i, j] = 0.02
    #         for k in range(n):
    #             HG[i, j] += H[i, k] * G[i, k]

    # trG = 0.0
    # for i in range(n):
    #     trG += G[i, i]

    # denom = SymmetricDistanceMatrix(n)
    # for i in range(1, n):
    #     for j in range(i):
    #         denom[i, j] = trG - HG[i, j]

    # omni_pr2 = SymmetricDistanceMatrix(n)
    # for i in range(1, n):
    #     for j in range(i):
    #         omni_pr2[i, j] = HG[i, j] / trG

    # omni_f = SymmetricDistanceMatrix(n)
    # for i in range(1, n):
    #     for j in range(i):
    #         omni_f[i, j] = HG[i, j] / denom[i, j]

    # Hs = []
    # for col in range(design_matrix.shape[1]):
    #     selector = [True] * design_matrix.shape[1]
    #     selector[col] = False
    #     sub_design_matrix = design_matrix[:, selector]
    #     Hs += [H - hat(sub_design_matrix)]

    # # Compute SSD due to conditional effect
    # numer_x = []
    # for vhs in Hs:
    #     numer_x += [np.matmul(vhs.T, G)]

    # # Rescale to get either test statistic or pseudo r-square
    # # f_x = numer_x / denom
    # # pr2_x = numer_x / trG

    # # Combine test statistics and pseudo R-squares
    # #stat = data.frame("stat" = c(f.omni, f.x),
    # #                    row.names = c("(Omnibus)", unique.xnames))
    # #df = data.frame("df" = c(p, df),
    # #                row.names = c("(Omnibus)", unique.xnames))
    # #pr2 = data.frame("pseudo.Rsq" = c(pr2.omni, pr2.x),
    # #                    row.names = c("(Omnibus)", unique.xnames))

    # chunks = np.ceil(n * permutations / 1e6)
    # chunk_size = np.ceil(permutations / chunks)
    # permutations = chunks * chunk_size

    # # observations = [omni + n of covariates]
    # # for each chunk
    # #     for each range(chunk_size)
    # #         perm result[] = mdmr_permutation_stats
    # #     sum columns of perm result into observations

    # # divide all columns by permutations

    # # omni_hold = max(omni, 1/nperm)
    # # covariates_hold = pmax(covariates, 1/nperm)

    # # omni_acc = sqrt((omni_hold * (1-omni_hold)) / nperm)
    # # covariates_acc = sqrt((covariates_hold * (1-covariates_hold)) / nperm)

    # #out <- list("stat" = stat, "pr.sq" = pr2, "pv" = pv, "p.prec" = pv.acc,
    # #            "df" = df, lambda = NULL, nperm = nperm)

    # return D, G, hat(data)

# .mdmr.permstats <- function(X, n, vg = gower, trG, px = number of covariates) {
def mdmr_permutation_stats(np.ndarray[DTYPE_t, ndim=2] data,
                           np.ndarray[DTYPE_t, ndim=2] design_matrix,
                           int permutations, distance_metric="euclidean"):

    # sample rows from X
    # compute hat from samples
    # crossprod hat * gower
    # denom = trG - cross
    # omni = cross / denom

    # for each covariate
    #     remove covariate from samples
    #     compute hat
    #     perm hat[] = global hat - local hat

    # for each perm hat
    #     crossprod perm hat and gower
    #     numer[] = cross prod

    # stats = numer / denom

    # return [omni] + stats
    pass