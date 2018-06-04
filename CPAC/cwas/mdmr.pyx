cimport cython
cimport numpy as np
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

# cimport data
# import data

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


@cython.boundscheck(False)
@cython.wraparound(False)
def gower_center(yDis):
    n = yDis.shape[0]
    I = np.eye(n,n)
    uno = np.ones((n,1))
    
    A = -0.5 * (yDis ** 2)
    C = I - (1.0 / n) * uno.dot(uno.T)
    G = C.dot(A).dot(C)
    
    return G


@cython.boundscheck(False)
@cython.wraparound(False)
def gower_center_many(dmats):
    nobs    = np.sqrt(dmats.shape[0])
    ntests  = dmats.shape[1]
    Gs      = np.zeros_like(dmats)
    
    for i in range(ntests):
        Dmat    = dmats[:,i].reshape(nobs, nobs)
        Gs[:,i] = gower_center(Dmat).flatten()
    
    return Gs


@cython.boundscheck(False)
@cython.wraparound(False)
def gen_h(x, cols=None, perms_ix=None):
    if perms_ix is not None:
        x = x.copy()
        return hatify(x[perms_ix][:, np.array(cols)])
    return hatify(x)


@cython.boundscheck(False)
@cython.wraparound(False)
def gen_h2_perms(x, cols, perms):
    permutations = perms.shape[0]
    observations = perms.shape[1]
    variables = x.shape[1]
    
    H2_perms = np.zeros((observations ** 2, permutations))
    for i in range(permutations):
        H = gen_h(x, cols, perms[i, :])
        other_cols = [
            i for i in range(variables) if i not in cols
        ]
        H2 = H - hatify(x[:, other_cols])

        H2_perms[:, i] = H2.flatten()
    
    return H2_perms


@cython.boundscheck(False)
@cython.wraparound(False)
def gen_ih_perms(x, cols, perms):
    permutations = perms.shape[0]
    observations = perms.shape[1]
    I            = np.eye(observations, observations)
    
    IH_perms = np.zeros((observations ** 2, permutations))
    for i in range(permutations):
        IH = I - gen_h(x, cols, perms[i,:])
        IH_perms[:,i] = IH.flatten()
    
    return IH_perms


@cython.wraparound(False)
def calc_ftest(Hs, IHs, Gs, m2, nm):
    N = Hs.T.dot(Gs)
    D = IHs.T.dot(Gs)
    F = (N / m2) / (D / nm)
    return F


@cython.boundscheck(False)
@cython.wraparound(False)
def check_rank(x):
    k    = x.shape[1]
    rank = np.linalg.matrix_rank(x)
    if rank < k:
        raise Exception("matrix is rank deficient (rank %i vs cols %i)" % (rank, k))


@cython.boundscheck(False)
@cython.wraparound(False)
def hatify(x):
    Q1, _ = np.linalg.qr(x)
    H = Q1.dot(Q1.T)
    return H


@cython.boundscheck(False)
@cython.wraparound(False)
def fperms_to_pvals(F_perms):
    nperms, ntests = F_perms.shape
    pvals = np.zeros(ntests)
    for i in range(ntests):
        j        = (F_perms[:, i] >= F_perms[0, i]).sum().astype('float')
        pvals[i] = j / nperms
    return pvals


@cython.boundscheck(False)
@cython.wraparound(False)
def mdmr(np.ndarray[DTYPE_t, ndim=2] ys,
         np.ndarray[DTYPE_t, ndim=2] x,
         cols,
         int permutations):

    check_rank(x)
    
    # Add intercept
    x = np.hstack((np.ones((x.shape[0], 1)), x))

    ntests = ys.shape[1]
    nobs   = x.shape[0]
    if nobs != np.sqrt(ys.shape[0]):
        raise Exception("# of observations incompatible between x and ys")
    
    Gs = gower_center_many(ys)

    df_among = len(cols)
    df_resid = nobs - x.shape[1]

    permutation_indexes = np.zeros((permutations + 1, nobs), dtype=np.int)
    permutation_indexes[0, :] = range(nobs)  # Omni
    for i in range(1, permutations + 1):
        permutation_indexes[i,:] = np.random.permutation(nobs)

    H2perms = gen_h2_perms(x, cols, permutation_indexes)
    IHperms = gen_ih_perms(x, cols, permutation_indexes)

    F_perms = calc_ftest(H2perms, IHperms, Gs,
                        df_among, df_resid)

    p_vals = fperms_to_pvals(F_perms)
    return p_vals, F_perms[0, :]
