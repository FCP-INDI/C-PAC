cimport cython
cimport numpy as np
import numpy as np
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
def check_rank(X):
    k    = X.shape[1]
    rank = np.linalg.matrix_rank(X)
    if rank < k:
        raise Exception("matrix is rank deficient (rank %i vs cols %i)" % (rank, k))


@cython.boundscheck(False)
@cython.wraparound(False)
def hatify(X):
    Q1, _ = np.linalg.qr(X)
    H = Q1.dot(Q1.T)
    return H


@cython.boundscheck(False)
@cython.wraparound(False)
def gower_center(Y):
    n = Y.shape[0]
    I = np.eye(n,n)
    uno = np.ones((n, 1))
    
    A = -0.5 * (Y ** 2)
    C = I - (1.0 / n) * uno.dot(uno.T)
    G = C.dot(A).dot(C)
    
    return G


@cython.boundscheck(False)
@cython.wraparound(False)
def gower_center_many(Ys):
    observations = np.sqrt(Ys.shape[0])
    tests        = Ys.shape[1]
    Gs           = np.zeros_like(Ys)
    
    for i in range(tests):
        D        = Ys[:, i].reshape(observations, observations)
        Gs[:, i] = gower_center(D).flatten()
    
    return Gs


@cython.boundscheck(False)
@cython.wraparound(False)
def gen_h(x, columns, permutations):
    return hatify(x[permutations][:, np.array(columns)])


@cython.boundscheck(False)
@cython.wraparound(False)
def gen_H2_perms(X, columns, permutation_indexes):
    permutations, observations = permutation_indexes.shape
    variables = X.shape[1]
    
    H2_permutations = np.zeros((observations ** 2, permutations))
    for i in range(permutations):
        perm_X = X[permutation_indexes[i, :]]
        cols_X = perm_X[:, columns]
        H = hatify(cols_X)
        other_columns = [i for i in range(variables) if i not in columns]
        H2 = H - hatify(X[:, other_columns])
        H2_permutations[:, i] = H2.flatten()
    
    return H2_permutations


@cython.boundscheck(False)
@cython.wraparound(False)
def gen_IH_perms(X, columns, permutation_indexes):
    permutations, observations = permutation_indexes.shape
    I            = np.eye(observations, observations)
    
    IH_permutations = np.zeros((observations ** 2, permutations))
    for i in range(permutations):
        IH = I - hatify(X[permutation_indexes[i, :]][:, columns])
        IH_permutations[:,i] = IH.flatten()
    
    return IH_permutations


@cython.wraparound(False)
def calc_ftest(Hs, IHs, Gs, m2, nm):
    N = Hs.T.dot(Gs)
    D = IHs.T.dot(Gs)
    F = (N / m2) / (D / nm)
    return F


@cython.boundscheck(False)
@cython.wraparound(False)
def fperms_to_pvals(F_perms):
    permutations, tests = F_perms.shape
    pvals = np.zeros(tests)
    for i in range(tests):
        j        = (F_perms[:, i] >= F_perms[0, i]).sum().astype('float')
        pvals[i] = j / permutations
    return pvals


@cython.boundscheck(False)
@cython.wraparound(False)
def mdmr(np.ndarray[DTYPE_t, ndim=2] D,
         np.ndarray[DTYPE_t, ndim=2] X,
         np.ndarray[ITYPE_t, ndim=1] columns,
         int permutations):

    check_rank(X)

    subjects = X.shape[0]
    if subjects != np.sqrt(D.shape[0]):
        raise Exception("# of subjects incompatible between X and D")
    
    X = np.hstack((np.ones((X.shape[0], 1)), X))
    columns = columns.copy()
    columns += 1

    Gs = gower_center_many(D)

    df_among = float(columns.shape[0])
    df_resid = float(subjects - X.shape[1])

    permutation_indexes = np.zeros((permutations + 1, subjects), dtype=np.int)
    permutation_indexes[0, :] = range(subjects)
    for i in range(1, permutations + 1):
        permutation_indexes[i,:] = np.random.permutation(subjects)

    H2perms = gen_H2_perms(X, columns, permutation_indexes)
    IHperms = gen_IH_perms(X, columns, permutation_indexes)

    F_perms = calc_ftest(H2perms, IHperms, Gs,
                        df_among, df_resid)

    p_vals = fperms_to_pvals(F_perms)
    
    return F_perms[0, :], p_vals 
