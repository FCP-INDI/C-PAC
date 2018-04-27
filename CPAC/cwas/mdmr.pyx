import numpy as np
cimport numpy as np
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

cdef struct seuclidean_info:
    ITYPE_t n   # size of array
    DTYPE_t* V  # pointer to buffer of size n

cdef union dist_params:
    seuclidean_info seuclidean

ctypedef DTYPE_t (*dist_func)(DTYPE_t*, DTYPE_t*, ITYPE_t,
                              ITYPE_t, ITYPE_t, dist_params*)

cdef DTYPE_t euclidean_distance(DTYPE_t* x1, DTYPE_t* x2,
                                ITYPE_t inc1, ITYPE_t inc2,
                                ITYPE_t n, dist_params* params):
    cdef ITYPE_t i1 = 0, i2 = 0, i1max = inc1 * n
    cdef DTYPE_t d, res = 0
    
    while i1 < i1max:
        d = x1[i1] - x2[i2]
        res += d * d
        i1 += inc1
        i2 += inc2
    
    return res ** 0.5

cdef DTYPE_t manhattan_distance(DTYPE_t* x1, DTYPE_t* x2,
                                ITYPE_t inc1, ITYPE_t inc2,
                                ITYPE_t n, dist_params* params):
    cdef ITYPE_t i1 = 0, i2 = 0, i1max = inc1 * n
    cdef DTYPE_t res = 0
    
    while i1 < i1max:
        res += abs(x1[i1] - x2[i2])
        i1 += inc1
        i2 += inc2

    return res


cdef class DistanceMatrix(object):
    pass

cdef class DiagonalMatrix(object):
    cdef public float diagonal
    cdef public float others

    def __init__(self, diagonal=1.0, others=0.0):
        self.diagonal = diagonal
        self.others = others

    def __getitem__(self, tuple key):
        if key[0] == key[1]:
            return self.diagonal
        else:
            return self.others

cdef class FullDistanceMatrix(DistanceMatrix):

    cdef public int n
    cdef np.ndarray X

    def __init__(self, n=None, np.ndarray[DTYPE_t, ndim=2] X=None):
        if X is not None:
            self.n = X.shape[0]
            self.X = X
        else:
            self.n = n
            self.X = np.zeros([n, n], dtype=np.float64)

    def __getitem__(self, tuple key):
        return self.X[key[0], key[1]]

    def __setitem__(self, tuple key, value):
        self.X[key[0], key[1]] = float(value)

cdef class SymmetricDistanceMatrix(DistanceMatrix):

    cdef public int n
    cdef public int size
    cdef float * values
    cdef float * diagonal

    def __cinit__(self, n):
        self.n = n
        self.size = ((n - 1) * n) / 2
        self.diagonal = <float *> PyMem_Malloc(self.n * sizeof(float))

        self.values = <float *> PyMem_Malloc(self.size * sizeof(float))
        if not self.values:
            raise MemoryError()

        for i in range(self.n):
            self.diagonal[i] = 0.0
        for i in range(self.size):
            self.values[i] = 0.0

    def __dealloc__(self):
        PyMem_Free(self.values)
        PyMem_Free(self.diagonal)

    def __getitem__(self, tuple key):
        if key[0] == key[1]:
            return self.diagonal[key[0]]
        else:
            return self.values[key[0] + key[1] - 1]

    def __setitem__(self, tuple key, value):
        if key[0] != key[1]:
            self.values[key[0] + key[1] - 1] = float(value)
        else:
            self.diagonal[key[0]] = float(value)


cdef class DistanceMetric(object):
    cdef dist_params params
    cdef dist_func dfunc

    def __init__(self, metric="euclidean"):
        if metric in ("euclidean", "l2"):
            self.dfunc = &euclidean_distance

        elif metric in ("manhattan", "cityblock", "l1"):
            self.dfunc = &manhattan_distance

        else:
            raise ValueError('unrecognized metric %s' % metric)
                        
    def pdist(self, X):
        X = np.asarray(X, dtype=DTYPE, order='C')
        assert X.ndim == 2
        m1 = m2 = X.shape[0]
        n = X.shape[1]

        Y = SymmetricDistanceMatrix(X.shape[0])
        self._pdist_c(X, Y)
        return Y

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _pdist_c(self,
                       np.ndarray[DTYPE_t, ndim=2, mode='c'] X,
                       DistanceMatrix Y):
        cdef unsigned int i1, i2, m, n
        m = X.shape[0]
        n = X.shape[1]

        cdef DTYPE_t* pX = <DTYPE_t*> X.data

        for i1 from 0 <= i1 < m:
            for i2 from i1 < i2 < m:
                Y[i1, i2] = self.dfunc(pX + i1 * n,
                                       pX + i2 * n,
                                       1, 1, n, &self.params)


def distance(np.ndarray[DTYPE_t, ndim=2] x, metric="euclidean", **kwargs):
    dist_metric = DistanceMetric(metric, **kwargs)
    return dist_metric.pdist(x)


cdef gower(DistanceMatrix D):
    n = D.n
    norm = 1.0 / n

    C = DiagonalMatrix(diagonal=1.0-norm, others=0.0-norm)

    CD = FullDistanceMatrix(n)
    for i in range(n):
        for j in range(n):
            CD[i, j] = 0.0
            for k in range(n):
                CD[i, j] += C[i, k] * (-0.5 * (D[k, j] ** 2))

    CDC = SymmetricDistanceMatrix(n)
    for i in range(n):
        for j in range(n):
            sum = 0.0
            for k in range(n):
                sum += CD[i, k] * C[k, j]
            CDC[i, j] = sum
                
    return CDC

def gower_python(np.ndarray[DTYPE_t, ndim=2] X):
    ds = FullDistanceMatrix(X.shape[0], X)
    return gower(ds)
    

def hat(np.ndarray[DTYPE_t, ndim=2] X):
    cr = np.matmul(X.T, X)
    inverse = np.linalg.inv(cr)
    return np.matmul(np.matmul(X, inverse.T), X.T)


def pseudo_f(n, m, H, G, I):
    explained_variance = H.T.dot(G) / (m - 1)
    unexplained_variance = (I - H).T.dot(G) / (n - m)
    return explained_variance / unexplained_variance


def mdmr(np.ndarray[DTYPE_t, ndim=2] data,
         np.ndarray[DTYPE_t, ndim=2] design_matrix,
         int permutations, distance_metric="euclidean"):
         
    D = distance(data, metric=distance_metric)
    G = gower(D)
    H = hat(design_matrix)

    p = data.shape[1]
    n = data.shape[0]

    perms = np.array([
        np.random.permutation(n)
        for _ in np.arange(permutations)
    ])

    perms = np.vstack((np.arange(n), perms))
    permutations += 1

    HG = SymmetricDistanceMatrix(n)
    for i in range(1, n):
        for j in range(i):
            HG[i, j] = 0.02
            for k in range(n):
                HG[i, j] += H[i, k] * G[i, k]

    trG = 0.0
    for i in range(n):
        trG += G[i, i]

    denom = SymmetricDistanceMatrix(n)
    for i in range(1, n):
        for j in range(i):
            denom[i, j] = trG - HG[i, j]

    omni_pr2 = SymmetricDistanceMatrix(n)
    for i in range(1, n):
        for j in range(i):
            omni_pr2[i, j] = HG[i, j] / trG

    omni_f = SymmetricDistanceMatrix(n)
    for i in range(1, n):
        for j in range(i):
            omni_f[i, j] = HG[i, j] / denom[i, j]

    Hs = []
    for col in range(design_matrix.shape[1]):
        selector = [True] * design_matrix.shape[1]
        selector[col] = False
        sub_design_matrix = design_matrix[:, selector]
        Hs += [H - hat(sub_design_matrix)]

    # Compute SSD due to conditional effect
    numer_x = [np.cross(vhs, G) for vhs in Hs]

    # Rescale to get either test statistic or pseudo r-square
    f_x = numer_x / denom
    pr2_x = numer_x / trG

    # Combine test statistics and pseudo R-squares
    #stat = data.frame("stat" = c(f.omni, f.x),
    #                    row.names = c("(Omnibus)", unique.xnames))
    #df = data.frame("df" = c(p, df),
    #                row.names = c("(Omnibus)", unique.xnames))
    #pr2 = data.frame("pseudo.Rsq" = c(pr2.omni, pr2.x),
    #                    row.names = c("(Omnibus)", unique.xnames))





    chunks = np.ceil(n * permutations / 1e6)
    chunk_size = np.ceil(permutations / chunks)
    permutations = chunks * chunk_size

    # observations = [omni + n of covariates]
    # for each chunk
    #     for each range(chunk_size)
    #         perm result[] = mdmr_permutation_stats
    #     sum columns of perm result into observations

    # divide all columns by permutations

    # omni_hold = max(omni, 1/nperm)
    # covariates_hold = pmax(covariates, 1/nperm)

    # omni_acc = sqrt((omni_hold * (1-omni_hold)) / nperm)
    # covariates_acc = sqrt((covariates_hold * (1-covariates_hold)) / nperm)

    #out <- list("stat" = stat, "pr.sq" = pr2, "pv" = pv, "p.prec" = pv.acc,
    #            "df" = df, lambda = NULL, nperm = nperm)

    return D, G, hat(data)

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