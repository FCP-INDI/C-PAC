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

    def to_numpy(self):
        c = np.zeros((self.n, self.n))
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    c[i][j] = self.diagonal[i]
                else:
                    c[i][j] = c[j][i] = self.values[i + j - 1]
        return c



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
    
