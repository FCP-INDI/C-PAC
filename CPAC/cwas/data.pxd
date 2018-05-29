import numpy as np
cimport numpy as np
cimport cython

ctypedef DTYPE_t (*dist_func)(DTYPE_t*, DTYPE_t*, ITYPE_t,
                              ITYPE_t, ITYPE_t, dist_params*)

cdef DTYPE_t euclidean_distance(DTYPE_t* x1, DTYPE_t* x2,
                                ITYPE_t inc1, ITYPE_t inc2,
                                ITYPE_t n, dist_params* params)

cdef DTYPE_t manhattan_distance(DTYPE_t* x1, DTYPE_t* x2,
                                ITYPE_t inc1, ITYPE_t inc2,
                                ITYPE_t n, dist_params* params)


cdef class DistanceMatrix

cdef class DiagonalMatrix(object):
    cdef public float diagonal
    cdef public float others

cdef class FullDistanceMatrix(DistanceMatrix):
    cdef public int n
    cdef public np.ndarray X

cdef class SymmetricDistanceMatrix(DistanceMatrix):
    cdef public int n
    cdef public int size
    cdef float * values
    cdef float * diagonal

cdef class DistanceMetric(object):
    cdef dist_params params
    cdef dist_func dfunc

cdef gower(DistanceMatrix D)