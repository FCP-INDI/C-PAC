#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import pearsonr

def isfc(subjects):

    # Network-based ISFC

    # We defined the ISFC correlation matrix between all nodes (voxels or ROIs)
    # across brains in the following manner. The neural signals Xi measured from
    # subject i, i = 1...k, are in the form of a p x n matrix that contains 
    # signals from p neural sources over n time points.
    dims = subjects.shape
    k, p, n = dims[0], np.prod(dims[1:-1]), dims[-1]
    subjects = subjects.reshape((k, p, n))

    # All time courses were z-scored within subjects to zero mean and unit
    # variance. 
    subjects -= subjects.mean(axis=(1, 2), keepdims=True)
    subjects /= subjects.std(axis=(1, 2), keepdims=True)
    C_hat = np.zeros((p, p))

    for i in range(k):
        not_i = list(set(range(k)) - {i})
        other_subjects = subjects[not_i]

        # Thus, the subject-based ISFC was calculated by the Pearson correlation
        # between single subject and the average of all other subjects.
        C = np.dot(subjects[i], other_subjects.mean(axis=0).T) / n

        # Fisher's r-to-z transformation was applied to each correlation
        # coefficient before averaging, in order to increase normality of the
        # distribution of correlation values,
        C_z_transformed = np.arctanh(C)
        C_hat += C_z_transformed

    # and averaged z values were then inverse transformed (z-to-r) to produce
    # average r-values.
    C_hat = np.tanh(C_hat / k)

    # The final ISFC matrix is given by:
    C_hat = (C_hat + C_hat.T) / 2
    # We imposed this symmetry since we consider the correlation between two
    # brain regions as unidirectional, as in FC

    return C_hat


def isfc_window(subjects, start, step):

    # Network correlation patterns over time.

    assert step > 0
    assert start < subjects.shape[-1] >= start + step

    # Network-based ISFC was calculated over a sliding window t_win, C_hat_t,
    # within each time interval (t, t + t_win).
    isfc_t = isfc(subjects[..., start:start+step])

    # We defined the network correlation pattern of a group, at time t,
    # as the lower off-diagonal terms of the symmetric correlation matrix.
    # @ASH: lower off-diagonal terms == lower-triangle of matrix
    isfc_t = isfc_t[np.tril_indices(isfc_t.shape[0], -1)]

    return isfc_t

    # The mean network correlation, at time interval (t, t + t_win) was 
    # defined as the mean of the lower off-diagonal terms of the
    # correlation matrix.
    # return isfc_t.mean()


def isfc_window_reliability_permutation(subjects, window_step=30):

    n = subjects.shape[-1]

    # beginning each time with a new random partition,
    ix = np.arange(subjects.shape[0])
    np.random.shuffle(ix)
    group_1, group_2 = np.array_split(ix, 2)

    # and computed the ISFC patterns in the DMN nodes at each
    # sliding window in each group.
    for t in range(0, n, window_step):

        # prevent invalid window
        if t + window_step > n:
            break

        # Network-based ISFC was calculated over a sliding window
        # t_win, C_hat_t, within each time interval (t, t + t_win).
        r, _ = pearsonr(
            isfc_window(subjects[group_1], t, window_step),
            isfc_window(subjects[group_2], t, window_step)
        )

        yield r


def isfc_window_reliability(subjects, permutations=100, window_step=30):

    # Reliability of network correlation patterns ('Network States') over time.

    # To assess the reliability ofthe network's state ('fingerprints') we
    # correlated the ISFC patterns across the two independent groups at each
    # sliding window.

    # We repeated this procedure 100 times,
    isfc_corrs = []

    for _ in range(permutations):
        isfc_corrs += list(isfc_window_reliability_permutation(subjects))

    isfc_corrs = np.array(isfc_corrs)

    # and calculated the mean correlation and the s.d. of
    # the mean across the 100 iterations.
    return isfc_corrs.mean(), isfc_corrs.std()


def test_isfc():

    subjects = np.random.uniform(0.0, 10.0, (
        36,  # subjects
        4,   # scans
        2,   # x
        3,   # y
        5,   # z
        200  # t
    ))

    for scan in range(subjects.shape[1]):
        u, o = isfc_window_reliability(subjects[:, scan])
        assert 0.0 < np.abs(u) < 0.05
        assert o < 0.1