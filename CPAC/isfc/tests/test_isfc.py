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

    # The mean network correlation, at time interval (t, t + t_win) was 
    # defined as the mean of the lower off-diagonal terms of the
    # correlation matrix.
    return isfc_t.mean()

def isfc_window_reliability(subjects):

    # Reliability of network correlation patterns ('Network States') over time.

    # To assess the reliability ofthe network's state ('fingerprints') we
    # correlated the ISFC patterns across the two independent groups at each
    # sliding window.

    # For each condition (rest, word scramble, intact story), we randomly
    # partitioned the group of 36 subjects into two independent groups of
    # 18 subjects,
    k = subjects.shape[0]
    n = subjects.shape[-1]

    window_step = 30

    isfc_means = []

    # We repeated this procedure 100 times,
    permutations = 100
    for _ in range(permutations):

        # beginning each time with a new random partition,
        ix = np.arange(k)
        np.random.shuffle(ix)
        splits = np.array_split(ix, 2)

        # and computed the ISFC patterns in the DMN nodes at each
        # sliding window in each group.
        for group in splits:
            for t in range(0, n, window_step):

                # prevent invalid window
                if t + window_step > n:
                    break

                # Network-based ISFC was calculated over a sliding window
                # t_win, C_hat_t, within each time interval (t, t + t_win).
                isfc_means += [isfc_window(subjects, t, window_step)]

    # and calculated the mean correlation and the s.d. of
    # the mean across the 100 iterations.
    isfc_means = np.array(isfc_means)
    return isfc_means.mean(), isfc_means.std()


def isfc_reliability(subjects):

    k = subjects.shape[0]
    n = subjects.shape[-1]

    isfc_ref = isfc(subjects)
    isfc_means = []

    permutations = 100
    for _ in range(permutations):
        ix = np.arange(k)
        np.random.shuffle(ix)
        splits = np.array_split(ix, 2)

        for group in splits:
            isfc_means += [isfc(subjects)]

    # and calculated the mean correlation and the s.d. of
    # the mean across the 100 iterations.
    isfc_means = np.array(isfc_means)
    return isfc_means.mean(), isfc_means.std()


def test_isfc():

    subjects = np.random.uniform(0.0, 10.0, (
        36,
        2,
        3,
        5,
        200
    ))

    u, o = isfc_reliability(subjects)

    assert np.isclose(u, 0.0, atol=1e-03)
    assert np.isclose(o, 0.0, atol=1e-03)