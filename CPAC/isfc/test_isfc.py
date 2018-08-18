import numpy as np
from scipy.stats import pearsonr


def test_isfc():

    subjects = np.random.uniform(0.0, 1.0, (
        20,
        2,
        3,
        5,
        200
    ))

    # We defined the ISFC correlation matrix between all nodes (voxels or ROIs)
    # across brains in the following manner. The neural signals Xi measured from
    # subject i, i = 1...k, are in the form of a p x n matrix that contains 
    # signals from p neural sources over n time points.
    k, x, y, z, n = subjects.shape
    p = x * y * z
    subjects = subjects.reshape((k, p, n))

    # All time courses were z-scored within subjects to zero mean and unit variance. 
    subjects -= subjects.mean(axis=(1, 2), keepdims=True)
    subjects /= subjects.std(axis=(1, 2), keepdims=True)
    C_hat = np.zeros((p, p))

    for i in range(k):
        not_i = list(set(range(k)) - {i})
        other_subjects = subjects[not_i]

        # Thus, the subject-based ISFC was calculated by the Pearson correlation between
        # single subject and the average of all other subjects.
        C = np.dot(subjects[i], other_subjects.mean(axis=0).T) / n

        # Fisher's r-to-z transformation was applied to each correlation coefficient
        # before averaging, in order to increase normality of the distribution of 
        # correlation values,
        C_z_transformed = np.arctanh(C)
        C_hat += C_z_transformed

    # and averaged z values were then inverse transformed (z-to-r) to produce
    # average r-values.
    C_hat = np.tanh(C_hat / k)

    # The final ISFC matrix is given by:
    C_hat = (C_hat + C_hat.T) / 2
    # We imposed this symmetry since we consider the correlation between two
    # brain regions as unidirectional, as in FC