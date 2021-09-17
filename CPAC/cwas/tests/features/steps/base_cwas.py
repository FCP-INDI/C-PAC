import numpy as np


def custom_corrcoef(X, Y=None):
    """Each of the columns in X will be correlated with each of the
    columns in Y. Each column represents a variable, with the rows
    containing the observations."""
    if Y is None:
        Y = X

    if X.shape[0] != Y.shape[0]:
        raise Exception("X and Y must have the same number of rows.")

    X = X.astype(float)
    Y = Y.astype(float)

    X -= X.mean(axis=0)[np.newaxis, ...]
    Y -= Y.mean(axis=0)

    xx = np.sum(X**2, axis=0)
    yy = np.sum(Y**2, axis=0)

    r = np.dot(X.T, Y)/np.sqrt(np.multiply.outer(xx, yy))

    return r
