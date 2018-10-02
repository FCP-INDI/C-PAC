import numpy as np
from sklearn.utils import check_random_state
from scipy.fftpack import fft, ifft
from scipy.stats import pearsonr


def ecdf(x):
    xp = np.sort(x)
    yp = np.arange(len(xp) + 1) / len(xp)
    return lambda q: yp[np.searchsorted(xp, q, side="right")]


def zscore(data, axis):
    data = data.copy()
    data -= data.mean(axis=axis, keepdims=True)
    data /= data.std(axis=axis, keepdims=True)
    np.copyto(data, 0.0, where=np.isnan(data))
    return data


def correlation(matrix1, matrix2,
                match_rows=False, z_scored=False, symmetric=False):
    d1 = matrix1.shape[-1]
    d2 = matrix2.shape[-1]

    assert d1 == d2
    assert matrix1.ndim <= 2
    assert matrix2.ndim <= 2
    if match_rows:
        assert matrix1.shape == matrix2.shape

    var = np.sqrt(d1 * d2)
    
    if not z_scored:
        matrix1 = zscore(matrix1, matrix1.ndim - 1)
        matrix2 = zscore(matrix2, matrix2.ndim - 1)

    if match_rows:
        return np.einsum('...i,...i', matrix1, matrix2) / var
    
    if matrix1.ndim >= matrix2.ndim:
        r = np.dot(matrix1, matrix2.T) / var
    else:
        r = np.dot(matrix2, matrix1.T) / var

    if symmetric:
        return (r + r.T) / 2
    
    return r


def phase_randomize(D, random_state=0):
    random_state = check_random_state(random_state)

    F = fft(D, axis=1)
    if D.shape[1] % 2 == 0:
        pos_freq = np.arange(1, D.shape[1] // 2)
        neg_freq = np.arange(D.shape[1] - 1, D.shape[1] // 2, -1)
    else:
        pos_freq = np.arange(1, (D.shape[1] - 1) // 2 + 1)
        neg_freq = np.arange(D.shape[1] - 1, (D.shape[1] - 1) // 2, -1)

    shift = random_state.rand(D.shape[0], len(pos_freq),
                              D.shape[2]) * 2 * np.pi

    F[:, pos_freq, :] *= np.exp(1j * shift)
    F[:, neg_freq, :] *= np.exp(-1j * shift)

    return np.real(ifft(F, axis=1))


def p_from_null(X, 
                max_null, min_null,
                two_sided=False):
    max_null_ecdf = ecdf(max_null)
    if two_sided:
        min_null_ecdf = ecdf(min_null)
        p = 2 * np.minimum(1 - max_null_ecdf(X), min_null_ecdf(X))
        p = np.minimum(p, 1)
    else:
        p = 1 - max_null_ecdf(X)
    return p
