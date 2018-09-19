import six
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

    F = fft(D, axis=2)
    if D.shape[2] % 2 == 0:
        pos_freq = np.arange(1, D.shape[2] // 2)
        neg_freq = np.arange(D.shape[2] - 1, D.shape[2] // 2, -1)
    else:
        pos_freq = np.arange(1, (D.shape[2] - 1) // 2 + 1)
        neg_freq = np.arange(D.shape[2] - 1, (D.shape[2] - 1) // 2, -1)

    shift = random_state.rand(D.shape[0], D.shape[1], len(pos_freq)) * 2 * np.pi

    F[:, :, pos_freq] *= np.exp(1j * shift)
    F[:, :, neg_freq] *= np.exp(-1j * shift)

    return np.real(ifft(F, axis=2))


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


def isc(D, collapse_subj=True):

    if isinstance(D, six.string_types):
        D = np.load(D)

    assert D.ndim == 3

    n_subj = D.shape[0]
    n_vox = D.shape[1]

    if collapse_subj:
        ISC = np.zeros(n_vox)
        for loo_subj in range(n_subj):
            ISC += correlation(
                D[loo_subj],
                np.mean(D[np.arange(n_subj) != loo_subj], axis=0),
                match_rows=True
            )
        ISC /= n_subj

    else:
        ISC = np.zeros((n_subj, n_vox))
        for loo_subj in range(n_subj):
            ISC[loo_subj] = correlation(
                D[loo_subj],
                np.mean(D[np.arange(n_subj) != loo_subj], axis=0),
                match_rows=True
            )

    return ISC


def isc_significance(ISC, min_null, max_null, two_sided=False):
    p = p_from_null(ISC,
                    max_null=max_null,
                    min_null=min_null,
                    two_sided=two_sided)
    return p


def isc_permutation(permutation, D, collapse_subj=True, random_state=0):

    if isinstance(D, six.string_types):
        D = np.load(D)

    min_null = 1
    max_null = -1

    n_subj, n_vox, _ = D.shape
    D = phase_randomize(D, random_state)

    if collapse_subj:
        ISC_null = np.zeros(n_vox)

    for loo_subj in range(n_subj):

        ISC_subj = \
            correlation(
                D[loo_subj],
                np.mean(D[np.arange(n_subj) != loo_subj], axis=0),
                match_rows=True
            )

        if collapse_subj:
            ISC_null += ISC_subj
        else:
            max_null = max(np.max(ISC_subj), max_null)
            min_null = min(np.min(ISC_subj), min_null)

    if collapse_subj:
        ISC_null /= n_subj
        max_null = np.max(ISC_null)
        min_null = np.min(ISC_null)

    return permutation, min_null, max_null


def isfc(D, collapse_subj=True):

    if isinstance(D, six.string_types):
        D = np.load(D)

    assert D.ndim == 3

    n_subj = D.shape[0]
    n_vox = D.shape[1]

    if collapse_subj:
        ISFC = np.zeros((n_vox, n_vox))
        for loo_subj in range(n_subj):
            ISFC += correlation(
                D[loo_subj],
                np.mean(D[np.arange(n_subj) != loo_subj], axis=0),
                symmetric=True
            )
        ISFC /= n_subj

    else:
        ISFC = np.zeros((n_subj, n_vox, n_vox))
        for loo_subj in range(n_subj):
            ISFC[loo_subj] = correlation(
                D[loo_subj],
                np.mean(D[np.arange(n_subj) != loo_subj], axis=0),
                symmetric=True
            )

    return ISFC


def isfc_significance(ISFC, min_null, max_null, two_sided=False):
    p = p_from_null(ISFC,
                    max_null=max_null,
                    min_null=min_null,
                    two_sided=two_sided)
    return p


def isfc_permutation(permutation, D, collapse_subj=True, random_state=0):

    if isinstance(D, six.string_types):
        D = np.load(D)

    min_null = 1
    max_null = -1

    n_subj, n_vox, _ = D.shape
    D = phase_randomize(D, random_state)

    if collapse_subj:
        ISFC_null = np.zeros((n_vox, n_vox))

    for loo_subj in range(n_subj):
        ISFC_subj = \
            correlation(
                D[loo_subj],
                np.mean(D[np.arange(n_subj) != loo_subj], axis=0),
                symmetric=True
            )

        if collapse_subj:
            ISFC_null += ISFC_subj
        else:
            max_null = max(np.max(ISFC_subj), max_null)
            min_null = min(np.min(ISFC_subj), min_null)
    
    if collapse_subj:
        ISFC_null /= n_subj
        max_null = np.max(ISFC_null)
        min_null = np.min(ISFC_null)

    return permutation, min_null, max_null