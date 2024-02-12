import numpy as np
from scipy.fftpack import fft, ifft

from CPAC.utils import check_random_state


def ecdf(x):
    xp = np.sort(x)
    yp = np.arange(len(xp) + 1) / len(xp)
    return lambda q: yp[np.searchsorted(xp, q, side="right")]


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
