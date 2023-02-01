import numpy as np
from CPAC.utils import correlation

from .utils import p_from_null, phase_randomize


def isfc(D, std=None, collapse_subj=True):

    assert D.ndim == 3

    n_vox, _, n_subj = D.shape
    n_subj_loo = n_subj - 1

    group_sum = np.add.reduce(D, axis=2)
    masked = None

    if collapse_subj:
        ISFC = np.zeros((n_vox, n_vox))
        for loo_subj in range(n_subj):
            loo_subj_ts = D[:, :, loo_subj]
            ISFC += correlation(
                loo_subj_ts,
                (group_sum - loo_subj_ts) / n_subj_loo,
                symmetric=True
            )
        ISFC /= n_subj

        if std:
            ISFC_avg = ISFC.mean()
            ISFC_std = ISFC.std()
            masked = (ISFC <= ISFC_avg + ISFC_std) | (ISFC >= ISFC_avg - ISFC_std)

    else:
        ISFC = np.zeros((n_vox, n_vox, n_subj))
        for loo_subj in range(n_subj):
            loo_subj_ts = D[:, :, loo_subj]
            ISFC[:, :, loo_subj] = correlation(
                loo_subj_ts,
                (group_sum - loo_subj_ts) / n_subj_loo,
                symmetric=True
            )

    if masked is not None:
        masked = np.all(masked, axis=1)
    else:
        masked = np.array([True] * n_vox)

    return ISFC, masked


def isfc_significance(ISFC, min_null, max_null, two_sided=False):
    p = p_from_null(ISFC,
                    max_null=max_null,
                    min_null=min_null,
                    two_sided=two_sided)
    return p


def isfc_permutation(permutation, D, masked, collapse_subj=True, random_state=0):

    print("Permutation", permutation)

    min_null = 1
    max_null = -1

    D = D[masked]

    n_vox, _, n_subj = D.shape
    n_subj_loo = n_subj - 1

    D = phase_randomize(D, random_state)

    if collapse_subj:
        ISFC_null = np.zeros((n_vox, n_vox))

    group_sum = np.add.reduce(D, axis=2)

    for loo_subj in range(n_subj):
        loo_subj_ts = D[:, :, loo_subj]
        ISFC_subj = \
            correlation(
                loo_subj_ts,
                (group_sum - loo_subj_ts) / n_subj_loo,
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
