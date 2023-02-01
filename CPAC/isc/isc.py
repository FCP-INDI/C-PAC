import numpy as np

from CPAC.utils import correlation

from .utils import p_from_null, phase_randomize


def isc(D, std=None, collapse_subj=True):

    assert D.ndim == 3

    n_vox, _, n_subj = D.shape
    n_subj_loo = n_subj - 1

    group_sum = np.add.reduce(D, axis=2)

    if collapse_subj:
        ISC = np.zeros(n_vox)
        for loo_subj in range(n_subj):
            loo_subj_ts = D[:, :, loo_subj]
            ISC += correlation(
                loo_subj_ts,
                (group_sum - loo_subj_ts) / n_subj_loo,
                match_rows=True
            )
        ISC /= n_subj

        if std:
            ISC_avg = ISC.mean()
            ISC_std = ISC.std()
            masked = (ISC <= ISC_avg + ISC_std * std) & (ISC >= ISC_avg - ISC_std * std)
        else:
            masked = np.array([True] * n_vox)

    else:
        ISC = np.zeros((n_subj, n_vox))
        for loo_subj in range(n_subj):
            loo_subj_ts = D[:, :, loo_subj]
            ISC[loo_subj] = correlation(
                loo_subj_ts,
                (group_sum - loo_subj_ts) / n_subj_loo,
                match_rows=True
            )
        
        masked = np.array([True] * n_vox)

    return ISC, masked


def isc_significance(ISC, min_null, max_null, two_sided=False):
    p = p_from_null(ISC,
                    max_null=max_null,
                    min_null=min_null,
                    two_sided=two_sided)
    return p


def isc_permutation(permutation, D, masked, collapse_subj=True, random_state=0):

    print("Permutation", permutation)

    min_null = 1
    max_null = -1

    D = D[masked]

    n_vox, _, n_subj = D.shape
    n_subj_loo = n_subj - 1
    D = phase_randomize(D, random_state)

    if collapse_subj:
        ISC_null = np.zeros(n_vox)

    group_sum = np.add.reduce(D, axis=2)

    for loo_subj in range(n_subj):
        loo_subj_ts = D[:, :, loo_subj]
        ISC_subj = \
            correlation(
                loo_subj_ts,
                (group_sum - loo_subj_ts) / n_subj_loo,
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
