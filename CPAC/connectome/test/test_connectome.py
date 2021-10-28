"""Tests for CPAC.connectome"""
import os
import random
import re
import numpy as np
import pytest
from CPAC.connectome.connectivity_matrix import compute_correlation
from CPAC.pipeline.schema import valid_options


@pytest.mark.parametrize('desc', [True, False])
@pytest.mark.parametrize('method', [
    option for option in valid_options['timeseries']['roi_paths'] if
    option.endswith('Corr') or option.endswith('Embed')
])
def test_compute_correlation(tmp_path, desc, method):
    """Test CPAC.connectome.connectivity_matrix.compute_correlation

    Parameters
    ----------
    desc : bool
        timeseries file already has a description?

    method : str
        Any correlation or embedding timeseries analysis option
    """
    desc = '_desc-existingDescription' if desc else ''

    timeseries_filepath = os.path.join(
        tmp_path,
        f'sub-fake_ses-fake{desc}_timeseries.1D'
    )

    # random numbers of ROIs & TRs
    n_roi = random.randint(1, 200)
    n_tr = random.randint(32, 4800)

    # save fake timeseries
    np.savetxt(
        timeseries_filepath,
        np.random.random_integers(0, 1, n_tr * n_roi).reshape(
            (n_roi, n_tr)
        )
    )

    # create the matrix
    matrix_filepath = compute_correlation(timeseries_filepath, method)
    matrix = np.load(matrix_filepath)

    try:
        # check for correct size output
        assert matrix.shape == (n_roi, n_roi)

        # check for correct filename
        method_string = method.replace(" ", r"\+")
        correct_suffix = re.compile(f'[+-]{method_string}_connectome.npy$')
        assert correct_suffix.search(matrix_filepath) is not None

    finally:
        # remove test files
        map(os.remove, [matrix_filepath, timeseries_filepath])
