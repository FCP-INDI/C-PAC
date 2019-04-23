import os
import tempfile
import pkg_resources as p
import numpy as np
from CPAC.nuisance.utils import find_offending_time_points


mocked_outputs = \
    p.resource_filename(
        "CPAC",
        os.path.join(
            'nuisance',
            'tests',
            'motion_statistics'
        )
    )


def test_find_offending_time_points():

    dl_dir = tempfile.mkdtemp()
    os.chdir(dl_dir)

    censored = find_offending_time_points(
        os.path.join(mocked_outputs, 'FD_J.1D'),
        os.path.join(mocked_outputs, 'FD_P.1D'),
        os.path.join(mocked_outputs, 'DVARS.1D'),
        2.0,
        2.0,
        '1.5SD'
    )

    censored = np.loadtxt(censored).astype(bool)

    assert set(np.where(np.logical_not(censored))[0].tolist()) == set([1, 3, 7])