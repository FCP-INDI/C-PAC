# Copyright (C) 2019-2025  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Test nuisance utilities."""

from importlib.resources import as_file, files
import os
import tempfile

import numpy as np
import pytest

from CPAC.nuisance.utils import calc_compcor_components, find_offending_time_points
from CPAC.utils.monitoring.custom_logging import getLogger

logger = getLogger("CPAC.nuisance.tests")

_mocked_outputs = files("CPAC").joinpath("nuisance/tests/motion_statistics")


@pytest.mark.skip(reason="needs refactoring")
def test_find_offending_time_points():
    dl_dir = tempfile.mkdtemp()
    os.chdir(dl_dir)
    with as_file(_mocked_outputs) as mocked_outputs:
        censored = find_offending_time_points(
            str(mocked_outputs / "FD_J.1D"),
            str(mocked_outputs / "FD_P.1D"),
            str(mocked_outputs / "DVARS.1D"),
            2.0,
            2.0,
            "1.5SD",
        )

    censored = np.loadtxt(censored).astype(bool)

    assert set(np.where(np.logical_not(censored))[0].tolist()) == {1, 3, 7}


@pytest.mark.skip(reason="needs local files not included in package")
def test_calc_compcor_components():
    data_filename = "/cc_dev/cpac_working/old_compcor/nuisance_0_0/_scan_test/_selector_CSF-2mmE-M_aC-WM-2mm-DPC5_G-M_M-SDB_P-2_BP-B0.01-T0.1/Functional_2mm_flirt/sub-M10978008_ses-NFB3_task-test_bold_calc_tshift_resample_volreg_calc_maths_flirt.nii.gz"
    mask_filename = "/cc_dev/cpac_working/old_compcor/nuisance_0_0/_scan_test/_selector_CSF-2mmE-M_aC-WM-2mm-DPC5_G-M_M-SDB_P-2_BP-B0.01-T0.1/aCompCor_union_masks/segment_seg_2_maths_flirt_mask.nii.gz"

    compcor_filename = calc_compcor_components(data_filename, 5, mask_filename)
    logger.info("compcor components written to %s", compcor_filename)

