# Copyright (C) 2021-2024  C-PAC Developers

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
"""Test the CWAS pipeline."""
from logging import basicConfig, INFO
import os
from urllib.error import URLError

import numpy as np
import pandas as pd
import pytest
import nibabel as nib
import nilearn.datasets

from CPAC.cwas.pipeline import create_cwas
from CPAC.utils.monitoring.custom_logging import getLogger

logger = getLogger("CPAC.cwas.tests")
basicConfig(format="%(message)s", level=INFO)


@pytest.mark.parametrize("z_score", [[0], [1], [0, 1], []])
def test_pipeline(z_score):
    """Test the CWAS pipeline with z-score forking options."""
    try:
        # pylint: disable=invalid-name
        cc = nilearn.datasets.fetch_atlas_craddock_2012()
    except URLError:
        logger.info("Could not fetch atlas, skipping test")
        return
    try:
        os.mkdir("/tmp/cwas")
    except:  # noqa: E722
        pass

    abide_data = nilearn.datasets.fetch_abide_pcp(n_subjects=10)
    pheno = pd.DataFrame.from_records(abide_data["phenotypic"])
    images = abide_data["func_preproc"]

    pheno = pheno[["FILE_ID", "AGE_AT_SCAN", "FIQ"]]

    pheno.to_csv("/tmp/cwas/pheno.csv")

    # Sanity ordering check
    assert all(FID in images[i] for i, FID in enumerate(pheno.FILE_ID))
    img = nib.load(cc["scorr_mean"])
    img_data = np.copy(img.get_fdata()[:, :, :, 10])
    img_data[img_data != 2] = 0.0  # noqa: PLR2004
    img = nib.Nifti1Image(img_data, img.affine)
    nib.save(img, "/tmp/cwas/roi.nii.gz")

    workflow = create_cwas("cwas", "/tmp/cwas", "/tmp/cwas")

    subjects = {FID: images[i] for i, FID in enumerate(pheno.FILE_ID)}

    roi = "/tmp/cwas/roi.nii.gz"
    regressor_file = "/tmp/cwas/pheno.csv"
    participant_column = "FILE_ID"
    columns = "AGE_AT_SCAN"
    permutations = 50
    parallel_nodes = 3

    workflow.inputs.inputspec.roi = roi
    workflow.inputs.inputspec.subjects = subjects
    workflow.inputs.inputspec.regressor = regressor_file
    workflow.inputs.inputspec.participant_column = participant_column
    workflow.inputs.inputspec.columns = columns
    workflow.inputs.inputspec.permutations = permutations
    workflow.inputs.inputspec.parallel_nodes = parallel_nodes
    workflow.inputs.inputspec.z_score = z_score

    workflow.run(plugin="Linear")
