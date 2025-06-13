# Copyright (C) 2025  C-PAC Developers

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
"""FSL registration connector tests."""

import pytest

from CPAC.registration.registration import FSL_registration_connector
from CPAC.utils.configuration import Configuration
from CPAC.utils.tests.test_utils import check_expected_keys


@pytest.mark.parametrize("sink_native_transforms", [True, False])
def test_fsl_registration_connector(sink_native_transforms: bool) -> None:
    """Test FSL registration connector with various configurations."""
    wf_name = "test_fsl_registration_connector"
    cfg = Configuration(
        {"registration_workflows": {"sink_native_transforms": sink_native_transforms}}
    )
    _, outputs = FSL_registration_connector(wf_name, cfg, opt="FSL")
    expected_keys = {
        "from-T1w_to-template_mode-image_desc-flirt_xfm",
        "from-template_to-T1w_mode-image_desc-flirt_xfm",
    }
    check_expected_keys(sink_native_transforms, outputs, expected_keys)
