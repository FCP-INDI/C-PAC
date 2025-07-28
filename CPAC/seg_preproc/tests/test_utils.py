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
"""Tests for segmentation utilities."""

import subprocess


def test_ants_joint_label_fusion_script() -> None:
    """Test antsJointLabelFusion.sh script can run in this environment."""
    try:
        subprocess.run(
            ["antsJointLabelFusion.sh"],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        # There's no explicit 'help' option, but if the script can run,
        # the error message does not contain the string "Error".
        if "Error" in e.stderr.decode():
            raise e
