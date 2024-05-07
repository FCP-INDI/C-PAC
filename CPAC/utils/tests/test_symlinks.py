# Copyright (C) 2019-2024  C-PAC Developers

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
from logging import basicConfig, INFO
import os
import tempfile

import pkg_resources as p

from CPAC.utils.monitoring.custom_logging import getLogger
from CPAC.utils.symlinks import create_symlinks

logger = getLogger("CPAC.utils.tests")
basicConfig(format="%(message)s", level=INFO)

mocked_outputs = p.resource_filename(
    "CPAC", os.path.join("utils", "tests", "test_symlinks-outputs.txt")
)


def test_symlinks():
    temp_dir = tempfile.mkdtemp(suffix="test_symlinks")

    paths = []
    with open(mocked_outputs, "r") as f:
        for _path in f.readlines():
            path = _path
            path = path.strip()
            if path:
                paths += [path]

    create_symlinks(
        temp_dir, "sym_links", "pipeline_benchmark-FNIRT", "1019436_1", paths
    )

    logger.info("Links created at %s", temp_dir)

    # TODO test the generated links

    # Normal resource case
    # Several resources within same key case
    # QC case
