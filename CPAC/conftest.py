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
"""Global fixtures for C-PAC tests."""

from pathlib import Path

from _pytest.tmpdir import TempPathFactory
from git import Repo
import pytest


@pytest.fixture(scope="session")
def bids_examples(tmp_path_factory: TempPathFactory) -> Path:
    """Get the BIDS examples dataset."""
    example_dir = tmp_path_factory.mktemp("bids-examples")
    if not example_dir.exists():
        Repo.clone_from(
            "https://github.com/bids-standard/bids-examples.git", str(example_dir)
        )
    return example_dir
