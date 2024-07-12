# Copyright (C) 2024  C-PAC Developers

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
"""Global pytest configuration."""

from pathlib import Path

import pytest


@pytest.fixture
def bids_examples(cache: pytest.Cache) -> Path:
    """Get cached example BIDS directories."""
    bids_dir = cache.mkdir("bids-examples").absolute()
    if not bids_dir.exists():
        from git import Repo

        Repo.clone_from("https://github.com/bids-standard/bids-examples.git", bids_dir)
    return bids_dir
