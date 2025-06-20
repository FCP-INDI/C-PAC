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
"""Open Science Framework testing utilities."""

import os
from pathlib import Path

import requests

FILES = {"residuals.nii.gz": "kyqad", "regressors.1D": "xzuyf"}


def download_file(file: str, destination: Path | str) -> Path:
    """Download a file from the Open Science Framework."""
    url = f"https://osf.io/download/{FILES[file]}"
    response = requests.get(
        url,
        headers={"Authorization": f"Bearer {os.getenv('OSF_DATA')}"},
        allow_redirects=True,
    )
    if not isinstance(destination, Path):
        destination = Path(destination)
    destination = destination / file if destination.is_dir() else destination
    if destination.exists():
        msg = f"File {destination} already exists. Please remove it before downloading."
        raise FileExistsError(msg)
    response.raise_for_status()
    with open(destination, "wb") as f:
        f.write(response.content)
    return destination
