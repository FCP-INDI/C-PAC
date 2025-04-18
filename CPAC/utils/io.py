# Copyright (C) 2012-2024  C-PAC Developers

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
"""Utilities for inputs and outputs."""

from pathlib import Path

from yaml import safe_load, YAMLError


def load_yaml(
    path: Path | str, desc: str = "YAML file", encoding="utf8"
) -> dict | list | str:
    """Try to load a YAML file to a Python object."""
    path = Path(path).absolute()
    try:
        with path.open("r", encoding=encoding) as _yaml:
            result = safe_load(_yaml)
    except FileNotFoundError as error:
        raise error
    except Exception as error:
        msg = f"{desc} is not in proper YAML format. Please check {path}"
        raise YAMLError(msg) from error
    return result
