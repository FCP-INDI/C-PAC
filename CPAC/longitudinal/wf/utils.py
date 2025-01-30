# Copyright (C) 2020-2025  C-PAC Developers

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
"""Utilities for longitudinal workflows."""


def select_session(
    session: str, output_brains: list[str], warps: list[str]
) -> tuple[str, str]:
    """Select output brain image and warp for given session."""
    try:
        return next(
            iter(brain_path for brain_path in output_brains if session in brain_path)
        ), next(iter(warp_path for warp_path in warps if session in warp_path))
    except StopIteration as stop_iteration:
        brain_paths_found = [
            brain_path for brain_path in output_brains if session in brain_path
        ]
        warps_found = [warp_path for warp_path in warps if session in warp_path]
        msg = ""
        if not brain_paths_found:
            msg += f"{session} not found in {output_brains}.\n"
        if not warps_found:
            msg += f"{session} not found in {warps}.\n"
        raise FileNotFoundError(msg) from stop_iteration
