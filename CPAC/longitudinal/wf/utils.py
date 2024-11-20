# -*- coding: utf-8 -*-
# Copyright (C) 2020-2024  C-PAC Developers

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

from pathlib import Path
from typing import Optional

from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.utils.interfaces.function import Function


def check_creds_path(creds_path: Optional[str], subject_id: str) -> Optional[str]:
    """Check credentials path."""
    if creds_path and "none" not in creds_path.lower():
        _creds_path = Path(creds_path)
        if _creds_path.exists():
            return str(_creds_path.absolute())
        err_msg = (
            'Credentials path: "%s" for subject "%s" was not '
            "found. Check this path and try again." % (creds_path, subject_id)
        )
        raise FileNotFoundError(err_msg)
    return None


def select_session(
    session: str, output_brains: list[str], warps: list[str]
) -> tuple[Optional[str], Optional[str]]:
    """Select output brain image and warp for given session."""
    brain_path = None
    warp_path = None
    for brain_path in output_brains:
        if f"{session}_" in brain_path:
            break
    for warp_path in warps:
        if f"{session}_" in warp_path:
            break
    return brain_path, warp_path


def select_session_node(unique_id: str, suffix: str = "") -> pe.Node:
    """Create a Node to select a single subject's output image and transform.

    Note
    ----
    FSL is the only currenlty implemented registration tool for longitudinal template
    generation, so it's hardcoded into the name of this node for
    feeding :py:meth:`~CPAC.utils.utils.check_prov_for_regtool`.
    """
    if suffix:
        suffix = f"_{suffix.lstrip('_')}"
    select_sess = pe.Node(
        Function(
            input_names=["session", "output_brains", "warps"],
            output_names=["brain_path", "warp_path"],
            function=select_session,
        ),
        name=f"longitudinal_select_FSL_{unique_id}{suffix}",
    )
    select_sess.inputs.session = unique_id
    return select_sess
