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
from typing import Any, cast, Optional

from networkx.classes.digraph import DiGraph
from nipype.interfaces.utility import IdentityInterface

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


@Function.sig_imports(
    [
        "from networkx.classes.digraph import DiGraph",
        "from CPAC.pipeline import nipype_pipeline_engine as pe",
        "from CPAC.longitudinal.wf.utils import get_output_from_graph",
    ]
)
def cross_graph_connections(
    wf1: DiGraph | pe.Workflow,
    wf2: pe.Workflow,
    node1: pe.Node | pe.Workflow,
    node2: pe.Node | pe.Workflow,
    output_name: str,
    input_name: str,
) -> None:
    """Make cross-graph connections appropriate to dry-run status.

    Parameters
    ----------
    wf1
        The results of the graph that already ran

    wf2
        The graph that runs second

    node1
        The node from ``wf1``

    node2
        The node from ``wf2``

    output_name
        The output name from ``node1``

    input_name
        The input name from ``node2``
    """
    if isinstance(wf1, pe.Workflow):  # dry run
        if isinstance(node1, pe.Workflow):
            sub_node_name, output_name = output_name.rsplit(".", 1)
            node1 = cast(pe.Node, node1.get_node(sub_node_name))
        wf2.connect(node1, output_name, node2, input_name)
    else:
        setattr(
            node2.inputs, input_name, get_output_from_graph(wf1, node1, output_name)
        )


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
    select_sess.set_input("session", f"{unique_id}_")
    return select_sess


def cross_pool_resources(name: str) -> pe.Node:
    """Return an IdentityInterface for cross-pool resources."""
    return pe.Node(
        IdentityInterface(fields=["from-longitudinal_to-template_mode-image_xfm"]),
        name=name,
    )


def get_output_from_graph(
    graph: DiGraph, node: pe.Node | pe.Workflow, output_name: str
) -> Any:
    """Get an output from a graph that has been run."""
    nodename = str(node.fullname)
    if isinstance(node, pe.Workflow):
        sub_node_name, output_name = output_name.rsplit(".", 1)
        nodename = f"{nodename}.{sub_node_name}"
    try:
        output = getattr(
            next(
                iter(
                    _node
                    for _node in graph
                    if _node.fullname.endswith(nodename)
                    or _node.fullname.endswith(f"{nodename}_")
                )
            ).result.outputs,
            output_name,
        )
    except StopIteration as stop_iteration:
        msg = f"{nodename} not found in completed workflow."
        raise FileNotFoundError(msg) from stop_iteration
    return output
