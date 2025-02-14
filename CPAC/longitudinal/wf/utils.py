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
from CPAC.pipeline.utils import get_edges_with_node
from CPAC.utils.interfaces.function import Function

_TEMPLATE_PATTERN = r"((sym)?template|longitudinal)"
LONGITUDINAL_TEMPLATE_PATTERN = f"from-{_TEMPLATE_PATTERN}_to-{_TEMPLATE_PATTERN}"


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
    wf: pe.Workflow,
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
    wf
        The graph that already ran

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
            node2.inputs, input_name, get_output_from_graph(wf, wf1, node1, output_name)
        )


def select_session(session: str, outputs: list[str]) -> str:
    """Select output brain image and warp for given session."""
    try:
        return next(iter(path for path in outputs if session in path))
    except StopIteration as stop_iteration:
        msg = f"{session} not found in {outputs}.\n"
        raise FileExistsError(msg) from stop_iteration


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
            input_names=["session", "outputs"],
            output_names=["path"],
            function=select_session,
        ),
        name=f"longitudinalSelect{suffix.title()}_{unique_id}",
    )
    select_sess.set_input("session", f"{unique_id}_")
    return select_sess


def cross_graph_identity(
    wf: pe.Workflow, graph: DiGraph, node: pe.Node | pe.Workflow, output_name: str
) -> pe.Node:
    """Sever connection to ``node_1``'s workflow while maintaining value."""
    identity_interface = pe.Node(
        IdentityInterface([output_name]),
        name="_".join(
            [str(getattr(node, "fullname", getattr(node, "name", ""))), output_name]
        ).replace(".", "_"),
    )
    identity_interface.set_input(
        output_name, get_output_from_graph(wf, graph, node, output_name)
    )
    return identity_interface


def get_output_from_graph(
    wf: pe.Workflow, graph: DiGraph, node: pe.Node | pe.Workflow, output_name: str
) -> Any:
    """Get an output from a graph that has been run."""
    nodename = str(node.fullname)
    if isinstance(node, pe.Workflow):
        sub_node_name, output_name = output_name.rsplit(".", 1)
        nodename = f"{nodename}.{sub_node_name}"
        edges = get_edges_with_node(node, output_name)
        for edge in reversed(edges):
            try:
                return get_output_from_graph(
                    wf,
                    graph,
                    edge[0],
                    next(
                        iter(
                            connection
                            for connection in edge[2]["connect"]
                            if connection[1] == output_name
                        )
                    )[0],
                )
            except StopIteration:
                continue
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
