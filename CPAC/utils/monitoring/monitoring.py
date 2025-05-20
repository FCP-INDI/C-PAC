# Copyright (C) 2018-2025  C-PAC Developers

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
# pylint: disable=too-many-lines,ungrouped-imports,wrong-import-order
"""Monitoring utilities for C-PAC."""

from datetime import datetime, timedelta
import glob
import json
import math
import os
import socketserver
import threading
from typing import Optional

import networkx as nx
from traits.trait_base import Undefined
from nipype.utils.profiler import log_nodes_cb as _nipype_log_nodes_cb

from CPAC.pipeline import nipype_pipeline_engine as pe
from .custom_logging import getLogger


def _safe_none_diff(
    self: "DatetimeWithSafeNone | _NoTime", other: "DatetimeWithSafeNone | _NoTime"
) -> datetime | timedelta:
    """Subtract between a datetime or timedelta or None."""
    if isinstance(self, _NoTime):
        return timedelta(0)
    if isinstance(other, DatetimeWithSafeNone):
        if isinstance(other, _NoTime):
            return timedelta(0)
        return self - other
    if isinstance(other, (datetime, timedelta)):
        return self._dt - other
    msg = f"Cannot subtract {type(other)} from {type(self)}"
    raise NotImplementedError(msg)


class _NoTime:
    """A wrapper for None values that can be used in place of a datetime object."""

    def __bool__(self) -> bool:
        """Return False for _NoTime."""
        return False

    def __int__(self) -> int:
        """Return 0 for _NoTime."""
        return 0

    def __repr__(self) -> str:
        """Return 'NoTime' for _NoTime."""
        return "NoTime"

    def __str__(self) -> str:
        """Return 'NoTime' for _NoTime."""
        return "NoTime"

    def __sub__(self, other: "DatetimeWithSafeNone | _NoTime") -> datetime | timedelta:
        """Subtract between None and a datetime or timedelta or None."""
        return _safe_none_diff(self, other)


NoTime = _NoTime()
"""A singleton None that can be used in place of a datetime object."""


class DatetimeWithSafeNone(datetime, _NoTime):
    """Time class that can be None or a time value."""

    def __new__(cls, dt: Optional[datetime]) -> "DatetimeWithSafeNone | _NoTime":
        """Create a new instance of the class."""
        return (
            NoTime
            if dt is None
            else datetime.__new__(
                cls,
                dt.year,
                dt.month,
                dt.day,
                dt.hour,
                dt.minute,
                dt.second,
                dt.microsecond,
                dt.tzinfo,
            )
        )

    def __bool__(self) -> bool:
        """Return True if not NoTime."""
        return self is not NoTime

    def __sub__(self, other: "DatetimeWithSafeNone | _NoTime") -> datetime | timedelta:
        """Subtract between a datetime or timedelta or None."""
        return _safe_none_diff(self, other)

    def __repr__(self) -> str:
        """Return the string representation of the datetime or NoTime."""
        if self:
            return datetime.__repr__(self)
        return "NoTime"

    def __str__(self) -> str:
        """Return the string representation of the datetime or NoTime."""
        return super().__str__()


def recurse_nodes(workflow, prefix=""):
    """Log initial information from all the nodes."""
    for node in nx.topological_sort(workflow._graph):
        if isinstance(node, pe.Workflow):
            for subnode in recurse_nodes(node, prefix + workflow.name + "."):
                yield subnode
        else:
            yield {
                "id": prefix + workflow.name + "." + node.name,
                "hash": node.inputs.get_hashval()[1],
            }


def log_nodes_initial(workflow):
    logger = getLogger("callback")
    for node in recurse_nodes(workflow):
        logger.debug(json.dumps(node))


def log_nodes_cb(node, status):
    # STATEMENT OF CHANGES:
    #     This function is derived from sources licensed under the Apache-2.0 terms,
    #     and this function has been changed.

    # CHANGES:
    #     * Skips logging MapNodes (since the sub-Nodes are logged)
    #     * Adds hash and input_data_shape to status dict
    #     * Drops duration from status dict
    #     * Sets number of threads used to math.ceil(cpu_percent/100)
    #     * Skips logging not-found Nodes
    #     * Sets `None` default for start and finish
    #     * Uses a MockLogger for the callback logger
    #     * Modified docstring to reflect local changes
    #     * Updated style to match C-PAC codebase

    # ORIGINAL WORK'S ATTRIBUTION NOTICE:
    #    Copyright (c) 2016, the CRN developers team.
    #    All rights reserved.

    #    Redistribution and use in source and binary forms, with or without
    #    modification, are permitted provided that the following conditions are met:

    #    * Redistributions of source code must retain the above copyright notice, this
    #      list of conditions and the following disclaimer.

    #    * Redistributions in binary form must reproduce the above copyright notice,
    #      this list of conditions and the following disclaimer in the documentation
    #      and/or other materials provided with the distribution.

    #   * Neither the name of niworkflows nor the names of its
    #      contributors may be used to endorse or promote products derived from
    #      this software without specific prior written permission.

    #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    #    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    #    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    #    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    #    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    #    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    #    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    #    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    #    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    #    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    #    Licensed under the Apache License, Version 2.0 (the "License");
    #    you may not use this file except in compliance with the License.
    #    You may obtain a copy of the License at

    #        http://www.apache.org/licenses/LICENSE-2.0

    #    Unless required by applicable law or agreed to in writing, software
    #    distributed under the License is distributed on an "AS IS" BASIS,
    #    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    #    See the License for the specific language governing permissions and
    #    limitations under the License.

    # Modifications copyright (C) 2019 - 2024  C-PAC Developers
    if status != "end":
        return

    from nipype.pipeline.engine import nodes

    logger = getLogger("callback")

    if isinstance(node, nodes.MapNode):
        return

    try:
        runtime = node.result.runtime
    except FileNotFoundError:
        runtime = {}
    runtime_threads = getattr(runtime, "cpu_percent", "N/A")
    if runtime_threads != "N/A":
        runtime_threads = math.ceil(runtime_threads / 100)

    status_dict = {
        "id": str(node),
        "hash": node.inputs.get_hashval()[1],
        "start": DatetimeWithSafeNone(getattr(runtime, "startTime", None)),
        "finish": DatetimeWithSafeNone(getattr(runtime, "endTime", None)),
        "runtime_threads": runtime_threads,
        "runtime_memory_gb": getattr(runtime, "mem_peak_gb", "N/A"),
        "estimated_memory_gb": node.mem_gb,
        "num_threads": node.n_procs,
    }

    if hasattr(node, "input_data_shape") and node.input_data_shape is not Undefined:
        status_dict["input_data_shape"] = node.input_data_shape

    if any(
        not isinstance(status_dict[label], datetime) for label in ["start", "finish"]
    ):
        status_dict["error"] = True

    logger.debug(json.dumps(status_dict))


log_nodes_cb.__doc__ = f"""{_nipype_log_nodes_cb.__doc__}

Modified from https://github.com/nipy/nipype/blob/5ab2fa0/nipype/utils/profiler.py#L112-L156
"""


class LoggingRequestHandler(socketserver.BaseRequestHandler):
    def handle(self):
        tree = {}

        logs = glob.glob(
            os.path.join(
                self.server.logging_dir, "pipeline_" + self.server.pipeline_name, "*"
            )
        )

        for log in logs:
            subject = log.split("/")[-1]
            tree[subject] = {}

            callback_file = os.path.join(log, "callback.log")

            if not os.path.exists(callback_file):
                continue

            with open(callback_file, "rb") as lf:
                for l in lf.readlines():  # noqa: E741
                    l = l.strip()  # noqa: E741,PLW2901
                    try:
                        node = json.loads(l)
                        if node["id"] not in tree[subject]:
                            tree[subject][node["id"]] = {"hash": node["hash"]}
                            if "start" in node and "finish" in node:
                                tree[subject][node["id"]]["start"] = node["start"]
                                tree[subject][node["id"]]["finish"] = node["finish"]

                        elif "start" in node and "finish" in node:
                            if tree[subject][node["id"]]["hash"] == node["hash"]:
                                tree[subject][node["id"]]["cached"] = {
                                    "start": node["start"],
                                    "finish": node["finish"],
                                }

                            # pipeline was changed, and we have a new hash
                            else:
                                tree[subject][node["id"]]["start"] = node["start"]
                                tree[subject][node["id"]]["finish"] = node["finish"]

                    except:
                        break

                tree = {s: t for s, t in tree.items() if t}

        headers = "HTTP/1.1 200 OK\nConnection: close\n\n"
        self.request.sendall(headers + json.dumps(tree) + "\n")


class LoggingHTTPServer(socketserver.ThreadingTCPServer, object):
    def __init__(
        self,
        pipeline_name,
        logging_dir="",
        host="",
        port=8080,
        request=LoggingRequestHandler,
    ):
        super(LoggingHTTPServer, self).__init__((host, port), request)

        if not logging_dir:
            logging_dir = os.getcwd()

        self.logging_dir = logging_dir
        self.pipeline_name = pipeline_name


def monitor_server(pipeline_name, logging_dir, host="0.0.0.0", port=8080):
    httpd = LoggingHTTPServer(
        pipeline_name, logging_dir, host, port, LoggingRequestHandler
    )

    server_thread = threading.Thread(target=httpd.serve_forever)
    server_thread.isDaemon = True
    server_thread.start()

    return server_thread
