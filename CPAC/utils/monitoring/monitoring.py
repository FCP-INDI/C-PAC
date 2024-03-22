import glob
import json
import math
import os
import socketserver
import threading

import networkx as nx
from traits.trait_base import Undefined
from nipype.utils.profiler import log_nodes_cb as _nipype_log_nodes_cb

from CPAC.pipeline import nipype_pipeline_engine as pe
from .custom_logging import getLogger


# Log initial information from all the nodes
def recurse_nodes(workflow, prefix=""):
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
        "start": getattr(runtime, "startTime", None),
        "finish": getattr(runtime, "endTime", None),
        "runtime_threads": runtime_threads,
        "runtime_memory_gb": getattr(runtime, "mem_peak_gb", "N/A"),
        "estimated_memory_gb": node.mem_gb,
        "num_threads": node.n_procs,
    }

    if hasattr(node, "input_data_shape") and node.input_data_shape is not Undefined:
        status_dict["input_data_shape"] = node.input_data_shape

    if status_dict["start"] is None or status_dict["finish"] is None:
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
                    l = l.strip()  # noqa: E741
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
