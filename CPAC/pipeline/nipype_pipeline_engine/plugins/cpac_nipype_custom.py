"""
Custom methods for Nipype pipeline plugins:
* _prerun_check method to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.

STATEMENT OF CHANGES:
    This file is derived from sources licensed under the Apache-2.0 terms,
    and this file has been changed.

CHANGES:
    * Supports just-in-time dynamic memory allocation
    * Supports overriding memory estimates via a log file and a buffer

ORIGINAL WORK'S ATTRIBUTION NOTICE:
    Copyright (c) 2009-2016, Nipype developers

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    Prior to release 0.12, Nipype was licensed under a BSD license.

Modifications Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.
"""
import gc
import json
import platform
import resource
import sys
from copy import deepcopy
from logging import INFO
from textwrap import indent
from traceback import format_exception
from nipype.pipeline.plugins.multiproc import logger
from numpy import flatnonzero
from CPAC.pipeline.nipype_pipeline_engine import MapNode, UNDEFINED_SIZE
from CPAC.utils.monitoring import log_nodes_cb


OVERHEAD_MEMORY_ESTIMATE: float = 1  # estimate of C-PAC + Nipype overhead (GB)


def get_peak_usage():
    """Function to return peak usage in GB.

    Parameters
    ----------
    None

    Returns
    -------
    float
    """
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    proc_peak = self_usage.ru_maxrss

    # getrusage.ru_maxrss in bytes on Macs
    if platform.system() == 'Darwin':
        proc_peak /= 1024.

    return proc_peak / 1024. / 1024.


def parse_previously_observed_mem_gb(callback_log_path):
    """Function to parse the previously observed memory usage.

    Parameters
    ----------
    callback_log_path : str
        Path to the callback.log file.

    Returns
    -------
    dict
        Dictionary of per-node memory usage.
    """
    with open(callback_log_path, 'r') as cbl:
        return {line['id'].split('.', 1)[-1]: line['runtime_memory_gb'] for
                line in [json.loads(line) for line in cbl.readlines()] if
                'runtime_memory_gb' in line}


# pylint: disable=too-few-public-methods, missing-class-docstring
class CpacNipypeCustomPluginMixin():
    def __init__(self, plugin_args=None):
        if not isinstance(plugin_args, dict):
            plugin_args = {}
        if 'status_callback' not in plugin_args:
            plugin_args['status_callback'] = log_nodes_cb
        if 'runtime' in plugin_args:
            self.runtime = {node_key: observation * (
                1 + plugin_args['runtime']['buffer'] / 100
            ) for node_key, observation in parse_previously_observed_mem_gb(
                plugin_args['runtime']['usage']).items()}
        super().__init__(plugin_args=plugin_args)
        self.peak = 0
        self._stats = None

    def _check_resources_(self, running_tasks):
        """
        Make sure there are resources available
        """
        free_memory_gb = self.memory_gb
        free_processors = self.processors
        for _, jobid in running_tasks:
            free_memory_gb -= min(self.procs[jobid].mem_gb, free_memory_gb)
            free_processors -= min(self.procs[jobid].n_procs, free_processors)

        return free_memory_gb, free_processors

    def _check_resources(self, running_tasks):
        """Make sure there are resources available, accounting for
        Nipype memory usage"""
        free_memory_gb, free_processors = self._check_resources_(running_tasks)

        # Nipype memory usage
        self.peak = get_peak_usage()
        free_memory_gb -= self.peak

        return free_memory_gb, free_processors

    def _clean_exception(self, jobid, graph):
        traceback = format_exception(*sys.exc_info())
        self._clean_queue(
            jobid, graph, result={"result": None,
                                  "traceback": traceback}
        )

    def _override_memory_estimate(self, node):
        """
        Override node memory estimate with provided runtime memory
        usage, buffered

        Parameters
        ----------
        node : nipype.pipeline.engine.nodes.Node

        Returns
        -------
        None
        """
        if hasattr(node, 'list_node_names'):
            for node_id in node.list_node_names():
                # drop top-level node name
                node_id = node_id.split('.', 1)[-1]
        else:
            node_id = node.fullname.split('.', 1)[-1]
        if self._match_for_overrides(node, node_id):
            return
        while '.' in node_id:  # iterate through levels of specificity
            node_id = node_id.rsplit('.', 1)[0]
            if self._match_for_overrides(node, node_id):
                return

    def _match_for_overrides(self, node, node_id):
        """Match node memory estimate with provided runtime memory usage key

        Parameters
        ----------
        node : nipype.pipeline.engine.nodes.Node

        node_id : str

        Returns
        -------
        bool : updated?
        """
        if node_id in self.runtime:
            node.override_mem_gb(self.runtime[node_id])
            return True
        partial_matches = [nid for nid in self.runtime if node_id in nid]
        if any(partial_matches):
            node.override_mem_gb(max(
                self.runtime[partial_match] for
                partial_match in partial_matches))
            return True
        return False

    def _prerun_check(self, graph):
        """Check if any node exeeds the available resources"""
        tasks_mem_gb = []
        tasks_num_th = []
        overrun_message_mem = None
        overrun_message_th = None
        for node in graph.nodes():
            if hasattr(self, 'runtime'):
                self._override_memory_estimate(node)
            elif hasattr(node, "throttle"):
                # for a throttled node without an observation run,
                # assume all available memory will be needed
                node._mem_gb = self.memory_gb - OVERHEAD_MEMORY_ESTIMATE
            try:
                node_memory_estimate = node.mem_gb
            except FileNotFoundError:
                # pylint: disable=protected-access
                node_memory_estimate = node._apply_mem_x(UNDEFINED_SIZE)
            node_memory_estimate += OVERHEAD_MEMORY_ESTIMATE
            if node_memory_estimate > self.memory_gb:
                tasks_mem_gb.append((node.name, node_memory_estimate))
            if node.n_procs > self.processors:
                tasks_num_th.append((node.name, node.n_procs))

        if tasks_mem_gb:
            overrun_message_mem = '\n'.join([
                f'\t{overrun[0]}: {overrun[1]} GB' for overrun in tasks_mem_gb
            ])
            logger.warning(
                "The following nodes are estimated to exceed the total amount "
                f"of memory available (%0.2fGB): \n{overrun_message_mem}",
                self.memory_gb,
            )

        if tasks_num_th:
            overrun_message_th = '\n'.join([
                f'\t{overrun[0]}: {overrun[1]} threads' for overrun in
                tasks_num_th])
            logger.warning(
                "Some nodes demand for more threads than available (%d): "
                f"\n{overrun_message_th}",
                self.processors,
            )

        if self.raise_insufficient and (tasks_mem_gb or tasks_num_th):
            raise RuntimeError("\n".join([msg for msg in [
                "Insufficient resources available for job:",
                overrun_message_mem, overrun_message_th
            ] if msg is not None]))

    def _send_procs_to_workers(self, updatehash=False, graph=None):
        """
        Sends jobs to workers when system resources are available.
        Customized from https://github.com/nipy/nipype/blob/79e2fdfc/nipype/pipeline/plugins/legacymultiproc.py#L311-L462
        to catch overhead deadlocks
        """  # noqa: E501  # pylint: disable=line-too-long
        # pylint: disable=too-many-branches, too-many-statements
        # Check to see if a job is available (jobs with all dependencies run)
        # See https://github.com/nipy/nipype/pull/2200#discussion_r141605722
        # See also https://github.com/nipy/nipype/issues/2372
        jobids = flatnonzero(
            ~self.proc_done & (self.depidx.sum(axis=0) == 0).__array__()
        )

        # Check available resources by summing all threads and memory used
        free_memory_gb, free_processors = self._check_resources(
            self.pending_tasks)

        num_pending = len(self.pending_tasks)
        num_ready = len(jobids)
        stats = (
            num_pending,
            num_ready,
            free_memory_gb,
            self.memory_gb,
            free_processors,
            self.processors,
        )
        if self._stats != stats:
            tasks_list_msg = ""

            if logger.level <= INFO:
                running_tasks = [
                    "  * %s" % self.procs[jobid].fullname
                    for _, jobid in self.pending_tasks
                ]
                if running_tasks:
                    tasks_list_msg = "\nCurrently running:\n"
                    tasks_list_msg += "\n".join(running_tasks)
                    tasks_list_msg = indent(tasks_list_msg, " " * 21)
            logger.info(
                "[%s] Running %d tasks, and %d jobs ready. Free "
                "memory (GB): %0.2f/%0.2f, Free processors: %d/%d.%s",
                type(self).__name__[:-len('Plugin')],
                num_pending,
                num_ready,
                free_memory_gb,
                self.memory_gb,
                free_processors,
                self.processors,
                tasks_list_msg,
            )
            self._stats = stats

        if self.raise_insufficient:
            if free_memory_gb < self.peak or free_processors == 0:
                logger.info("No resources available. Potential deadlock")
                return

            if num_ready + num_pending == 0:
                logger.info(
                    "No tasks are being run, and no jobs can "
                    "be submitted to the queue. Potential deadlock"
                )
                return

        jobids = self._sort_jobs(jobids,
                                 scheduler=self.plugin_args.get("scheduler"))

        # Run garbage collector before potentially submitting jobs
        gc.collect()

        # Submit jobs
        for jobid in jobids:
            force_allocate_job = False
            # First expand mapnodes
            if isinstance(self.procs[jobid], MapNode):
                try:
                    num_subnodes = self.procs[jobid].num_subnodes()
                except Exception:  # pylint: disable=broad-except
                    self._clean_exception(jobid, graph)
                    self.proc_pending[jobid] = False
                    continue
                if num_subnodes > 1:
                    submit = self._submit_mapnode(jobid)
                    if not submit:
                        continue

            # Check requirements of this job
            next_job_gb = min(self.procs[jobid].mem_gb, self.memory_gb)
            next_job_th = min(self.procs[jobid].n_procs, self.processors)

            # If node does not fit, skip at this moment
            if not self.raise_insufficient and (
                num_pending == 0 and num_ready > 0
            ):
                force_allocate_job = True
                free_processors -= 1
            if not force_allocate_job and (
                next_job_th > free_processors or next_job_gb > free_memory_gb
            ):
                logger.debug(
                    "Cannot allocate job %s ID=%d (%0.2fGB, %d threads).",
                    self.procs[jobid].fullname,
                    jobid,
                    next_job_gb,
                    next_job_th,
                )
                continue

            free_memory_gb -= next_job_gb
            free_processors -= next_job_th
            logger.debug(
                "Allocating %s ID=%d (%0.2fGB, %d threads). Free: "
                "%0.2fGB, %d threads.",
                self.procs[jobid].fullname,
                jobid,
                next_job_gb,
                next_job_th,
                free_memory_gb,
                free_processors,
            )

            # change job status in appropriate queues
            self.proc_done[jobid] = True
            self.proc_pending[jobid] = True

            # If cached and up-to-date just retrieve it, don't run
            if self._local_hash_check(jobid, graph):
                continue

            # updatehash and run_without_submitting are also run locally
            if updatehash or self.procs[jobid].run_without_submitting:
                logger.debug("Running node %s on master thread",
                             self.procs[jobid])
                try:
                    self.procs[jobid].run(updatehash=updatehash)
                except Exception:  # pylint: disable=broad-except
                    self._clean_exception(jobid, graph)

                # Release resources
                self._task_finished_cb(jobid)
                self._remove_node_dirs()
                free_memory_gb += next_job_gb
                free_processors += next_job_th
                # Display stats next loop
                self._stats = None

                # Clean up any debris from running node in main process
                gc.collect()
                continue

            # Task should be submitted to workers
            # Send job to task manager and add to pending tasks
            if self._status_callback:
                self._status_callback(self.procs[jobid], "start")
            tid = self._submit_job(deepcopy(self.procs[jobid]),
                                   updatehash=updatehash)
            if tid is None:
                self.proc_done[jobid] = False
                self.proc_pending[jobid] = False
            else:
                self.pending_tasks.insert(0, (tid, jobid))
            # Display stats next loop
            self._stats = None
