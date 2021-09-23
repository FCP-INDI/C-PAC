"""
Custom methods for Nipype pipeline plugins:
* _prerun_check method to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.
"""
import gc
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


# pylint: disable=too-few-public-methods, missing-class-docstring
class CpacNipypeCustomPluginMixin():
    def __init__(self, plugin_args=None):
        if not isinstance(plugin_args, dict):
            plugin_args = {}
        if 'status_callback' not in plugin_args:
            plugin_args['status_callback'] = log_nodes_cb
        self.peak = 0
        self._stats = None
        super().__init__(plugin_args)

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

    def _prerun_check(self, graph):
        """Check if any node exeeds the available resources"""
        tasks_mem_gb = []
        tasks_num_th = []
        overrun_message_mem = None
        overrun_message_th = None
        # estimate of C-PAC + Nipype overhead (GB):
        overhead_memory_estimate = 1
        for node in graph.nodes():
            try:
                node_memory_estimate = node.mem_gb
            except FileNotFoundError:
                # pylint: disable=protected-access
                node_memory_estimate = node._apply_mem_x(UNDEFINED_SIZE)
            node_memory_estimate += overhead_memory_estimate
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
        Customized from https://github.com/nipy/nipype/commit/79e2fdfc38759bc0853e4051b99ba4c37587d65f
        to catch overhead deadlocks
        """  # noqa E501  # pylint: disable=line-too-long
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

        stats = (
            len(self.pending_tasks),
            len(jobids),
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
                len(self.pending_tasks),
                len(jobids),
                free_memory_gb,
                self.memory_gb,
                free_processors,
                self.processors,
                tasks_list_msg,
            )
            self._stats = stats

        if free_memory_gb < self.peak or free_processors == 0:
            logger.debug("No resources available")
            return

        if len(jobids) + len(self.pending_tasks) == 0:
            logger.debug(
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
            if next_job_th > free_processors or next_job_gb > free_memory_gb:
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
