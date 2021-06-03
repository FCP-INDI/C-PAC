"""

Override Nipype's LegacyMultiproc:
* _prerun_check to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.

"""

import platform
import resource
from nipype.pipeline.plugins.legacymultiproc import (
    LegacyMultiProcPlugin as LegacyMultiProc)
from CPAC.pipeline.nipype_pipeline_engine import UNDEFINED_SIZE
from .base import logger, run_node


def get_peak_usage():
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    proc_peak = self_usage.ru_maxrss

    # getrusage.ru_maxrss in bytes on Macs
    if platform.system() == 'Darwin':
        proc_peak /= 1024.

    return proc_peak / 1024. / 1024.


class LegacyMultiProcPlugin(LegacyMultiProc):
    def __init__(self, plugin_args=None):
        super().__init__(plugin_args)

    def _prerun_check(self, graph):
        """Check if any node exeeds the available resources"""
        tasks_mem_gb = []
        tasks_num_th = []
        overrun_message_mem = None
        overrun_message_th = None
        for node in graph.nodes():
            try:
                node_memory_estimate = node.mem_gb
            except FileNotFoundError:
                node_memory_estimate = node._apply_mem_x(UNDEFINED_SIZE)
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

    def _check_resources(self, running_tasks):
        """
        Make sure there are resources available, accounting for Nipype
        memory usage
        """

        free_memory_gb, free_processors = super()._check_resources(
            running_tasks)

        # Nipype memory usage
        peak = get_peak_usage()
        free_memory_gb -= peak

        return free_memory_gb, free_processors

    def _submit_job(self, node, updatehash=False):
        self._taskid += 1

        # Don't allow streaming outputs
        if getattr(node.interface, "terminal_output", "") == "stream":
            node.interface.terminal_output = "allatonce"

        self._task_obj[self._taskid] = self.pool.apply_async(
            run_node, (node, updatehash, self._taskid),
            callback=self._async_callback
        )

        logger.debug(
            "[LegacyMultiProc] Submitted task %s (taskid=%d).",
            node.fullname,
            self._taskid,
        )
        return self._taskid
