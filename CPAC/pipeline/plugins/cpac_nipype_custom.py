"""
Custom methods for Nipype pipeline plugins:
* _prerun_check method to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.
"""
import platform
import resource
from nipype.pipeline.plugins.multiproc import logger


def get_peak_usage():
    self_usage = resource.getrusage(resource.RUSAGE_SELF)
    proc_peak = self_usage.ru_maxrss

    # getrusage.ru_maxrss in bytes on Macs
    if platform.system() == 'Darwin':
        proc_peak /= 1024.

    return proc_peak / 1024. / 1024.


class CpacNipypeCustomPlugin():
    def _check_resources(self, running_tasks):
        """Make sure there are resources available, accounting for
        Nipype memory usage"""
        free_memory_gb, free_processors = super(
        )._check_resources(running_tasks)

        # Nipype memory usage
        peak = get_peak_usage()
        free_memory_gb -= peak

        return free_memory_gb, free_processors

    def _prerun_check(self, graph):
        """Check if any node exeeds the available resources"""
        tasks_mem_gb = []
        tasks_num_th = []
        overrun_message_mem = None
        overrun_message_th = None
        for node in graph.nodes():
            if node.mem_gb > self.memory_gb:
                tasks_mem_gb.append((node.name, node.mem_gb))
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
