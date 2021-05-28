"""
Override Nipype's LegacyMultiproc:
* _prerun_check to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.
"""
import platform
import resource
from nipype.pipeline.plugins.legacymultiproc import \
    LegacyMultiProcPlugin as LegacyMultiProc
from .verbose_prerun import VerbosePrerun

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

    _prerun_check = VerbosePrerun._prerun_check

    def _check_resources(self, running_tasks):
        """
        Make sure there are resources available, accounting for Nipype memory usage
        """

        free_memory_gb, free_processors = super()._check_resources(running_tasks)

        # Nipype memory usage
        peak = get_peak_usage()
        free_memory_gb -= peak

        return free_memory_gb, free_processors

