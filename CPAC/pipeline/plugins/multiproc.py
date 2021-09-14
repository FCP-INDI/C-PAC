"""
Override Nipype's MultiProc:
* _prerun_check to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.
"""
from nipype.pipeline.plugins.multiproc import MultiProcPlugin as MultiProc
from .cpac_nipype_custom import CpacNipypeCustomPlugin


class MultiProcPlugin(MultiProc, CpacNipypeCustomPlugin):
    def __init__(self, plugin_args=None):
        super().__init__(plugin_args)

    _check_resources = CpacNipypeCustomPlugin._check_resources
    _prerun_check = CpacNipypeCustomPlugin._prerun_check
