"""Override Nipype's MultiProc's _prerun_check to tell which
Nodes use too many resources."""
from nipype.pipeline.plugins.multiproc import MultiProcPlugin as MultiProc
from .verbose_prerun import VerbosePrerun


class MultiProcPlugin(MultiProc, VerbosePrerun):
    def __init__(self, plugin_args=None):
        super().__init__(plugin_args)

    _prerun_check = VerbosePrerun._prerun_check
