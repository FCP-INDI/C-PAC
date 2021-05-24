"""Override Nipype's LegacyMultiProc's _prerun_check to tell which
Nodes use too many resources."""
from nipype.pipeline.plugins.legacymultiproc import \
    LegacyMultiProcPlugin as LegacyMultiProc
from .verbose_prerun import VerbosePrerun


class LegacyMultiProcPlugin(LegacyMultiProc):
    def __init__(self, plugin_args=None):
        super().__init__(plugin_args)

    _prerun_check = VerbosePrerun._prerun_check
