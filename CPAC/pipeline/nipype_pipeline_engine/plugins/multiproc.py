"""
Override Nipype's LegacyMultiProc:
* _prerun_check to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.
"""
from nipype.pipeline.plugins.multiproc import \
    LegacyMultiProcPlugin as LegacyMultiProc
from .cpac_nipype_custom import CpacNipypeCustomPluginMixin


class LegacyMultiProcPlugin(CpacNipypeCustomPluginMixin, LegacyMultiProc):
    # pylint: disable=too-few-public-methods
    __doc__ = LegacyMultiProc.__doc__
    _check_resources = CpacNipypeCustomPluginMixin._check_resources
    _prerun_check = CpacNipypeCustomPluginMixin._prerun_check
