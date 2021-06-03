"""Import Nipype's pipeline plugins and selectively override"""
from nipype.pipeline.plugins import *  # noqa F401,F403
from .base import logfile, run_node  # noqa F401
# Override LegacyMultiProc
from .legacymultiproc import LegacyMultiProcPlugin  # noqa F401