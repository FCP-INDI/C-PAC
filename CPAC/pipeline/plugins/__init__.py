"""Import Nipype's pipeline plugins and selectively override"""
from nipype.pipeline.plugins import *  # noqa F401,F403
# Override LegacyMultiProc
from .legacymultiproc import LegacyMultiProcPlugin  # noqa F401
