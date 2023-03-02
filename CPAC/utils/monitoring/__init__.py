'''Module to customize Nipype's process monitoring for use in C-PAC

See https://fcp-indi.github.io/docs/developer/nodes for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.utils.profiler.html for Nipype's documentation.'''  # noqa: E501  # pylint: disable=line-too-long
from .config import LOGTAIL, WARNING_FREESURFER_OFF_WITH_DATA
from .custom_logging import failed_to_start, getLogger, set_up_logger
from .monitoring import LoggingHTTPServer, LoggingRequestHandler, \
                        log_nodes_cb, log_nodes_initial, monitor_server, \
                        recurse_nodes

__all__ = ['failed_to_start', 'getLogger', 'LoggingHTTPServer',
           'LoggingRequestHandler', 'log_nodes_cb', 'log_nodes_initial',
           'LOGTAIL', 'monitor_server', 'recurse_nodes', 'set_up_logger',
           'WARNING_FREESURFER_OFF_WITH_DATA']
