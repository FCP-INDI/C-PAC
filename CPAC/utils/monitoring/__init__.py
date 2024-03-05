# Copyright (C) 2021-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Customize Nipype's process monitoring for use in C-PAC.

See https://fcp-indi.github.io/docs/developer/nodes for C-PAC-specific documentation.
See https://nipype.readthedocs.io/en/latest/api/generated/nipype.utils.profiler.html for Nipype's documentation.
"""  # pylint: disable=line-too-long
from .config import LOGTAIL, WARNING_FREESURFER_OFF_WITH_DATA
from .custom_logging import (
    failed_to_start,
    FMLOGGER,
    getLogger,
    IFLOGGER,
    set_up_logger,
    UTLOGGER,
    WFLOGGER,
)
from .monitoring import (
    log_nodes_cb,
    log_nodes_initial,
    LoggingHTTPServer,
    LoggingRequestHandler,
    monitor_server,
    recurse_nodes,
)

__all__ = [
    "failed_to_start",
    "FMLOGGER",
    "getLogger",
    "IFLOGGER",
    "LoggingHTTPServer",
    "LoggingRequestHandler",
    "log_nodes_cb",
    "log_nodes_initial",
    "LOGTAIL",
    "monitor_server",
    "recurse_nodes",
    "set_up_logger",
    "UTLOGGER",
    "WARNING_FREESURFER_OFF_WITH_DATA",
    "WFLOGGER",
]
