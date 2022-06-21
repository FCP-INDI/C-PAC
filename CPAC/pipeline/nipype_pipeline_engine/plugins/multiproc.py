"""
Override Nipype's MultiProc:
* _prerun_check to tell which Nodes use too many resources.
* _check_resources to account for the main process' memory usage.

Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""
from nipype.pipeline.plugins.multiproc import MultiProcPlugin as MultiProc
from .cpac_nipype_custom import CpacNipypeCustomPluginMixin


class MultiProcPlugin(CpacNipypeCustomPluginMixin, MultiProc):
    # pylint: disable=too-few-public-methods
    __doc__ = MultiProc.__doc__
    _check_resources = CpacNipypeCustomPluginMixin._check_resources
    _prerun_check = CpacNipypeCustomPluginMixin._prerun_check
