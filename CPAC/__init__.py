"""
Configurable Pipeline for the Analysis of Connectomes
=====================================================

CPAC is a configurable, open-source, Nipype-based, automated processing
pipeline for resting state functional MRI (R-fMRI) data, for use by
both novice and expert users.

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

from .info import __version__
version = __version__
_url_version = 'nightly' if __version__.endswith(
    '-dev') or __version__.endswith('.dev') else f'v{__version__.lstrip("v")}'
docs_prefix = f'https://fcp-indi.github.io/docs/{_url_version}'
license_notice = f"""Copyright (C) 2022  C-PAC Developers.

This program comes with ABSOLUTELY NO WARRANTY. This is free software,
and you are welcome to redistribute it under certain conditions. For
details, see {docs_prefix}/license or the COPYING and
COPYING.LESSER files included in the source code."""
