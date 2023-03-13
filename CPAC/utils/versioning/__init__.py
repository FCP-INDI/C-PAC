# Copyright (C) 2022  C-PAC Developers

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
"""Gather and report dependency versions alphabetically"""
from subprocess import PIPE, Popen
import sys
from pkg_resources import working_set


def cli_version(command, dependency=None, in_result=True, delimiter=' ',
                formatting=None):
    """Collect a version from a CLI

    Parameters
    ----------
    command : str

    dependency : str, optional
        software name to report if not included in result

    in_result : boolean, optional
        parse software name from result?

    delimiter : str, optional
        if parsing software name from result, what's the delimiter?

    formatting : func, optional
        if result needs any formatting, function to do so

    Returns
    -------
    dict
        {software: version}
    """
    with Popen(command, stdout=PIPE, shell=True) as _command:
        _version = _command.stdout.read().decode('utf-8')
    if formatting is not None:
        _version = formatting(_version)
    if in_result:
        return dict([tuple(_version.split(delimiter, 1))])
    return {dependency: _version}


def _version_sort(_version_item):
    """Key to report by case-insensitive dependecy name"""
    return _version_item[0].lower()

def first_line(stdout):
    """Return first line of stdout"""
    if '\n' in stdout:
        return stdout.split('\n', 1)[0]
    return stdout


PYTHON_PACKAGES = dict(sorted({
  d.key: d.version for d in list(working_set)}.items(), key=_version_sort))
REPORTED = dict(sorted({
  **cli_version('ldd --version', formatting=first_line),
  'Python': sys.version.replace('\n', ' ').replace('  ', ' ')
}.items(), key=_version_sort))

__all__ = ['PYTHON_PACKAGES', 'REPORTED']
