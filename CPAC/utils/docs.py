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
"""Utilties for documentation."""
import ast
from urllib import request
from urllib.error import ContentTooShortError, HTTPError, URLError
from CPAC import __version__
from CPAC.utils import versioning


def docstring_parameter(*args, **kwargs):
    """Decorator to parameterize docstrings.
    Use double-curly-braces ({{}}) for literal curly braces.

    Examples
    --------
    >>> @docstring_parameter('test', answer='Yes it does.')
    ... def do_nothing():
    ...     '''Does this {} do anything? {answer}'''
    ...     pass
    >>> print(do_nothing.__doc__)
    Does this test do anything? Yes it does.
    >>> @docstring_parameter('test', answer='It should not.')
    ... def how_about_now():
    ...     '''How about {{ this }}?'''
    ...     pass
    >>> print(how_about_now.__doc__)
    How about { this }?
    """
    def dec(obj):
        if obj.__doc__ is None:
            obj.__doc__ = ''
        obj.__doc__ = obj.__doc__.format(*args, **kwargs)
        return obj
    return dec


def _docs_url_prefix():
    """Function to determine the URL prefix for this version of C-PAC"""
    def _url(url_version):
        return f'https://fcp-indi.github.io/docs/{url_version}'
    url_version = f'v{__version__}'
    try:
        request.urlopen(  # pylint: disable=consider-using-with
                        _url(url_version))
    except (ContentTooShortError, HTTPError, URLError):
        if 'dev' in url_version:
            url_version = 'nightly'
        else:
            url_version = 'latest'
    return _url(url_version)


def version_report() -> str:
    """A formatted block of versions included in CPAC's environment"""
    version_list = []
    for pkg, version in versioning.REPORTED.items():
        version_list.append(f'{pkg}: {version}')
        if pkg == 'Python':
            version_list.append('  Python packages')
            version_list.append('  ---------------')
            for ppkg, pversion in versioning.PYTHON_PACKAGES.items():
                version_list.append(f'  {ppkg}: {pversion}')
    return '\n'.join(version_list)


DOCS_URL_PREFIX = _docs_url_prefix()
