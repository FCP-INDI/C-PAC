# Copyright (C) 2022-2024  C-PAC Developers

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
"""Utilties for C-PAC documentation."""
from functools import wraps
from typing import Callable, Optional
from urllib import request
from urllib.error import ContentTooShortError, HTTPError, URLError
from warnings import warn

from CPAC import __version__
from CPAC.utils import versioning


def deprecated(
    version: Optional[str] = None, explanation: Optional[str] = None
) -> Callable:
    """Mark a function as deprecated.

    Parameters
    ----------
    version : str, optional
        The version in which the function was deprecated.

    explanation : str, optional
        An explanation of why the function was deprecated.

    Returns
    -------
    Callable
        The decorated function.
    """

    def decorator(func: Callable) -> Callable:
        if func.__doc__ is None:
            func.__doc__ = ""

        note = ".. deprecated::"
        if version:
            note += f" {version}"
        if explanation:
            note += f"\n   {explanation}\n"
        func.__doc__ = note + "\n" + func.__doc__

        @wraps(func)
        def new_func(*args, **kwargs) -> Callable:
            """Warn that the function is deprecated."""
            _warning = f"Call to deprecated function '{func.__qualname__}'."
            if explanation:
                _warning += f" {explanation}\n"
            warn(
                _warning,
                category=DeprecationWarning,
                stacklevel=2,
            )
            return func(*args, **kwargs)

        return new_func

    return decorator


def docstring_parameter(*args, **kwargs) -> Callable:
    """Parameterize docstrings.

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

    def dec(obj: Callable) -> Callable:
        if obj.__doc__ is None:
            obj.__doc__ = ""
        obj.__doc__ = obj.__doc__.format(*args, **kwargs)
        return obj

    return dec


def _docs_url_prefix() -> str:
    """Determine the URL prefix for this version of C-PAC."""

    def _url(url_version: str) -> str:
        return f"https://fcp-indi.github.io/docs/{url_version}"

    url_version = f"v{__version__}"
    try:
        request.urlopen(  # pylint: disable=consider-using-with
            _url(url_version)
        )
    except (ContentTooShortError, HTTPError, URLError):
        if "dev" in url_version:
            url_version = "nightly"
        else:
            url_version = "latest"
    return _url(url_version)


def version_report() -> str:
    """Return a formatted block of versions included in CPAC's environment."""
    version_list = []
    for pkg, version in versioning.REPORTED.items():
        version_list.append(f"{pkg}: {version}")
        if pkg == "Python":
            version_list.append("  Python packages")
            version_list.append("  ---------------")
            for ppkg, pversion in versioning.PYTHON_PACKAGES.items():
                version_list.append(f"  {ppkg}: {pversion}")
    return "\n".join(version_list)


def outdent_lines(docstring: str, spaces: int = 4) -> str:
    """Outdent lines in a string by specified number of spaces.

    Only outdents lines that are at least that indented.
    Useful for combining docstrings.

    Examples
    --------
    >>> import re
    >>> re.findall(r'^    Only.*$', outdent_lines.__doc__, flags=re.MULTILINE)
    ['    Only outdents lines that are at least that indented.']
    >>> re.findall(r'^Only.*$', outdent_lines.__doc__, flags=re.MULTILINE)
    []
    >>> re.findall(r'^    Only.*$', outdent_lines(outdent_lines.__doc__),
    ...     flags=re.MULTILINE)
    []
    >>> re.findall(r'^Only.*$', outdent_lines(outdent_lines.__doc__),
    ...     flags=re.MULTILINE)
    ['Only outdents lines that are at least that indented.']
    >>> re.findall(r'^ Only.*$', outdent_lines(outdent_lines.__doc__, 3),
    ...     flags=re.MULTILINE)
    [' Only outdents lines that are at least that indented.']
    """
    new_docstring = []
    for line in docstring.split("\n"):
        if line.startswith(" " * spaces):
            new_docstring.append(line[spaces:])
        else:
            new_docstring.append(line)
    return "\n".join(new_docstring)


DOCS_URL_PREFIX = _docs_url_prefix()
