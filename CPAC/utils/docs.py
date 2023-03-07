"""Utilties for documentation."""
import ast
from urllib import request
from urllib.error import ContentTooShortError, HTTPError, URLError
from CPAC import __version__


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


DOCS_URL_PREFIX = _docs_url_prefix()
