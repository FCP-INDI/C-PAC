"""Utilties for documentation."""
from urllib import request
from urllib.error import HTTPError
from CPAC import __version__


def docstring_parameter(*args, **kwargs):
    """Decorator to parameterize docstrings.

    Examples
    --------
    >>> @docstring_parameter('test', answer='Yes it does.')
    ... def do_nothing():
    ...     '''Does this {} do anything? {answer}'''
    ...     pass
    >>> print(do_nothing.__doc__)
    Does this test do anything? Yes it does.
    """
    def dec(obj):
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
    except HTTPError:
        if 'dev' in url_version:
            url_version = 'nightly'
        else:
            url_version = 'latest'
    return _url(url_version)


DOCS_URL_PREFIX = _docs_url_prefix()
