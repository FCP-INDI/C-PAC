"""Utilties for documentation."""
import ast
from urllib import request
from urllib.error import ContentTooShortError, HTTPError, URLError
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


def grab_docstring_dct(fn):
    """Function to grab a NodeBlock dictionary from a docstring.

    Parameters
    ----------
    fn : function
        The NodeBlock function with the docstring to be parsed.

    Returns
    -------
    dct : dict
        A NodeBlock configuration dictionary.
    """
    fn_docstring = fn.__doc__
    init_dct_schema = ['name', 'config', 'switch', 'option_key',
                       'option_val', 'inputs', 'outputs']
    if 'Node Block:' in fn_docstring:
        fn_docstring = fn_docstring.split('Node Block:')[1]
    fn_docstring = fn_docstring.lstrip().replace('\n', '')
    dct = ast.literal_eval(fn_docstring)
    for key in init_dct_schema:
        if key not in dct.keys():
            raise Exception('\n[!] Developer info: At least one of the '
                            'required docstring keys in your node block '
                            'is missing.\n\nNode block docstring keys:\n'
                            f'{init_dct_schema}\n\nYou provided:\n'
                            f'{dct.keys()}\n\nDocstring:\n{fn_docstring}\n\n')
    return dct


def retry_docstring(orig):
    """Decorator to autodocument retries.

    Examples
    --------
    >>> @retry_docstring(grab_docstring_dct)
    ... def do_nothing():
    ...     '''Does this do anything?'''
    ...     pass
    >>> print(do_nothing.__doc__)
    Does this do anything?
    Retries the following after a failed QC check:
    Function to grab a NodeBlock dictionary from a docstring.
    <BLANKLINE>
        Parameters
        ----------
        fn : function
            The NodeBlock function with the docstring to be parsed.
    <BLANKLINE>
        Returns
        -------
        dct : dict
            A NodeBlock configuration dictionary.
    <BLANKLINE>
    """
    def retry(obj):
        if obj.__doc__ is None:
            obj.__doc__ = ''
        origdoc = (f'{orig.__module__}.{orig.__name__}' if
                   orig.__doc__ is None else orig.__doc__)
        obj.__doc__ = '\n'.join([
            obj.__doc__, 'Retries the following after a failed QC check:',
            origdoc])
        return obj
    return retry


DOCS_URL_PREFIX = _docs_url_prefix()
