"""Utilties for documentation."""
import ast
from urllib import request
from urllib.error import ContentTooShortError, HTTPError, URLError
from CPAC import __version__
from CPAC.pipeline.nodeblock import NodeBlockData


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

    # Decorator node blocks

    if isinstance(fn, NodeBlockData):
        return fn.legacy_nodeblock_dict()

    # Legacy (docstring) node blocks

    fn_docstring = fn.__doc__
    init_dct_schema = ['name', 'config', 'switch', 'option_key',
                       'option_val', 'inputs', 'outputs']
    if 'Node Block:' in fn_docstring:
        fn_docstring = fn_docstring.split('Node Block:')[1]
    fn_docstring = fn_docstring.lstrip().replace('\n', '').rstrip()
    dct = ast.literal_eval(fn_docstring)
    for key in init_dct_schema:
        if key not in dct.keys():
            raise Exception('\n[!] Developer info: At least one of the '
                            'required docstring keys in your node block '
                            'is missing.\n\nNode block docstring keys:\n'
                            f'{init_dct_schema}\n\nYou provided:\n'
                            f'{dct.keys()}\n\nDocstring:\n{fn_docstring}\n\n')
    return dct


DOCS_URL_PREFIX = _docs_url_prefix()
