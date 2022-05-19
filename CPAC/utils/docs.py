"""Utilties for documentation."""
import ast


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


def grab_docstring_dct(fn_docstring):
    """Function to grab a NodeBlock dictionary from a docstring.

    Parameters
    ----------
    fn_docstring : str
        The docstring of the function to be parsed.

    Returns
    -------
    dct : dict
        A NodeBlock configuration dictionary.
    """
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
