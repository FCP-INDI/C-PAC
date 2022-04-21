"""Utilties for documentation."""
import re


def camel_case(snake_case):
    """
    Function to convert snake_case to camelCase. For each term after
    the first, uppercases the first letter. For each term before the
    last, lowercases the last letter. Otherwise keeps casing.

    Parameters
    ----------
    snake_case : str

    Returns
    -------
    str

    Examples
    --------
    >>> camel_case('snake_case')
    'snakeCase'
    >>> camel_case('with_ALLCAPS')
    'withALLCAPS'
    >>> camel_case('with_interNaL-caPs')
    'withInterNalCaPs'
    >>> camel_case('-'.join([camel_case(f'{key}-{value}') for key, value
    ...  in {'regressor': 'aCompCor', 'brain_extraction': 'FSL'}.items()]))
    'regressorACompCorBrainExtractionFSL'
    """
    pieces = re.split(r'[\-_]', snake_case)
    if len(pieces) > 1:
        pieces = [*[case_terminal_letter(piece, 'lower', 'last') for
                    piece in pieces[:-1]], pieces[-1]]
        return ''.join([pieces[0],
                        *[case_terminal_letter(piece, 'upper', 'first') for
                          piece in pieces[1:]]])
    return snake_case


def case_terminal_letter(cased_str, case, terminal):
    """Function to capitalize or lowercase the first or last letter of
    a cased string.

    Parameters
    ----------
    cased_str : str

    case : str
        'upper' or 'lower'

    terminal : str
        'first' or 'last'

    Examples
    --------
    >>> case_terminal_letter('BIDS', 'lower', 'first')
    'bIDS'
    >>> case_terminal_letter('BIDS', 'lower', 'last')
    'BIDs'
    >>> case_terminal_letter('bids', 'upper', 'first')
    'Bids'
    >>> case_terminal_letter('bids', 'upper', 'last')
    'bidS'
    """
    if not cased_str:
        return cased_str
    if not isinstance(cased_str, str):
        cased_str = str(cased_str)
    if terminal == 'first':
        return getattr(cased_str, case)()[0] + cased_str[1:]
    if terminal == 'last':
        return cased_str[:-1] + getattr(cased_str[-1], case)()[-1]


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
