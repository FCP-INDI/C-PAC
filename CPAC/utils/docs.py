"""Utilties for documentation."""


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
        try:
            obj.__doc__ = obj.__doc__.format(*args, **kwargs)
        except KeyError as key_error:
            print(key_error)
        return obj
    return dec
