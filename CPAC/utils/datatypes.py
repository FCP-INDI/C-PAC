"""Custom datatypes for C-PAC"""


class ListFromItem(list):
    """Subclass of list to coerce non-lists into lists

    Examples
    --------
    >>> list('one')
    ['o', 'n', 'e']
    >>> ListFromItem('one')
    ['one']
    >>> list(['one'])
    ['one']
    >>> ListFromItem(['one'])
    ['one']
    >>> list()
    []
    >>> ListFromItem()
    []
    >>> list(None)
    Traceback (most recent call last):
    ...
    TypeError: 'NoneType' object is not iterable
    >>> ListFromItem(None)
    []
    """
    def __init__(self, /, *args, **kwargs):
        """Initialize ListFromItem"""
        if len(args) == 1 and not isinstance(args[0], (list, tuple)):
            if args[0] is None:
                args = ()
            else:
                args = ([args[0]],)
        list.__init__(self, *args, **kwargs)
