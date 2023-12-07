"""Custom datatypes for C-PAC"""
class ItemFromList:  # pylint: disable=too-few-public-methods
    """Coerce single-item lists into just the only item in the list.
    Returns item if item is not a list, set, or tuple.
    Raises CoerceInvalid if impossible.

    Examples
    --------
    >>> ItemFromList(['seagull'])
    'seagull'
    >>> ItemFromList(['two', 'seagulls'])
    Traceback (most recent call last):
        ...
    voluptuous.error.CoerceInvalid: Cannot coerce list of length 2 to item
    >>> ItemFromList('string')
    'string'
    """
    def __new__(cls, list_of_one, msg=None):
        """Initialize item from list"""
        from voluptuous import CoerceInvalid, Length, LengthInvalid
        if not isinstance(list_of_one, (list, set, tuple)):
            return list_of_one
        try:
            Length(max=1)(list_of_one)
        except LengthInvalid as length_invalid:
            raise CoerceInvalid(
                f'Cannot coerce list of length {len(list_of_one)} to item'
                if msg is None else msg) from length_invalid
        return list_of_one[0]


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
    def __init__(self, *args, **kwargs):
        """Initialize ListFromItem"""
        if len(args) == 1 and not isinstance(args[0], (list, tuple)):
            if args[0] is None:
                args = ()
            else:
                args = ([args[0]],)
        list.__init__(self, *args, **kwargs)
