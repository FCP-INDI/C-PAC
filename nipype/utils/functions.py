# -*- coding: utf-8 -*-
"""
Handles custom functions used in Function interface. Future imports
are avoided to keep namespace as clear as possible.
"""
import inspect
from textwrap import dedent


def getsource(function):
    """Returns the source code of a function"""
    return dedent(inspect.getsource(function))
