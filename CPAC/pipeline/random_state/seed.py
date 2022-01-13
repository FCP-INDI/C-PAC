'''Functions to set, check, and log random seed'''
import random

import numpy as np


def random_random_seed():
    '''Returns a random postive integer up to 2147483647

    Parameters
    ----------
    None

    Returns
    -------
    random_seed : int

    Examples
    --------
    >>> 0 < random_random_seed() <= np.iinfo(np.int32).max
    True
    '''
    return random.randint(1, np.iinfo(np.int32).max)
