'''Functions to set, check, and log random seed'''
import os
import random

import numpy as np

from CPAC.utils.monitoring.custom_logging import set_up_logger

_seed = {'seed': None}


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


def random_seed():
    '''Function to access current random seed

    Parameters
    ----------
    None

    Returns
    -------
    seed : int or None
    '''
    if _seed['seed'] == 'random':
        return random_random_seed()
    return _seed['seed']


def set_up_random_state(seed, log_dir=None):
    '''Prepare C-PAC for random seed setting and logging

    Parameters
    ----------
    seed : int, 'random', or None

    log_dir : str, optional

    Returns
    -------
    workflow

    Examples
    --------
    >>> from CPAC.pipeline.random_state import random_seed
    >>> set_up_random_state('random')
    >>> isinstance(random_seed(), int)
    True
    >>> set_up_random_state('rando')
    Traceback (most recent call last):
    ValueError: Valid random seeds are positive integers up to 2147483647, "random", or None, not rando
    >>> set_up_random_state(100)
    >>> random_seed()
    100
    >>> set_up_random_state(0)
    Traceback (most recent call last):
    ValueError: Valid random seeds are positive integers up to 2147483647, "random", or None, not 0
    >>> set_up_random_state(None)
    >>> random_seed()

    '''  # noqa E501  # pylint: disable=line-too-long
    if seed is not None:
        if (seed != 'random' and not (
            isinstance(seed, int) and
            (0 < int(seed) <= np.iinfo(np.int32).max)
        )):
            raise ValueError('Valid random seeds are positive integers up to '
                             f'2147483647, "random", or None, not {seed}')
    try:
        _seed['seed'] = int(seed)
    except (TypeError, ValueError):
        _seed['seed'] = seed
    if log_dir is None:
        log_dir = os.getcwd()
    set_up_logger('random', level='info', log_dir=log_dir)
