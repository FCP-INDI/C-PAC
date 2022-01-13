'''Functions to set, check, and log random seed'''
import os
import random

import numpy as np

from CPAC.utils.monitoring import set_up_logger


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


def set_up_random_state(workflow, seed, log_dir=None):
    '''Prepare C-PAC for random seed setting and logging

    Parameters
    ----------
    workflow : nipype.pipeline.engine.Workflow

    seed : int, 'random', or None

    log_dir : str, optional

    Returns
    -------
    workflow

    Examples
    --------
    >>> from CPAC.pipeline.nipype_pipeline_engine import Workflow
    >>> wf = set_up_random_state(Workflow('random_workflow'), 'random')
    >>> wf.seed
    'random'
    >>> wf = set_up_random_state(Workflow('random_workflow'), 'rando')
    Traceback (most recent call last):
    ValueError: Valid random seeds are positive integers up to 2147483647, "random", or None, not rando
    >>> wf = set_up_random_state(Workflow('random_workflow'), 100)
    >>> wf.seed
    100
    >>> wf = set_up_random_state(Workflow('random_workflow'), 0)
    Traceback (most recent call last):
    ValueError: Valid random seeds are positive integers up to 2147483647, "random", or None, not 0
    >>> wf = set_up_random_state(Workflow('random_workflow'), None)
    >>> wf.seed
    Traceback (most recent call last):
    AttributeError: 'Workflow' object has no attribute 'seed'
    '''  # noqa E501  # pylint: disable=line-too-long
    if seed is not None:
        if (seed != 'random' and not (
            isinstance(seed, int) and (0 < seed <= np.iinfo(np.int32).max)
        )):
            raise ValueError('Valid random seeds are positive integers up to '
                             f'2147483647, "random", or None, not {seed}')
        workflow.seed = seed
    if log_dir is None:
        log_dir = os.getcwd()
    set_up_logger('random', level='info', log_dir=log_dir)
    return workflow
