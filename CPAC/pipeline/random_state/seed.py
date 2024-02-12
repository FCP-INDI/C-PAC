# Copyright (C) 2022-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
'''Functions to set, check, and log random seed'''
import os
import random

import numpy as np
from nipype.interfaces.ants.registration import Registration
from nipype.interfaces.ants.segmentation import Atropos
from nipype.interfaces.freesurfer.preprocess import ApplyVolTransform, ReconAll
from nipype.interfaces.fsl.maths import MathsCommand
from nipype.interfaces.fsl.utils import ImageMaths

from CPAC.utils.interfaces.ants import AI
from CPAC.utils.monitoring.custom_logging import getLogger, set_up_logger

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
        _seed['seed'] = random_random_seed()
    return _seed['seed']


def random_seed_flags():
    '''Function to return dictionary of flags with current random seed.

    Parameters
    ----------
    None

    Returns
    -------
    dict of {'functions': {function: function}, 'interfaces': {list or tuple}}
        Developer note: sequence matters here! Only the first match
        will be applied!

        In the 'functions' sub-dictionary, each key should be a function
        used in a :py:class:`CPAC.utils.interfaces.function.Function`,
        and each value should be a function that takes a single
        parameter (a function) to apply to the key-function such that,
        once the value-function is applied, the random seed is applied
        when the key-function is run.

        In the 'interfaces' sub-dictionary, each key should be an
        :py:class:`nipype.interfaces.base.core.Interface`, and each
        value should be either a list of strings to add to the
        interface's flags/args or a tuple of (list of strings to add to,
        list of strings to remove from) the interface's flags/args.

    Examples
    --------
    >>> list(random_seed_flags().keys())
    ['functions', 'interfaces']
    >>> all([isinstance(random_seed_flags()[key], dict) for key in [
    ...     'functions', 'interfaces']])
    True
    >>> rs = set_up_random_state('random')
    >>> list(random_seed_flags().keys())
    ['functions', 'interfaces']
    >>> all([isinstance(random_seed_flags()[key], dict) for key in [
    ...     'functions', 'interfaces']])
    True
    '''
    from CPAC.registration.utils import hardcoded_reg
    seed = random_seed()
    if seed is None:
        return {'functions': {}, 'interfaces': {}}
    return {
        'functions': {
            # function: lambda function to apply to function source
            hardcoded_reg: lambda fn_string: fn_string.replace(
                'regcmd = ["antsRegistration"]',
                f'regcmd = ["antsRegistration", "--random-seed", \"{seed}\"]'
            )
        },
        'interfaces': {
            # interface: [flags to apply]
            # OR
            # interface: ([flags to apply], [flags to remove])
            #
            # ANTs
            # NOTE: Atropos gives the option "Initialize internal random number
            #       generator with a random seed. Otherwise, initialize with a
            #       constant seed number," so for Atropos nodes, the built-in
            #       Atropos constant seed is used if a seed is specified for
            #       C-PAC
            AI: _reusable_flags()['ANTs'],
            Registration: _reusable_flags()['ANTs'],
            Atropos: (['--use-random-seed 0'],
                      [flag for one in ['', ' 1'] for flag in
                      [f'--use-random-seed{one}', f'-r{one}']]),
            # FreeSurfer
            ReconAll: ['-norandomness', f'-rng-seed {seed}'],
            ApplyVolTransform: [f'--seed {seed}'],
            # FSL
            ImageMaths: _reusable_flags()['FSL'],
            MathsCommand: _reusable_flags()['FSL']
        }
    }


def _reusable_flags():
    seed = random_seed()
    return {
        'ANTs': [f'--random-seed {seed}'],
        'FSL': [f'-seed {seed}']
    }


def set_up_random_state(seed):
    '''Set global random seed

    Parameters
    ----------
    seed : int, 'random', or None

    Returns
    -------
    seed: int or None

    Examples
    --------
    >>> isinstance(set_up_random_state('random'), int)
    True
    >>> set_up_random_state('rando')
    Traceback (most recent call last):
    ValueError: Valid random seeds are positive integers up to 2147483647, "random", or None, not rando
    >>> set_up_random_state(100)
    100
    >>> set_up_random_state(0)
    Traceback (most recent call last):
    ValueError: Valid random seeds are positive integers up to 2147483647, "random", or None, not 0
    >>> set_up_random_state(None)

    '''  # noqa: E501  # pylint: disable=line-too-long
    if seed is not None:
        if seed == 'random':
            seed = random_random_seed()
        else:
            try:
                seed = int(seed)
                assert 0 < seed <= np.iinfo(np.int32).max
            except(ValueError, TypeError, AssertionError):
                raise ValueError('Valid random seeds are positive integers up to '
                                    f'2147483647, "random", or None, not {seed}')
    
    _seed['seed'] = seed
    return random_seed()


def set_up_random_state_logger(log_dir):
    '''Prepare C-PAC for logging random seed use.

    Parameters
    ----------
    log_dir : str
    '''
    set_up_logger('random', level='info', log_dir=log_dir, mock=True)
    getLogger('random').info('seed: %s', random_seed())
