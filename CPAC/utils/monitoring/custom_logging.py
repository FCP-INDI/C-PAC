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
"""Funtions for logging."""
import logging
import os
import subprocess
from sys import exc_info as sys_exc_info
from traceback import print_exception
from nipype import logging as nipype_logging
from CPAC.utils.docs import docstring_parameter
from CPAC.utils.monitoring.config import MOCK_LOGGERS


def failed_to_start(log_dir, exception):
    """Launch a failed-to-start logger for a run that failed to start.
    Must be called from within an ``except`` block.

    Parameters
    ----------
    log_dir : str
        path to logging directory

    exception : Exception
    """
    logger = set_up_logger('failedToStart', 'failedToStart.log', 'error',
                           log_dir, True)
    logger.exception('C-PAC failed to start')
    logger.exception(exception)


def getLogger(name):  # pylint: disable=invalid-name
    """Function to get a mock logger if one exists, falling back on
    real loggers.

    Parameters
    ----------
    name : str

    Returns
    -------
    logger : CPAC.utils.monitoring.custom_logging.MockLogger or logging.Logger
    """
    if name in MOCK_LOGGERS:
        return MOCK_LOGGERS[name]
    logger = nipype_logging.getLogger(name)
    return logging.getLogger(name) if logger is None else logger


def log_failed_subprocess(cpe):
    """Pass STDERR from a subprocess to the interface's logger

    Parameters
    ----------
    cpe : subprocess.CalledProcessError
    """
    logger = getLogger('nipype.interface')
    logger.error("%s\nExit code %s", cpe.output, cpe.returncode)


def log_subprocess(cmd, *args, raise_error=True, **kwargs):
    """Pass STDERR and STDOUT from subprocess to interface's logger.
    This function is nearly a drop-in replacement for
    `subprocess.check_output`.

    Caveat: if you're assigning to a variable (like

    >>> output = subprocess.check_output(cmd)  # doctest: +SKIP

    ), the new function also returns the command's exit code, so you can just
    assign that to a throwaway variable if you don't want it

    >>> output, _ = log_subprocess(cmd)  # doctest: +SKIP

    or subscript the command like

    >>> output = log_subprocess(cmd)[0]  # doctest: +SKIP

    . If you're not assigning to a variable, it doesn't matter and just

    >>> log_subprocess(cmd)  # doctest: +SKIP

    should work just like

    >>> subprocess.check_output(cmd)  # doctest: +SKIP

    Parameters
    ----------
    cmd : str
        command to run with `subprocess.check_output`

    raise_error : boolean
        raise any exception after logging

    args, kwargs : any
        pass-through arguments for subprocess.check_output

    Returns
    -------
    output : str

    exit_code : int
    """
    logger = getLogger('nipype.interface')
    try:
        output = subprocess.check_output(cmd, *args, stderr=subprocess.STDOUT,
                                         universal_newlines=True, **kwargs)
        logger.info(output)
    except subprocess.CalledProcessError as cpe:
        log_failed_subprocess(cpe)
        if raise_error:
            raise
        return cpe.output, cpe.returncode
    return output, 0


# pylint: disable=too-few-public-methods
class MockHandler:
    """Handler for MockLogger."""
    def __init__(self, filename):
        self.baseFilename = filename  # pylint: disable=invalid-name


# pylint: disable=too-few-public-methods
class MockLogger:
    """Mock logging.Logger to provide the same API without keeping the
    logger in memory."""
    def __init__(self, name, filename, level, log_dir):
        self.name = name
        self.level = level
        self.handlers = [MockHandler(os.path.join(log_dir, filename))]
        MOCK_LOGGERS[name] = self
        for loglevel in ['debug', 'info', 'warning', 'error', 'critical']:
            # set up log methods for all built-in levels
            setattr(self, loglevel, self._factory_log(loglevel))

    def exception(self, msg, *args, exc_info=True, **kwargs):
        # pylint: disable=missing-function-docstring,no-member
        return self.error(msg, *args, exc_info=exc_info, **kwargs)

    exception.__doc__ = logging.exception.__doc__

    def _factory_log(self, level):
        r"""Generate a log method like `self.log(message)` for a given
        built-in level."""
        @docstring_parameter(level=level)
        def _log(message, *items, exc_info=False):
            """Log a message if logging level >= {level}. See `Logging Levels <https://docs.python.org/3/library/logging.html#levels>`_ for a list of levels."""
            if self.level == 0 or self.level >= getattr(logging, level.upper(),
                                                        logging.NOTSET):
                with open(self.handlers[0].baseFilename, 'a',
                          encoding='utf-8') as log_file:
                    if exc_info and isinstance(message, Exception):
                        value, traceback = sys_exc_info()[1:]
                        print_exception(_lazy_sub(message, *items),
                                        value=value, tb=traceback,
                                        file=log_file)
                    else:
                        print(_lazy_sub(message, *items), file=log_file)
        return _log

    def delete(self):
        """Delete the mock logger from memory."""
        del MOCK_LOGGERS[self.name]


def _lazy_sub(message, *items):
    """Given lazy-logging syntax, return string with substitutions

    Parameters
    ----------
    message : str

    items : tuple

    Returns
    -------
    str

    Examples
    --------
    >>> _lazy_sub('no substitution')
    'no substitution'
    >>> _lazy_sub('%s substitution', 'yes')
    'yes substitution'
    >>> _lazy_sub('%s substitution %s', 'yes', 'again')
    'yes substitution again'
    """
    try:
        return str(message) % items
    except (AttributeError, TypeError):
        return str([message, *items])


def set_up_logger(name, filename=None, level=None, log_dir=None, mock=False,
                  overwrite_existing=False):
    r"""Function to initialize a logger

    Parameters
    ----------
    name : str
        logger name (for subsequent calls to ``logging.getLogger``) to
        write to the same log file)

    filename : str, optional
        filename to write log to. If not specified, filename will be
        the same as ``name`` with the extension ``log``

    level : str, optional
        one of ``{critical, error, warning, info, debug, notset}``,
        case-insensitive

    log_dir : str, optional

    mock : bool, optional
        if ``True``, return a ``CPAC.utils.monitoring.MockLogger``
        instead of a ``logging.Logger``

    Returns
    -------
    logger : logging.Handler
        initialized logging Handler

    Examples
    --------
    >>> lg = set_up_logger('test')
    >>> lg.handlers[0].baseFilename.split('/')[-1]
    'test.log'
    >>> lg.level
    0
    >>> lg = set_up_logger('second_test', 'specific_filename.custom', 'debug')
    >>> lg.handlers[0].baseFilename.split('/')[-1]
    'specific_filename.custom'
    >>> lg.level
    10
    >>> lg = set_up_logger('third_test', mock=True)
    >>> getLogger('third_test') == lg
    True
    >>> 'third_test' in MOCK_LOGGERS
    True
    >>> lg.delete()
    >>> 'third_test' in MOCK_LOGGERS
    False
    """
    if filename is None:
        filename = f'{name}.log'
    try:
        level = getattr(logging, level.upper())
    except AttributeError:
        level = logging.NOTSET
    if log_dir is None:
        log_dir = os.getcwd()
    filepath = os.path.join(log_dir, filename)
    if overwrite_existing and os.path.exists(filepath):
        with open(filepath, 'w') as log_file:
            log_file.write('')
    if mock:
        return MockLogger(name, filename, level, log_dir)
    logger = getLogger(name)
    logger.setLevel(level)
    handler = logging.FileHandler(filepath)
    logger.addHandler(handler)
    return logger
