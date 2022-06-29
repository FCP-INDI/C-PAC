'''Funtions for logging.'''
import logging
import os

from CPAC.utils.docs import docstring_parameter
from CPAC.utils.monitoring.config import MOCK_LOGGERS


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
    return logging.getLogger(name)


# pylint: disable=too-few-public-methods
class MockHandler:
    '''Handler for MockLogger.'''
    def __init__(self, baseFilename):
        self.baseFilename = baseFilename  # pylint: disable=invalid-name


# pylint: disable=too-few-public-methods
class MockLogger:
    '''Mock logging.Logger to provide the same API without keeping the
    logger in memory.'''
    def __init__(self, name, filename, level, log_dir):
        self.name = name
        self.level = level
        self.handlers = [MockHandler(os.path.join(log_dir, filename))]
        MOCK_LOGGERS[name] = self
        for loglevel in ['debug', 'info', 'warning', 'error', 'critical']:
            # set up log methods for all built-in levels
            setattr(self, loglevel, self._factory_log(loglevel))

    def _factory_log(self, level):
        r"""Generate a log method like `self.log(message)` for a given
        built-in level."""
        @docstring_parameter(level=level)
        def _log(message):
            """Log a message if logging level >= {level}. See `Logging Levels <https://docs.python.org/3/library/logging.html#levels>`_ for a list of levels."""
            if self.level == 0 or self.level >= getattr(logging, level.upper(),
                                                        logging.NOTSET):
                with open(self.handlers[0].baseFilename, 'a') as log_file:
                    print(message, file=log_file)
        return _log

    def delete(self):
        '''Delete the mock logger from memory.'''
        del MOCK_LOGGERS[self.name]


def set_up_logger(name, filename=None, level=None, log_dir=None, mock=False,
                  overwrite_existing=False):
    r'''Function to initialize a logger

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
    '''
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
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.FileHandler(filepath)
    logger.addHandler(handler)
    return logger
