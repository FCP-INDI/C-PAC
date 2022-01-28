'''Funtions for logging.'''
import logging
import os


def set_up_logger(name, filename=None, level=None, log_dir=None):
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
    '''
    if filename is None:
        filename = f'{name}.log'
    try:
        level = getattr(logging, level.upper())
    except AttributeError:
        level = logging.NOTSET
    if log_dir is None:
        log_dir = os.getcwd()
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.FileHandler(os.path.join(log_dir, filename))
    logger.addHandler(handler)
    return logger
