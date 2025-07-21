# Copyright (C) 2022-2025  C-PAC Developers

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
from pathlib import Path
import subprocess
from sys import exc_info as sys_exc_info
from traceback import print_exception
from typing import Literal, Optional, Sequence, TYPE_CHECKING, TypeAlias

import yaml
from nipype import config as nipype_config, logging as nipype_logging

from CPAC.utils.docs import docstring_parameter
from CPAC.utils.monitoring.config import MOCK_LOGGERS

if TYPE_CHECKING:
    from CPAC.utils.configuration import Configuration
LogLevel: TypeAlias = (
    Literal[
        "CRITICAL",
        "critical",
        "Critical",
        "DEBUG",
        "debug",
        "Debug",
        "ERROR",
        "error",
        "Error",
        "INFO",
        "info",
        "Info",
        "NOTSET",
        "notset",
        "Notset",
        "NotSet",
        "notSet",
        "WARNING",
        "warning",
        "Warning",
    ]
    | int
)


def failed_to_start(log_dir, exception):
    """Launch a failed-to-start logger for a run that failed to start.

    Must be called from within an ``except`` block.

    Parameters
    ----------
    log_dir : str
        path to logging directory

    exception : Exception
    """
    logger = set_up_logger("failedToStart", "failedToStart.log", "error", log_dir, True)
    logger.exception("C-PAC failed to start")
    logger.exception(exception)


def getLogger(name: str) -> "logging.Logger | MockLogger":  # pylint: disable=invalid-name
    """Get a mock logger if one exists, falling back on real loggers."""
    if name in MOCK_LOGGERS:
        return MOCK_LOGGERS[name]
    logger = nipype_logging.getLogger(name)
    if logger is None:
        logger = logging.getLogger(name)
        if not logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(message)s"))
            logger.setLevel(logging.INFO)
            logger.addHandler(handler)
    return logger


# Nipype built-in loggers
IFLOGGER = getLogger("nipype.interface")
FMLOGGER = getLogger("nipype.filemanip")
UTLOGGER = getLogger("nipype.utils")
WFLOGGER = getLogger("nipype.workflow")


def log_failed_subprocess(cpe):
    """Pass STDERR from a subprocess to the interface's logger.

    Parameters
    ----------
    cpe : subprocess.CalledProcessError
    """
    IFLOGGER.error("%s\nExit code %s", cpe.output, cpe.returncode)


def log_subprocess(cmd, *args, raise_error=True, **kwargs):
    """Pass STDERR and STDOUT from subprocess to interface's logger.

    This function is nearly a drop-in replacement for `subprocess.check_output`.

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
    try:
        output = subprocess.check_output(
            cmd, *args, stderr=subprocess.STDOUT, universal_newlines=True, **kwargs
        )
        IFLOGGER.info(output)
    except subprocess.CalledProcessError as cpe:
        log_failed_subprocess(cpe)
        if raise_error:
            raise
        return cpe.output, cpe.returncode
    return output, 0


class ListToSetYamlLoader(yaml.Loader):
    """Custom YAML loader to convert lists to sets."""

    def construct_sequence(  # pyright: ignore[reportIncompatibleMethodOverride]
        self, node, deep=False
    ) -> set[str]:
        """Convert YAML sequence to a set."""
        return set(super().construct_sequence(node, deep))


ListToSetYamlLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_SEQUENCE_TAG,
    ListToSetYamlLoader.construct_sequence,
)


# pylint: disable=too-few-public-methods
class MockHandler:
    """Handler for MockLogger."""

    def __init__(self, filename):
        self.baseFilename = filename  # pylint: disable=invalid-name


# pylint: disable=too-few-public-methods
class MockLogger:
    """Mock logging.Logger to provide API without keeping the logger in memory."""

    def __init__(
        self, name: str, filename: str, level: LogLevel, log_dir: Path | str
    ) -> None:
        """Initialize a mock logger."""
        self.name = name
        self.level = level
        self.handlers = [MockHandler(os.path.join(log_dir, filename))]
        MOCK_LOGGERS[name] = self
        for loglevel in ["debug", "info", "warning", "error", "critical"]:
            # set up log methods for all built-in levels
            setattr(self, loglevel, self._factory_log(loglevel))

    def exception(self, msg, *args, exc_info=True, **kwargs):  # noqa: D102
        return self.error(msg, *args, exc_info=exc_info, **kwargs)

    exception.__doc__ = logging.exception.__doc__

    def _factory_log(self, level):
        r"""Generate a log method like `self.log(message)` for a given built-in level."""

        @docstring_parameter(level=level)
        def _log(message, *items, exc_info=False):
            """Log a message if logging level >= {level}. See `Logging Levels <https://docs.python.org/3/library/logging.html#levels>`_ for a list of levels."""
            if self.level == 0 or self.level >= getattr(
                logging, level.upper(), logging.NOTSET
            ):
                with open(
                    MockLogger._get_first_file_handler(self.handlers).baseFilename,
                    "a",
                    encoding="utf-8",
                ) as log_file:
                    if exc_info and isinstance(message, Exception):
                        value, traceback = sys_exc_info()[1:]
                        print_exception(
                            _lazy_sub(message, *items),
                            value=value,
                            tb=traceback,
                            file=log_file,
                        )
                    else:
                        print(_lazy_sub(message, *items), file=log_file)

        return _log

    def delete(self):
        """Delete the mock logger from memory."""
        del MOCK_LOGGERS[self.name]

    @staticmethod
    def _get_first_file_handler(
        handlers: Sequence[logging.Handler | MockHandler],
    ) -> Optional[logging.FileHandler | MockHandler]:
        """Given a list of Handlers, return the first FileHandler found or return None."""
        for handler in handlers:
            if isinstance(handler, (logging.FileHandler, MockHandler)):
                return handler
        return None

    def yaml_contents(self) -> dict:
        """If the logger's first handler is a YAML file, return the contents and delete them from the logger."""
        file = self._get_first_file_handler(self.handlers)
        if hasattr(file, "baseFilename"):
            file = Path(getattr(file, "baseFilename"))
            if file.suffix == ".yml":
                with file.open("r", encoding="utf-8") as f:
                    contents = yaml.load(f.read(), Loader=ListToSetYamlLoader)
                with file.open("w", encoding="utf-8") as f:
                    f.write("")
                return contents
            error = TypeError
            msg = f"Could not load YAML contents from {file}"
        else:
            error = FileNotFoundError
            msg = f"Could not find file handler for {self.name}"
        raise error(msg)


def _lazy_sub(message, *items):
    """Given lazy-logging syntax, return string with substitutions.

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


def set_up_logger(
    name: str,
    filename: Optional[str] = None,
    level: Optional[LogLevel] = None,
    log_dir: Optional[Path | str] = None,
    mock: bool = False,
) -> logging.Logger | MockLogger:
    r"""Initialize a logger.

    Parameters
    ----------
    name
        logger name (for subsequent calls to ``logging.getLogger``) to
        write to the same log file)

    filename
        filename to write log to. If not specified, filename will be
        the same as ``name`` with the extension ``log``

    level
        https://docs.python.org/3/library/logging.html#levels

    log_dir

    mock
        if ``True``, return a ``CPAC.utils.monitoring.MockLogger``
        instead of a ``logging.Logger``

    Returns
    -------
    logger
        initialized logger

    Examples
    --------
    >>> lg = set_up_logger('test')
    >>> MockLogger._get_first_file_handler(lg.handlers).baseFilename.split('/')[-1]
    'test.log'
    >>> lg.level
    0
    >>> lg = set_up_logger('second_test', 'specific_filename.custom', 'debug')
    >>> MockLogger._get_first_file_handler(lg.handlers).baseFilename.split('/')[-1]
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
        filename = f"{name}.log"
    if isinstance(level, str):
        try:
            level = getattr(logging, level.upper())
        except AttributeError:
            pass
    if not level:
        level = logging.NOTSET
    log_dir = Path(log_dir) if log_dir else Path.cwd()
    filepath = log_dir / filename
    if not filepath.exists():
        filepath.parent.mkdir(parents=True, exist_ok=True)
    if mock:
        return MockLogger(name, filename, level, log_dir)
    logger = getLogger(name)
    if isinstance(logger, MockLogger):
        return logger
    logger.setLevel(level)
    handler = logging.FileHandler(filepath)
    logger.addHandler(handler)
    return logger


def init_loggers(
    subject_id: str,
    cpac_config: "Configuration",
    log_dir: str,
    mock: bool = True,
    longitudinal: bool = False,
) -> None:
    """Set up and configure loggers."""
    from CPAC.utils.datasource import bidsier_prefix

    if "subject_id" not in cpac_config:
        cpac_config["subject_id"] = subject_id
    set_up_logger(
        f"{cpac_config['subject_id']}_expectedOutputs",
        filename=f"{bidsier_prefix(cpac_config['subject_id'])}_expectedOutputs.yml",
        level="info",
        log_dir=log_dir,
        mock=mock,
    )

    if cpac_config["pipeline_setup", "Debugging", "verbose"]:
        set_up_logger("CPAC.engine", level="debug", log_dir=log_dir, mock=True)

    nipype_config.update_config(
        {
            "logging": {
                "log_directory": log_dir,
                "log_to_file": bool(
                    getattr(
                        cpac_config["pipeline_setup", "log_directory"],
                        "run_logging",
                        True,
                    )
                ),
            },
            "execution": {
                "crashfile_format": "txt",
                "resource_monitor_frequency": 0.2,
                "stop_on_first_crash": cpac_config[
                    "pipeline_setup", "system_config", "fail_fast"
                ],
            },
        }
    )

    nipype_config.enable_resource_monitor()

    nipype_logging.update_logging(nipype_config)
