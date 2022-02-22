"""Test to check if all expected outputs were generated"""
import os
from logging import getLogger, Logger

import yaml

from CPAC.pipeline.engine import ExpectedOutputs
from CPAC.utils.monitoring.custom_logging import set_up_logger


def check_outputs(output_dir):
    """Check if all expected outputs were generated

    Parameters
    ----------
    output_dir : str
        Path to the output directory

    Returns
    -------
    message :str
    """
    outputs_log = getLogger('expected_outputs')
    missing_outputs = ExpectedOutputs()
    if isinstance(outputs_log, Logger) and len(outputs_log.handlers):
        outputs_log = getattr(outputs_log.handlers[0], 'baseFilename', None)
    else:
        outputs_log = None
    if outputs_log is None:
        message = 'Could not find expected outputs log file'
    else:
        with open(outputs_log, 'r') as f:
            expected_outputs = yaml.safe_load(f)
        for subdir in expected_outputs:
            for filename in expected_outputs['subdir']:
                if not os.path.exists(os.path.join(output_dir, subdir,
                                      filename)):
                    missing_outputs += (subdir, filename)
        if len(missing_outputs):
            missing_yaml = yaml.dump(missing_outputs)
            missing.log.error(missing_yaml)
            missing_log = set_up_logger('missing_outputs')
            try:
                log_note = 'Missing outputs have been logged in ' \
                           f'{missing_log.handlers[0].baseFilename}'
            except (AttributeError, IndexError):
                log_note = ''
            message = '\n'.join([string for string in [
                'Missing expected outputs:', missing_yaml, log_note
            ] if len(string)])
        else:
            message = 'All expected outputs were generated'
    return message
