"""Test to check if all expected outputs were generated"""
import fnmatch
import os
from logging import Logger

import yaml

from CPAC.utils.monitoring.custom_logging import getLogger, MockLogger, \
                                                 set_up_logger


def check_outputs(output_dir, log_dir, pipe_name, unique_id):
    """Check if all expected outputs were generated

    Parameters
    ----------
    output_dir : str
        Path to the output directory for the participant pipeline

    log_dir : str
        Path to the log directory for the participant pipeline

    pipe_name : str

    unique_id : str

    Returns
    -------
    message :str
    """
    outputs_logger = getLogger(f'expected_outputs_{unique_id}')
    missing_outputs = ExpectedOutputs()
    container = os.path.join(f'cpac_{pipe_name}', unique_id)
    if (
        isinstance(outputs_logger, (Logger, MockLogger)) and
        len(outputs_logger.handlers)
    ):
        outputs_log = getattr(outputs_logger.handlers[0], 'baseFilename', None)
    else:
        outputs_log = None
    if outputs_log is None:
        message = 'Could not find expected outputs log file'
    else:
        with open(outputs_log, 'r') as expected_outputs_file:
            expected_outputs = yaml.safe_load(expected_outputs_file.read())
        for subdir, filenames in expected_outputs.items():
            observed_outputs = os.listdir(
                os.path.join(output_dir, container, subdir))
            for filename in filenames:
                if not fnmatch.filter(observed_outputs,
                                      f'*{unique_id}*'
                                      f'{filename.replace(unique_id, "")}*'):
                    missing_outputs += (subdir, filename)
        if missing_outputs:
            missing_log = set_up_logger(f'missing_outputs_{unique_id}',
                                        filename='missing_outputs.yml',
                                        level='info', log_dir=log_dir,
                                        mock=True)
            missing_log.info(missing_outputs)
            try:
                log_note = 'Missing outputs have been logged in ' \
                           f'{missing_log.handlers[0].baseFilename}'
            except (AttributeError, IndexError):
                log_note = ''
            message = '\n'.join([string for string in [
                'Missing expected outputs:', str(missing_outputs), log_note
            ] if string])
            missing_log.delete()
        else:
            message = 'All expected outputs were generated'
    outputs_logger.delete()
    return message


class ExpectedOutputs:
    r'''Class to hold expected outputs for a pipeline

    Attributes
    ----------
    expected_outputs : dict
        dictionary of expected output subdirectories and files therein

    Methods
    -------
    add(subdir, output)
        Add an expected output to the expected outputs dictionary

    Examples
    --------
    >>> expected_outputs = ExpectedOutputs()
    >>> expected_outputs.add('anat', 'T1w')
    >>> expected_outputs.add('anat', 'T1w')  # shouldn't be added again
    >>> expected_outputs.add('func', 'task-rest_bold.nii.gz')
    >>> expected_outputs.add('func', 'desc-preproc_bold.json')
    >>> dict(expected_outputs)['anat']
    ['T1w']
    >>> dict(expected_outputs)['func']
    ['desc-preproc_bold.json', 'task-rest_bold.nii.gz']
    >>> str(expected_outputs)
    'anat:\n- T1w\nfunc:\n- desc-preproc_bold.json\n- task-rest_bold.nii.gz\n'
    >>> expected_outputs
    anat:
    - T1w
    func:
    - desc-preproc_bold.json
    - task-rest_bold.nii.gz
    >>> len(expected_outputs)
    3
    '''
    def __init__(self):
        self.expected_outputs = {}

    def __bool__(self):
        return bool(len(self))

    def __iter__(self):
        yield from {subdir: sorted(list(filename)) for
                    subdir, filename in self.expected_outputs.items()}.items()

    def __iadd__(self, other):
        if not isinstance(other, tuple) or not len(other) == 2:
            raise TypeError(
                f'{self.__module__}.{self.__class__.__name__} requires a '
                "tuple of ('subdir', 'output') for addition")
        self.add(*other)
        return self

    def __len__(self):
        return len([filepath for subdir, filepaths in
                    self.expected_outputs.items() for filepath in filepaths])

    def __repr__(self):
        return str(self).rstrip()

    def __str__(self):
        return yaml.dump(dict(self))

    def add(self, subdir, output):
        '''Add an expected output to the expected outputs dictionary

        Parameters
        ----------
        subdir : str
            subdirectory of expected output

        output : str
            filename of expected output
        '''
        if subdir in self.expected_outputs:
            self.expected_outputs[subdir].add(output)
        else:
            self.expected_outputs[subdir] = {output}
