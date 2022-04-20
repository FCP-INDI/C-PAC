'''Module to support backwards compaitibility in schema validation.'''
from nipype import logging
from voluptuous import Any, Required, Schema
from voluptuous.validators import ExactSequence, Maybe
from CPAC.utils.docs import docstring_parameter
from CPAC.utils.utils import dct_diff, delete_nested_value, \
                             lookup_nested_value, set_nested_value
from .constants import Number, OVERRIDABLE_DEFAULTS

logger = logging.getLogger('nipype.workflow')


class UpgradingDict(dict):
    '''Custom dictionary class to track if a dictionary has been updated.'''
    def __init__(self, *args):
        self.updated = False
        if args:
            self.data = args[0]
        dict.__init__(self, *args)

    def __setitem__(self, key, val):
        if dct_diff(self.data.get(key), val):
            self.updated = True
        dict.__setitem__(self, key, val)


def backwards_compatible(config_dict):
    '''Fucntion to apply backwards compatibility patches to a config
       dict for validation.

    Parameters
    ----------
    config_dict : dict

    Returns
    -------
    dict
    '''
    return kick_forward(config_dict)


def _kick_1_8_0_to_1_8_1(config_dict):
    '''Update config from 1.8.0 to 1.8.1'''
    for key_sequence in {
        ('anatomical_preproc', 'non_local_means_filtering'),
        ('anatomical_preproc', 'n4_bias_field_correction')
    }:
        config_dict = _now_runswitch(config_dict, key_sequence)
        for combiners in {
            ((
                ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                 'CSF_label'),
            ), ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                'CSF_label')),
            ((
                ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                 'left_GM_label'),
                ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                 'right_GM_label')
            ), ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                'GM_label')),
            ((
                ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                 'left_WM_label'),
                ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                 'right_WM_label')
            ), ('segmentation', 'tissue_segmentation', 'ANTs_Prior_Based',
                'WM_label'))
        }:
            config_dict = _combine_labels(config_dict, *combiners)
        try:
            calculate_motion_first = lookup_nested_value(
                config_dict,
                ['functional_preproc', 'motion_estimates_and_correction',
                 'calculate_motion_first']
            )
        except KeyError:
            calculate_motion_first = None
        if calculate_motion_first is not None:
            del config_dict['functional_preproc'][
                'motion_estimates_and_correction'][
                'calculate_motion_first']
            config_dict = set_nested_value(config_dict, [
                'functional_preproc', 'motion_estimates_and_correction',
                'motion_estimates', 'calculate_motion_first'
            ], calculate_motion_first)
    return config_dict


def _kick_1_8_4_to_1_8_5(config_dict):
    '''Update config from 1.8.4 to 1.8.5'''
    def _set_name(motion_estimate_filter):
        if isinstance(motion_estimate_filter, list):
            return [_set_name(m_e_filter) for m_e_filter in
                    motion_estimate_filter]
        if (isinstance(motion_estimate_filter, dict) and
                ('name' not in motion_estimate_filter or
                 motion_estimate_filter['name'] is None)):
            motion_estimate_filter['name'] = 'motion_estimate_filter'
        return motion_estimate_filter
    motion_estimate_filter = {
        'keys': ['functional_preproc', 'motion_estimates_and_correction',
                 'motion_estimate_filter', 'filters']}
    try:
        motion_estimate_filter['value'] = lookup_nested_value(
            config_dict, motion_estimate_filter['keys'])
    except KeyError:
        motion_estimate_filter['value'] = None
    if motion_estimate_filter['value'] is not None:
        if ('name' not in motion_estimate_filter['value'] or
                motion_estimate_filter['value']['name'] is None):
            config_dict = set_nested_value(
                config_dict, motion_estimate_filter['keys'],
                _set_name(motion_estimate_filter['value']))
    return config_dict


_kick = {
    '1.8.0 to 1.8.1': _kick_1_8_0_to_1_8_1,
    '1.8.4 to 1.8.5': _kick_1_8_4_to_1_8_5
}


@docstring_parameter(available_versions=str(list(_kick.keys()) + ['all']))
def kick_forward(config_dict, versions='all'):
    '''Function to apply backwards compatibility patches to a config

    Parameters
    ----------
    config_dict : dict

    versions : str
        versions of C-PAC to apply the patches for.
        Available options: {available_versions}
    '''
    if versions == 'all':
        for version, fxn in _kick.items():
            config_dict = fxn(UpgradingDict(config_dict))
            if config_dict.updated:
                logger.info('Upgraded config from C-PAC %s', version)
    else:
        config_dict = _kick[versions](config_dict)
        logger.info('Upgraded config from C-PAC %s', versions)
    return dict(config_dict)

# Adding examples here to avoid having to escape all the curly braces for the
# docstring_parameter decorator.
kick_forward.__doc__ += r'''
    Examples
    --------
    Starting with 1.8.0
    >>> zero = {'anatomical_preproc': {
    ...     'non_local_means_filtering': True,
    ...     'n4_bias_field_correction': True
    ... }, 'functional_preproc': {
    ...     'motion_estimates_and_correction': {
    ...         'calculate_motion_first': False
    ...     }
    ... }, 'segmentation': {
    ...     'tissue_segmentation': {
    ...         'ANTs_Prior_Based': {
    ...             'CSF_label': 0,
    ...             'left_GM_label': 1,
    ...             'right_GM_label': 2,
    ...             'left_WM_label': 3,
    ...             'right_WM_label': 4}}}}
    >>> updated_apb = kick_forward(zero)[
    ...     'segmentation']['tissue_segmentation']['ANTs_Prior_Based']
    >>> updated_apb['CSF_label']
    [0]
    >>> updated_apb['GM_label']
    [1, 2]
    >>> updated_apb['WM_label']
    [3, 4]

    Starting with 1.8.1
    >>> one = {'anatomical_preproc': {
    ...     'non_local_means_filtering': True,
    ...     'n4_bias_field_correction': True
    ... }, 'functional_preproc': {
    ...     'motion_estimates_and_correction': {
    ...         'calculate_motion_first': False
    ...     }
    ... }, 'segmentation': {
    ...     'tissue_segmentation': {
    ...         'ANTs_Prior_Based': {
    ...             'CSF_label': [0],
    ...             'GM_label': [1, 2],
    ...             'WM_label': [3, 4]}}}}
    >>> updated_apb = kick_forward(one)[
    ...     'segmentation']['tissue_segmentation']['ANTs_Prior_Based']
    >>> updated_apb['CSF_label']
    [0]
    >>> updated_apb['GM_label']
    [1, 2]
    >>> updated_apb['WM_label']
    [3, 4]
    '''


def _combine_labels(config_dict, list_to_combine, new_key):
    '''
    Helper function to combine formerly separate keys into a
    combined key.

    Parameters
    ----------
    config_dict: dict

    key_sequence: iterable of lists or tuples

    new_key: list or tuple

    Returns
    -------
    updated_config_dict: dict
    '''
    new_value = []
    any_old_values = False
    for _to_combine in list_to_combine:
        try:
            old_value = lookup_nested_value(config_dict, _to_combine)
        except KeyError:
            old_value = None
        if old_value is not None:
            any_old_values = True
            if isinstance(old_value, (list, set, tuple)):
                for value in old_value:
                    new_value.append(value)
            else:
                new_value.append(old_value)
            config_dict = delete_nested_value(config_dict, _to_combine)
    if any_old_values:
        return set_nested_value(config_dict, new_key, new_value)
    return config_dict


deprecated_syntaxes = {
    'motion_estimate_filter': {
        'pre-1.8.5': Schema(Required(
            Any(
              {  # no motion estimate filter
               **OVERRIDABLE_DEFAULTS['motion_estimate_filter'],
               'run': Maybe(Any(ExactSequence([False]), ExactSequence([]),
                            False))},
              {  # notch filter with breathing_rate_*
               **OVERRIDABLE_DEFAULTS['motion_estimate_filter'],
               'filter_type': 'notch',
               'filter_order': int,
               'breathing_rate_min': Number,
               'breathing_rate_max': Number},
              {  # notch filter with manual parameters
               **OVERRIDABLE_DEFAULTS['motion_estimate_filter'],
               'filter_type': 'notch',
               'filter_order': int,
               'breathing_rate_min': None,
               'breathing_rate_max': None,
               'center_frequency': Number,
               'filter_bandwidth': Number},
              {  # lowpass filter with breathing_rate_min
               **OVERRIDABLE_DEFAULTS['motion_estimate_filter'],
               'filter_type': 'lowpass',
               'filter_order': int,
               'breathing_rate_min': Number},
              {  # lowpass filter with lowpass_cutoff
               'filter_type': 'lowpass',
               'filter_order': int,
               'breathing_rate_min': None,
               'lowpass_cutoff': Number,
              })
        ))
    }
}


def _now_runswitch(config_dict, key_sequence):
    '''
    Helper function to convert a formerly forkable value to a
    runswitch.

    Parameters
    ----------
    config_dict : dict

    key_sequence : list or tuple

    Returns
    -------
    updated_config_dict : dict
    '''
    try:
        old_forkable = lookup_nested_value(config_dict, key_sequence)
    except KeyError:
        return config_dict
    if isinstance(old_forkable, (bool, list)):
        return set_nested_value(
            config_dict, key_sequence, {'run': old_forkable})
    return config_dict
