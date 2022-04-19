'''Module to support backwards compaitibility in schema validation.'''
from voluptuous import Any, In, Required, Schema
from voluptuous.validators import ExactSequence, Maybe
from CPAC.utils.utils import delete_nested_value, lookup_nested_value, \
                             set_nested_value
from .constants import forkable, Number


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
    return _changes_1_8_0_to_1_8_1(config_dict)


def _changes_1_8_0_to_1_8_1(config_dict):
    '''
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
    >>> updated_apb = _changes_1_8_0_to_1_8_1(zero)[
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
    >>> updated_apb = _changes_1_8_0_to_1_8_1(one)[
    ...     'segmentation']['tissue_segmentation']['ANTs_Prior_Based']
    >>> updated_apb['CSF_label']
    [0]
    >>> updated_apb['GM_label']
    [1, 2]
    >>> updated_apb['WM_label']
    [3, 4]
    '''
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
            'motion_estimates_and_correction']['calculate_motion_first']
        config_dict = set_nested_value(config_dict, [
            'functional_preproc', 'motion_estimates_and_correction',
            'motion_estimates', 'calculate_motion_first'
        ], calculate_motion_first)

    return config_dict


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


OVERRIDABLE_DEFAULTS = {
    'motion_estimate_filter': {
        {'run': forkable,
         'filter_type': Maybe(In({'notch', 'lowpass'})),
         'filter_order': Maybe(int),
         'breathing_rate_min': Maybe(Number),
         'breathing_rate_max': Maybe(Number),
         'center_frequency': Maybe(Number),
         'filter_bandwidth': Maybe(Number),
         'lowpass_cutoff': Maybe(Number)}}
}


deprecated_syntaxes = {
    'motion_estimate_filter': {
        'pre-1.8.5': {Schema(Required(
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
        ))}
    }
}


def _now_runswitch(config_dict, key_sequence):
    '''
    Helper function to convert a formerly forkable value to a
    runswitch.

    Parameters
    ----------
    config_dict: dict

    key_sequence: list or tuple

    Returns
    -------
    updated_config_dict: dict
    '''
    try:
        old_forkable = lookup_nested_value(config_dict, key_sequence)
    except KeyError:
        return config_dict
    if isinstance(old_forkable, (bool, list)):
        return set_nested_value(
            config_dict, key_sequence, {'run': old_forkable})
    return config_dict
