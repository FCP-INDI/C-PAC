'''Custom validation error messages'''
from nipype import logging
from voluptuous.error import NotInInvalid, RequiredFieldInvalid
from CPAC.utils.utils import lookup_nested_value, set_nested_value
from .constants import VALID_OPTIONS

logger = logging.getLogger('nipype.workflow')


def _remove_from_dict_in_list(data, subdata, path, item, keys_to_ignore,
                              message=None):
    '''Remove keys from a dict in a list

    Parameters
    ----------
    data : dict
        The data that was being validated

    subdata : dict
        The dict in the list in which the keys are being removed

    path : list
        The list of keys to the list in data

    item : int
        The index of subdata in the list

    keys_to_ignore : list
        The keys to be removed from subdata

    message : str or None
        The message to be displayed if a key is removed

    Returns
    -------
    dict
        Updated data
    '''
    if not isinstance(item, int):
        item = int(item)
    ignored = False
    for key in keys_to_ignore:
        if subdata.get(key) is not None:
            if message:
                logger.warning(message.format(key=key,
                                              path=str(path + [item])))
            ignored = True
            del subdata[key]
        elif key in subdata:
            del subdata[key]
    if ignored:
        new_data = lookup_nested_value(data, path)
        new_data[item] = subdata
        return set_nested_value(data, path, subdata)
    return data


def _handle_motion_estimate_filter(error, data, item):
    path = ['functional_preproc', 'motion_estimates_and_correction',
            'motion_estimate_filter', 'filters']
    filter_types = {filter_option['filter_type'] for option in
                    VALID_OPTIONS['motion_estimate_filter'].validators for
                    filter_option in option}
    filters = lookup_nested_value(data, path)
    m_e_f = filters[int(item)]
    m_e_filter = {key: m_e_f.get(key) for key in [
        'filter_type', 'filter_order', 'breathing_rate_min',
        'breathing_rate_max', 'lowpass_cutoff', 'center_frequency',
        'filter_bandwidth']}
    if m_e_filter['filter_type'] not in filter_types:
        raise NotInInvalid('Valid options for \'filter_type\' are ' +
                           str(filter_types) + ", but '" +
                           str(m_e_filter['filter_type']) + "' was "
                           'provided', path=path)
    if m_e_f['filter_type'] == 'lowpass':
        if m_e_filter['breathing_rate_min'] is not None:
            pass
        elif m_e_filter['lowpass_cutoff'] is not None:
            pass
        else:
            raise RequiredFieldInvalid('Lowpass filter requires either '
                                       '\'breathing_rate_min\' or '
                                       '\'lowpass_cutoff\' to be set.',
                                       path=path)
    if m_e_f['filter_type'] == 'notch':
        if (
            m_e_filter['breathing_rate_min'] is not None and
            m_e_filter['breathing_rate_max'] is not None
        ):
            return _remove_from_dict_in_list(data, m_e_filter, path, item, [
                'center_frequency', 'filter_bandwidth', 'lowpass_cutoff'],
                message='\'{key}\' is not used for notch filter when '
                        '\'breathing_rate_min\' and \'breathing_rate_max\' '
                        'are provided. Dropping from {path}.')
        elif (
            m_e_filter['center_frequency'] is not None and
            m_e_filter['filter_bandwidth'] is not None
        ):
            pass
        else:
            raise RequiredFieldInvalid('Notch filter requires either '
                                       '\'breathing_rate_min\' and '
                                       '\'breathing_rate_max\' to be set, or '
                                       '\'center_frequency\' and '
                                       '\'filter_bandwidth\' to be set.',
                                       path=path)
    raise error


HANDLED_KEYS = {
    "'functional_preproc', 'motion_estimates_and_correction', "
    "'motion_estimate_filter', 'filters'": _handle_motion_estimate_filter
}


def handle_custom_error(error, data):
    '''Handle custom error messages

    Parameters
    ----------
    error : voluptuous.error.Error

    data : dict
        The data that was being validated

    Raises
    ------
    Exception
    '''
    if hasattr(error, 'path'):
        for handled_key in HANDLED_KEYS:
            error_path_str = str(error.path)
            starting_key = f'[{handled_key}'
            if error_path_str.startswith(starting_key):
                item = error_path_str.split(f'{starting_key}, ', 1)[
                    1].split(',', 1)[0]
                return HANDLED_KEYS[handled_key](error, data, item)
    raise error
