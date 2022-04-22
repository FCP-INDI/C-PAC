'''Validation schema for C-PAC pipeline configurations'''
# pylint: disable=too-many-lines
from voluptuous.error import Error as VoluptuousError

from CPAC.utils import coerce_to_list
from CPAC.utils.utils import lookup_nested_value, set_nested_value
from .backwards_compatibility import backwards_compatible
from .config import ALWAYS_LISTS, latest_schema, TOGGLED_OPTIONS
from .exceptions import handle_custom_error


def schema(config_dict):
    '''Function to test the schema validity of a given config dict.

    Parameters
    ----------
    config_dict : dict

    Returns
    -------
    voluptuous.schema_builder.Schema
        validated configuration schema
    '''
    for keys in ALWAYS_LISTS:
        config_dict = coerce_to_list(config_dict, keys)
    config_dict = backwards_compatible(config_dict)
    for option in TOGGLED_OPTIONS:
        # Allow options to be mutually incompatible if those options are off
        switch = lookup_nested_value(config_dict, option['switch'])
        if True not in switch:
            latest_schema.schema = set_nested_value(latest_schema.schema,
                                                    option['key'],
                                                    option['Off'])
    try:
        return latest_schema(config_dict)
    except VoluptuousError as voluptuous_error:
        return latest_schema(
            handle_custom_error(voluptuous_error, config_dict))


schema.schema = latest_schema.schema
