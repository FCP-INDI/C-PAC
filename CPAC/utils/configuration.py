import os
import warnings

from CPAC.pipeline.schema import schema


class Configuration(object):
    """Class to set dictionary keys as map attributes.

    Given a Configuration `c`, and a list or tuple of an attribute name
    and nested keys `keys = ['attribute', 'key0', 'key1']` or
    `keys = ('attribute', 'key0', 'key1')`, the value 'value' nested in

    c.attribute = {'key0': {'key1': 'value'}}

    can be accessed (get and set) in any of the following ways (and
    more):

    c.attribute['key0']['key1']
    c['attribute']['key0']['key1']
    c.attribute['key0', 'key1']
    c['attribute', 'key0', 'key1']
    c.attribute[keys[1:]]
    c[keys]
    """
    def __init__(self, config_map):
        config_map = schema(self.nonestr_to_None(config_map))
        for key in config_map:
            # set FSLDIR to the environment $FSLDIR if the user sets it to
            # 'FSLDIR' in the pipeline config file
            if key == 'FSLDIR':
                if config_map[key] == 'FSLDIR':
                    if os.environ.get('FSLDIR'):
                        config_map[key] = os.environ['FSLDIR']
            setattr(self, key, config_map[key])
        self.__update_attr()

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)
        elif isinstance(key, tuple) or isinstance(key, list):
            return self.get_nested(self, key)
        else:
            self.key_type_error(key)

    def __setitem__(self, key, value):
        if isinstance(key, str):
            setattr(self, key, value)
        elif isinstance(key, tuple) or isinstance(key, list):
            self.set_nested(self, key, value)
        else:
            self.key_type_error(key)

    def nonestr_to_None(self, d):
        """
        recursive method to type convert 'None' to None in nested
        config

        Parameters
        ----------
        d: any
            config item to check

        Returns
        -------
        d: any
            same item, same type, but with 'none' strings converted to
            Nonetypes
        """
        if isinstance(d, str) and d.lower() == 'none':
            return None
        elif isinstance(d, list):
            return [self.nonestr_to_None(i) for i in d]
        elif isinstance(d, set):
            return {self.nonestr_to_None(i) for i in d}
        elif isinstance(d, dict):
            return {i: self.nonestr_to_None(d[i]) for i in d}
        else:
            return d

    def return_config_elements(self):
        # this returns a list of tuples
        # each tuple contains the name of the element in the yaml config file
        # and its value
        attributes = [
            (attr, getattr(self, attr))
            for attr in dir(self)
            if not callable(attr) and not attr.startswith("__")
        ]
        return attributes

    # method to find any pattern ($) in the configuration
    # and update the attributes with its pattern value
    def update_attr(self):
        from string import Template

        # TODO remove function from here
        def check_pattern(orig_key, tags=None):
            temp = Template(orig_key)
            if type(temp.template) is str:
                patterns = temp.pattern.findall(temp.template)
                pattern_map = {}
                keep_tags_map = {}
                for pattern in patterns:
                    for val in [_f for _f in pattern if _f]:
                        if val:
                            if not pattern_map.get(val):
                                if tags is not None and val not in tags:
                                    keep_tags_map[val] = '${' + val + '}'
                                    continue
                                if getattr(self, val):
                                    pattern_map[val] = getattr(self, val)
                                else:
                                    raise ValueError(
                                        f"No value found for attribute %s. "
                                        "Please check the configuration "
                                        "file" % val)
                if pattern_map:
                    pattern_map.update(keep_tags_map)
                    return check_pattern(
                        Template(orig_key).substitute(pattern_map), tags=tags)
                else:
                    return orig_key
            else:
                return orig_key

        def check_path(key):
            if type(key) is str and '/' in key:
                if not os.path.exists(key):
                    warnings.warn(
                        "Invalid path- %s. Please check your configuration "
                        "file" % key)

        attributes = [(attr, getattr(self, attr)) for attr in dir(self)
                      if not callable(attr) and not attr.startswith("__")]

        template_list = ['template_brain_only_for_anat',
                         'template_skull_for_anat',
                         'ref_mask',
                         'template_brain_only_for_func',
                         'template_skull_for_func',
                         'template_symmetric_brain_only',
                         'template_symmetric_skull',
                         'dilated_symmetric_brain_mask']

        for attr_key, attr_value in attributes:

            if attr_key in template_list:
                new_key = check_pattern(attr_value, 'FSLDIR')
            else:
                new_key = check_pattern(attr_value)
            setattr(self, attr_key, new_key)

    __update_attr = update_attr

    def update(self, key, val):
        setattr(self, key, val)

    def __copy__(self):
        newone = type(self)({})
        newone.__dict__.update(self.__dict__)
        newone.__update_attr()
        return newone

    def get_nested(self, d, keys):
        if isinstance(keys, str):
            return d[keys]
        elif isinstance(keys, tuple) or isinstance(keys, list):
            if len(keys) > 1:
                return self.get_nested(d[keys[0]], keys[1:])
            else:
                return d[keys[0]]

    def set_nested(self, d, keys, value):
        if isinstance(keys, str):
            d[keys] = value
        elif isinstance(keys, tuple) or isinstance(keys, list):
            if len(keys) > 1:
                d[keys[0]] = self.set_nested(d[keys[0]], keys[1:], value)
            else:
                d[keys[0]] = value
        return d

    def key_type_error(self, key):
        raise KeyError(' '.join([
                'Configuration key must be a string, list, or tuple;',
                type(key).__name__,
                f'`{str(key)}`',
                'was given.'
            ]))
