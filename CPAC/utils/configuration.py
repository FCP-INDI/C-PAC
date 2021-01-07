import os
import warnings

"""
Class to set dictionary keys as map attributes
"""
class Configuration(object):

    def __init__(self, config_map):
        for key in config_map:
            if config_map[key] == 'None':
                config_map[key] = None
            if isinstance(config_map[key], dict):
                for subkey in config_map[key]:
                    if config_map[key][subkey] == 'None':
                        config_map[key][subkey] = None
            # set FSLDIR to the environment $FSLDIR if the user sets it to
            # 'FSLDIR' in the pipeline config file
            if key == 'FSLDIR':
                if config_map[key] == 'FSLDIR':
                    if os.environ.get('FSLDIR'):
                        config_map[key] = os.environ['FSLDIR']
            setattr(self, key, config_map[key])
        self.__update_attr()

    def return_config_elements(self):
        # this returns a list of tuples
        # each tuple contains the name of the element in the yaml config file
        # and its value

        # TODO ASH use __dict__ to retrieve elements
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
        def check_pattern(orig_key, tags = None):
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
                                    raise ValueError("No value found for attribute %s. "\
                                                    "Please check the configuration file" %val)
                if pattern_map:
                    pattern_map.update(keep_tags_map)
                    return check_pattern(Template(orig_key).substitute(pattern_map), tags=tags)
                else:
                    return orig_key
            else:
                return orig_key
                
        def check_path(key):
            if type(key) is str and '/' in key:
                if not os.path.exists(key):
                    warnings.warn("Invalid path- %s. Please check your configuration file"%key)
            
        attributes = [(attr, getattr(self, attr)) for attr in dir(self) \
                      if not callable(attr) and not attr.startswith("__")]     
        
        template_list = ['template_brain_only_for_anat',
                        'template_skull_for_anat',
                        'template_brain_mask_for_anat',
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
                #check_path(new_key)
            setattr(self, attr_key, new_key)

    __update_attr = update_attr
    
    def update(self, key, val):
        setattr(self, key, val)
        
    def __copy__(self):
        newone = type(self)({})
        newone.__dict__.update(self.__dict__)
        newone.__update_attr()
        return newone
