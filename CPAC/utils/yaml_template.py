import re
import yaml
from CPAC.utils.configuration import DEFAULT_PIPELINE_FILE
from CPAC.utils.utils import update_nested_dict


def create_yaml_from_template(d, template=DEFAULT_PIPELINE_FILE):
    """Save dictionary to a YAML file, keeping the structure
    (such as first level comments and ordering) from the template

    It may not be fully robust to YAML structures, but it works
    for C-PAC config files!

    Parameters
    ----------
    d: dict

    template: str
        path to template

    Examples
    --------
    >>> import json, yaml
    >>> from CPAC.utils.configuration import DEFAULT_PIPELINE_FILE
    >>> json.dumps(yaml.safe_load(
    ...     create_yaml_from_template({})), sort_keys=True
    ... ) == json.dumps(yaml.safe_load(
    ...     open(DEFAULT_PIPELINE_FILE, 'r')), sort_keys=True)
    True
    """
    def _count_indent(line):
        '''Helper method to determine indentation level

        Parameters
        ----------
        line: str

        Returns
        -------
        number_of_indents: int

        Examples
        --------
        >>> _count_indent('No indent')
        0
        >>> _count_indent('    Four spaces')
        2
        '''
        return (len(line) - len(line.lstrip())) // 2

    def _format_key(key, level):
        '''Helper method to format YAML keys

        Parameters
        ----------
        key: str
        level: int

        Returns
        -------
        yaml: str

        Examples
        --------
        >>> _format_key('base', 0)
        '\nbase: '
        >>> _format_key('indented', 2)
        '\n    indented:'
        '''
        return f'\n{" " * level * 2}{key}: '

    def _format_list_items(l, line_level):
        '''Helper method to handle lists in the YAML

        Parameters
        ----------
        l: list
        
        line_level: int
        
        Returns
        -------
        yaml: str
        
        Examples
        --------
        >>> _format_list_items([1, 2, {'nested': 3}], 0)
        '  - 1\n  - 2\n  - nested: 3'
        >>> _format_list_items([1, 2, {'nested': [3, {'deep': [4]}]}], 1)
        '    - 1\n    - 2\n    - nested:\n      - 3\n      - deep:\n        - 4'
        '''  # noqa
        indent = " " * (2 * line_level + 2)
        return '\n'.join([
            f'{indent}{li}' for li in yaml.dump(l).split('\n')
        ]).rstrip()

    def _lookup_value(d, keys):
        '''Helper method to look up nested values

        Paramters
        ---------
        d: dict
        keys: list or tuple

        Returns
        -------
        yaml: str or dict

        Examples
        --------
        >>> _lookup_value({'nested': {'True': True}}, ['nested', 'True'])
        'On'
        >>> _lookup_value({'nested': {'None': None}}, ['nested', 'None'])
        ''
        '''
        if len(keys) == 1:
            value = d.get(keys[0])
            if value is None:
                return ''
            if isinstance(value, bool):
                if value == False:  # noqa E712
                    return 'Off'
                return 'On'
            return d.get(keys[0])
        else:
            return _lookup_value(d.get(keys[0], {}), keys[1:])

    # set starting values
    output = ''
    comment = ''
    space_match = r'^\s+.*'
    level = 0
    nest = []
    list_item = False
    list_level = 0
    line_level = 0

    # load default values
    d_default = yaml.safe_load(open(template, 'r'))

    # update values
    d = update_nested_dict(d_default, d)

    # generate YAML from template with updated values
    with open(template, 'r') as f:
        for line in f:

            # persist comments and frontmatter
            if line.startswith('%') or line.startswith('---') or re.match(
                r'^\s*#.*$', line
            ):
                list_item = False
                line = line.strip('\n')
                comment += f'\n{line}'
            elif len(line.strip()):
                if re.match(space_match, line):
                    line_level = _count_indent(line)
                else:
                    line_level = 0

                # handle lists as a unit
                if list_item:
                    if line_level < list_level - 1:
                        list_item = False
                        level = list_level
                        list_level = 0
                elif line.lstrip().startswith('-'):
                    list_item = True
                    list_level = line_level - 1

                else:

                    # extract dict key
                    key_group = re.match(
                        r'^\s*([a-z0-9A-Z_/][\sa-z0-9A-Z_/\.-]+)\s*:', line)
                    if key_group:
                        key = key_group.group(1).strip()

                        # calculate key depth
                        if line_level == level:
                            if level > 0:
                                nest = nest[:-1] + [key]
                            else:
                                nest = [key]
                        elif line_level == level + 1:
                            nest += [key]
                        elif line_level < level:
                            nest = nest[:line_level] + [key]

                        # get updated value for key
                        value = _lookup_value(d, nest)

                        # prepend comment from template
                        if len(comment.strip()):
                            output += comment

                        # write YAML
                        output += _format_key(key, line_level)
                        if isinstance(value, list):
                            output += '\n' + _format_list_items(
                                value, line_level)
                        elif not isinstance(value, dict):
                            output += str(value)

                        # reset variables for loop
                        comment = '\n'
                        level = line_level

    return output.lstrip('\n')
