import os
import re
import yaml
from optparse import OptionError
from CPAC.utils.configuration import Configuration, DEFAULT_PIPELINE_FILE
from CPAC.utils.utils import dct_diff, load_preconfig, lookup_nested_value


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
    >>> import yaml
    >>> from CPAC.utils.configuration import Configuration
    >>> Configuration(yaml.safe_load(create_yaml_from_template({}))).dict(
    ...    ) == Configuration({}).dict()
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

    def _create_import_dict(diff):
        '''Method to return a dict of only changes given a nested dict
        of (dict1_value, dict2_value) tuples

        Parameters
        ----------
        diff: dict
            output of `dct_diff`

        Returns
        -------
        dict
            dict of only changed values

        Examples
        --------
        >>> _create_import_dict({'anatomical_preproc': {'brain_extraction': {'extraction': {'using': (['3dSkullStrip'], ['niworkflows-ants'])}}}})
        {'anatomical_preproc': {'brain_extraction': {'extraction': {'using': ['niworkflows-ants']}}}}
        '''  # noqa
        if isinstance(diff, tuple) and len(diff) == 2:
            return diff[1]
        if isinstance(diff, dict):
            i = {}
            for k in diff:
                try:
                    i[k] = _create_import_dict(diff[k])
                except KeyError:
                    continue
            return {k: i[k] for k in i if i[k]}
        return diff

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

    # set starting values
    output = ''
    comment = ''
    space_match = r'^\s+.*'
    level = 0
    nest = []
    list_item = False
    list_level = 0
    line_level = 0
    template_name = template
    try:
        template = load_preconfig(template)
    except OptionError:
        if 'default' in template.lower():
            template = DEFAULT_PIPELINE_FILE
        assert os.path.exists(template) or os.path.islink(template), \
            f'{template_name} is not a defined preconfig or a valid path.'
    template_included = False

    # load default values
    d_default = Configuration(yaml.safe_load(open(template, 'r'))).dict()

    if (
        template == DEFAULT_PIPELINE_FILE or
        not dct_diff(
            yaml.safe_load(open(DEFAULT_PIPELINE_FILE, 'r')), d_default)
    ):
        template_name = 'default'

    # update values
    d = _create_import_dict(dct_diff(d_default, d))

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
                        if not template_included:
                            # prepend comment from template
                            if len(comment.strip()):
                                comment = re.sub(
                                    r'(?<=# based on )(.* pipeline)',
                                    f'{template_name} pipeline',
                                    comment
                                )
                                output += comment
                                output += f'\nFROM: {template_name}\n'
                                comment = ''
                            template_included = True
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

                        # only include updated and new values
                        try:
                            # get updated value for key
                            value = lookup_nested_value(d, nest)

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
                        except KeyError:
                            # clear comment for excluded key
                            comment = '\n'

                        # reset variables for loop
                        comment = '\n'
                        level = line_level
            else:
                comment += '\n'

    # update to include new values
    new_keys = _create_import_dict(dct_diff(yaml.safe_load(output), d))
    if new_keys:
        output += '\n'
        output += yaml.dump(new_keys).strip()

    return output.lstrip('\n')
