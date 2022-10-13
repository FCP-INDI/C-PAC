"""Copyright (C) 2022  C-PAC Developers

This file is part of C-PAC.

C-PAC is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

C-PAC is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>."""
import os
import re
from datetime import datetime
from hashlib import sha1
import yaml
from CPAC.utils.configuration import Configuration, DEFAULT_PIPELINE_FILE
from CPAC.utils.utils import update_config_dict, update_pipeline_values_1_8, \
                             YAML_BOOLS

YAML_LOOKUP = {yaml_str: key for key, value in YAML_BOOLS.items() for
               yaml_str in value}


class YamlTemplate():  # pylint: disable=too-few-public-methods
    """A class to link YAML comments to the contents of a YAML file"""
    def __init__(self, original_yaml):
        """
        Parameters
        ----------
        original_yaml : str
            raw YAML or path to YAML file
        """
        if os.path.exists(original_yaml):
            with open(original_yaml, 'r', encoding='utf-8') as _f:
                original_yaml = _f.read()
        self.comments = {}
        self.original = original_yaml
        self._dict = yaml.safe_load(self.original)
        self._parse_comments()

    get_nested = Configuration.get_nested

    def dump(self, new_dict, parents=None):
        """Dump a YAML file from a new dictionary with the comments from
        the template dictionary

        Parameters
        ----------
        new_dict : dict

        parents : list of str

        Returns
        -------
        str
        """
        parents = parents if parents is not None else []
        line_level = len(parents)
        _dump = []
        for key in (self.get_nested(new_dict, parents) if
                    parents else new_dict):
            keys = [*parents, key]
            comment = self.comments.get('.'.join(keys))
            try:
                value = self.get_nested(new_dict, keys)
            except KeyError:  # exclude unincluded keys
                continue
            indented_key = f'{indent(line_level, 0)}{key}:'
            if comment:
                _dump += ['']
                _dump += [indent(line_level, 0) + line for line in comment]
            else:
                while _dump and _dump[-1].strip() == '':
                    _dump = _dump[:-1]
            if value is not None:
                if isinstance(value, dict):
                    _dump += [indented_key, self.dump(new_dict, keys)]
                elif isinstance(value, list):
                    list_line = _format_list_items(value, line_level + 1)
                    if '\n' in list_line:
                        _dump += [indented_key, *list_line.split('\n')]
                    else:
                        _dump += [f'{indented_key} {list_line}']
                elif isinstance(value, bool) or (isinstance(value, str) and
                                                 value.lower() in YAML_LOOKUP):
                    if isinstance(value, str):
                        value = YAML_LOOKUP[value.lower()]
                    value = 'On' if value is True else 'Off'
                    _dump += [f'{indented_key} {value}']
                else:
                    _dump += [f'{indented_key} {value}']
            else:
                _dump += [indented_key]
        while _dump and _dump[0] == '\n':
            _dump = _dump[1:]
        return re.sub('\n{3,}', '\n\n', '\n'.join(_dump)).rstrip() + '\n'

    def _parse_comments(self):
        yaml_lines = self.original.split('\n')
        comment = []
        key = []
        for line in yaml_lines:
            line_level = _count_indent(line)
            stripped_line = line.strip()
            if stripped_line.startswith('#'):
                comment.append(stripped_line)
            elif not any(stripped_line.startswith(seq) for
                         seq in ('%YAML', '---')):
                if ':' in stripped_line:
                    line_key = stripped_line.split(':', 1)[0]
                    if line_level == 0:
                        key = [line_key]
                    else:
                        key = [*key[:line_level], line_key]
                    self.comments['.'.join(key)] = comment
                    comment = []


def _count_indent(line):
    '''Helper method to determine indentation level

    Parameters
    ----------
    line : str

    Returns
    -------
    number_of_indents : int

    Examples
    --------
    >>> _count_indent('No indent')
    0
    >>> _count_indent('    Four spaces')
    2
    '''
    return (len(line) - len(line.lstrip())) // 2


def create_yaml_from_template(d,  # pylint: disable=invalid-name
                              template=DEFAULT_PIPELINE_FILE,
                              include_all=False):
    """Save dictionary to a YAML file, keeping the structure
    (such as first level comments and ordering) from the template

    It may not be fully robust to YAML structures, but it works
    for C-PAC config files!

    Parameters
    ----------
    d : dict or Configuration

    template : str
        path to template or YAML as a string

        .. versionchanged:: 1.8.5
           used to take path to template or name of preconfig

    include_all : bool
        include every key, even those that are unchanged

    Examples
    --------
    >>> import yaml
    >>> from CPAC.utils.configuration import Configuration, Preconfiguration
    >>> Configuration(yaml.safe_load(create_yaml_from_template({}))).dict(
    ...     ) == Configuration({}).dict()
    True
    >>> fmriprep_options = Preconfiguration('fmriprep-options')
    Loading the 'fmriprep-options' pre-configured pipeline.
    >>> fmriprep_options - Configuration({}) != {}
    True
    >>> fmriprep_options - fmriprep_options
    {}
    >>> fmriprep_options - Preconfiguration('fmriprep-options')
    Loading the 'fmriprep-options' pre-configured pipeline.
    {}
    >>> fmriprep_options - Configuration({'FROM': 'fmriprep-options'})
    Loading the 'fmriprep-options' pre-configured pipeline.
    {}
    >>> fmriprep_options - Configuration(yaml.safe_load(
    ...     create_yaml_from_template(fmriprep_options, include_all=True)))
    {}
    >>> fmriprep_options - Configuration(yaml.safe_load(
    ...     create_yaml_from_template(fmriprep_options, include_all=False)))
    {}
    >>> different_sca = Configuration({'pipeline_setup': {
    ...     'pipeline_name': 'different_SCA'},
    ...     'seed_based_correlation_analysis': {'run': 'y',
    ...     'norm_timeseries_for_DR': 'Off'}})
    >>> (Configuration(yaml.safe_load(create_yaml_from_template(
    ...     different_sca))) - Configuration()).get(
    ...     'seed_based_correlation_analysis') not in (None, {})
    True
    """
    yaml_template = YamlTemplate(template)
    d = d.dict() if isinstance(d, Configuration) else d
    return yaml_template.dump(new_dict=d)


def _format_list_items(l,  # noqa: E741  # pylint:disable=invalid-name
                       line_level):
    '''Helper method to handle lists in the YAML

    Parameters
    ----------
    l : list

    line_level : int

    Returns
    -------
    yaml : str

    Examples
    --------
    # >>> _format_list_items([1, 2, {'nested': 3}], 0
    # ... ) == r'\n  - 1\n  - 2\n  - nested: 3'
    # True
    # >>> _format_list_items([1, 2, {'nested': [3, {'deep': [4]}]}], 1
    # ... ) == (r'\n    - 1\n    - 2\n    - nested:\n      - 3\n'
    # ...       '      - deep:\n        - 4')
    # True
    '''
    # keep short, simple lists in square brackets
    if all(isinstance(item, (str, bool, int, float)) for item in l):
        preformat = str([yaml_bool(item) for item in l])
        if len(preformat) < 50:
            return preformat.replace("'", '').replace('"', '')
    # list long or complex lists on lines with indented '-' lead-ins
    return '\n' + '\n'.join([
        f'{indent(line_level)}{li}' for li in yaml.dump(
            yaml_bool(l)
        ).replace("'On'", 'On').replace("'Off'", 'Off').split('\n')
    ]).rstrip()


def hash_data_config(sub_list):
    '''Function to generate a short SHA1 hash from a data config
    subject list of dicts

    Parameters
    ----------
    sub_list : list of dicts

    Returns
    -------
    data_config_hash : str, len(8)

    Examples
    --------
    >>> sub_list = [{'site_id': f'site{i}', 'subject_id': f'sub{i}',
    ...              'unique_id': f'uid{i}'} for i in range(1, 4)]
    >>> sub_list[0]
    {'site_id': 'site1', 'subject_id': 'sub1', 'unique_id': 'uid1'}
    >>> hash_data_config(sub_list)
    '6f49a278'
    '''
    return sha1('_'.join([','.join([run.get(key, '') for run in sub_list]) for
                key in ['site_id', 'subject_id',
                        'unique_id']]).encode('utf-8')).hexdigest()[:8]


def indent(line_level, plus=2):
    '''Function to return an indent string for a given level

    Parameters
    ----------
    line_level : int
        The level of indentation to return

    Returns
    -------
    str
        The string of spaces to use for indentation
    '''
    return " " * (2 * line_level + plus)


def yaml_bool(value):
    '''Helper function to give On/Off value to bools

    Parameters
    ----------
    value : any

    Returns
    -------
    value : any

    Examples
    --------
    >>> yaml_bool(True)
    'On'
    >>> yaml_bool([False, 'On', True])
    ['Off', 'On', 'On']
    '''
    if isinstance(value, str):
        lookup_value = value.lower()
        if lookup_value in YAML_LOOKUP:
            value = YAML_LOOKUP[lookup_value]
    elif isinstance(value, list):
        return [yaml_bool(item) for item in value]
    elif isinstance(value, dict):
        return {k: yaml_bool(value[k]) for k in value}
    if isinstance(value, bool):
        if value is True:
            return 'On'
        return 'Off'
    return value


def upgrade_pipeline_to_1_8(path):
    '''Function to upgrade a C-PAC 1.7 pipeline config to C-PAC 1.8

    Parameters
    ----------
    path : str

    Returns
    -------
    None

    Outputs
    -------
    {path}.{now}.bak
        original file

    path
        upgraded file
    '''
    # back up original config
    now = datetime.isoformat(datetime.now()).replace(':', '_')
    backup = f'{path}.{now}.bak'
    print(f'Backing up {path} to {backup} and upgrading to C-PAC 1.8')
    original = open(path, 'r').read()
    open(backup, 'w').write(original)
    # upgrade and overwrite
    orig_dict = yaml.safe_load(original)
    # set Regressor 'Name's if not provided
    regressors = orig_dict.get('Regressors')
    if isinstance(regressors, list):
        for i, regressor in enumerate(regressors):
            if 'Name' not in regressor:
                regressor['Name'] = f'Regressor-{str(i + 1)}'
    if 'pipelineName' in orig_dict and len(original.strip()):
        middle_dict, leftovers_dict, complete_dict = update_config_dict(
            orig_dict)
        open(path, 'w').write(create_yaml_from_template(
            update_pipeline_values_1_8(middle_dict))
        )
        if leftovers_dict:
            open(f'{path}.rem', 'w').write(yaml.dump(leftovers_dict))
