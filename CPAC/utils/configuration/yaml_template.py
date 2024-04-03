#!/usr/bin/env python3
# Copyright (C) 2022  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Functions to create YAML configuration files from templates."""
from copy import deepcopy
import os
import re
from typing import Optional, Union
from datetime import datetime
from hashlib import sha1
from click import BadParameter
import yaml

from CPAC.utils.configuration import Configuration, Preconfiguration, \
                                     preconfig_yaml
from CPAC.utils.utils import update_config_dict, update_pipeline_values_1_8, \
                             YAML_BOOLS

YAML_LOOKUP = {yaml_str: key for key, value in YAML_BOOLS.items() for
               yaml_str in value}


class YamlTemplate():  # pylint: disable=too-few-public-methods
    """A class to link YAML comments to the contents of a YAML file

    Attributes
    ----------
    comments : dict
        Flat dictionary with ``'.'``-delimited pseudo-nested structure.
        E.g., comments for ``{'pipeline_setup': {'pipeline_name': value}}``
        would be keyed
        ``{'pipeline_setup': comment0, 'pipeline_setup.pipeline_name: comment1}`` to
        allow comments at each level of depth.

    dump : method

    get_nested : method

    original : str
    """
    def __init__(self, original_yaml, base_config=None):
        """
        Parameters
        ----------
        original_yaml : str
            raw YAML or path to YAML file

        base_config : Configuration, optional
        """
        try:
            original_yaml = preconfig_yaml(original_yaml)
        except BadParameter:
            pass
        if os.path.exists(original_yaml):
            with open(original_yaml, 'r', encoding='utf-8') as _f:
                original_yaml = _f.read()
        self.comments = {}
        self.template = original_yaml
        if base_config is None:
            if isinstance(self.template, dict):
                self._dict = self.template
            if isinstance(self.template, str):
                self._dict = yaml.safe_load(self.template)
        else:
            self._dict = base_config.dict()
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
    
        # SSOT FSLDIR
        try:  # Get from current config
            fsldir = self.get_nested(new_dict,
                                     ['pipeline_setup', 'system_config',
                                      'FSLDIR'])
        except KeyError:  # Get from imported base
            fsldir = self.get_nested(self._dict,
                                     ['pipeline_setup', 'system_config',
                                      'FSLDIR'])

        # Add YAML version directive to top of document and ensure
        # C-PAC version comment and 'FROM' are at the top of the YAML
        # output
        if parents is None:
            parents = []
            _dump = ['%YAML 1.1', '---']
            if 'pipeline_setup' not in new_dict:
                new_dict['pipeline_setup'] = None
        else:
            _dump = []
        # Prepare for indentation
        line_level = len(parents)
        # Get a safely mutable copy of the dict
        loop_dict = deepcopy(self.get_nested(new_dict, parents) if
                             parents else new_dict)
        # Grab special key to print first
        import_from = loop_dict.pop('FROM', None)
        # Iterate through mutated dict
        for key in loop_dict:
            # List of progressively-indented key strings
            keys = [*parents, key]
            # Comments are stored in a flat dictionary with
            # '.'-delimited pseudonested keys
            comment = self.comments.get('.'.join(keys))
            # This exception should only happen from mutations
            # introduced this function
            try:
                value = self.get_nested(new_dict, keys)
            except KeyError:  # exclude unincluded keys
                continue

            # Print comment if there's one above this key in the template
            if comment:
                if key != 'pipeline_setup':
                    _dump += ['']  # Add a blank line above the comment
                _dump += [indent(line_level, 0) + line for line in comment]
            # Print 'FROM' between preamble comment and rest of config
            # if applicable
            if key == 'pipeline_setup' and import_from is not None:
                _dump += [f'FROM: {import_from}', '']
            # Apply indentation to key
            indented_key = f'{indent(line_level, 0)}{key}:'
            # Print YAML-formatted value
            if value is not None:
                # SSOT FSLDIR
                if (isinstance(value, str) and fsldir in value and
                        key != 'FSLDIR'):
                    value = re.sub(r'\$*FSLDIR', '$FSLDIR',
                                   value.replace(fsldir, '$FSLDIR'))
                if isinstance(value, dict):
                    _dump += [indented_key, self.dump(new_dict, keys)]
                elif isinstance(value, list):
                    list_line = _format_list_items(value, line_level)
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
            elif key != 'pipeline_setup':
                _dump += [indented_key]
        # Normalize line spacing and return YAML string
        return re.sub('\n{3,}', '\n\n', '\n'.join(_dump)).rstrip() + '\n'

    def _parse_comments(self):
        # Split YAML into lines
        yaml_lines = self.template.split('\n')
        # Initialize comment and key
        comment = []
        key = []
        for line in yaml_lines:
            # Calculate indentation
            line_level = _count_indent(line)
            # Remove indentation and trailing whitespace
            stripped_line = line.strip()
            # Collect a line of a comment
            if stripped_line.startswith('#'):
                comment.append(stripped_line)
            # If a line is not a comment line:
            elif not any(stripped_line.startswith(seq) for
                         seq in ('%YAML', '---')):
                # If the line is a key
                if ':' in stripped_line:
                    # Set the key for the comments dictionary
                    line_key = stripped_line.split(':', 1)[0].strip()
                    if line_level == 0:
                        key = [line_key]
                    else:
                        key = [*key[:line_level], line_key]
                    # Store the full list of comment lines
                    self.comments['.'.join(key)] = comment
                    # Reset the comment variable to collect the next comment
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


def create_yaml_from_template(
        d: Union[Configuration, dict],  # pylint: disable=invalid-name
        template: str = 'default', import_from: Optional[str] = None,
        skip_env_check: Optional[bool] = False) -> str:
    """Save dictionary to a YAML file, keeping the structure
    (such as first level comments and ordering) from the template

    It may not be fully robust to YAML structures, but it works
    for C-PAC config files!

    Parameters
    ----------
    d : dict or Configuration

    template : str
        path to template, name of preconfig, or YAML as a string

    import_from : str, optional
        name of a preconfig. Full config is generated if omitted

    skip_env_check : bool, optional
        skip environment check (for validating a config without running)

    Examples
    --------
    >>> import yaml
    >>> from CPAC.utils.configuration import Configuration, Preconfiguration
    >>> Configuration(yaml.safe_load(create_yaml_from_template({}))).dict(
    ...     ) == Configuration({}).dict()
    True
    >>> fmriprep_options = Preconfiguration('fmriprep-options')
    >>> fmriprep_options - Configuration({}) != {}
    True
    >>> fmriprep_options - fmriprep_options
    {}
    >>> fmriprep_options - Preconfiguration('fmriprep-options')
    {}
    >>> fmriprep_options - Configuration({'FROM': 'fmriprep-options'})
    {}
    >>> fmriprep_options - Configuration(yaml.safe_load(
    ...     create_yaml_from_template(fmriprep_options, import_from=None)))
    {}
    >>> fmriprep_options - Configuration(yaml.safe_load(
    ...     create_yaml_from_template(fmriprep_options,
    ...                               import_from='default')))
    {}
    >>> fmriprep_options - Configuration(yaml.safe_load(
    ...     create_yaml_from_template(fmriprep_options, import_from='blank')))
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
    if import_from is None:  # full config
        d = d.dict() if isinstance(d, Configuration) else d
        base_config = None
    else:  # config based on preconfig
        d = Configuration(d) if not isinstance(d, Configuration) else d
        base_config = Preconfiguration(import_from, skip_env_check=skip_env_check)
        d = (d - base_config).left
        d.update({'FROM': import_from})
    yaml_template = YamlTemplate(template, base_config)
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
    >>> print(_format_list_items([1, 2, {'nested': 3}], 0))
      - 1
      - 2
      - nested: 3
    >>> print(
    ...     _format_list_items([1, 2, {'nested': [3, {'deep': [4]}]}], 1))
        - 1
        - 2
        - nested:
          - 3
          - deep:
            - 4
    '''
    # keep short, simple lists in square brackets
    if all(isinstance(item, (str, bool, int, float)) for item in l):
        preformat = str([yaml_bool(item) for item in l])
        if len(preformat) < 50:
            return preformat.replace("'", '').replace('"', '')
    # list long or complex lists on lines with indented '-' lead-ins
    return '\n'.join([
        f'{indent(line_level)}{li}' for li in yaml.dump(
            yaml_bool(l), sort_keys=False
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
        # if 'Name' is a key, promote that item to the top
        return {**({'Name': value['Name']} if 'Name' in value else {}),
                **{k: yaml_bool(value[k]) for k in value if k != 'Name'}}
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
    with open(path, 'r', encoding='utf-8') as _f:
        original = _f.read()
    with open(backup, 'w', encoding='utf-8') as _f:
        _f.write(original)
    # upgrade and overwrite
    orig_dict = yaml.safe_load(original)
    # set Regressor 'Name's if not provided
    regressors = orig_dict.get('Regressors')
    if isinstance(regressors, list):
        for i, regressor in enumerate(regressors):
            if 'Name' not in regressor:
                regressor['Name'] = f'Regressor-{str(i + 1)}'
    if 'pipelineName' in orig_dict and len(original.strip()):
        middle_dict, leftovers_dict, _complete_dict = update_config_dict(
            orig_dict)
        with open(path, 'w', encoding='utf-8') as _f:
            _f.write(create_yaml_from_template(
                update_pipeline_values_1_8(middle_dict)))
        if leftovers_dict:
            with open(f'{path}.rem', 'w', encoding='utf-8') as _f:
                _f.write(yaml.dump(leftovers_dict))


def update_a_preconfig(preconfig, import_from):
    """
    Parameters
    ----------
    preconfig : str

    import_from : str
    """
    import sys
    print(f'Updating {preconfig} preconfig…', file=sys.stderr)
    updated = create_yaml_from_template(Preconfiguration(preconfig,
                                                         skip_env_check=True),
                                        import_from=import_from, skip_env_check=True)
    with open(preconfig_yaml(preconfig), 'w', encoding='utf-8') as _f:
        _f.write(updated)


def update_all_preconfigs():
    """Update all other preconfigs with comments from default"""
    from CPAC.pipeline import ALL_PIPELINE_CONFIGS
    not_from_blank = ('anat-only', 'blank', 'default', 'fx-options',
                      'nhp-macaque', 'preproc', 'rbc-options')
    update_a_preconfig('blank', None)
    for preconfig in ('anat-only', 'preproc'):
        update_a_preconfig(preconfig, 'default')
    for preconfig in ('fx-options', 'rbc-options'):
        update_a_preconfig(preconfig, 'fmriprep-options')
    update_a_preconfig('nhp-macaque', 'monkey')
    for preconfig in (_ for _ in ALL_PIPELINE_CONFIGS if
                      _ not in not_from_blank):
        update_a_preconfig(preconfig, 'blank')


if __name__ == '__main__':
    update_all_preconfigs()
