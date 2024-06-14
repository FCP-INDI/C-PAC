# Copyright (C) 2022-2024  C-PAC Developers

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
"""C-PAC Configuration class and related functions."""

import os
import re
from typing import Optional
from warnings import warn

from click import BadParameter
import pkg_resources as p
import yaml

from .diff import dct_diff

CONFIG_KEY_TYPE = str | list[str]
SPECIAL_REPLACEMENT_STRINGS = {r"${resolution_for_anat}", r"${func_resolution}"}


class ConfigurationDictUpdateConflation(SyntaxError):
    """Custom exception to clarify similar methods."""

    def __init__(self):
        self.msg = (
            "`Configuration().update` requires a key and a value. "
            "Perhaps you meant `Configuration().dict().update`?"
        )
        super().__init__()


class Configuration:
    """
    Class to set dictionary keys as map attributes.

    If the given dictionary includes the key ``FROM``, that key's value
    will form the base of the Configuration object with the values in
    the given dictionary overriding matching keys in the base at any
    depth. If no ``FROM`` key is included, the base Configuration is
    the default Configuration.

    ``FROM`` accepts either the name of a preconfigured pipleine or a
    path to a YAML file.

    Given a Configuration ``c``, and a list or tuple of an attribute name
    and nested keys ``keys = ['attribute', 'key0', 'key1']`` or
    ``keys = ('attribute', 'key0', 'key1')``, the value 'value' nested in

    .. code-block:: python

        c.attribute = {'key0': {'key1': 'value'}}

    can be accessed (get and set) in any of the following ways (and
    more):

    .. code-block:: python

        c.attribute['key0']['key1']
        c['attribute']['key0']['key1']
        c['attribute', 'key0', 'key1']
        c[keys]

    Examples
    --------
    >>> c = Configuration({})
    >>> c['pipeline_setup', 'pipeline_name']
    'cpac-blank-template'
    >>> c = Configuration({'pipeline_setup': {
    ...     'pipeline_name': 'example_pipeline'}})
    >>> c['pipeline_setup', 'pipeline_name']
    'example_pipeline'
    >>> c['pipeline_setup', 'pipeline_name'] = 'new_pipeline2'
    >>> c['pipeline_setup', 'pipeline_name']
    'new_pipeline2'

    >>> from CPAC.utils.tests.configs import SLACK_420349

    # test "FROM: /path/to/file"
    >>> slack_420349_filepath = Configuration(
    ...     yaml.safe_load(SLACK_420349['filepath']))
    >>> slack_420349_filepath['pipeline_setup', 'pipeline_name']
    'slack_420349_filepath'

    # test "FROM: preconfig"
    >>> slack_420349_preconfig = Configuration(
    ...    yaml.safe_load(SLACK_420349['preconfig']))
    >>> slack_420349_preconfig['pipeline_setup', 'pipeline_name']
    'slack_420349_preconfig'
    """

    def __init__(
        self, config_map: Optional[dict] = None, skip_env_check: bool = False
    ) -> None:
        """Initialize a Configuration instance.

        Parameters
        ----------
        config_map : dict, optional

        skip_env_check : bool, optional
        """
        from CPAC.pipeline.schema import schema
        from CPAC.utils.utils import lookup_nested_value, update_nested_dict

        if config_map is None:
            config_map = {}
        if skip_env_check:
            config_map["skip env check"] = True

        base_config = config_map.pop("FROM", None)
        if base_config:
            if base_config.lower() in ["default", "default_pipeline"]:
                base_config = "default"
            # import another config (specified with 'FROM' key)
            try:
                base_config = Preconfiguration(
                    base_config, skip_env_check=skip_env_check
                )
            except BadParameter:
                base_config = configuration_from_file(base_config)
            config_map = update_nested_dict(base_config.dict(), config_map)
        else:
            # base everything on blank pipeline for unspecified keys
            config_map = update_nested_dict(
                preconfig_yaml("blank", load=True), config_map
            )

        config_map = self._nonestr_to_None(config_map)

        try:
            regressors = lookup_nested_value(
                config_map,
                ["nuisance_corrections", "2-nuisance_regression", "Regressors"],
            )
        except KeyError:
            regressors = []
        if isinstance(regressors, list):
            for i, regressor in enumerate(regressors):
                # set Regressor 'Name's if not provided
                if "Name" not in regressor:
                    regressor["Name"] = f"Regressor-{i + 1!s}"
                # make Regressor 'Name's Nipype-friendly
                regressor["Name"] = nipype_friendly_name(regressor["Name"])

        config_map = schema(config_map)

        # remove 'skip env check' now that the config is validated
        if "skip env check" in config_map:
            del config_map["skip env check"]
        # remove 'FROM' before setting attributes now that it's imported
        if "FROM" in config_map:
            del config_map["FROM"]

        if skip_env_check:
            for key in config_map:
                # set attribute
                setattr(self, key, self.set_without_ENV(config_map[key]))
        else:
            # set FSLDIR to the environment $FSLDIR if the user sets it to
            # 'FSLDIR' in the pipeline config file
            _FSLDIR = config_map.get("FSLDIR")
            if _FSLDIR and bool(re.match(r"^[\$\{]{0,2}?FSLDIR[\}]?$", _FSLDIR)):
                config_map["FSLDIR"] = os.environ["FSLDIR"]
            for key in config_map:
                # set attribute
                setattr(self, key, self.set_from_ENV(config_map[key]))
        self._update_attr()

        # set working directory as an environment variable
        os.environ["CPAC_WORKDIR"] = self["pipeline_setup", "working_directory", "path"]

    def __str__(self):
        return f"C-PAC Configuration ('{self['pipeline_setup', 'pipeline_name']}')"

    def __repr__(self):
        """Show Configuration as a dict when accessed directly."""
        return str(self.dict())

    def __copy__(self):
        newone = type(self)({})
        newone.__dict__.update(self.__dict__)
        newone._update_attr()
        return newone

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)
        if isinstance(key, (list, tuple)):
            return self.get_nested(self, key)
        self.key_type_error(key)
        return None

    def __setitem__(self, key, value):
        if isinstance(key, str):
            setattr(self, key, value)
        elif isinstance(key, (list, tuple)):
            self.set_nested(self, key, value)
        else:
            self.key_type_error(key)

    def __sub__(self: "Configuration", other: "Configuration"):
        """Return the set difference between two Configurations.

        Examples
        --------
        >>> diff = (Preconfiguration('fmriprep-options')
        ...         - Preconfiguration('default'))
        >>> diff['pipeline_setup']['pipeline_name']
        ('cpac_fmriprep-options', 'cpac-default-pipeline')
        >>> diff['pipeline_setup']['pipeline_name'].s_value
        'cpac_fmriprep-options'
        >>> diff['pipeline_setup']['pipeline_name'].t_value
        'cpac-default-pipeline'
        >>> diff.s_value['pipeline_setup']['pipeline_name']
        'cpac_fmriprep-options'
        >>> diff.t_value['pipeline_setup']['pipeline_name']
        'cpac-default-pipeline'
        >>> diff['pipeline_setup']['pipeline_name'].left
        'cpac_fmriprep-options'
        >>> diff.left['pipeline_setup']['pipeline_name']
        'cpac_fmriprep-options'
        >>> diff['pipeline_setup']['pipeline_name'].minuend
        'cpac_fmriprep-options'
        >>> diff.minuend['pipeline_setup']['pipeline_name']
        'cpac_fmriprep-options'
        >>> diff['pipeline_setup']['pipeline_name'].right
        'cpac-default-pipeline'
        >>> diff.right['pipeline_setup']['pipeline_name']
        'cpac-default-pipeline'
        >>> diff['pipeline_setup']['pipeline_name'].subtrahend
        'cpac-default-pipeline'
        >>> diff.subtrahend['pipeline_setup']['pipeline_name']
        'cpac-default-pipeline'
        """
        return dct_diff(self.dict(), other.dict())

    def dict(self):
        """Show contents of a C-PAC configuration as a dict."""
        return {k: v for k, v in self.__dict__.items() if not callable(v)}

    def keys(self):
        """Show toplevel keys of a C-PAC configuration dict."""
        return self.dict().keys()

    def _nonestr_to_None(self, d):
        """Recursive method to type convert 'None' to None in nested config.

        Parameters
        ----------
        d : any
            config item to check

        Returns
        -------
        d : any
            same item, same type, but with 'none' strings converted to
            Nonetypes
        """
        if isinstance(d, str) and d.lower() == "none":
            return None
        if isinstance(d, list):
            return [self._nonestr_to_None(i) for i in d]
        if isinstance(d, set):
            return {self._nonestr_to_None(i) for i in d}
        if isinstance(d, dict):
            return {i: self._nonestr_to_None(d[i]) for i in d}
        return d

    def set_from_ENV(self, conf):  # pylint: disable=invalid-name
        """Replace strings like $VAR and ${VAR} with environment variable values.

        Parameters
        ----------
        conf : any

        Returns
        -------
        conf : any

        Examples
        --------
        >>> import os
        >>> os.environ['SAMPLE_VALUE_SFE'] = '/example/path'
        >>> c = Configuration()
        >>> c.set_from_ENV({'key': {'nested_list': [
        ...     1, '1', '$SAMPLE_VALUE_SFE/extended']}})
        {'key': {'nested_list': [1, '1', '/example/path/extended']}}
        >>> c.set_from_ENV(['${SAMPLE_VALUE_SFE}', 'SAMPLE_VALUE_SFE'])
        ['/example/path', 'SAMPLE_VALUE_SFE']
        >>> del os.environ['SAMPLE_VALUE_SFE']
        """
        if isinstance(conf, list):
            return [self.set_from_ENV(item) for item in conf]
        if isinstance(conf, dict):
            return {key: self.set_from_ENV(conf[key]) for key in conf}
        if isinstance(conf, str):
            # set any specified environment variables
            # (only matching all-caps plus `-` and `_`)
            # like `${VAR}`
            _pattern1 = r"\${[A-Z\-_]*}"
            # like `$VAR`
            _pattern2 = r"\$[A-Z\-_]*(?=/|$)"
            # replace with environment variables if they exist
            for _pattern in [_pattern1, _pattern2]:
                _match = re.search(_pattern, conf)
                if _match:
                    _match = _match.group().lstrip("${").rstrip("}")
                    conf = re.sub(_pattern, os.environ.get(_match, f"${_match}"), conf)
        return conf

    def set_without_ENV(self, conf):  # pylint: disable=invalid-name
        """Retain strings like $VAR and ${VAR} when setting attributes.

        Parameters
        ----------
        conf : any

        Returns
        -------
        conf : any

        Examples
        --------
        >>> import os
        >>> os.environ['SAMPLE_VALUE_SFE'] = '/example/path'
        >>> c = Configuration()
        >>> c.set_without_ENV({'key': {'nested_list': [
        ...     1, '1', '$SAMPLE_VALUE_SFE/extended']}})
        {'key': {'nested_list': [1, '1', '$SAMPLE_VALUE_SFE/extended']}}
        >>> c.set_without_ENV(['${SAMPLE_VALUE_SFE}', 'SAMPLE_VALUE_SFE'])
        ['${SAMPLE_VALUE_SFE}', 'SAMPLE_VALUE_SFE']
        >>> del os.environ['SAMPLE_VALUE_SFE']
        """
        if isinstance(conf, list):
            return [self.set_without_ENV(item) for item in conf]
        if isinstance(conf, dict):
            return {key: self.set_without_ENV(conf[key]) for key in conf}
        return conf

    def sub_pattern(self, pattern, orig_key):
        return orig_key.replace(pattern, self[pattern[2:-1].split(".")])

    def check_pattern(self, orig_key, tags=None):
        if tags is None:
            tags = []
        if isinstance(orig_key, dict):
            return {k: self.check_pattern(orig_key[k], tags) for k in orig_key}
        if isinstance(orig_key, list):
            return [self.check_pattern(item) for item in orig_key]
        if not isinstance(orig_key, str):
            return orig_key
        template_pattern = r"\${.*}"
        r = re.finditer(template_pattern, orig_key)
        for i in r:
            pattern = i.group(0)
            if isinstance(pattern, str) and len(pattern) and pattern not in tags:
                try:
                    orig_key = self.sub_pattern(pattern, orig_key)
                except AttributeError as ae:
                    if pattern not in SPECIAL_REPLACEMENT_STRINGS:
                        warn(str(ae), category=SyntaxWarning)
        return orig_key

    # method to find any pattern ($) in the configuration
    # and update the attributes with its pattern value
    def _update_attr(self):
        def check_path(key):
            if isinstance(key, str) and "/" in key:
                if not os.path.exists(key):
                    warn(f"Invalid path- {key}. Please check your configuration file")

        attributes = [
            (attr, getattr(self, attr))
            for attr in dir(self)
            if not callable(attr) and not attr.startswith("__")
        ]

        template_list = [
            "template_brain_only_for_anat",
            "template_skull_for_anat",
            "ref_mask",
            "template_brain_only_for_func",
            "template_skull_for_func",
            "template_symmetric_brain_only",
            "template_symmetric_skull",
            "dilated_symmetric_brain_mask",
        ]

        for attr_key, attr_value in attributes:
            if attr_key in template_list:
                new_key = self.check_pattern(attr_value, "FSLDIR")
            else:
                new_key = self.check_pattern(attr_value)
            setattr(self, attr_key, new_key)

    def update(self, key, val=ConfigurationDictUpdateConflation()):
        if isinstance(key, dict):
            raise ConfigurationDictUpdateConflation
        if isinstance(val, Exception):
            raise val
        setattr(self, key, val)

    def get_nested(self, _d, keys):
        if _d is None:
            _d = {}
        if isinstance(keys, str):
            return _d[keys]
        if isinstance(keys, (list, tuple)):
            if len(keys) > 1:
                return self.get_nested(_d[keys[0]], keys[1:])
            return _d[keys[0]]
        return _d

    def set_nested(self, d, keys, value):  # pylint: disable=invalid-name
        if isinstance(keys, str):
            d[keys] = value
        elif isinstance(keys, (list, tuple)):
            if len(keys) > 1:
                d[keys[0]] = self.set_nested(d[keys[0]], keys[1:], value)
            else:
                d[keys[0]] = value
        return d

    def _check_if_switch(self, key: CONFIG_KEY_TYPE, error: bool = False) -> bool:
        """Check if a given entity is a switch.

        Parameters
        ----------
        key : str or list of str
            key to check

        error : bool
            raise a TypeError if not a switch

        Returns
        -------
        bool
            True if the given key is a switch, False otherwise

        Examples
        --------
        >>> c = Configuration()
        >>> c._check_if_switch('anatomical_preproc')
        False
        >>> c._check_if_switch(['anatomical_preproc'])
        False
        >>> c._check_if_switch(['anatomical_preproc', 'run'])
        True
        """
        _maybe_switch = self[key]
        if isinstance(_maybe_switch, bool):
            return True
        if isinstance(_maybe_switch, list):
            _answer = all(isinstance(_, bool) for _ in _maybe_switch)
            if _answer:
                return _answer
        if error:
            msg = f"`{key}` is not a switch in {self!s}."
            raise TypeError(msg)
        return False

    def _switch_bool(self, key: CONFIG_KEY_TYPE, value: bool, exclusive: bool) -> bool:
        """Return True if the key is set to the given value or False otherwise.

        Parameters
        ----------
        key : str or list of str
            key to check

        value : bool
            value to check for

        exclusive : bool
            return False if forking (both True and False)

        Returns
        -------
        bool
            True if the given key is set to the given value or False
            otherwise. If exclusive is True, return False if the key
            is set to both True and False.
        """
        if not (exclusive and self.switch_is_on_off(key)):
            if isinstance(self[key], bool):
                return self[key] is value
            if isinstance(self[key], list):
                return value in self[key]
        return False

    def switch_is_off(self, key: CONFIG_KEY_TYPE, exclusive: bool = False) -> bool:
        """Return True if the key is set to 'off' OR 'on' and 'off' or False otherwise.

        Used for tracking forking.

        Parameters
        ----------
        key : str or list of str
            key to check

        exclusive : bool, optional, default: False
            return False if the key is set to 'on' and 'off'

        Returns
        -------
        bool
            True if key is set to 'off', False if not set to 'off'.
            If exclusive is set to True, return False if the key is
            set to 'on' and 'off'.

        Examples
        --------
        >>> c = Configuration()
        >>> c.switch_is_off(['nuisance_corrections', '2-nuisance_regression',
        ...                  'run'])
        True
        >>> c = Configuration({'nuisance_corrections': {
        ...     '2-nuisance_regression': {'run': [True, False]}}})
        >>> c.switch_is_off(['nuisance_corrections', '2-nuisance_regression',
        ...                  'run'])
        True
        >>> c.switch_is_off(['nuisance_corrections', '2-nuisance_regression',
        ...                  'run'], exclusive=True)
        False
        """
        self._check_if_switch(key, True)
        return self._switch_bool(key, False, exclusive)

    def switch_is_on(self, key: CONFIG_KEY_TYPE, exclusive: bool = False) -> bool:
        """Return True if the key is set to 'on' OR 'on' and 'off' or False otherwise.

        Used for tracking forking.

        Parameters
        ----------
        key : str or list of str
            key to check

        exclusive : bool, optional, default: False
            return False if the key is set to 'on' and 'off'

        Returns
        -------
        bool
            True if key is set to 'on', False if not set to 'on'.
            If exclusive is set to True, return False if the key is
            set to 'on' and 'off'.

        Examples
        --------
        >>> c = Configuration()
        >>> c.switch_is_on(['nuisance_corrections', '2-nuisance_regression',
        ...                 'run'])
        False
        >>> c = Configuration({'nuisance_corrections': {
        ...     '2-nuisance_regression': {'run': [True, False]}}})
        >>> c.switch_is_on(['nuisance_corrections', '2-nuisance_regression',
        ...                 'run'])
        True
        >>> c.switch_is_on(['nuisance_corrections', '2-nuisance_regression',
        ...                 'run'], exclusive=True)
        False
        """
        self._check_if_switch(key, True)
        return self._switch_bool(key, True, exclusive)

    def switch_is_on_off(self, key: CONFIG_KEY_TYPE) -> bool:
        """Return True if the key is set to both 'on' and 'off' or False otherwise.

        Used for tracking forking.

        Parameters
        ----------
        key : str or list of str
            key to check

        Returns
        -------
        bool
            True if key is set to 'on' and 'off', False otherwise

        Examples
        --------
        >>> c = Configuration()
        >>> c.switch_is_on_off(['nuisance_corrections',
        ...                     '2-nuisance_regression', 'run'])
        False
        >>> c = Configuration({'nuisance_corrections': {
        ...     '2-nuisance_regression': {'run': [True, False]}}})
        >>> c.switch_is_on_off(['nuisance_corrections',
        ...                     '2-nuisance_regression', 'run'])
        True
        """
        self._check_if_switch(key, True)
        if isinstance(self[key], list):
            return True in self[key] and False in self[key]
        return False

    def key_type_error(self, key):
        """Raise a KeyError if an inappropriate type of key is attempted."""
        raise KeyError(
            " ".join(
                [
                    "Configuration key must be a string, list, or tuple;",
                    type(key).__name__,
                    f"`{key!s}`",
                    "was given.",
                ]
            )
        )


def check_pname(p_name: str, pipe_config: Configuration) -> str:
    """Check / set `p_name`, the str representation of a pipeline for use in filetrees.

    Parameters
    ----------
    p_name : str or None

    pipe_config : Configuration

    Returns
    -------
    p_name

    Examples
    --------
    >>> c = Configuration()
    >>> check_pname(None, c)
    'pipeline_cpac-blank-template'
    >>> check_pname('cpac-default-pipeline', c)
    'pipeline_cpac-default-pipeline'
    >>> check_pname('pipeline_cpac-default-pipeline', c)
    'pipeline_cpac-default-pipeline'
    >>> check_pname('different-name', Configuration())
    'pipeline_different-name'
    >>> p_name = check_pname(None, Preconfiguration('blank'))
    >>> p_name
    'pipeline_cpac-blank-template'
    >>> p_name = check_pname(None, Preconfiguration('default'))
    >>> p_name
    'pipeline_cpac-default-pipeline'
    """
    if p_name is None:
        p_name = f'pipeline_{pipe_config["pipeline_setup", "pipeline_name"]}'
    elif not p_name.startswith("pipeline_"):
        p_name = f"pipeline_{p_name}"
    return p_name


def collect_key_list(config_dict):
    """Return a list of lists of keys for a nested dictionary.

    Parameters
    ----------
    config_dict : dict

    Returns
    -------
    key_list : list

    Examples
    --------
    >>> collect_key_list({'test': {'nested': 1, 'dict': 2}})
    [['test', 'nested'], ['test', 'dict']]
    """
    key_list = []
    for key in config_dict:
        if isinstance(config_dict[key], dict):
            for inner_key_list in collect_key_list(config_dict[key]):
                key_list.append([key, *inner_key_list])
        else:
            key_list.append([key])
    return key_list


def configuration_from_file(config_file):
    """Load a Configuration from a pipeline config file.

    Parameters
    ----------
    config_file : str
        path to configuration file

    Returns
    -------
    Configuration
    """
    with open(config_file, "r", encoding="utf-8") as config:
        return Configuration(yaml.safe_load(config))


def preconfig_yaml(preconfig_name="default", load=False):
    """Get the path to a preconfigured pipeline's YAML file.

    Raises BadParameter if an invalid preconfig name is given.

    Parameters
    ----------
    preconfig_name : str

    load : boolean
        return dict if True, str if False

    Returns
    -------
    str or dict
        path to YAML file or dict loaded from YAML
    """
    from CPAC.pipeline import ALL_PIPELINE_CONFIGS, AVAILABLE_PIPELINE_CONFIGS

    if preconfig_name not in ALL_PIPELINE_CONFIGS:
        msg = (
            f"The pre-configured pipeline name '{preconfig_name}' you "
            "provided is not one of the available pipelines.\n\nAvailable "
            f"pipelines:\n{AVAILABLE_PIPELINE_CONFIGS!s}\n"
        )
        raise BadParameter(
            msg,
            param="preconfig",
        )
    if load:
        with open(preconfig_yaml(preconfig_name), "r", encoding="utf-8") as _f:
            return yaml.safe_load(_f)
    return p.resource_filename(
        "CPAC",
        os.path.join("resources", "configs", f"pipeline_config_{preconfig_name}.yml"),
    )


class Preconfiguration(Configuration):
    """A preconfigured Configuration.

    Parameters
    ----------
    preconfig : str
        The canonical name of the preconfig to load
    """

    def __init__(self, preconfig, skip_env_check=False):
        super().__init__(
            config_map=preconfig_yaml(preconfig, True), skip_env_check=skip_env_check
        )


def set_subject(
    sub_group, pipe_config: "Configuration", p_name: Optional[str] = None
) -> TUPLE[str, str, str]:
    """Set pipeline name and log directory path for a given sub_dict.

    Parameters
    ----------
    sub_dict : dict

    pipe_config : CPAC.utils.configuration.Configuration

    p_name : str, optional
        pipeline name string

    Returns
    -------
    subject_id : str

    p_name : str
        pipeline name string

    log_dir : str
        path to subject log directory

    Examples
    --------
    >>> from tempfile import TemporaryDirectory
    >>> from CPAC.utils.configuration import Configuration
    >>> sub_dict = {'site_id': 'site1', 'subject_id': 'sub1',
    ...             'unique_id': 'uid1'}
    >>> with TemporaryDirectory() as tmpdir:
    ...     subject_id, p_name, log_dir = set_subject(
    ...         sub_dict, Configuration({'pipeline_setup': {'log_directory':
    ...             {'path': tmpdir}}}))
    >>> subject_id
    'sub1_uid1'
    >>> p_name
    'pipeline_cpac-blank-template'
    >>> log_dir.endswith(f'{p_name}/{subject_id}')
    True
    """
    subject_id = sub_group[0][0]
    if sub_group[0][1]:
        subject_id += f"_{sub_group[0][1]}"
    p_name = check_pname(p_name, pipe_config)
    log_dir = os.path.join(
        pipe_config.pipeline_setup["log_directory"]["path"], p_name, subject_id
    )
    if not os.path.exists(log_dir):
        os.makedirs(os.path.join(log_dir))
    return subject_id, p_name, log_dir


def nipype_friendly_name(name: str) -> str:
    """Replace each sequence of non-alphanumeric characters...

    ...with a single underscore and remove any leading underscores.

    Parameters
    ----------
    name : str

    Returns
    -------
    str
    """
    return re.sub(r"[^a-zA-Z0-9]+", "_", name).lstrip("_")
