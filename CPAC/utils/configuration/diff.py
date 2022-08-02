'''Tools for configuration differences

Copyright (C) 2022  C-PAC Developers

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
License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.'''


def dct_diff(dct1, dct2):
    '''Function to compare 2 nested dicts, dropping values unspecified
    in the second. Adapted from https://github.com/FCP-INDI/CPAC_regtest_pack/blob/9056ef63/cpac_pipe_diff.py#L31-L78

    Parameters
    ----------
    dct1 : dict or CPAC.utils.Configuration

    dct2 : dict or CPAC.utils.Configuration

    Returns
    -------
    diff : dict
       each value is a DiffValue of values from dct1, dct2 for each
       differing key between original dicts or a subdictionary thereof

    Example
    -------
    >>> import yaml
    >>> from . import DEFAULT_PIPELINE_FILE, Preconfiguration
    >>> def read_yaml_file(yaml_file):
    ...     return yaml.safe_load(open(yaml_file, 'r'))
    >>> pipeline = read_yaml_file(DEFAULT_PIPELINE_FILE)
    >>> dct_diff(pipeline, pipeline)
    {}
    >>> pipeline2 = Preconfiguration('fmriprep-options')  # doctest: +NORMALIZE_WHITESPACE
    Loading the 'fmriprep-options' pre-configured pipeline.
    >>> dct_diff(pipeline, pipeline2)['pipeline_setup']['pipeline_name']
    ('cpac-default-pipeline', 'cpac_fmriprep-options')
    '''  # pylint: disable=line-too-long  # noqa: E501
    diff = DiffDict()
    for key in dct1:
        if isinstance(dct1[key], dict):
            if not isinstance(dct2, dict):
                try:
                    dct2 = dct2.dict()
                except AttributeError:
                    raise TypeError(f'{dct2} is not a dict.')
            diff[key] = dct_diff(dct1[key], dct2.get(key, {}))
            if diff[key] == {}:
                del diff[key]
        else:
            dct1_val = dct1.get(key)
            dct2_val = dct2.get(key) if isinstance(dct2, dict) else None

            if dct1_val != dct2_val:
                diff[key] = DiffValue(dct1_val, dct2_val)

    # add any new keys
    if isinstance(dct2, dict):
        for key in dct2:
            if key not in dct1:
                diff[key] = dct2[key]

        # only return non-empty diffs
        return DiffDict({k: v for k, v in diff.items() if k in dct2})

    return DiffDict()


def diff_dict(diff):
    '''Method to return a dict of only changes given a nested dict
    of (dict1_value, dict2_value) tuples

    Parameters
    ----------
    diff : dict
        output of `dct_diff`

    Returns
    -------
    dict
        dict of only changed values

    Examples
    --------
    >>> diff_dict({'anatomical_preproc': {
    ...     'brain_extraction': {'extraction': {
    ...         'run': DiffValue([True], False),
    ...         'using': DiffValue(['3dSkullStrip'],
    ...                            ['niworkflows-ants'])}}}})
    {'anatomical_preproc': {'brain_extraction': {'extraction': {'run': False, 'using': ['niworkflows-ants']}}}}
    '''  # noqa: E501  # pylint: disable=line-too-long
    if isinstance(diff, DiffValue):
        return diff.t_value
    if isinstance(diff, dict):
        i = DiffDict()
        for k in diff:
            try:
                j = diff_dict(diff[k])
                if j != {}:
                    i[k] = j
            except KeyError:
                continue
        return i
    return diff


class DiffDict(dict):
    '''Class to sematically store a dictionary of set differences from
    Configuration(S) - Configuration(T)

    Attributes
    ----------
    left : dict
        dictionary of differing values from Configuration(S)
        (alias for s_value)

    minuend : dict
        dictionary of differing values from Configuration(S)
        (alias for s_value)

    right : dict
        dictionary of differing values from Configuration(T)
        (alias for t_value)

    subtrahend : dict
        dictionary of differing values from Configuration(T)
        (alias for t_value)

    s_value : dict
        dictionary of differing values from Configuration(S)

    t_value : dict
        dictionary of differing values from Configuration(T)
    '''
    def __init__(self, *args, **kwargs):
        '''Dictionary of difference Configuration(S) - Configuration(T).

        Each value in a DiffDict should be either a DiffDict or a DiffValue.
        '''
        super().__init__(*args, **kwargs)
        self.left = self.minuend = self.s_value = self._s_value()
        self.right = self.subtrahend = self.t_value = self._t_value()

    def _return_one_value(self, which_value):
        return_dict = {}
        for k, v in self.items():
            if isinstance(v, (DiffDict, DiffValue)):
                return_dict[k] = getattr(v, which_value)
            else:
                return_dict[k] = v
        return return_dict

    def _s_value(self):
        '''Get a dictionary of only the differing 'S' values that differ
        in S - T'''
        return self._return_one_value('s_value')

    def _t_value(self):
        '''Get a dictionary of only the differing 'T' values that differ
        in S - T'''
        return self._return_one_value('t_value')


class DiffValue:
    '''Class to semantically store values of set difference from
    Configuration(S) - Configuration(T)

    Attributes
    ----------
    left : any
        value from Configuration(S) (alias for s_value)

    minuend : dict
        value from Configuration(S) (alias for s_value)

    right : dict
        value from Configuration(T) (alias for t_value)

    subtrahend : dict
        value from Configuration(T) (alias for t_value)

    s_value : any
        value from Configuration(S)

    t_value : any
        value from Configuration(T)
    '''
    def __init__(self, s_value, t_value):
        '''Different values from Configuration(S) - Configuration(T)

        Parameters
        ----------
        s_value : any
           value from Configuration(S)

        t_value : any
           value from Configuration(T)
        '''
        self.left = self.minuend = self.s_value = s_value
        self.right = self.subtrahend = self.t_value = t_value

    def __len__(self):
        return 2  # self.__repr__ should always be a 2-tuple

    def __repr__(self):
        return str(tuple((self.s_value, self.t_value)))
