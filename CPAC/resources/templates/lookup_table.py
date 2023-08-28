# Copyright (C) 2022-2023  C-PAC Developers

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
"""Utilities for determining BIDS standard template identifiers
(https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#standard-template-identifiers)
from in-container template paths"""
from os import environ, path as op
from re import findall, search
from typing import Optional
from bids.layout import parse_file_entities
from numpy import loadtxt
from CPAC.utils.typing import TUPLE

LOOKUP_TABLE = {row[0].replace(r'$FSLDIR', environ['FSLDIR']): (
    row[1], str(row[2]) if row[2] else None) for row in
    loadtxt(op.join(op.dirname(__file__), 'BIDS_identifiers.tsv'),
            dtype='str', delimiter='\t')}


def format_identifier(identifier: str, desc: Optional[str] = None) -> str:
    '''Function to create an identifier string from a name and description

    Parameters
    ----------
    identifier : str

    desc : str, optional

    Returns
    -------
    str

    Examples
    --------
    >>> format_identifier('CC', '200')
    'CC_desc-200'
    >>> format_identifier('AAL')
    'AAL'
    '''
    if desc:
        return f'{identifier}_desc-{desc}'
    return identifier


def lookup_identifier(template_path: str) -> TUPLE[str, None]:
    '''Function to return a standard template identifier for a packaged
    template, if known. Otherwise, returns the literal string
    'template'

    Parameters
    ----------
    template_path : str

    Returns
    -------
    identifier : str

    desc : str or None

    Examples
    --------
    >>> lookup_identifier('$FSLDIR/data/standard/'
    ...                   'MNI152_T1_1mm_brain.nii.gz')
    ('MNI152NLin6ASym', None)
    >>> lookup_identifier('/code/CPAC/resources/templates/'
    ...                   'tpl-MNI152NLin2009cAsym_res-01_label-brain_'
    ...                   'probseg.nii.gz')
    ('MNI152NLin2009cAsym', None)
    >>> lookup_identifier('/cpac_templates/chd8_functional_template_sk.nii')
    ('template', None)
    >>> lookup_identifier('/cpac_templates/CC200.nii.gz')
    ('CC', '200')
    >>> lookup_identifier('s3://templateflow/tpl-MNI152NLin2009cAsym/'
    ...                   'tpl-MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_'
    ...                   'desc-200Parcels17Networks_dseg.nii.gz')
    ('Schaefer2018', '200Parcels17Networks')
    '''
    if r'$' in template_path:
        bash_var_pattern = r'(\$[\w]+(?=/|\s)|\${\w+})'
        matches = findall(bash_var_pattern, template_path)
        if matches is not None:
            for match in matches:
                bash_var = match.lstrip('${').rstrip('}')
                if bash_var in environ:
                    template_path = template_path.replace(match,
                                                          environ[bash_var])
    for key, value in LOOKUP_TABLE.items():
        if search(key, template_path) is not None:
            return value
    _bidsy = parse_file_entities(template_path)
    if 'atlas' in _bidsy:
        return _bidsy['atlas'], _bidsy.get('desc')
    return 'template', None
