'''Utilities for determining BIDS standard template identifiers
(https://bids-specification.readthedocs.io/en/stable/99-appendices/08-coordinate-systems.html#standard-template-identifiers)
from in-container template paths'''
import os
import re
import numpy as np

LOOKUP_TABLE = {row[0]: row[1] for row in
                np.loadtxt(os.path.join(os.path.dirname(__file__),
                                        'BIDS_identifiers.tsv'),
                dtype='str', delimiter='\t')}


def lookup_identifier(template_path: str) -> str:
    '''Function to return a standard template identifier for a packaged
    template, if known. Otherwise, returns the literal string
    'template'

    Parameters
    ----------
    template_path : str

    Returns
    -------
    identifier : str

    Examples
    --------
    >>> lookup_identifier('/usr/share/fsl/5.0/data/standard/'
    ...                   'MNI152_T1_1mm_brain.nii.gz')
    'MNI152NLin6ASym'
    >>> lookup_identifier('/code/CPAC/resources/templates/'
    ...                   'tpl-MNI152NLin2009cAsym_res-01_label-brain_'
    ...                   'probseg.nii.gz')
    'MNI152NLin2009cAsym'
    >>> lookup_identifier('/cpac_templates/chd8_functional_template_sk.nii')
    'template'
    '''
    for key, value in LOOKUP_TABLE.items():
        if re.search(key, template_path) is not None:
            return value
    return 'template'
