"""Configs for testing"""
import os
from pkg_resources import resource_filename
with open(resource_filename("CPAC",
                            os.path.join('utils', 'tests', 'configs',
                                         'neurostars_23786.yml')),
          'r', encoding='utf-8') as _f:
    NEUROSTARS_23786 = _f.read()
with open(resource_filename("CPAC",
                            os.path.join('utils', 'tests', 'configs',
                                         'neurostars_24035.yml')),
          'r', encoding='utf-8') as _f:
    NEUROSTARS_24035 = _f.read()
__all__ = ['NEUROSTARS_23786', 'NEUROSTARS_24035']
