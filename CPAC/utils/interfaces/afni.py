# Copyright (C) 2023  C-PAC Developers

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
"""AFNI Nipype interfaces with customized functionality"""
# pylint: disable=too-few-public-methods
from nipype.interfaces.afni.preprocess import ECM as _ECM, \
    ECMInputSpec as _ECMInputSpec
import semver
import traits.api as traits
from CPAC.utils.versioning import REPORTED

AFNI_SEMVER = REPORTED.get('AFNI', '0.0.0').split(':')[0]
"""Validated AFNI Semver string"""
try:
    AFNI_SEMVER = str(semver.Version.parse(AFNI_SEMVER))
except ValueError:
    _major, _minor, _patch = [int(part) for part in AFNI_SEMVER.split('.')]
    AFNI_SEMVER = str(semver.Version.parse(f'{_major}.{_minor}.{_patch}'))
    del _major, _minor, _patch
AFNI_GTE_21_1_1 = semver.compare(AFNI_SEMVER, '21.1.1') >= 0
"""AFNI version >= 21.1.1?"""


class ECMInputSpec(_ECMInputSpec):
    """ECMInputSpec + 'do_binary' option"""
    do_binary = traits.Bool(desc='whether to perform the ECM calculation on a '
                            'binarized version of the connectivity matrix',
                            argstr='-do_binary')


class ECM(_ECM):
    """ECM + 'do_binary' option"""
    input_spec = ECMInputSpec


__all__ = ['AFNI_GTE_21_1_1', 'ECM']
