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
"""Global registration guardrail values"""
from traits.trait_base import Undefined
_REGISTRATION_GUARDRAILS = {}


def __getattr__(name):
    """Get global values"""
    if name == 'retry_on_first_failure':
        return (_REGISTRATION_GUARDRAILS.get('best_of') == 1 and
                _REGISTRATION_GUARDRAILS.get(name) is True)
    return _REGISTRATION_GUARDRAILS.get(name, Undefined)


def update(_dict) -> None:
    """Set registration guardrails

    Parameters
    ----------
    _dict : dict
        keys and values to update

    Returns
    -------
    None
    """
    _REGISTRATION_GUARDRAILS.update(_dict)
