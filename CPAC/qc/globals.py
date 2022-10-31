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
"""Global QC values"""
_REGISTRATION_GUARDRAIL_THRESHOLDS = {'thresholds': {}}


def registration_guardrail_thresholds() -> dict:
    """Get registration guardrail thresholds

    Returns
    -------
    dict
    """
    return _REGISTRATION_GUARDRAIL_THRESHOLDS['thresholds']


def update_thresholds(thresholds) -> None:
    """Set a registration guardrail threshold

    Parameters
    ----------
    thresholds : dict of {str: float or int}

    Returns
    -------
    None
    """
    _REGISTRATION_GUARDRAIL_THRESHOLDS['thresholds'].update(thresholds)
