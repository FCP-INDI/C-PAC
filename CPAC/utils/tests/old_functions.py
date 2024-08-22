# Copyright (C) 2012-2024  C-PAC Developers

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
"""Functions from before refactoring."""


def check(params_dct, subject_id, scan_id, val_to_check, throw_exception):
    """https://github.com/FCP-INDI/C-PAC/blob/96db8b0b65ab1d5f55fb3b895855af34d72c17e4/CPAC/utils/utils.py#L630-L653"""
    if val_to_check not in params_dct:
        if throw_exception:
            raise Exception(
                f"Missing Value for {val_to_check} for participant " f"{subject_id}"
            )
        return None
    if isinstance(params_dct[val_to_check], dict):
        ret_val = params_dct[val_to_check][scan_id]
    else:
        ret_val = params_dct[val_to_check]
    if ret_val == "None":
        if throw_exception:
            raise Exception(
                f"'None' Parameter Value for {val_to_check} for participant "
                f"{subject_id}"
            )
        else:
            ret_val = None
    if ret_val == "" and throw_exception:
        raise Exception(
            f"Missing Value for {val_to_check} for participant " f"{subject_id}"
        )
    return ret_val


def check2(val):
    """https://github.com/FCP-INDI/C-PAC/blob/96db8b0b65ab1d5f55fb3b895855af34d72c17e4/CPAC/utils/utils.py#L745-L746"""
    return val if val == None or val == "" or isinstance(val, str) else int(val)


def try_fetch_parameter(scan_parameters, subject, scan, keys):
    """https://github.com/FCP-INDI/C-PAC/blob/96db8b0b65ab1d5f55fb3b895855af34d72c17e4/CPAC/utils/utils.py#L679-L703"""
    scan_parameters = dict((k.lower(), v) for k, v in scan_parameters.items())
    for key in keys:
        key = key.lower()
        if key not in scan_parameters:
            continue
        if isinstance(scan_parameters[key], dict):
            value = scan_parameters[key][scan]
        else:
            value = scan_parameters[key]
        if value == "None":
            return None
        if value is not None:
            return value
    return None
