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
"""Custom registration exceptions"""


class BadRegistrationError(ValueError):
    """Exception for when a QC measure for a registration falls below a
    specified threshold"""
    def __init__(self, *args, metric=None, value=None, threshold=None,
                 **kwargs):
        """
        Parameters
        ----------
        metric : str
            QC metric

        value : float
            calculated QC value

        threshold : float
            specified threshold
        """
        msg = "Registration failed quality control"
        if all(arg is not None for arg in (metric, value, threshold)):
            msg += f" ({metric}: {value} < {threshold})"
        msg += "."
        super().__init__(msg, *args, **kwargs)
