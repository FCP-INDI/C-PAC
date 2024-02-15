# Copyright (C) 2024  C-PAC Developers

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
"""Test Function interface."""
from pytest import raises

from CPAC.utils.interfaces.function.function import Function


def test_autologger():
    """Test autoassignment of `logger` and `iflogger`."""

    def faux_fxn():
        return luigi_mario  # noqa: F821

    interface = Function(function=faux_fxn)
    with raises(NameError) as name_error:
        interface.run()
        assert "name 'luigi_mario' is not defined" in str(name_error.value)

    def faux_fxn():
        return logger, iflogger

    interface = Function(function=faux_fxn, outputs=["logger", "iflogger"])
    res = interface.run()
    logger, iflogger = res.outputs.out
    assert logger.name == "nipype.workflow"
    assert iflogger.name == "nipype.interface"
