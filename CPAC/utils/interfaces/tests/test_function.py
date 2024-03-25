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

from pytest import mark, raises

from CPAC.utils.interfaces.function.function import Function


def faux_fxn(_loggers: bool = True):
    """Require autoassignment (for testing)."""
    if _loggers:
        return WFLOGGER, IFLOGGER  # noqa: F821
    return luigi_mario  # noqa: F821


@mark.parametrize("as_module", [True, False])
def test_autologger(as_module: bool) -> None:
    """Test autoassignment of global Nipype loggers`."""
    interface = Function(
        function=faux_fxn, input_names=["_loggers"], as_module=as_module
    )
    interface.inputs._loggers = False
    with raises(NameError) as name_error:
        interface.run()
        assert "name 'luigi_mario' is not defined" in str(name_error.value)

    interface = Function(
        function=faux_fxn,
        input_names=["_loggers"],
        output_names=["logger", "iflogger"],
        as_module=as_module,
    )
    interface.inputs._loggers = True
    res = interface.run()
    assert res.outputs.logger.name == "nipype.workflow"
    assert res.outputs.iflogger.name == "nipype.interface"
