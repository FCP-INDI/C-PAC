# from https://github.com/nipy/nipype/blob/0.13.1/nipype/interfaces/utility/wrappers.py

# CHANGES:
#     * Removes Python 2 imports
#     * Adds `as_module` argument and property
#     * Adds `sig_imports` decorator
#     * Automatically imports global Nipype loggers in function nodes

# ORIGINAL WORK'S ATTRIBUTION NOTICE:
#     Copyright (c) 2009-2016, Nipype developers

#     Licensed under the Apache License, Version 2.0 (the "License");
#     you may not use this file except in compliance with the License.
#     You may obtain a copy of the License at

#         http://www.apache.org/licenses/LICENSE-2.0

#     Unless required by applicable law or agreed to in writing, software
#     distributed under the License is distributed on an "AS IS" BASIS,
#     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#     See the License for the specific language governing permissions and
#     limitations under the License.

#     Prior to release 0.12, Nipype was licensed under a BSD license.

# Modifications Copyright (C) 2018-2024 C-PAC Developers

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
"""Interface for wrapping Python functions.

Like the built-in nipype Function interace, except includes
- `as_module` to allow module.function name
- `sig_imports` to set necessary imports on function nodes with a decorator
"""

from ast import FunctionDef, parse
from importlib import import_module
import inspect
from typing import Callable, Optional

from nipype.interfaces.base import (
    isdefined,
)
from nipype.interfaces.io import add_traits, IOBase
from nipype.interfaces.utility.wrappers import Function as NipypeFunction
from nipype.utils.filemanip import ensure_list
from nipype.utils.functions import getsource

from CPAC.utils.docs import outdent_lines

_AUTOLOGGING_IMPORTS = [
    "from CPAC.utils.monitoring.custom_logging import FMLOGGER, IFLOGGER, UTLOGGER,"
    " WFLOGGER"
]


def _as_module(fxn: str, ns: dict) -> tuple[str, dict]:
    """Get full module name and namespace."""
    module = inspect.getmodule(fxn).__name__
    return f"{module}.{fxn.__name__}", _module_imports(module, ns, fxn.__name__)


def get_function_name_from_source(function_source: str) -> str:
    r"""Get the name of a function from its source code.

    Parameters
    ----------
    function_source: str
        The source code of the function.

    Returns
    -------
    str
        The name of the function.

    Examples
    --------
    >>> get_function_name_from_source("def fake_function():\n    return")
    'fake_function'
    >>> get_function_name_from_source("not a def")
    Traceback (most recent call last):
       ...
    ValueError: No function definition found in the provided source.
    >>> get_function_name_from_source("class FakeClass:\n    pass")
    Traceback (most recent call last):
       ...
    ValueError: No function definition found in the provided source.
    """
    value_error = ValueError("No function definition found in the provided source.")
    try:
        for node in parse(function_source).body:
            if isinstance(node, FunctionDef):
                return node.name
    except SyntaxError as syntax_error:
        raise value_error from syntax_error
    raise value_error


def create_function_from_source(
    function_source: str, imports: Optional[list[str]] = None, ns: Optional[dict] = None
):
    """Return a function object from a function source.

    Parameters
    ----------
    function_source : unicode string
        unicode string defining a function
    imports : list of strings
        list of import statements in string form that allow the function
        to be executed in an otherwise empty namespace
    ns : dict
        namespace dictionary
    """
    if ns is None:
        ns = {}
    import_keys = []
    try:
        if imports is not None:
            for statement in imports:
                exec(statement, ns)
            import_keys = list(ns.keys())
        exec(function_source, ns)

    except Exception as e:
        msg = f"Error executing function\n{function_source}\n"
        msg += (
            "Functions in connection strings have to be standalone. "
            "They cannot be declared either interactively or inside "
            "another function or inline in the connect string. Any "
            "imports should be done inside the function."
        )
        raise RuntimeError(msg) from e
    ns_funcs = list(set(ns) - {*import_keys, "__builtins__"})
    assert len(ns_funcs) == 1, "Function or inputs are ill-defined"
    return ns[ns_funcs[0]]


class Function(NipypeFunction):
    """Can automatically set a module name on the interface.

    Automatically imports global Nipype loggers.
    """

    def __init__(
        self,
        input_names: Optional[str | list[str]] = None,
        output_names: Optional[str | list[str]] = "out",
        function: Optional[Callable] = None,
        imports: Optional[list[str]] = None,
        as_module: bool = False,
        **inputs,
    ):
        """Initialize a :py:func:`~CPAC.utils.interfaces.function.Function` interface.

        Parameters
        ----------
        input_names
            names corresponding to function inputs
            if ``None``, derive input names from function argument names
        output_names
            names corresponding to function outputs (default: 'out').
            if list of length > 1, has to match the number of outputs
        function
            callable python object. must be able to execute in an
            isolated namespace (possibly in concert with the `imports`
            parameter)
        imports
            list of import statements that allow the function to execute
            in an otherwise empty namespace. If these collide with
            imports defined via the :py:meth:`Function.sig_imports`
            decorator, the imports given as a parameter here will take
            precedence over those from the decorator.
        """
        super(IOBase, self).__init__(**inputs)
        ns = {}
        if imports is None:
            imports = []
        if function:
            if hasattr(function, "ns_imports"):
                # prepend the ns_imports from the decorator to
                # the paramater imports
                _ns_imports = [
                    "from CPAC.utils.interfaces.function import Function",
                    *function.ns_imports,
                ]
                imports = _ns_imports if not imports else [*_ns_imports, *imports]
            if as_module:
                self.inputs.function_str, ns = _as_module(function, ns)
            elif hasattr(function, "__call__"):
                try:
                    self.inputs.function_str = getsource(function)
                except IOError as os_error:
                    msg = (
                        "Interface Function does not accept "
                        "function objects defined interactively "
                        "in a python session"
                    )
                    raise ValueError(msg) from os_error
                else:
                    if input_names is None:
                        fninfo = function.__code__
            elif isinstance(function, (str, bytes)):
                self.inputs.function_str = function
                if input_names is None:
                    fninfo = create_function_from_source(function, imports, ns).__code__
            else:
                msg = "Unknown type of function"
                raise TypeError(msg)
            if input_names is None:
                try:
                    input_names = fninfo.co_varnames[: fninfo.co_argcount]
                except NameError:
                    input_names = []

        self.as_module = as_module
        self.inputs.on_trait_change(self._set_function_string, "function_str")
        self._input_names = ensure_list(input_names)
        self._output_names = ensure_list(output_names)
        add_traits(self.inputs, list(self._input_names))
        self.imports = [*imports, *_AUTOLOGGING_IMPORTS]
        self._out = {}
        for name in self._output_names:
            self._out[name] = None

    @staticmethod
    def sig_imports(imports: list[str]) -> Callable:
        """Set an ``ns_imports`` attribute on a function for Function-node functions.

        This can be useful for classes needed for decorators, typehints
        and for avoiding redefinitions.

        Parameters
        ----------
        imports
            import statements to import the function in an otherwise empty
            namespace. If these collide with imports defined via the
            :py:meth:`Function.__init__` method, the imports given as a parameter here
            will be overridden by those from the initializer.

        Returns
        -------
        func : function

        Examples
        --------
        See the defintion of
        :py:func:`~CPAC.generate_motion_statistics.calculate_FD_J` to see the
        decorator tested here being applied.
        >>> from CPAC.generate_motion_statistics import calculate_FD_J
        >>> calc_fdj = Function(input_names=['in_file', 'calc_from', 'center'],
        ...                     output_names=['out_file'],
        ...                     function=calculate_FD_J,
        ...                     as_module=True)
        >>> calc_fdj.imports == ["from CPAC.utils.interfaces.function import Function",
        ...     *calculate_FD_J.ns_imports, *_AUTOLOGGING_IMPORTS]
        True
        >>> from inspect import signature
        >>> from nipype.utils.functions import (getsource,
        ...     create_function_from_source)
        >>> f = create_function_from_source(getsource(calculate_FD_J),
        ...                                 calc_fdj.imports)
        >>> inspect.signature(calculate_FD_J) == inspect.signature(f)
        True
        """

        def _imports(func: Callable) -> Callable:
            setattr(func, "ns_imports", imports)
            return func

        return _imports

    def _set_function_string(self, obj, name, old, new):
        if name == "function_str":
            ns = {}
            if self.as_module:
                self.inputs.function_str, ns = _as_module(new, ns)
            elif hasattr(new, "__call__"):
                function_source = getsource(new)
                fninfo = new.__code__
            elif isinstance(new, (str, bytes)):
                function_source = new
                fninfo = create_function_from_source(new, self.imports, ns).__code__
            self.inputs.trait_set(
                trait_change_notify=False, **{"%s" % name: function_source}
            )
            # Update input traits
            input_names = fninfo.co_varnames[: fninfo.co_argcount]
            new_names = set(input_names) - set(self._input_names)
            add_traits(self.inputs, list(new_names))
            self._input_names.extend(new_names)

    def _run_interface(self, runtime):
        # Create function handle
        ns = {}
        if self.as_module:
            pieces = self.inputs.function_str.split(".")
            module = ".".join(pieces[:-1])
            function = pieces[-1]
            ns = _module_imports(module, ns, function)
            try:
                function_str = inspect.getsource(
                    getattr(import_module(module), function)
                )
            except ImportError as import_error:
                msg = f"Could not import module: {self.inputs.function_str}"
                raise RuntimeError(msg) from import_error
        else:
            function_str = self.inputs.function_str
        function_handle = create_function_from_source(function_str, self.imports, ns)

        # Get function args
        args = {}
        for name in self._input_names:
            value = getattr(self.inputs, name)
            if isdefined(value):
                args[name] = value
        out = function_handle(**args)
        if len(self._output_names) == 1:
            self._out[self._output_names[0]] = out
        else:
            if isinstance(out, tuple) and (len(out) != len(self._output_names)):
                msg = "Mismatch in number of expected outputs"
                raise RuntimeError(msg)

            for idx, name in enumerate(self._output_names):
                self._out[name] = out[idx]

        return runtime


Function.__doc__ = "\n\n".join(
    [NipypeFunction.__doc__.rstrip(), outdent_lines(Function.__doc__)]
)


def _module_imports(module: str, ns: dict, fxn: str) -> dict:
    """Import module-level imports to a namespace."""
    exec(f"from {module} import *", ns)
    try:
        exec(f"del {fxn}", ns)  # We'll redefine the function itself...
    except NameError:
        pass  # ...unless the function isn't defined in a module
    return ns
