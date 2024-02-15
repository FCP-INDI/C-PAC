# from https://github.com/nipy/nipype/blob/0.13.1/nipype/interfaces/utility/wrappers.py

# CHANGES:
#     * Adds `as_module` argument and property
#     * Adds `sig_imports` decorator
#     * Automatically assigns logger and iflogger variables in function nodes

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
from builtins import bytes, str
import inspect
from typing import Callable, List

from nipype import logging
from nipype.interfaces.base import (
    isdefined,
)
from nipype.interfaces.io import add_traits, IOBase
from nipype.interfaces.utility.wrappers import Function as NipypeFunction
from nipype.utils.filemanip import ensure_list
from nipype.utils.functions import create_function_from_source, getsource

from CPAC.utils.docs import outdent_lines

iflogger = logging.getLogger("nipype.interface")


class Function(NipypeFunction):
    """Can automatically set a module name on the interface.

    Automatically sets ``iflogger`` and ``logger`` variables to "nipype.interface" and
    "nipype.workflow" respectively.
    """

    def __init__(
        self,
        input_names=None,
        output_names="out",
        function=None,
        imports=None,
        as_module=False,
        **inputs,
    ):
        """Initialize a :py:func`~CPAC.utils.interfaces.function.Function` interface.

        Parameters
        ----------
        input_names : single str or list or None
            names corresponding to function inputs
            if ``None``, derive input names from function argument names
        output_names : single str or list
            names corresponding to function outputs (default: 'out').
            if list of length > 1, has to match the number of outputs
        function : callable
            callable python object. must be able to execute in an
            isolated namespace (possibly in concert with the ``imports``
            parameter)
        imports : list of strings
            list of import statements that allow the function to execute
            in an otherwise empty namespace. If these collide with
            imports defined via the :py:meth:`Function.sig_imports`
            decorator, the imports given as a parameter here will take
            precedence over those from the decorator.
        """
        super(IOBase, self).__init__(**inputs)
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
                module = inspect.getmodule(function).__name__
                full_name = "%s.%s" % (module, function.__name__)
                self.inputs.function_str = full_name
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
                    fninfo = create_function_from_source(function, imports).__code__
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
        self.imports = [
            *imports,
            "from CPAC.utils.monitoring.custom_logging import getLogger",
            "iflogger = getLogger('nipype.interface')",
            "logger = getLogger('nipype.workflow')",
        ]
        self._out = {}
        for name in self._output_names:
            self._out[name] = None

    @staticmethod
    def sig_imports(imports: List[str]) -> Callable:
        """Set an ``ns_imports`` attribute on a function for Function-node functions.

        This can be useful for classes needed for decorators, typehints
        and for avoiding redefinitions.

        Parameters
        ----------
        imports : list of str
            import statements to import the function in an otherwise empty
            namespace. If these collide with imports defined via the
            :py:meth:`Function.__init__` initialization method, the
            imports given as a parameter here will be overridden by
            those from the initializer.

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
        >>> calc_fdj.imports  # doctest: +NORMALIZE_WHITESPACE
        ['from CPAC.utils.interfaces.function import Function',
         'import os',
         'import sys',
         'from typing import Optional',
         'import numpy as np',
         'from CPAC.utils.pytest import skipif',
         'from CPAC.utils.typing import LITERAL, TUPLE']
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
            if self.as_module:
                module = inspect.getmodule(new).__name__
                full_name = "%s.%s" % (module, new.__name__)
                self.inputs.function_str = full_name
            elif hasattr(new, "__call__"):
                function_source = getsource(new)
                fninfo = new.__code__
            elif isinstance(new, (str, bytes)):
                function_source = new
                fninfo = create_function_from_source(new, self.imports).__code__
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
        if self.as_module:
            import importlib

            pieces = self.inputs.function_str.split(".")
            module = ".".join(pieces[:-1])
            function = pieces[-1]
            try:
                function_str = inspect.getsource(
                    getattr(importlib.import_module(module), function)
                )
            except ImportError as import_error:
                msg = f"Could not import module: {self.inputs.function_str}"
                raise RuntimeError(msg) from import_error
        else:
            function_str = self.inputs.function_str
        function_handle = create_function_from_source(function_str, self.imports)

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
