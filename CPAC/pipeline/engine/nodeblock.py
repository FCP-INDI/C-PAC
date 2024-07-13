# Copyright (C) 2023-2024  C-PAC Developers

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
"""Class and decorator for NodeBlock functions."""

from typing import Any, Callable, Optional

NODEBLOCK_INPUTS = list[str | list | tuple]


class NodeBlockFunction:
    """Store a reference to the nodeblock function and all of its meta-data."""

    def __init__(
        self,
        func: Callable,
        name: Optional[str] = None,
        config: Optional[list[str]] = None,
        switch: Optional[list[str] | list[list[str]]] = None,
        option_key: Optional[str | list[str]] = None,
        option_val: Optional[str | list[str]] = None,
        inputs: Optional[NODEBLOCK_INPUTS] = None,
        outputs: Optional[list[str] | dict[str, Any]] = None,
    ) -> None:
        self.func = func
        """Nodeblock function reference."""
        self.name: Optional[str] = name
        """Used in the graph and logging to identify the NodeBlock and its component nodes."""
        self.config: Optional[list[str]] = config
        """
        Indicates the nested keys in a C-PAC pipeline configuration should configure a NodeBlock built from this
        function. If config is set to ``None``, then all other configuration-related entities must be specified from the
        root of the configuration.
        """
        self.switch: Optional[list[str] | list[list[str]]] = switch
        """
        Indicates any keys that should evaluate to True for this NodeBlock to be active. A list of lists of strings
        indicates multiple switches that must all be True to run, and is currently only an option if config is set to
        ``None``.
        """
        self.option_key: Optional[str | list[str]] = option_key
        """
        Indicates the nested keys (starting at the nested key indicated by config) that should configure this NodeBlock.
        """
        self.option_val: Optional[str | list[str]] = option_val
        """Indicates values for which this NodeBlock should be active."""
        if inputs is None:
            inputs = []
        self.inputs: list[str | list | tuple] = inputs
        """ResourcePool keys indicating resources needed for the NodeBlock's functionality."""
        self.outputs: Optional[list[str] | dict[str, Any]] = outputs
        """
        ResourcePool keys indicating resources generated or updated by the NodeBlock, optionally including metadata
        for the outputs' respective sidecars.
        """

        # Forward function attributes similar to functools.update_wrapper:
        # https://docs.python.org/3/library/functools.html#functools.update_wrapper
        self.__module__ = func.__module__
        self.__name__ = func.__name__
        self.__qualname__ = func.__qualname__
        self.__annotations__ = func.__annotations__
        self.__doc__ = "".join(
            [
                _.replace("        ", "")
                for _ in [func.__doc__, "", "", NodeBlockFunction.__call__.__doc__]
                if _ is not None
            ]
        ).rstrip()

    # all node block functions have this signature
    def __call__(self, wf, cfg, strat_pool, pipe_num, opt=None):
        """

        Parameters
        ----------
        wf : ~nipype.pipeline.engine.workflows.Workflow

        cfg : ~CPAC.utils.configuration.Configuration

        strat_pool

        pipe_num : int

        opt : str, optional

        Returns
        -------
        wf : ~nipype.pipeline.engine.workflows.Workflow

        out : dict
        """
        return self.func(wf, cfg, strat_pool, pipe_num, opt)

    def legacy_nodeblock_dict(self):
        """Return nodeblock metadata as a dictionary.

        Helper for compatibility reasons.
        """
        return {
            "name": self.name,
            "config": self.config,
            "switch": self.switch,
            "option_key": self.option_key,
            "option_val": self.option_val,
            "inputs": self.inputs,
            "outputs": self.outputs,
        }

    def __repr__(self) -> str:
        """Return reproducible string representation of a NodeBlockFunction."""
        return (
            f"NodeBlockFunction({self.func.__module__}."
            f'{self.func.__name__}, "{self.name}", '
            f"config={self.config}, switch={self.switch}, "
            f"option_key={self.option_key}, option_val="
            f"{self.option_val}, inputs={self.inputs}, "
            f"outputs={self.outputs})"
        )

    def __str__(self) -> str:
        """Return string representation of a NodeBlockFunction."""
        return f"NodeBlockFunction({self.name})"


def nodeblock(
    name: Optional[str] = None,
    config: Optional[list[str]] = None,
    switch: Optional[list[str] | list[list[str]]] = None,
    option_key: Optional[str | list[str]] = None,
    option_val: Optional[str | list[str]] = None,
    inputs: Optional[NODEBLOCK_INPUTS] = None,
    outputs: Optional[list[str] | dict[str, Any]] = None,
):
    """
    Define a node block.

    Connections to the pipeline configuration and to other node blocks.

    Parameters
    ----------
    name
        Used in the graph and logging to identify the NodeBlock and its component nodes.
    config
        Indicates the nested keys in a C-PAC pipeline configuration should configure a NodeBlock built from this
        function. If config is set to ``None``, then all other configuration-related entities must be specified from the
        root of the configuration.
    switch
        Indicates any keys that should evaluate to True for this NodeBlock to be active. A list of lists of strings
        indicates multiple switches that must all be True to run, and is currently only an option if config is set to
        ``None``.
    option_key
        Indicates the nested keys (starting at the nested key indicated by config) that should configure this NodeBlock.
    option_val
        Indicates values for which this NodeBlock should be active.
    inputs
        ResourcePool keys indicating files needed for the NodeBlock's functionality.
    outputs
        ResourcePool keys indicating files generated or updated by the NodeBlock, optionally including metadata
        for the outputs' respective sidecars.
    """
    return lambda func: NodeBlockFunction(
        func,
        name if name is not None else func.__name__,
        config,
        switch,
        option_key,
        option_val,
        inputs,
        outputs,
    )
