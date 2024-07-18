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
"""Classes and decorator for :py:class:`NodeBlock`\u200bs and :py:class:`NodeBlockFunction`\u200bs."""

from typing import Any, Callable, Optional, TYPE_CHECKING

import yaml
from nipype import config, logging  # type: ignore [import-untyped]
from nipype.pipeline.engine import Workflow  # type: ignore[import-untyped]

from CPAC.utils.configuration.configuration import Configuration
from CPAC.utils.monitoring import (
    WFLOGGER,
)

if TYPE_CHECKING:
    from CPAC.pipeline.engine.resource import ResourceData, StratPool

NODEBLOCK_INPUTS = list[str | list | tuple]
NODEBLOCK_OUTPUTS = list[str] | dict[str, Any]
PIPELINE_BLOCKS = list["NodeBlockFunction | PIPELINE_BLOCKS"]


class NodeBlockFunction:
    """Store a reference to the nodeblock function and all of its meta-data."""

    def __init__(
        self,
        func: Callable,
        name: str,
        config: Optional[list[str]] = None,
        switch: Optional[list[str] | list[list[str]]] = None,
        option_key: Optional[str | list[str]] = None,
        option_val: Optional[str | list[str]] = None,
        inputs: Optional[NODEBLOCK_INPUTS] = None,
        outputs: Optional[NODEBLOCK_OUTPUTS] = None,
    ) -> None:
        self.func = func
        """`Nodeblock` function reference."""
        self.name: str = name
        """Used in the graph and logging to identify the :py:class:`NodeBlock` and its component :py:class:`~nipype.pipeline.engine.Node`\u200bs."""
        self.config: Optional[list[str]] = config
        """
        Indicates the nested keys in a C-PAC pipeline :py:class:`Configuration`
        should configure a `NodeBlock` built from this function. If `config` is set to
        `None`, then all other configuration-related entities must be specified from the
        root of the configuration.
        """
        self.switch: Optional[list[str] | list[list[str]]] = switch
        """
        Indicates any keys that should evaluate to `True` for this :py:class:`NodeBlock`
        to be active. A list of lists of strings indicates multiple `switch`\u200bes
        that must all be `True` to run, and is currently only an option if `config` is
        set to `None`.
        """
        self.option_key: Optional[str | list[str]] = option_key
        """
        Indicates the nested keys (starting at the nested key indicated by `config`)
        that should configure this :py:class:`NodeBlock`.
        """
        self.option_val: Optional[str | list[str]] = option_val
        """Indicates values for which this :py:class:`NodeBlock` should be active."""
        self.inputs: list[str | list | tuple] = inputs if inputs else []
        """:py:class:`~CPAC.pipeline.engine.resource.ResourcePool` keys indicating
        resources needed for the :py:class:`NodeBlock`\u200b's functionality."""
        self.outputs: list[str] | dict[str, Any] = outputs if outputs else []
        """
        :py:class:`~CPAC.pipeline.engine.resource.ResourcePool` keys indicating
        resources generated or updated by the `NodeBlock`, optionally including metadata
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

    def __call__(
        self,
        wf: Workflow,
        cfg: Configuration,
        strat_pool: "StratPool",
        pipe_num: Optional[int | str],
        opt: Optional[str] = None,
    ) -> tuple[Workflow, dict[str, "ResourceData"]]:
        """Call a `NodeBlockFunction`.

        All `NodeBlockFunction`\u200bs have the same signature.
        """
        return self.func(wf, cfg, strat_pool, pipe_num, opt)

    def legacy_nodeblock_dict(self):
        """Return :py:class:`NodeBlock` metadata as a dictionary.

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
        """Return reproducible string representation of a `NodeBlockFunction`."""
        return (
            f"NodeBlockFunction({self.func.__module__}."
            f'{self.func.__name__}, "{self.name}", '
            f"config={self.config}, switch={self.switch}, "
            f"option_key={self.option_key}, option_val="
            f"{self.option_val}, inputs={self.inputs}, "
            f"outputs={self.outputs})"
        )

    def __str__(self) -> str:
        """Return string representation of a `NodeBlockFunction`."""
        return f"NodeBlockFunction({self.name})"


class NodeBlock:
    """A :py:class:`Workflow` subgraph composed of :py:class:`NodeBlockFunction`\u200bs."""

    def __init__(
        self,
        node_block_functions: NodeBlockFunction | PIPELINE_BLOCKS,
        debug: bool = False,
    ) -> None:
        """Create a `NodeBlock` from a list of py:class:`NodeBlockFunction`\u200bs."""
        if not isinstance(node_block_functions, list):
            node_block_functions = [node_block_functions]

        self.node_blocks: dict[str, Any] = {}

        for node_block_function in node_block_functions:  # <---- sets up the NodeBlock object in case you gave it a list of node blocks instead of a single one - for option forking.
            self.input_interface = []
            if isinstance(node_block_function, tuple):
                self.input_interface = node_block_function[1]
                node_block_function = node_block_function[0]  # noqa: PLW2901
                if not isinstance(self.input_interface, list):
                    self.input_interface = [self.input_interface]

            if not isinstance(node_block_function, NodeBlockFunction):
                # If the object is a plain function `__name__` will be more useful than `str()`
                obj_str = (
                    node_block_function.__name__  # type: ignore [attr-defined]
                    if hasattr(node_block_function, "__name__")
                    else str(node_block_function)
                )
                msg = f'Object is not a nodeblock: "{obj_str}"'
                raise TypeError(msg)

            name = node_block_function.name
            self.name = name
            self.node_blocks[name] = {}

            if self.input_interface:
                for interface in self.input_interface:
                    for orig_input in node_block_function.inputs:
                        if isinstance(orig_input, tuple):
                            list_tup = list(orig_input)
                            if interface[0] in list_tup:
                                list_tup.remove(interface[0])
                                list_tup.append(interface[1])
                                node_block_function.inputs.remove(orig_input)
                                node_block_function.inputs.append(tuple(list_tup))
                        elif orig_input == interface[0]:
                            node_block_function.inputs.remove(interface[0])
                            node_block_function.inputs.append(interface[1])

            for key, val in node_block_function.legacy_nodeblock_dict().items():
                self.node_blocks[name][key] = val

            self.node_blocks[name]["block_function"] = node_block_function

            # TODO: fix/replace below
            self.outputs: dict[str, Optional[str]] = {}
            for out in node_block_function.outputs:
                self.outputs[out] = None

            self.options: list[str] | dict[str, Any] = ["base"]
            if node_block_function.outputs is not None:
                self.options = node_block_function.outputs

            WFLOGGER.info("Connecting %s...", name)
            if debug:
                config.update_config({"logging": {"workflow_level": "DEBUG"}})
                logging.update_logging(config)
                WFLOGGER.debug(
                    '"inputs": %s\n\t "outputs": %s%s',
                    node_block_function.inputs,
                    list(self.outputs.keys()),
                    f'\n\t"options": {self.options}'
                    if self.options != ["base"]
                    else "",
                )
                config.update_config({"logging": {"workflow_level": "INFO"}})
                logging.update_logging(config)

    def check_output(self, outputs: NODEBLOCK_OUTPUTS, label: str, name: str) -> None:
        """Check if a label is listed in a `NodeBlock`\u200b's `outputs`.

        Raises
        ------
        NameError
            If a mismatch is found.
        """
        if label not in outputs:
            msg = (
                f'\n[!] Output name "{label}" in the block '
                "function does not match the outputs list "
                f'{outputs} in Node Block "{name}"\n'
            )
            raise NameError(msg)

    @staticmethod
    def list_blocks(
        pipeline_blocks: PIPELINE_BLOCKS, indent: Optional[int] = None
    ) -> str:
        """List :py:class:`NodeBlockFunction`\u200bs line by line.

        Parameters
        ----------
        pipeline_blocks
            list of :py:class:`NodeBlockFunction`\u200bs

        indent
            number of spaces after a tab indent

        Returns
        -------
        str
            formatted list of :py:class:`NodeBlockFunction`\u200bs
        """
        blockstring = yaml.dump(
            [
                getattr(
                    block,
                    "__name__",
                    getattr(
                        block,
                        "name",
                        yaml.safe_load(NodeBlock.list_blocks(list(block)))
                        if isinstance(block, (tuple, list, set))
                        else str(block),
                    ),
                )
                for block in pipeline_blocks
            ]
        )
        if isinstance(indent, int):
            blockstring = "\n".join(
                [
                    "\t" + " " * indent + line.replace("- - ", "- ")
                    for line in blockstring.split("\n")
                ]
            )
        return blockstring


def nodeblock(
    name: Optional[str] = None,
    config: Optional[list[str]] = None,
    switch: Optional[list[str] | list[list[str]]] = None,
    option_key: Optional[str | list[str]] = None,
    option_val: Optional[str | list[str]] = None,
    inputs: Optional[NODEBLOCK_INPUTS] = None,
    outputs: Optional[list[str] | dict[str, Any]] = None,
):
    """Define a :py:class:`NodeBlockFunction`\u200b.

    Connections to the pipeline :py:class:`Configuration` and to other :py:class:`NodeBlockFunction`\u200bs.

    Parameters
    ----------
    name
        Used in the graph and logging to identify the :py:class:`NodeBlock` and its
        component :py:class:`~nipype.pipeline.engine.Node`\u200bs.
        The :py:class:`NodeBlockFunction`\u200b's `.__name__` is used if `name` is not
        provided.

    config
        Indicates the nested keys in a C-PAC pipeline :py:class:`Configuration` should
        configure a :py:class:`NodeBlock` built from this
        :py:class:`NodeBlockFunction`\u200b. If `config` is set to `None`, then all other
        :py:class:`Configuration`\u200b-related entities must be specified from the root
        of the :py:class:`Configuration`\u200b.

    switch
        Indicates any keys that should evaluate to `True` for this :py:class:`NodeBlock`
        to be active. A list of lists of strings indicates multiple switches that must
        all be `True` to run, and is currently only an option if config is set to
        `None`.

    option_key
        Indicates the nested keys (starting at the nested key indicated by `config`)
        that should configure this :py:class:`NodeBlock`\u200b.

    option_val
        Indicates values for which this :py:class:`NodeBlock` should be active.

    inputs
        ResourcePool keys indicating files needed for the :py:class:`NodeBlock`\u200b's
        functionality.

    outputs
        :py:class:`~CPAC.pipeline.engine.resource.ResourcePool` keys indicating files
        generated or updated by the :py:class:`NodeBlock`, optionally including metadata
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
