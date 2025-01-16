#!/usr/bin/env python
# Copyright (C) 2025  C-PAC Developers

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
"""Inspect inputs and outputs for NodeBlockFunctions."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import ast
from dataclasses import dataclass, field
import importlib
from importlib.resources import files
import inspect
from itertools import chain, product
import os
from pathlib import Path
from typing import Any, cast, Iterable, Optional
from unittest.mock import patch

from traits.trait_errors import TraitError
import yaml

from CPAC.pipeline.engine import template_dataframe
from CPAC.pipeline.nodeblock import NodeBlockFunction
from CPAC.utils.monitoring import UTLOGGER
from CPAC.utils.outputs import Outputs


def import_nodeblock_functions(
    package_name: str, exclude: Optional[list[str]] = None
) -> list[NodeBlockFunction]:
    """
    Import all functions with the @nodeblock decorator from all modules and submodules in a package.

    Parameters
    ----------
    package_name
        The name of the package to import from.

    exclude
        A list of module names to exclude from the import.
    """
    if exclude is None:
        exclude = []
    functions: list[NodeBlockFunction] = []
    package = importlib.import_module(package_name)
    package_path = package.__path__[0]  # Path to the package directory

    for root, _, package_files in os.walk(package_path):
        for file in package_files:
            if file.endswith(".py") and file != "__init__.py":
                # Get the module path
                rel_path = os.path.relpath(os.path.join(root, file), package_path)
                module_name = f"{package_name}.{rel_path[:-3].replace(os.sep, '.')}"
                if module_name in exclude:
                    continue

                # Import the module
                try:
                    with patch.dict(
                        "sys.modules", {exclusion: None for exclusion in exclude}
                    ):
                        module = importlib.import_module(module_name)
                except (ImportError, TraitError, ValueError) as e:
                    UTLOGGER.debug(f"Failed to import {module_name}: {e}")
                    continue
                # Extract nodeblock-decorated functions from the module
                for _name, obj in inspect.getmembers(
                    module, predicate=lambda obj: isinstance(obj, NodeBlockFunction)
                ):
                    functions.append(obj)

    return functions


@dataclass
class ResourceSourceList:
    """A list of resource sources without duplicates."""

    sources: list[str] = field(default_factory=list)

    def __add__(self, other: "str | list[str] | ResourceSourceList") -> list[str]:
        """Add a list of sources to the list."""
        if isinstance(other, str):
            if not other:
                other = "(dummy node)"
            other = [other]
        new_set = {*self.sources, *other}
        return sorted(new_set, key=str.casefold)

    def __contains__(self, item: str) -> bool:
        """Check if a source is in the list."""
        return item in self.sources

    def __delitem__(self, key: int) -> None:
        """Delete a source by index."""
        del self.sources[key]

    def __eq__(self, value: Any) -> bool:
        """Check if the lists of sources are the same."""
        return set(self) == set(value)

    def __getitem__(self, item: int) -> str:
        """Get a source by index."""
        return self.sources[item]

    def __hash__(self) -> int:
        """Get the hash of the list of sources."""
        return hash(self.sources)

    def __iadd__(
        self, other: "str | list[str] | ResourceSourceList"
    ) -> "ResourceSourceList":
        """Add a list of sources to the list."""
        self.sources = self + other
        return self

    def __iter__(self):
        """Iterate over the sources."""
        return iter(self.sources)

    def __len__(self) -> int:
        """Get the number of sources."""
        return len(self.sources)

    def __repr__(self) -> str:
        """Get the reproducable string representation of the sources."""
        return f"ResourceSourceList({(self.sources)})"

    def __reversed__(self) -> list[str]:
        """Get the sources reversed."""
        return list(reversed(self.sources))

    def __setitem__(self, key: int, value: str) -> None:
        """Set a source by index."""
        self.sources[key] = value

    def __sorted__(self) -> list[str]:
        """Get the sources sorted."""
        return sorted(self.sources, key=str.casefold)

    def __str__(self) -> str:
        """Get the string representation of the sources."""
        return str(self.sources)


@dataclass
class ResourceIO:
    """NodeBlockFunctions that use a resource for IO."""

    name: str
    """The name of the resource."""
    output_from: ResourceSourceList | list[str] = field(
        default_factory=ResourceSourceList
    )
    """The functions that output the resource."""
    output_to: ResourceSourceList | list[str] = field(
        default_factory=ResourceSourceList
    )
    """The subdirectory the resource is output to."""
    input_for: ResourceSourceList | list[str] = field(
        default_factory=ResourceSourceList
    )
    """The functions that use the resource as input."""

    def __post_init__(self) -> None:
        """Handle optionals."""
        if isinstance(self.output_from, list):
            self.output_from = ResourceSourceList(self.output_from)
        if isinstance(self.output_to, list):
            self.output_to = ResourceSourceList(self.output_to)
        if isinstance(self.input_for, list):
            self.input_for = ResourceSourceList(self.input_for)

    def __str__(self) -> str:
        """Return string representation for ResourceIO instance."""
        return f"{{{self.name}: {{'input_for': {self.input_for!s}, 'output_from': {self.output_from!s}}}}})"

    def as_dict(self) -> dict[str, list[str]]:
        """Return the ResourceIO as a built-in dictionary type."""
        return {
            k: v
            for k, v in {
                "input_for": [str(source) for source in self.input_for],
                "output_from": [str(source) for source in self.output_from],
                "output_to": [str(source) for source in self.output_to],
            }.items()
            if v
        }


def cli_parser() -> Namespace:
    """Parse command line argument."""
    parser = ArgumentParser(
        description="Inventory resources for C-PAC NodeBlockFunctions.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-o",
        "--output",
        nargs="?",
        help="The output file to write the inventory to.",
        type=Path,
        default=Path("resource_inventory.yaml"),
    )
    return parser.parse_args()


def _flatten_io(io: list[Iterable]) -> list[str]:
    """Given a list of strings or iterables thereof, flatten the list to all strings."""
    if all(isinstance(resource, str) for resource in io):
        return cast(list[str], io)
    while not all(isinstance(resource, str) for resource in io):
        io = list(
            chain.from_iterable(
                [
                    resource if not isinstance(resource, str) else [resource]
                    for resource in io
                ]
            )
        )
    return cast(list[str], io)


class MultipleContext(list):
    """Subclass of list to store multilpe contexts."""


class DirectlySetResources(ast.NodeVisitor):
    """Class to track resources set directly, rather than through NodeBlocks."""

    def __init__(self) -> None:
        """Initialize the visitor."""
        super().__init__()
        self._context: dict[str, Any] = {}
        self.dynamic_resources: dict[str, ResourceSourceList] = {}
        self._history: dict[str, list[Any]] = {}
        self.resources: dict[str, ResourceSourceList] = {}

    def assign_resource(self, resource: str, value: str) -> None:
        """Assign a value to a resource."""
        target = self.dynamic_resources if r".*" in value else self.resources
        resource = str(resource)
        if resource not in target:
            target[resource] = ResourceSourceList()
        target[resource] += value

    @property
    def context(self) -> dict[str, Any]:
        """Return the context."""
        return self._context

    @context.setter
    def context(self, value: tuple[Iterable, Any]) -> None:
        """Set the context."""
        key, _value = value
        if not isinstance(key, str):
            for subkey in key:
                self.context = subkey, _value
        else:
            self._context[key] = _value
            if key not in self._history:
                self._history[key] = []
            self._history[key].append(_value)

    def lookup_context(self, variable: str) -> str | MultipleContext:
        """Plug in variable."""

        def lookup() -> str | list[str]:
            """Look up context."""
            if variable in self.context:
                if self.context[variable] == variable:
                    history = list(self._history[variable])
                    while history and history[-1] == variable:
                        history.pop()
                    if history:
                        return history[-1]
                return self.context[variable]
            return ".*"

        context = lookup()
        if isinstance(context, list):
            context = MultipleContext(context)
        return context

    @staticmethod
    def handle_multiple_contexts(contexts: list[str | list[str]]) -> list[str]:
        """Parse multiple contexts."""
        if isinstance(contexts, list):
            return MultipleContext(
                [
                    "".join(list(ctx))
                    for ctx in product(
                        *[
                            context if isinstance(context, list) else [context]
                            for context in contexts
                        ]
                    )
                ]
            )
        return contexts

    def parse_ast(self, node: Any) -> Any:
        """Parse AST."""
        if not isinstance(node, ast.AST):
            if isinstance(node, str):
                return node
            if not isinstance(node, Iterable):
                return str(node)
            if isinstance(node, ast.Dict):
                return {
                    self.parse_ast(key): self.parse_ast(value)
                    for key, value in dict(zip(node.keys, node.values)).items()
                }
            if isinstance(node, (MultipleContext, list, set, tuple)):
                return type(node)(self.parse_ast(subnode) for subnode in node)
        if isinstance(node, ast.FormattedValue):
            if hasattr(node, "value") and hasattr(node.value, "id"):
                return self.lookup_context(getattr(node.value, "id"))
        if isinstance(node, ast.JoinedStr):
            node_values = [self.parse_ast(value) for value in node.values]
            if any(isinstance(value, MultipleContext) for value in node_values):
                return self.handle_multiple_contexts(node_values)
            return "".join(str(item) for item in node_values)
        if isinstance(node, ast.Dict):
            return {
                self.parse_ast(key): self.parse_ast(value)
                for key, value in dict(zip(node.keys, node.values)).items()
            }
        for attr in ["values", "elts"]:
            if hasattr(node, attr):
                return [self.parse_ast(subnode) for subnode in getattr(node, attr)]
        for attr in ["value", "id"]:
            if hasattr(node, attr):
                return self.parse_ast(getattr(node, attr))
        return r".*"  # wildcard for regex matching

    def visit_Assign(self, node: ast.Assign) -> None:
        """Visit an assignment."""
        value = self.parse_ast(node.value)
        for target in node.targets:
            resource = self.parse_ast(target)
            self.context = resource, value
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        """Visit a function call."""
        if isinstance(node.func, ast.Attribute) and node.func.attr == "set_data":
            value = self.parse_ast(node.args[5])
            if hasattr(node.args[0], "value"):
                resource: str = getattr(node.args[0], "value")
            elif hasattr(node.args[0], "id"):
                resource = self.lookup_context(getattr(node.args[0], "id"))
                if isinstance(resource, MultipleContext):
                    if len(resource) == len(value):
                        for k, v in zip(resource, value):
                            self.assign_resource(k, v)
                    else:
                        for resource_context in resource:
                            self.assign_resource(resource_context, value)
                    self.generic_visit(node)
                    return
            elif isinstance(node.args[0], ast.JoinedStr):
                resource = self.parse_ast(node.args[0])
            else:
                self.generic_visit(node)
                return
            self.assign_resource(resource, value)
        self.generic_visit(node)

    def visit_For(self, node: ast.For) -> None:
        """Vist for loop."""
        # This is probably too specific,
        # will need to be updated if we add more out-of-nodeblock settings.
        target = self.parse_ast(node.target)
        if (
            hasattr(node.iter, "func")
            and hasattr(node.iter.func, "value")
            and hasattr(node.iter.func.value, "id")
        ):
            context = self.context.get(self.parse_ast(node.iter.func.value.id), ".*")
            if isinstance(target, list) and isinstance(context, dict):
                self.context = target[0], list(context.keys())
        else:
            self.context = target, self.parse_ast(node.iter)
        self.generic_visit(node)


def find_directly_set_resources(package_name: str) -> dict[str, ResourceSourceList]:
    """Find all resources set explicitly via :pyy:method:`~CPAC.pipeline.engine.ResourcePool.set_data`.

    Parameters
    ----------
    package_name
        The name of the package to search for resources.

    Returns
    -------
    dict
        A dictionary containing the name of the resource and the name of the functions that set it.
    """
    resources: dict[str, ResourceSourceList] = {}
    dynamic_resources: dict[str, ResourceSourceList] = {}
    for dirpath, _, filenames in os.walk(str(files(package_name))):
        for filename in filenames:
            if filename.endswith(".py"):
                filepath = os.path.join(dirpath, filename)
                with open(filepath, "r", encoding="utf-8") as file:
                    tree = ast.parse(file.read(), filename=filepath)
                    directly_set = DirectlySetResources()
                    directly_set.visit(tree)
                    for resource in directly_set.resources:
                        if resource not in resources:
                            resources[resource] = ResourceSourceList()
                        resources[resource] += directly_set.resources[resource]
                    for resource in directly_set.dynamic_resources:
                        if resource not in dynamic_resources:
                            dynamic_resources[resource] = ResourceSourceList()
                        dynamic_resources[resource] += directly_set.dynamic_resources[
                            resource
                        ]
    # for dynamic_key, dynamic_value in dynamic_resources.items():
    #     dynamic_resource = re.compile(dynamic_key)
    #     for resource in resources.keys():
    #         if dynamic_resource.search(resource):
    #             resources[resource] += dynamic_value
    return resources


def resource_inventory(package: str = "CPAC") -> dict[str, ResourceIO]:
    """Gather all inputs and outputs for a list of NodeBlockFunctions."""
    resources: dict[str, ResourceIO] = {}
    # Node block function inputs and outputs
    for nbf in import_nodeblock_functions(
        package,
        [
            # No nodeblock functions in these modules that dynamically isntall torch
            "CPAC.unet.__init__",
            "CPAC.unet._torch",
        ],
    ):
        nbf_name = f"{nbf.__module__}.{nbf.__qualname__}"
        if hasattr(nbf, "inputs"):
            for nbf_input in _flatten_io(cast(list[Iterable], nbf.inputs)):
                if nbf_input:
                    if nbf_input not in resources:
                        resources[nbf_input] = ResourceIO(
                            nbf_input, input_for=[nbf_name]
                        )
                    else:
                        resources[nbf_input].input_for += nbf_name
        if hasattr(nbf, "outputs"):
            for nbf_output in _flatten_io(cast(list[Iterable], nbf.outputs)):
                if nbf_output:
                    if nbf_output not in resources:
                        resources[nbf_output] = ResourceIO(
                            nbf_output, output_from=[nbf_name]
                        )
                    else:
                        resources[nbf_output].output_from += nbf_name
    # Template resources set from pipeline config
    templates_from_config_df = template_dataframe()
    for _, row in templates_from_config_df.iterrows():
        output_from = f"pipeline configuration: {row.Pipeline_Config_Entry}"
        if row.Key not in resources:
            resources[row.Key] = ResourceIO(row.Key, output_from=[output_from])
        else:
            resources[row.Key].output_from += output_from
    # Hard-coded resources
    for resource, functions in find_directly_set_resources(package).items():
        if resource not in resources:
            resources[resource] = ResourceIO(resource, output_from=functions)
        else:
            resources[resource].output_from += functions
    # Outputs
    for _, row in Outputs.reference.iterrows():
        if row.Resource not in resources:
            resources[row.Resource] = ResourceIO(
                row.Resource, output_to=[row["Sub-Directory"]]
            )
        else:
            resources[row.Resource].output_to += row["Sub-Directory"]
    return dict(sorted(resources.items(), key=lambda item: item[0].casefold()))


def dump_inventory_to_yaml(inventory: dict[str, ResourceIO]) -> str:
    """Dump NodeBlock Interfaces to a YAML string."""
    return yaml.dump(
        {key: value.as_dict() for key, value in inventory.items()}, sort_keys=False
    )


def main() -> None:
    """Save the NodeBlock inventory to a file."""
    args = cli_parser()
    with args.output.open("w") as file:
        file.write(dump_inventory_to_yaml(resource_inventory("CPAC")))


if __name__ == "__main__":
    main()
