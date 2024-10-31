#!/usr/bin/env python3
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
"""Update issue templates."""

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
import re
import sys
from typing import Any, Optional

from git import Repo
from github import Github
from github.GithubException import UnknownObjectException
import yaml

C_PAC_ROOT = str(Path().parent.parent)
sys.path.append(C_PAC_ROOT)
from CPAC.pipeline import AVAILABLE_PIPELINE_CONFIGS  # noqa: E402


@dataclass
class TemplateItem:
    """Store index and contents together."""

    index: int
    contents: dict[str, Any]


@dataclass
class Template:
    """Store path and contents together."""

    path: Path
    contents: dict[str, Any]

    def get_template_item(self, key: str) -> TemplateItem:
        """Given a key, return the index and contents of that item in the body."""
        for i, item in enumerate(self.contents.get("body", [])):
            if item.get("id") == key:
                return TemplateItem(i, item)
        return TemplateItem(-1, {})

    def update_c_pac_version(self) -> None:
        """Update placeholder version of C-PAC in issue template.

        In priority order, tries:

        1. The remote ``origin`` repository on GitHub's latest tag
        2. ``FCP-INDI/C-PAC`` on GitHub's latest tag
        3. The current local version
        """
        item = self.get_template_item("c-pac-version")
        repo = "FCP-INDI/C-PAC"
        try:
            match: Optional[re.Match] = re.search(
                r"github\.com[:/](.+/.+?)(?:\.git)?$",
                Repo(C_PAC_ROOT).remotes.origin.url,
            )
        except AttributeError:
            match = None
        if match:
            repo = match.groups()[0]
        try:
            version: Optional[str] = (
                Github().get_repo(repo).get_latest_release().tag_name
            )
        except UnknownObjectException:
            version = None
        if not version:
            from CPAC import __version__ as version
        if "attributes" in item.contents:
            self.contents["body"][item.index]["attributes"]["placeholder"] = version

    def update_preconfigs(self) -> None:
        """Update list of available preconfigs."""
        options: list[dict[str, str]] = [
            {"label": preconfig}
            for preconfig in [
                "default",
                *[
                    preconfig
                    for preconfig in AVAILABLE_PIPELINE_CONFIGS
                    if preconfig not in ["blank", "default"]
                    and "deprecated" not in preconfig
                ],
            ]
        ]
        item = self.get_template_item("preconfig")
        if "attributes" in item.contents:
            self.contents["body"][item.index]["attributes"]["options"] = options

    def write(self) -> None:
        """Write contents back to file."""
        with self.path.open("w", encoding="utf-8") as _f:
            yaml.dump(self.contents, _f)


def load_issue_template(template_name: str) -> Template:
    """Get an issue template's filepath and contents."""
    template_path: Path = Path(str(files("CPAC"))).parent.joinpath(
        f".github/ISSUE_TEMPLATE/{template_name}.yaml"
    )
    with template_path.open("r", encoding="utf-8") as _f:
        contents = yaml.safe_load(_f)
    return Template(template_path, contents)


def update_all(template_name: str) -> None:
    """Load, update, and write updated issue template."""
    template = load_issue_template(template_name)
    template.update_preconfigs()
    template.update_c_pac_version()
    template.write()


if __name__ == "__main__":
    for template in ["1-bug_report"]:
        update_all(template)
