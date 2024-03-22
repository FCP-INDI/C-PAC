# STATEMENT OF CHANGES:
#     This file is derived from sources licensed under the Apache-2.0 terms,
#     and this file has been changed.

# CHANGES:
#     * Symlinks data instead of copying
#     * Uses pathlib instead of py.path
#     * Removes namespace configurations (we aren't using them here)
#     * Docstrings updated accordingly
#     * Style modifications

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

# Modifications Copyright (C) 2024  C-PAC Developers

# This file is part of C-PAC.
"""Configure Pytest environment for interfaces tests.

Modified from https://github.com/nipy/nipype/blob/a17de8e/nipype/conftest.py
"""

from contextlib import contextmanager
from os import chdir, getcwd
from pathlib import Path
from shutil import copytree, rmtree
from tempfile import mkdtemp

from pytest import fixture
import nipype

from CPAC.utils.typing import PATHSTR

NIPYPE_DATADIR = Path(nipype.__file__).parent / "testing/data"
TEMP_FOLDER = Path(mkdtemp())
DATA_DIR = TEMP_FOLDER / "data"
copytree(NIPYPE_DATADIR, DATA_DIR, symlinks=True)


@contextmanager
def cwd(path: PATHSTR):
    """Change the current working directory."""
    orig_wd = getcwd()
    try:
        chdir(path)
        yield
    finally:
        chdir(orig_wd)


@fixture(autouse=True)
def _docdir(request):
    """Change doctest working directory to one with test data available.

    Grabbed from https://stackoverflow.com/a/46991331
    """
    # Trigger ONLY for the doctests.
    doctest_plugin = request.config.pluginmanager.getplugin("doctest")
    if isinstance(request.node, doctest_plugin.DoctestItem):
        # Chdir only for the duration of the test.
        with cwd(DATA_DIR):
            yield
    else:
        # For normal tests, we have to yield, since this is a yield-fixture.
        yield


def pytest_unconfigure(config):
    """Delete temp folder after session is finished."""
    rmtree(TEMP_FOLDER)
