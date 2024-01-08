# Copyright (c) 2009-2013, NIPY Developers
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the NIPY Developers nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Modifications Copyright (C) 2022-2023  C-PAC Developers

# This file is part of C-PAC.
import toml

ga_tracker = 'UA-19224662-10'


# Below are legacy reconstructions of previously dynamically defined metadata
# that are now defined in pyproject.toml.
# These are kept for backwards compatibility.

def _read_pyproject_toml():
    """Read pyproject.toml file."""
    with open('pyproject.toml') as f:
        pyproject = toml.load(f)
    return pyproject


_pyproject_toml = _read_pyproject_toml()
__version__ = _pyproject_toml['tool']['poetry']['version']
REQUIREMENTS = list(_pyproject_toml['tool']['poetry']['dependencies'].keys()) + \
    list(_pyproject_toml['tool']['poetry']['dev-dependencies'].keys())
UNET_REQUIREMENTS = list(_pyproject_toml['tool']['poetry']['dev-dependencies'].keys())
