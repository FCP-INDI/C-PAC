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
"""This file contains defines parameters for CPAC that we use to fill
settings in setup.py, the CPAC top-level docstring, and for building the
docs.  In setup.py in particular, we exec this file, so it cannot import CPAC.
This script was borrowed from and inspired by nipype's info.py file
(https://github.com/nipy/nipype/blob/08391871/nipype/info.py)."""
# CPAC version information.  An empty _version_extra corresponds to a
# full release.  'dev' as a _version_extra string means this is a development
# version
_version_major = 1
_version_minor = 8
_version_micro = 7
_version_extra = 'dev1'


def get_cpac_gitversion():
    """CPAC version as reported by the last commit in git

    Returns
    -------
    None or str

        Version of C-PAC according to git.
    """
    import os
    import subprocess

    gitpath = os.path.realpath(os.path.join(os.path.dirname(__file__),
                                            os.path.pardir))

    gitpathgit = os.path.join(gitpath, '.git')
    if not os.path.exists(gitpathgit):
        return None

    ver = None

    try:
        o, _ = subprocess.Popen('git describe --always', shell=True,
                                cwd=gitpath, stdout=subprocess.PIPE
                                ).communicate()
    except Exception:
        pass
    else:
        ver = o.decode().strip().split('-')[-1]

    return ver


if 'dev' in _version_extra:
    gitversion = get_cpac_gitversion()
    if gitversion:
        _version_extra = f'{_version_extra}+{gitversion}'


__version__ = "%s.%s.%s" % (_version_major,
                            _version_minor,
                            _version_micro)

if _version_extra:
    __version__ += ".%s" % _version_extra

ga_tracker = 'UA-19224662-10'

CLASSIFIERS = ["Development Status :: 4 - Beta",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# pylint: disable=invalid-name
description = 'Configurable Pipeline for the Analysis of Connectomes'

# Note: this long_description is actually a copy/paste from the top-level
# README.md, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
# pylint: disable=invalid-name
long_description = """
============================================================
C-PAC: Configurable Pipeline for the Analysis of Connectomes
============================================================

A configurable, open-source, Nipype-based, automated processing pipeline for
resting state fMRI data.
Designed for use by both novice users and experts, C-PAC brings the power,
flexibility and elegance of Nipype to users in a plug-and-play fashion; no
programming required.

Website
-------

CPAC website is located here:  https://fcp-indi.github.io/


Documentation
-------------

User documentation can be found here: https://fcp-indi.github.io/docs/user

Developer documention can be found here: https://fcp-indi.github.io/docs/developer

Documentation pertaining to this latest release can be found here: https://fcp-indi.github.io/docs/latest


Dicussion Forum
---------------

CPAC Discussion forum is located here: https://neurostars.org/tag/cpac

Troubleshooting and Help
------------------------

This is a beta version of CPAC, which means that it is still under active
development. As such, although we have done our best to ensure a stable
pipeline, there will likely still be a few bugs that we did not catch. If you
find a bug, have a question that is not answered in the User Guide, or would
like to suggest a new feature, please create an issue on CPAC github issue
page: https://github.com/FCP-INDI/C-PAC/issues?state=open
"""  # noqa: E501
CYTHON_MIN_VERSION = '0.12.1'
NAME = 'CPAC'
MAINTAINER = "C-PAC developers"
MAINTAINER_EMAIL = "CNL@childmind.org"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "https://fcp-indi.github.io"
DOWNLOAD_URL = "https://github.com/FCP-INDI/C-PAC"
LICENSE = "LGPL-3.0-or-later"
AUTHOR = "C-PAC developers"
AUTHOR_EMAIL = "CNL@childmind.org"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
ISRELEASE = _version_extra == ''
VERSION = __version__
STATUS = 'stable'
REQUIREMENTS = [
    "boto3",
    "ciftify @ git+https://git@github.com/fcp-indi/ciftify#egg=ciftify",
    "click",
    "click-aliases",
    "configparser",
    "cython",
    "flowdump==0.1.2",
    "future",
    "INDI-Tools @ git+https://git@github.com/FCP-INDI/INDI-Tools.git#egg=INDI-Tools",
    "lockfile",
    "joblib",
    "matplotlib",
    "networkx",
    "nibabel",
    "nilearn",
    "nipype",
    "nose",
    "numpy",
    "pandas",
    "pathvalidate",
    "patsy",
    "prov",
    "psutil",
    "PyBASC",
    "pybids",
    "pygraphviz",
    "PyPEER @ git+https://git@github.com/ChildMindInstitute/PyPEER.git@6965d2b2bea0fef824e885fec33a8e0e6bd50a97#egg=PyPEER",
    "python-dateutil",
    "pyyaml",
    "scikit-learn",
    "scipy",
    "sdcflows",
    "semver",
    "traits",
    "voluptuous>=0.12.0",
    "xvfbwrapper"
]
UNET_REQUIREMENTS = [
    "torch==1.13.1",
    "torchvision==0.14.1"
]
