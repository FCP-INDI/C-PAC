""" This file contains defines parameters for CPAC that we use to fill
settings in setup.py, the CPAC top-level docstring, and for building the
docs.  In setup.py in particular, we exec this file, so it cannot import CPAC.
This script was borrowed from and inspired by nipype's info.py file.
"""


# CPAC version information.  An empty _version_extra corresponds to a
# full release.  '.dev' as a _version_extra string means this is a development
# version
_version_major = 0
_version_minor = 3
_version_micro = 8
_version_extra = '.dev'

def get_cpac_gitversion():
    """CPAC version as reported by the last commit in git

    Returns
    -------
    None or str
      Version of NiPype according to git.
    """
    import os
    import subprocess
    try:
        import CPAC
        gitpath = os.path.realpath(os.path.join(os.path.dirname(CPAC.__file__),
                                                os.path.pardir))
    except:
        gitpath = os.getcwd()
    gitpathgit = os.path.join(gitpath, '.git')
    if not os.path.exists(gitpathgit):
        return None
    ver = None
    try:
        o, _ = subprocess.Popen('git describe --always', shell=True, cwd=gitpath,
                                stdout=subprocess.PIPE).communicate()
    except Exception:
        try:
            o, _ = subprocess.Popen('git describe --always', shell=True, cwd=gitpath,
                                    stdout=subprocess.PIPE).communicate()
        except Exception:
            pass
        pass
    else:
        ver = o.strip().split('-')[-1]
    return ver

if '.dev' in _version_extra:
    gitversion = get_cpac_gitversion()
    if gitversion:
        _version_extra = '.' + gitversion + '-' + 'dev'

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
__version__ = "%s.%s.%s%s" % (_version_major,
                              _version_minor,
                              _version_micro,
                              _version_extra)

CLASSIFIERS = ["Development Status :: 3 - Alpha", # 3 - Alpha
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: BSD License", # TODO: check if this is true
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

description  = 'Configural Pipeline for the Analysis of Connectomes'

# Note: this long_description is actually a copy/paste from the top-level
# README.md, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
long_description = \
"""
============================================================
C-PAC: Configurable Pipeline for the Analysis of Connectomes
============================================================

A configurable, open-source, Nipype-based, automated processing pipeline for resting state fMRI data. 
Designed for use by both novice users and experts, C-PAC brings the power, flexibility and elegance 
of Nipype to users in a plug-and-play fashion; no programming required.

Website
-------

CPAC website is located here:  http://fcp-indi.github.com/


Documentation
-------------

User documentation can be found here: http://fcp-indi.github.com/docs/user/index.html

Developer documention can ne found here: http://fcp-indi.github.com/docs/developer/index.html


Dicussion Forum
---------------

CPAC Discussion forum is located here: http://www.nitrc.org/forum/forum.php?forum_id=3567

Troubleshooting and Help
------------------------

This is an alpha version of CPAC, which means that it is still under active development. As such, although we have done our best to ensure a stable pipeline, 
there will likely still be a few bugs that we did not catch. If you find a bug, have a question that is not answered in the User Guide, or would like to suggest a new feature, 
please create an issue on CPAC github issue page: https://github.com/FCP-INDI/C-PAC/issues?state=open
"""

# versions
CYTHON_MIN_VERSION      = '0.12.1'
MATPLOTLIB_MIN_VERSION  = '1.1'
JINJA_MIN_VERSION       = '2.6'
PYLOCKFILE_MIN_VERSION  = '0.9'
PYYAML_MIN_VERSION      = '2.0'

NAME                = 'CPAC'
MAINTAINER          = "cpac developers"
MAINTAINER_EMAIL    = "XXX"
DESCRIPTION         = description
LONG_DESCRIPTION    = long_description
URL                 = "http://fcp-indi.github.io"
DOWNLOAD_URL        = "https://github.com/FCP-INDI/C-PAC"
LICENSE             = "BSD license" # TODO: figure out if this is right
CLASSIFIERS         = CLASSIFIERS
AUTHOR              = "cpac developmers"
AUTHOR_EMAIL        = "XXX"
PLATFORMS           = "OS Independent"
MAJOR               = _version_major
MINOR               = _version_minor
MICRO               = _version_micro
ISRELEASE           = _version_extra == ''
VERSION             = __version__
REQUIRES            = ["matplotlib (>=1.2)", "Jinja2 (>=2.6)", "pylockfile (>=0.9)",
                       "pyyaml (>=3.0)"]
STATUS              = 'stable'
