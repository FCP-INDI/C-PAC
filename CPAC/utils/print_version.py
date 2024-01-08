"""Print C-PAC version for CI."""

import importlib.metadata


def print_cpac_version():
    """Print the version of CPAC"""
    print(importlib.metadata.version('CPAC'))
