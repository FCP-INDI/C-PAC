# -*- coding: utf-8 -*-
"""Utilities for ALFF."""

from pathlib import Path


def get_opt_string(mask: Path | str) -> str:
    """
    Return option string for 3dTstat.

    Parameters
    ----------
    mask : string
        Path to mask file

    Returns
    -------
    opt_str : string
        Command args

    """
    return f" -stdev -mask {mask}"
