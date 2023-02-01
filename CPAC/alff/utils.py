# -*- coding: utf-8 -*-

def get_opt_string(mask):
    """
    Method to return option string for 3dTstat
    
    Parameters
    ----------
    mask : string
        Path to mask file
    
    Returns
    -------
    opt_str : string
        Command args
    
    """

    opt_str = " -stdev -mask %s" % mask

    return opt_str
