import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def set_gauss(fwhm):

    """
    Compute the sigma value, given Full Width Half Max. 
    Further it builds an operand string and returns it

    Parameters
    ----------

    fwhm : float

    Returns
    -------

    op_string : string

    """

    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string


def get_img_nvols(in_files):

    """
    Calculates the number of volumes in the given nifti image

    Parameters
    ----------

    in_files : string (nifti file)

    Returns
    -------

    out : int
        number of volumes of input nifti file

    """

    out = None
    from nibabel import load
    img = load(in_files)
    hdr = img.get_header()
    nvols = None
    if len(hdr.get_data_shape()) > 3:
        nvols = int(hdr.get_data_shape()[3])
    else:
        nvols = 1
    out = nvols

    return out


def get_operand_expression(nvols):

    """
    Generates operand string

    Parameters
    ----------

    nvols : int

    Returns
    -------

    expr : string

    """

    expr = None
    vol = int(nvols)
    expr = ('a*sqrt(%d-3)' % vol)

    return expr

