import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


def set_gauss(fwhm):

    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string


def getImgNVols(in_files, stopIdx, startIdx):

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()

            if len(hdr.get_data_shape()) > 3:
                nvols = int(hdr.get_data_shape()[3])
            else:
                nvols = 1
            out.append(nvols)
        return out

    else:
        img = load(in_files)
        hdr = img.get_header()
        if len(hdr.get_data_shape()) > 3:
            nvols = int(hdr.get_data_shape()[3])
        else:
            nvols = 1
        return [nvols]


def getEXP(nvols):

    expr = []
    for vol in nvols:
        vol = int(vol)
        expr.append("'a*sqrt('%d'-3)'" % vol)

    return expr

