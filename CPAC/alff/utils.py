import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

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


def getImgTR(in_files, TRa):

    out = []
    from nibabel import load
    if(isinstance(in_files, list)):
        for in_file in in_files:
            img = load(in_file)
            hdr = img.get_header()
            tr = float(hdr.get_zooms()[3])
            if tr > 10:
                tr = float(float(tr) / 1000.0)
            if not (TRa == None):
                diff = None
                if TRa > tr:
                    diff = TRa - tr
                else:
                    diff = tr - TRa

                if (diff > 0.001):
                    print "Warning: specified TR  %f and TR in image header  %f do not match:" % (TRa, tr)
                out.append(TRa)
            else:
                out.append(tr)
        return out
    else:
        img = load(in_files)
        hdr = img.get_header()
        tr = float(hdr.get_zooms()[3])
        if tr > 10:
            tr = float(float(tr) / 1000.0)
        if not (TRa == None):
            diff = None
            if TRa > tr:
                diff = TRa - tr
            else:
                diff = tr - TRa

            if (diff > 0.001):
                print "Warning: specified TR  %f and TR in image header  %f do not match:" % (TRa, tr)
            return [TRa]
        else:
            return [tr]


def getN1(TR, nvols, HP):

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n1 = int("%1.0f" % (float(n_lp - 1.0)))

    return n1


def getN2(TR, nvols, LP, HP):

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n_hp = float(LP) * float(int(nvols)) * float(TR)
    n2 = int("%1.0f" % (float(n_hp - n_lp + 1.0)))

    return n2


def getOpString(mean, std_dev):

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))

    op_string = str1 + " -mas %s"

    return op_string


def set_op_str(n2):

    strs = []
    for n in n2:
        str = "-Tmean -mul %f" % (n)
        strs.append(str)
    return strs


def set_op1_str(nvols):

    strs = []
    for vol in nvols:
        str = '-Tmean -mul %d -div 2' % (int(vol))
        strs.append(str)

    return strs


def set_gauss(fwhm):

    op_string = ""

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string



def takemod(nvols):

    decisions = []
    for vol in nvols:
        mod = int(int(vol) % 2)

        if mod == 1:
            decisions.append([0])
        else:
            decisions.append([1])

    return decisions
