import os
import sys
import re
import commands
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util


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

    from nibabel import load

    nvols = None
    img = load(in_files)
    hdr = img.get_header()
    if len(hdr.get_data_shape()) > 3:
        nvols = int(hdr.get_data_shape()[3])
    else:
        nvols = 1
    return nvols



def get_img_tr(in_files, TRa):


    """
    Computes the Temporal Resolution parameter from the image header.
    If TRa is not None then it compares the computed TR and the supplied
    if the difference is greater than 0.001 then it returns supplied TRa as the
    TR of the image and flashes a warning message else returns TRa

    Parameters
    ----------

    in_files : string (nifti file)

    TRa : float
        user supplied Temporal Resolution

    Returns
    -------

    out : float
        TR for input file

    """

    from nibabel import load
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
        return TRa
    else:
        return tr



def get_N1(nvols, TR, HP):

    """
    Get the low frequency point

    Parameters
    ----------

    TR : float
        Temporal Resolution

    nvols : int
        Number of volumes

    HP : float
        HighPass Low Cutoff Frequency

    Returns
    -------

    n1 : float

    Low Frequency Point

    """

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n1 = int("%1.0f" % (float(n_lp - 1.0)))

    return n1



def get_N2(nvols, TR, LP, HP):

    """
    Get the high frequency point

    Parameters
    ----------

    TR : float
        Temporal Resolution

    nvols : int
        Number of volumes

    LP : float
        LowPass High Cutoff Frequency
    
    HP : float
        HighPass Low Cutoff Frequency

    Returns
    -------

    n2 : float

    High Frequency Point

    """

    n_lp = float(HP) * float(int(nvols)) * float(TR)
    n_hp = float(LP) * float(int(nvols)) * float(TR)
    n2 = int("%1.0f" % (float(n_hp - n_lp + 1.0)))

    return n2



def get_operand_string(mean, std_dev):

    """
    Generate the Operand String to be used in workflow nodes to supply 
    mean and std deviation to alff workflow nodes

    Parameters
    ----------

    mean: string
        mean value in string format
    
    std_dev : string
        std deviation value in string format


    Returns
    -------

    op_string : string


    """

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))

    op_string = str1 + " -mas %s"

    return op_string




def set_op_str(n2):

    """
    Build the operand string used by the workflow nodes

    Parameters
    ----------

    n2 : float

    Returns
    -------

     strs : string

        operand string

    """

    strs = "-Tmean -mul %f" % (n2)
    return strs


def set_op1_str(nvols):

    """
    Build the operand string used by the workflow nodes

    Parameters
    ----------

    nvols : int

    Returns
    -------

     strs : string

        operand string

    """

    strs = '-Tmean -mul %d -div 2' % (int(nvols))

    return strs


def get_opt_string(mask):
    """
    Method to return option string for 3dTstat
    
    Parameters
    ----------
    mask : string (file)
    
    Returns
    -------
    opt_str : string
    
    """
    opt_str = " -stdev -mask %s" %mask
    return opt_str



def takemod(nvols):

    """
    Determine if the input is odd or even values and 
    return a of 0 and 1 depending on the truth value

    Parameters
    ----------

    nvols : int

    Returns
    -------

    decisions : int
    """

    mod = int(int(nvols) % 2)

    if mod == 1:
        return 0
    else:
        return 1

