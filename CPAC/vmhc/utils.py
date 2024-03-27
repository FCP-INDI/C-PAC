def get_img_nvols(in_files):
    """
    Calculates the number of volumes in the given nifti image.

    Parameters
    ----------
    in_files : string (nifti file)

    Returns
    -------
    out : int
        number of volumes of input nifti file

    """
    from nibabel import load

    img = load(in_files)
    hdr = img.header
    nvols = None
    if len(hdr.get_data_shape()) > 3:
        nvols = int(hdr.get_data_shape()[3])
    else:
        nvols = 1
    return nvols


def get_operand_expression(nvols):
    """
    Generates operand string.

    Parameters
    ----------
    nvols : int

    Returns
    -------
    expr : string

    """
    vol = int(nvols)
    return "a*sqrt(%d-3)" % vol
