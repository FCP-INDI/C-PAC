import numpy as np
import nibabel as nb

"""
These set of functions calculate the memory required for different CWAS steps.
"""

b2gb = lambda x: float(x)/(1024**3)
gb2b = lambda x: float(x)*(1024**3)
b2mb = lambda x: float(x)/(1024**2)
mb2b = lambda x: float(x)*(1024**2)
b2kb = lambda x: float(x)/(1024**1)
kb2b = lambda x: float(x)*(1024**1)
b2b  = lambda x: float(x)*(1024**0)

def array_memory(shape, dtype, otype='gb'):
    """docstring for predict_memory_usage"""
    nitems      = np.prod(shape)
    item_size   = np.dtype(dtype).itemsize
    func_convert = {
        "gb": b2gb, 
        "mb": b2mb, 
        "kb": b2kb, 
        'b':  b2b
    }
    ret = func_convert[otype](nitems * item_size)
    return ret
    
    
###
# Memory for Functional (4D) Data
###

def func_memory(nvoxs, ntpts, dtype):
    return array_memory(nvoxs*ntpts, dtype)

def subjs_func_ntpts_from_list(func_file_list):
    func_ntpts = []
    for func_file in func_file_list:
        img = nb.load(func_file)
        if len(img.shape) != 4:
            raise ValueError("Shape of functional is not of length 4")
        func_ntpts.append(img.shape[-1])
    return func_ntpts

def subjs_func_memory(nvoxs, func_ntpts, dtype):
    mem = [ func_memory(nvoxs, ntpts, dtype) for ntpts in func_ntpts ]
    tot = np.array(mem).sum()
    return tot


###
# Memory for Computing Connectivity Maps
###

def conn_map_memory(nvoxs, dtype):
    """Calculate the memory for one whole-brain connectivity map"""
    return array_memory(nvoxs, dtype)

def nvoxs_with_conn_map(memlimit, nvoxs, nsubjs, dtype):
    """Determines the number of voxels to compute connectivity maps for at once"""
    x_voxels = float(memlimit)/(nsubjs * conn_map_memory(nvoxs, dtype))
    x_voxels = np.floor(x_voxels)
    x_voxels = min(x_voxels, nvoxs)
    return x_voxels


###
# Memory for Distance Matrices
###

def dmat_memory(nvoxs, nobs, dtype):
    return array_memory(nvoxs*(nobs**2), dtype)


###
# Memory for MDMR
###

def perm_mats_memory(nperms, nobs, dtype):
    """
    Calculate memory requirement for permutation indices 
    and permuted matrices for H2s and IHs.
    """
    perms   = array_memory((nobs**2)*nperms, dtype)
    H2perms = array_memory((nobs**2)*nperms, dtype)
    IHperms = array_memory((nobs**2)*nperms, dtype)
    return perms + H2perms + IHperms

def fperms_memory(nperms, nvoxs, dtype):
    Fperms = array_memory(nperms*nvoxs, dtype)
    return Fperms

def nvoxs_per_mdmr(memlimit, nperms, nobs, nvoxs, dtype):
    """Determine the number of voxels to compute with MDMR at once"""
    # divide by number of gower centered matrices and number of Fperms
    x_voxels = float(memlimit)/(array_memory(nobs**2, dtype) + array_memory(nperms, dtype))
    x_voxels = np.floor(x_voxels)
    x_voxels = min(x_voxels, nvoxs)
    return x_voxels

