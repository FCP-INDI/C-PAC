#in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
#mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')

def gen_voxel_timeseries(in_file):

    """
    Extracts voxelwise timeseries and return a np array with them.
    
    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    Returns
    -------

    data_array =  voxel(x,y,z) * timepoints

    """
    import numpy as np    
    import nibabel as nb
    
    img_data = nb.load(in_file).get_data()
    #TR = datafile.header['pixdim'][4]
    n_samples = img_data.shape[3]

    #print img_data.shape
    (n_x, n_y, n_z, n_t) = img_data.shape
    voxel_data_array = np.reshape(img_data, (n_x*n_y*n_z, n_t), order='F')
    
    
return  voxel_data_array

def gen_roi_timeseries(in_file, mask_file):

    """
    Extracts ROI timeseries and return a np array with them.
    

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 
        
    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)
        Must be 3D    

    Returns
    -------

    data_array = ROI_number * timepoints

    """
    import numpy as np    
    import nibabel as nb
    
    img_data = nb.load(in_file).get_data()
    #TR = datafile.header['pixdim'][4]
    #n_samples = img_data.shape[3]

    mask_data = nb.load(mask_file).get_data()
    # Cast as rounded-up integer
    mask_data = np.int64(np.ceil(mask_data)) #ROI numbers, int64 is enough

    if mask_data.shape != img_data.shape[:3]:
        raise Exception('Invalid Shape Error.'\
                        'Please check the voxel dimensions.'\
                        'Data and roi should have'\
                        'same shape')

    nodes = np.unique(mask_data).tolist()
    nodes.remove(0) #quits the ROI number '0'
    roi_data_array = np.zeros((len(nodes),n_samples)) #,dtype=float change?

    # Be carefull with number of ROIs and np-arrays
    nodes.sort()
    for n in nodes:
        node_array = img_data[mask_data == n]
        avg = np.mean(node_array, axis=0)
        roi_data_array[n-1] = np.round(avg, 6)
    
return  roi_data_array

def compute_corr(timeseries):

    """
    Computes the Network Correlation Matrix for a timeseries * timepoints
    nparray

    Parameters
    ----------

    in_file : nparray: timeseries * timepoints


    Returns
    -------

    out_file : correlation matrix

    """

    import numpy as np    
 
    corr_matrix=np.corrcoef(timeseries)
    
    
    ## IN CASE WE WOULD LIKE TO ADD A THRESHOLD FEATURE (set inside the function 
    # or from input)
    # C.corrcoef[C.corrcoef<0.7] = 0
    return  corr_matrix


def compute_te(in_file, mask_file):

    """
    Computes the Pairwise Transfer Entropy Matrix for ROIs in the mask
    For now, we are computing GC, since GC and TE are equivalent for Gaussian
    signals.
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.103.238701

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)
        Must be 3D

    Returns
    -------

    out_file : Transfer Entropy matrix

    """

    import numpy as np
    import nibabel as nb
    
    #NEED TO  WORK ON THIS
    
    return  g1
    

def entropy(*X):
    
    """
    Same function for entropy as for joint entropy
    """

    import numpy as np
    import itertools 
    
    n_insctances = len(X[0])
    H = 0
    for classes in itertools.product(*[set(x) for x in X]):
        v = np.array([True] * n_insctances)
        for predictions, c in zip(X, classes):
            v = np.logical_and(v, predictions == c)
        p = np.mean(v)
        H += -p * np.log2(p) if p > 0 else 0
    return H  
    
def mutual_information(X,Y):

    Hx = entropy(X)
    Hy = entropy(Y)
    Hxy = entropy(X,Y)
    MI = Hx + Hy - Hxy
    
    return MI

def cond_entropy(X,Y):
    
    """
    Conditional entropy H(X|Y) = H(Y,X) - H(Y). X conditioned on Y
    """

    Hy = entropy(Y)
    Hyx = entropy(Y,X)
    CE = Hyx - Hy    
    
    return CE
    
#def entropy(*X):
#    n_insctances = len(X[0])
#    H = 0
#    for classes in itertools.product(*[set(x) for x in X]):
#        v = reduce(np.logical_and, (predictions, c for predictions, c in zip(X, classes)))
#        p = np.mean(v)
#        H += -p * np.log2(p) if p > 0 else 0
#    return H   
#   
#def entropy(*X):
#    return = np.sum(-p * np.log2(p) if p > 0 else 0 for p in \
#    (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product (*[set(x) for x in X])))

    
