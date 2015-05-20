def compute_corr(in_file, mask_file):

    """
    Computes the Network Correlation Matrix for ROIs in the mask

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)
        Must be 3D

    Returns
    -------

    out_file : correlation matrix

    """

    import numpy as np
    from matplotlib.mlab import csv2rec
    
    from CPAC.timeseries import gen_roi_timeseries  
    
    #from http://nipy.org/nitime/examples/resting_state_fmri.html    
    
    #Import the time-series objects:
    from nitime.timeseries import TimeSeries
    #Import the analysis objects:
    from nitime.analysis import CorrelationAnalyzer
    #Import utility functions:
    from nitime.utils import percent_change
    
    import nibabel as nib
    
    #in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    n1_img = nib.load(in_file)
    TR = n1_img.header['pixdim'][4]
    
    
    output_type = [True,False] #list of boolean for csv and npz file formats
    
    data = gen_roi_timeseries(in_file, mask_file, output_type)
    #from this files gen_roi_timeseries
    #once we have the time series:
    data_rec = csv2rec(data[2])
    
    #Extract information:
    roi_names = np.array(data_rec.dtype.names)
    n_samples = data_rec.shape[0]
    
    #Make an empty container for the data
    data = np.zeros((len(roi_names), n_samples))
    
    for n_idx, roi in enumerate(roi_names):
        data[n_idx] = data_rec[roi]
    
    #Normalize the data:
    data = percent_change(data)
    
    T = TimeSeries(data, sampling_interval=TR)
    T.metadata['roi'] = roi_names    
    #Initialize the correlation analyzer
    C = CorrelationAnalyzer(T)    
    
    #fig01 = drawmatrix_channels(C.corrcoef, roi_names, size=[10., 10.], color_anchor=0)
    
    
    ## IN CASE WE WOULD LIKE TO ADD A THRESHOLD FEATURE (set inside the function 
    # or from input)
    # C.corrcoef[C.corrcoef<0.7] = 0
    return  C.corrcoef


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
    from matplotlib.mlab import csv2rec
    
    from CPAC.timeseries import gen_roi_timeseries  
    
    #from http://nipy.org/nitime/examples/granger_fmri.html 
    
    #Import the time-series objects:
    import nitime.analysis as nta
    import nitime.timeseries as ts
    import nitime.utils as tsu
    
    import nibabel as nib
    
    #in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
    #mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')
    n1_img = nib.load(in_file)
    TR = n1_img.header['pixdim'][4]
    
    
#    res_fname = (in_file)
#    res_mask_fname = (mask_file)
#
#    res_img = nb.load(res_fname)
#    res_mask_img = nb.load(res_mask_fname)
#
#    res_data = res_img.get_data()
#    res_mask_data = res_mask_img.get_data()
#
#    print res_data.shape
#    (n_x, n_y, n_z, n_t) = res_data.shape
#
#    res_data = np.reshape(res_data, (n_x*n_y*n_z, n_t), order='F').T
#
#    Ranks_res_data = np.tile((np.zeros((1, (res_data.shape)[1]))), [(res_data.shape)[0], 1])
    
    
        
    output_type = [True,False] #list of boolean for csv and npz file formats
    data = gen_roi_timeseries(in_file, mask_file, output_type)
    #from this files gen_roi_timeseries
    #once we have the time series:
    data_rec = csv2rec(data[2])
    
    #Extract information:
    roi_names = np.array(data_rec.dtype.names)
    nseq = len(roi_names)
    n_samples = data_rec.shape[0]
    
    #Make an empty container for the data
    data = np.zeros((nseq, n_samples))
    
    
    
    
    for n_idx, roi in enumerate(roi_names):
        data[n_idx] = data_rec[roi]
    
    #Normalize the data in each of the ROIs to be in units 
    #of % change and initialize the TimeSeries object:
    pdata = tsu.percent_change(data)
    time_series = ts.TimeSeries(pdata, sampling_interval=TR)
    #We initialize the GrangerAnalyzer object, while specifying 
    #the order of the autoregressive model to be 1 
    #(predict the current behavior of the time-series 
    #based on one time-point back).
    G = nta.GrangerAnalyzer(time_series, order=1)  
    
    #We are only interested in the physiologically relevant frequency band
    
    g1 = np.mean(G.causality_xy, -1) #is this what we are looking for
    
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
#                
    
    
#    def pearson_correlation(x, y):
#    '''Computes pearson correlations on matrices
#    Parameters
#    ----------
#    x: np.ndarray or Dataset
#        PxM array
#    y: np.ndarray or Dataset or None (the default).
#        PxN array. If None, then y=x.
#    Returns
#    -------
#    c: np.ndarray
#        MxN array with c[i,j]=r(x[:,i],y[:,j])
#    Notes
#    -----
#    Unlike numpy. this function behaves like matlab's 'corr' function.
#    Its numerical precision is slightly lower than numpy's correlate function.
#    Unlike scipy's 'pearsonr' function it does not return p values.
#    TODO integrate with CorrCoef
#    '''
#
#
#    xd = x - np.mean(x, axis=0)
#    yd = y - np.mean(y, axis=0)
#
#    if xd.shape[0] != yd.shape[0]:
#        raise ValueError("Shape mismatch: %s != %s" % (xd.shape, yd.shape))
#
#    # normalize
#    n = 1. / (x.shape[0] - 1) # normalize
#
#    # standard deviation
#    xs = (n * np.sum(xd * xd, axis=0)) ** -.5
#    ys = (n * np.sum(yd * yd, axis=0)) ** -.5
#
#    return n * np.dot(xd.T, yd) * np.tensordot(xs, ys, 0)