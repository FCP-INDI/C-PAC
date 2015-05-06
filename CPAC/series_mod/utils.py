def compute_corr(in_file, mask_file, TR):

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
    
    from cpac.timeseries import gen_roi_timeseries  
    
    #from http://nipy.org/nitime/examples/resting_state_fmri.html    
    
    #Import the time-series objects:
    from nitime.timeseries import TimeSeries
    #Import the analysis objects:
    from nitime.analysis import CorrelationAnalyzer
    #Import utility functions:
    from nitime.utils import percent_change
    
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


def compute_te(in_file, mask_file, TR):

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
    
    from cpac.timeseries import gen_roi_timeseries  
    
    #from http://nipy.org/nitime/examples/granger_fmri.html 
    
    #Import the time-series objects:
    import nitime.analysis as nta
    import nitime.timeseries as ts
    import nitime.utils as tsu
    
    #freq band of interest (these could be an input paramenter)
    f_ub = 0.15
    f_lb = 0.02
    
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
    freq_idx_G = np.where((G.frequencies > f_lb) * (G.frequencies < f_ub))[0]
    g1 = np.mean(G.causality_xy[:, :, freq_idx_G], -1) #is this what we are looking for? ?¿ ask
    
    return  g1