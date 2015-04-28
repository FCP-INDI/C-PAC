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

    return  C.corrcoef


