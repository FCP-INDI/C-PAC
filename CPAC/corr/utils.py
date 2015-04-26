

def compute_corr(in_file, mask_file, TR):

    """
    Computes the Network Correlation Matrix for ROIs in the mask

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)

    Returns
    -------

    out_file : nifti file
        ReHo map of the input EPI image

    """

    import nibabel as nb
    import numpy as np
    import os
    import sys
    
    from cpac.timeseries import gen_roi_timeseries  
    
    #from http://nipy.org/nitime/examples/resting_state_fmri.html    
    
    import nitime
    #Import the time-series objects:
    from nitime.timeseries import TimeSeries
    #Import the analysis objects:
    from nitime.analysis import CorrelationAnalyzer
    #Import utility functions:
    from nitime.utils import percent_change
    
    res_fname = (in_file)
    res_mask_fname = (mask_file) 
    out_file = None
    
    #from this files gen_roi_timeseries
    #once we have the time series:
    
    #data_path = os.path.join(nitime.__path__[0], 'data')
    #data_rec = csv2rec(os.path.join(data_path, 'fmri_timeseries.csv'))
    
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
    
    
    
    


    #create corr_file with C.Corrcoeff
    #ROI-names needed?

    out_file = corr_file

    return out_file


