import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def bandpass_voxels(realigned_file, bandpass_freqs, sample_period = None):
    """
    Performs ideal bandpass filtering on each voxel time-series.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    bandpass_freqs : tuple
        Tuple containing the bandpass frequencies.
    sample_period : float, optional
        Length of sampling period in seconds.  If not specified,
        this value is read from the nifti file provided.
        
    Returns
    -------
    bandpassed_file : string
        Path of filtered output (nifti file).
    
    """
    
    import os
    import nibabel as nb
    import numpy as np

    def ideal_bandpass(data, sample_period, bandpass_freqs):
        #Derived from YAN Chao-Gan 120504 based on REST.
        from scipy.fftpack import fft, ifft
        
    #    sample_period = T
    #    LowCutoff = 10.
    #    HighCutoff = 15.
    #    data = x
        
        def nextpow2(n):
            x = np.log2(n)
            return 2**np.ceil(x)
        
        sample_freq = 1./sample_period
        sample_length = data.shape[0]
        
        data_p = np.zeros(nextpow2(sample_length))
        data_p[:sample_length] = data
        
        LowCutoff, HighCutoff = bandpass_freqs
        
        if(LowCutoff is None): #No lower cutoff (low-pass filter)
            low_cutoff_i = 0
        elif(LowCutoff > sample_freq/2.): #Cutoff beyond fs/2 (all-stop filter)
            low_cutoff_i = int(data_p.shape[0]/2)
        else:
            low_cutoff_i = np.ceil(LowCutoff*data_p.shape[0]*sample_period).astype('int')
        
        if(HighCutoff > sample_freq/2. or HighCutoff is None): #Cutoff beyond fs/2 or unspecified (become a highpass filter)
            high_cutoff_i = int(data_p.shape[0]/2)
        else:
            high_cutoff_i = np.fix(HighCutoff*data_p.shape[0]*sample_period).astype('int')
        
        freq_mask = np.zeros_like(data_p, dtype='bool')
        freq_mask[low_cutoff_i:high_cutoff_i+1] = True
        freq_mask[data_p.shape[0]-high_cutoff_i:data_p.shape[0]+1-low_cutoff_i] = True
        
        
        f_data = fft(data_p)
        f_data[freq_mask != True] = 0.
        data_bp = np.real_if_close(ifft(f_data)[:sample_length])
        
        return data_bp

    nii = nb.load(realigned_file)
    data = nii.get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0
    Y = data[mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    
    if not sample_period:
        hdr = nii.get_header()
        sample_period = float(hdr.get_zooms()[3])
        # Sketchy check to convert TRs in millisecond units
        if sample_period > 20.0:
            sample_period /= 1000.0

    print 'Frequency filtering using sample period:', sample_period, 'sec'

    Y_bp = np.zeros_like(Y)
    for j in range(Y.shape[1]):
        Y_bp[:,j] = ideal_bandpass(Yc[:,j], sample_period, bandpass_freqs)
        
    data[mask] = Y_bp.T
    img = nb.Nifti1Image(data, header=nii.get_header(), affine=nii.get_affine())
    bandpassed_file = os.path.join(os.getcwd(), 'bandpassed_demeaned_filtered.nii.gz')
    img.to_filename(bandpassed_file)
    
    return bandpassed_file

def calc_residuals(subject,
                   selector,
                   wm_mask_file = None,
                   csf_mask_file = None,
                   gm_mask_file = None,
                   motion_file = None,
                   compcor_ncomponents = 0):
    """
    Calculates residuals of nuisance regressors for every voxel for a subject.
    
    Parameters
    ----------
    subject : string
        Path of a subject's realigned nifti file.
    selector : dictionary
        Dictionary of selected regressors.  Keys are  represented as a string of the regressor name and keys 
        are True/False.  See notes for an example.
    wm_mask_file : string, optional
        Path to subject's white matter mask (in the same space as the subject's functional file)
    csf_mask_file : string, optional
        Path to subject's cerebral spinal fluid mask (in the same space as the subject's functional file)
    gm_mask_file : string, optional
        Path to subject's grey matter mask (in the same space as the subject's functional file)
    compcor_ncomponents : integer, optional
        The first `n` principal of CompCor components to use as regressors.  Default is 0.
        
    Returns
    -------
    residual_file : string
        Path of residual file in nifti format
    regressors_file : string
        Path of matlab file of regressors used
        
    Notes
    -----
    
    Example of selector parameter:
    
    >>> selector = {'compcor' : True,
    >>> 'wm' : True,
    >>> 'csf' : True,
    >>> 'gm' : True,
    >>> 'global' : True,
    >>> 'pc1' : True,
    >>> 'motion' : True,
    >>> 'linear' : True,
    >>> 'quadratic' : True}
    
    
    """
    import numpy as np
    import nibabel as nb
    import os
    import numpy as np
    import scipy
    from CPAC.nuisance import calc_compcor_components
    
    nii = nb.load(subject)
    data = nii.get_data().astype(np.float64)
    global_mask = (data != 0).sum(-1) != 0
    
    regressor_map = {'constant' : np.ones((data.shape[3],1))}
    if(selector['compcor']):
        print 'compcor_ncomponents', compcor_ncomponents
        wm_mask = nb.load(wm_mask_file).get_data().astype('bool')
        csf_mask = nb.load(csf_mask_file).get_data().astype('bool')
        regressor_map['compcor'] = calc_compcor_components(data, compcor_ncomponents, wm_mask, csf_mask)
        
    if(selector['wm']):
        wm_mask = nb.load(wm_mask_file).get_data().astype('bool')
        regressor_map['wm'] = data[wm_mask].mean(0)
        
    if(selector['csf']):
        csf_mask = nb.load(csf_mask_file).get_data().astype('bool')
        regressor_map['csf'] = data[csf_mask].mean(0)
        
    if(selector['gm']):
        gm_mask = nb.load(gm_mask_file).get_data().astype('bool')
        regressor_map['gm'] = data[gm_mask].mean(0)
        
    if(selector['global']):
        regressor_map['global'] = data[global_mask].mean(0)
        
    if(selector['pc1']):
        bdata = data[global_mask].T
        bdatac = bdata - np.tile(bdata.mean(0), (bdata.shape[0], 1))
        U, S, Vh = np.linalg.svd(bdatac, full_matrices=False)
        regressor_map['pc1'] = U[:,0]
        
    if(selector['motion']):
        motion = np.genfromtxt(motion_file)
        regressor_map['motion'] = motion
        
    if(selector['linear']):
        regressor_map['linear'] = np.arange(0, data.shape[3])
    
    if(selector['quadratic']):
        regressor_map['quadratic'] = np.arange(0, data.shape[3])**2
    
    print 'Regressors include', regressor_map.keys()
    
    X = np.zeros((data.shape[3], 1))
    for rname, rval in regressor_map.items():
        X = np.hstack((X, rval.reshape(rval.shape[0],-1)))
    X = X[:,1:]
    
    print 'Regressors dim', X.shape, 'starting regression'
    
    Y = data[global_mask].T
    B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    Y_res = Y - X.dot(B)
    
    data[global_mask] = Y_res.T
    
    print 'Writing residual and regressors'
    img = nb.Nifti1Image(data, header=nii.get_header(), affine=nii.get_affine())
    residual_file = os.path.join(os.getcwd(), 'residual.nii.gz')
    img.to_filename(residual_file)
    
    regressors_file = os.path.join(os.getcwd(), 'nuisance_regressors.mat')
    
    scipy.io.savemat(regressors_file, regressor_map, oned_as='column')
    
    return residual_file, regressors_file

def create_nuisance(name='nuisance'):
    """
    Workflow for the removal of various signals considered to be noise in resting state
    fMRI data.  The residual signals for linear regression denoising is performed in a single
    model.  Therefore the residual time-series will be orthogonal to all signals.
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
    
    Returns
    -------
    nuisance : nipype.pipeline.engine.Workflow
        Nuisance workflow.
        
    Notes
    -----
    
    Workflow Inputs::

        inputspec.subject : string (nifti file)
            Path of the subject's realigned nifti file.
        inputspec.wm_mask : string (nifti file)
            Corresponding white matter mask.
        inputspec.csf_mask : string (nifti file)
            Corresponding cerebral spinal fluid mask.
        inputspec.gm_mask : string (nifti file)
            Corresponding grey matter mask.
        inputspec.motion_components : string (text file)
            Corresponding rigid-body motion parameters.  Matrix in the file should be of shape 
            (`T`, `R`), `T` timepoints and `R` motion parameters.
        inputspec.selector : dictionary
        inputspec.compcor_ncomponents : integer
        
    Workflow Outputs::

        outputspec.subject : string (nifti file)
        outputspec.regressors : string (mat file)
    
    Nuisance Procedure:
    
    1. Compute nuisance regressors based on input selections.
    2. Calculate residuals with respect to these nuisance regressors in a
       single model for every voxel.
    
    Workflow Graph:
    
    .. image:: ../images/nuisance.dot.png
        :width: 500
    
    Detailed Workflow Graph:
    
    .. image:: ../images/nuisance_detailed.dot.png
        :width: 500    
    
    """
    nuisance = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                       'wm_mask',
                                                       'csf_mask',
                                                       'gm_mask',
                                                       'motion_components',
                                                       'selector',
                                                       'compcor_ncomponents']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                        'regressors']),
                         name='outputspec')
    
    calc_r = pe.Node(util.Function(input_names=['subject',
                                                'selector',
                                                'wm_mask_file',
                                                'csf_mask_file',
                                                'gm_mask_file',
                                                'motion_file',
                                                'compcor_ncomponents'],
                                   output_names=['residual_file',
                                                'regressors_file'],
                                   function=calc_residuals),
                     name='residuals')
    
    nuisance.connect(inputspec, 'subject',
                     calc_r, 'subject')
    nuisance.connect(inputspec, 'wm_mask',
                     calc_r, 'wm_mask_file')
    nuisance.connect(inputspec, 'csf_mask',
                     calc_r, 'csf_mask_file')
    nuisance.connect(inputspec, 'gm_mask',
                     calc_r, 'gm_mask_file')
    nuisance.connect(inputspec, 'motion_components',
                     calc_r, 'motion_file')
    nuisance.connect(inputspec, 'selector',
                     calc_r, 'selector')
    nuisance.connect(inputspec, 'compcor_ncomponents',
                     calc_r, 'compcor_ncomponents')
    nuisance.connect(calc_r, 'residual_file',
                     outputspec, 'subject')
    nuisance.connect(calc_r, 'regressors_file',
                     outputspec, 'regressors')
    
    return nuisance