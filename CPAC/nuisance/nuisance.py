import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util

def bandpass_voxels(realigned_file, sample_period, bandpass_freqs):
    """
    Performs ideal bandpass filtering on each voxel time-series.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    sample_period : float
        Length of sampling period in seconds.
    bandpass_freqs : tuple
        Tuple containing the bandpass frequencies.
    
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
    
    Y_bp = np.zeros_like(Y)
    for j in range(Y.shape[1]):
        Y_bp[:,j] = ideal_bandpass(Yc[:,j], sample_period, bandpass_freqs)
        
    data[mask] = Y_bp.T
    img = nb.Nifti1Image(data, header=nii.get_header(), affine=nii.get_affine())
    bandpassed_file = os.path.join(os.getcwd(), 'bandpassed_demeaned_filtered.nii.gz')
    img.to_filename(bandpassed_file)
    
    return bandpassed_file


def linear_detrend_voxels(realigned_file):
    """
    Removes linear trend for each voxel time-series.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
        
    Returns
    -------
    detrended_file : string
        Path of detrended output (nifti file).
    
    """
    
    import os
    import nibabel as nb
    import numpy as np
    from scipy.signal import detrend
    
    nii = nb.load(realigned_file)
    data = nii.get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0
    Y = detrend(data[mask], axis=1, type='linear').T
    data[mask] = Y.T
    img = nb.Nifti1Image(data, header=nii.get_header(), affine=nii.get_affine())
    detrended_file = os.path.join(os.getcwd(), 'linear_detrended.nii.gz')
    img.to_filename(detrended_file)
    
    return detrended_file

def extract_compcor_components(nc,
                               realigned_file,
                               wm_mask,
                               csf_mask):
    """
    Extracts the largest principal components found in white matter (WM) and cerebral spinal fluid
    (CSF).  CompCor algorithm based on [1]_.  Note that WM and CSF masks are not eroded.
    
    Parameters
    ----------
    nc : integer
        Number of the largest principal components (in decsending order) to extract.
    realigned_file : string
        Path of a realigned nifti file.
    wm_mask : string
        Path of a white matter mask (nifti file).
    csf_mask : string
        Path of a cerebral spinal fluid mask (nifti file).
        
    Returns
    -------
    components_file : string
        Path of the text file containing the extracted principal components.  Each component is
        stored as a column in the matrix.
    
    References
    ----------
    .. [1] Y. Behzadi, K. Restom, J. Liau, and T. T. Liu, A component based noise correction method (CompCor) for BOLD and perfusion based fMRI., NeuroImage, vol. 37, no. 1, pp. 90-101, Aug. 2007.
     
    """
    
    import os
    import nibabel as nb
    import scipy as sp
    import numpy as np
    from scipy.signal import detrend

    data = nb.load(realigned_file).get_data().astype('float64')
    wm_mask = nb.load(wm_mask).get_data().astype('float64')
    csf_mask = nb.load(csf_mask).get_data().astype('float64')
    print 'Data and masks loaded.'
    wmcsf_mask = (csf_mask + wm_mask).astype('bool')

    print 'Detrending and centering data'
    Y = detrend(data[wmcsf_mask], axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))

    print 'Calculating SVD decomposition of Y*Y\''
    U, S, Vh = np.linalg.svd(np.dot(Yc, Yc.T))

    components_file = os.path.join(os.getcwd(), 'noise_components.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, U[:, :nc])
    return components_file

def extract_global_component(realigned_file):
    """
    Extracts the global signal, defined as the average time-series of all non-zero voxels.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    
    Returns
    -------
    components_file : string
        Path of the text file containing global signal component
    
    """
    
    import os
    import nibabel as nb
    import numpy as np
    from utils import mean_roi_signal

    data = nb.load(realigned_file).get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0  # Global Mask
    print 'Data loaded.'
#    Y = data[mask].T
#    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    glb_comp = mean_roi_signal(data, mask)

    components_file = os.path.join(os.getcwd(), 'global_component.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, glb_comp)

    return components_file

def extract_frommask_component(realigned_file, mask_file):
    """
    Extracts the average time-series of voxels inside a given mask.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    mask_file : string
        Path of a mask file
    
    Returns
    -------
    components_file : string
        Path of the text file containing extracted component
    
    """
    import os
    import nibabel as nb
    import numpy as np
    from utils import mean_roi_signal
    from nipype import logging
    iflogger = logging.getLogger('interface')

    data = nb.load(realigned_file).get_data().astype('float64')
    mask = nb.load(mask_file).get_data().astype('float64')
    iflogger.info('Data and mask loaded.')
    mask_comp = mean_roi_signal(data, mask.astype('bool'))

    components_file = os.path.join(os.getcwd(), 'mask_mean_component.txt')
    iflogger.info('Saving components file:' + components_file)
    np.savetxt(components_file, mask_comp)

    return components_file

def extract_firstprinc_component(realigned_file):
    """
    Extracts the first principal component of all non-zero voxels.  Note the data is centered in 
    voxel space only.  Without centering in time space, the first principal component will likely 
    be highly correlated with the global signal
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    
    Returns
    -------
    components_file : string
        Path of the text file containing the first principal component
    """
    import os
    import nibabel as nb
    import numpy as np

    data = nb.load(realigned_file).get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0  # Global Mask
    print 'Data loaded.'
    Y = data[mask].T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))

    print 'Calculating SVD decomposition of Y'
    U, S, Vh = np.linalg.svd(Yc, full_matrices=False)

    components_file = os.path.join(os.getcwd(), 'firstprinc_component.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, U[:, 0])

    return components_file

def extract_linear_component(realigned_file):
    """
    Creates a linear component of the same length as the input time-series.  This linear
    component starts at `0` and ends at `(timepoints - 1)`
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    
    Returns
    -------
    components_file : string
        Path of the text file containing the linear component
    """
    import os
    import nibabel as nb
    import numpy as np

    data = nb.load(realigned_file).get_data().astype('float64')
    lt = np.arange(0, data.shape[-1])
    
    components_file = os.path.join(os.getcwd(), 'linear_component.txt')
    print 'Saving components file:', components_file
    np.savetxt(components_file, lt)

    return components_file

def create_filter_matrix(global_component,
                         compcor_components,
                         wm_component,
                         csf_component,
                         gm_component,
                         firstprinc_component,
                         linear_component,
                         motion_components,
                         selector):
    """
    Combines the nuisance factors into one components file.  Nuisance selection structure 
    based on https://github.com/satra/BrainImagingPipelines/tree/master/fmri
    """
    import numpy as np
    import os

    def try_import(fname):
        try:
            a = np.genfromtxt(fname)
            return a
        except:
            return np.array([])

    options = np.array([global_component,
                         compcor_components,
                         wm_component,
                         csf_component,
                         gm_component,
                         firstprinc_component,
                         linear_component,
                         motion_components])
    fieldnames = np.array(['global',
                           'compcor',
                           'wm',
                           'csf',
                           'gm',
                           'firstprinc',
                           'lt',
                           'motion'])

    selector = np.array(selector)  # Use selector as an index mask
    #Grab component filenames of according to selector
    filenames = fieldnames[selector]
    filter_file = os.path.abspath('filter_%s.txt' % '_'.join(filenames))

    z = None
    for i, opt in enumerate(options[selector]):
        a = try_import(opt)
        if len(a.shape) == 1:
            a = np.array([a]).T
        print a.shape
        if i == 0:
            z = a
        else:
            z = np.hstack((z, a))

    print 'Writing filter design matrix of size', z.shape,\
         'to file', filter_file
    np.savetxt(filter_file, z)
    return filter_file

def median_angle_correct(target_angle_deg, realigned_file):
    """
    Performs median angle correction on fMRI data.  Median angle correction algorithm
    based on [1]_.
    
    Parameters
    ----------
    target_angle_deg : float
        Target median angle to adjust the time-series data.
    realigned_file : string
        Path of a realigned nifti file.
    
    Returns
    -------
    corrected_file : string
        Path of corrected file (nifti file).
    angles_file : string
        Path of numpy file (.npy file) containing the angles (in radians) of all voxels with 
        the 5 largest principal components.
    
    References
    ----------
    .. [1] H. He and T. T. Liu, "A geometric view of global signal confounds in resting-state functional MRI," NeuroImage, Sep. 2011.
    
    """
    import numpy as np
    import nibabel as nb
    import os
    from scipy.stats.stats import pearsonr

    def shiftCols(pc, A, dtheta):
        pcxA = np.dot(pc, A)
        x = A - np.dot(pc[:, np.newaxis], pcxA[np.newaxis, :])

        theta = np.arccos(np.dot(pc.T, A))
        theta_new = theta + dtheta

        x /= np.tile(np.sqrt((x * x).sum(0)), (x.shape[0], 1))
        v_new = np.dot(pc[:, np.newaxis],\
             np.cos(theta_new)[np.newaxis, :]) + (np.sin(theta_new) * x)

        return v_new

    def writeToFile(data, nii, fname):
        img_whole_y = nb.Nifti1Image(data,\
            header=nii.get_header(), affine=nii.get_affine())
        img_whole_y.to_filename(fname)

    nii = nb.load(os.path.join(realigned_file))
    data = nii.get_data().astype(np.float64)
    print realigned_file, "subject data dimensions:", data.shape

    mask = (data != 0).sum(-1) != 0

    Y = data[mask].T

    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yn = Yc / np.tile(np.sqrt((Yc * Yc).sum(0)), (Yc.shape[0], 1))
    U, S, Vh = np.linalg.svd(Yn, full_matrices=False)

    G = Yc.mean(1)
    corr_gu = pearsonr(G, U[:, 0])
    PC1 = U[:, 0] if corr_gu[0] >= 0 else -U[:, 0]
    print 'Correlation of Global and U:', corr_gu

    median_angle = np.median(np.arccos(np.dot(PC1.T, Yn)))
    print 'Median Angle:', (180.0 / np.pi) * median_angle,\
        'Target Angle:', target_angle_deg
    angle_shift = (np.pi / 180) * target_angle_deg - median_angle
    if(angle_shift > 0):
        print 'Shifting all vectors by',\
            (180.0 / np.pi) * angle_shift, 'degrees.'
        Ynf = shiftCols(PC1, Yn, angle_shift)
    else:
        print 'Median Angle >= Target Angle, skipping correction'
        Ynf = Yn

    corrected_file = os.path.join(os.getcwd(), 'median_angle_corrected.nii.gz')
    angles_file = os.path.join(os.getcwd(), 'angles_U5_Yn.npy')

    print 'Writing U[:,0:5] angles to file...', angles_file
    angles_U5_Yn = np.arccos(np.dot(U[:, 0:5].T, Yn))
    np.save(angles_file, angles_U5_Yn)

    print 'Writing correction to file...', corrected_file
    data = np.zeros_like(data)
    data[mask] = Ynf.T
    writeToFile(data, nii, corrected_file)

    return corrected_file, angles_file

def extract_residuals(realigned_file, regressors_file):
    """
    Extracts the residual time-series of all voxels with respect to the given regressors.  The
    output file will be in an orthogonal subspace to the regressors.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    regressors_file : string
        Path of text file of regressors.  Matrix in the file should be of shape (`T`, `R`), `T` 
        timepoints and `R` regressors.
        
    Returns
    -------
    residual_file : string
        Path of residual output file (nifti file).
    """
    import os
    import nibabel as nb
    import numpy as np

    nii = nb.load(realigned_file)
    data = nii.get_data().astype('float64')
    mask = (data != 0).sum(-1) != 0
    Y = data[mask].T

    X = np.genfromtxt(regressors_file)
    if len(X.shape) <= 1:
        X = X[:,np.newaxis]

    X = np.hstack((X,np.ones((X.shape[0],1)))) #Add constant regressor to model
    B = np.dot(np.dot(np.linalg.inv(np.dot(X.T,X)), X.T), Y)
    XB = np.dot(X,B)
    Y_res = Y - XB
    
    data[mask] = Y_res.T
    
    img = nb.Nifti1Image(data, header=nii.get_header(), affine=nii.get_affine())
    residual_file = os.path.join(os.getcwd(), 'residual.nii.gz')
    img.to_filename(residual_file)
    
    return residual_file

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

        inputspec.realigned_file : string (nifti file)
            Path of a realigned nifti file.
        inputspec.wm_mask : string (nifti file)
            Corresponding white matter mask.
        inputspec.csf_mask : string (nifti file)
            Corresponding cerebral spinal fluid mask.
        inputspec.gm_mask : string (nifti file)
            Corresponding grey matter mask.
        inputspec.motion_components : string (text file)
            Corresponding rigid-body motion parameters.  Matrix in the file should be of shape 
            (`T`, `R`), `T` timepoints and `R` motion parameters.
    
    
    Workflow Outputs::
    
        outputspec.residual_file : string (nifti file)
        outputspec.median_angle_corrected_file : string (nifti file)
        outputspec.angles : string (numpy file)
    
    Workflow Graph:
    
    .. image:: ../images/nuisance.dot.png
        :width: 500
        
    Detailed Workflow Graph:
    
    .. image:: ../images/nuisance_detailed.dot.png
        :width: 500
    
    """
    inputspec = pe.Node(util.IdentityInterface(fields=['realigned_file',
                                                       'wm_mask',
                                                       'csf_mask',
                                                       'gm_mask',
                                                       'motion_components'
                                                      ]),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=[
                                            'residual_file',
                                            'median_angle_corrected_file',
                                            'angles']),
                         name='outputspec')

    nuisance = pe.Workflow(name=name)

    inputnode_selector = pe.Node(util.IdentityInterface(fields=['selector']),
                             name='selector_input')

    inputnode_nc = pe.Node(util.IdentityInterface(fields=['nc']),
                             name='nc_input')

    inputnode_target_angle_deg = pe.Node(util.IdentityInterface(\
                             fields=['target_angle_deg']),
                             name='target_angle_deg_input')

    linear_detrend = pe.MapNode(util.Function(input_names=['realigned_file'],
                                              output_names=['detrended_file'],
                                              function=linear_detrend_voxels),
                                              name='linear_detrend',
                                              iterfield=['realigned_file'])

    median_angle = pe.MapNode(util.Function(input_names=['target_angle_deg',
                                                         'realigned_file'],
                                            output_names=['corrected_file',
                                                          'angles'],
                                            function=median_angle_correct),
                                            name='median_angle',
                                            iterfield=['realigned_file'])

    compcor = pe.MapNode(util.Function(input_names=['nc',
                                                    'realigned_file',
                                                    'wm_mask',
                                                    'csf_mask'],
                                       output_names=['noise_components'],
                                       function=extract_compcor_components),
                                       name='compcor',
                                       iterfield=['realigned_file',
                                                  'wm_mask',
                                                  'csf_mask'])

    glb_sig = pe.MapNode(util.Function(input_names=['realigned_file'],
                                       output_names=['global_component'],
                                       function=extract_global_component),
                                       name='glb_sig',
                                       iterfield=['realigned_file'])

    gm_sig = pe.MapNode(util.Function(input_names=['realigned_file',
                                                   'mask_file'],
                                      output_names=['mask_mean_component'],
                                      function=extract_frommask_component),
                                      name='gm_sig',
                                      iterfield=['realigned_file',
                                                 'mask_file'])
    wm_sig = gm_sig.clone(name='wm_sig')
    csf_sig = gm_sig.clone(name='csf_sig')

    fp1_sig = pe.MapNode(util.Function(input_names=['realigned_file'],
                                       output_names=['firstprinc_component'],
                                       function=extract_firstprinc_component),
                                       name='fp1_sig',
                                       iterfield=['realigned_file'])
    
    lt_sig = pe.MapNode(util.Function(input_names=['realigned_file'],
                                      output_names=['linear_component'],
                                      function=extract_linear_component),
                                      name='lt_sig',
                                      iterfield=['realigned_file'])

    addoutliers = pe.MapNode(util.Function(input_names=['global_component',
                                                        'compcor_components',
                                                        'wm_component',
                                                        'csf_component',
                                                        'gm_component',
                                                        'firstprinc_component',
                                                        'linear_component',
                                                        'motion_components',
                                                        'selector'],
                                           output_names=['filter_file'],
                                           function=create_filter_matrix),
                                           name='create_nuisance_filter',
                                           iterfield=['global_component',
                                                      'compcor_components',
                                                      'wm_component',
                                                      'csf_component',
                                                      'gm_component',
                                                      'firstprinc_component',
                                                      'linear_component',
                                                      'motion_components'])

#    remove_noise = pe.MapNode(fsl.FilterRegressor(filter_all=True),
#                              name='regress_nuisance',
#                              iterfield=['design_file', 'in_file'])
    remove_noise = pe.MapNode(util.Function(input_names=['realigned_file',
                                                         'regressors_file'],
                                            output_names=['residual_file'],
                                            function=extract_residuals),
                                            name='regress_nuisance',
                                            iterfield=['realigned_file', 'regressors_file'])

#    nuisance.connect(inputspec, 'realigned_file',
#                             linear_detrend, 'realigned_file')
    nuisance.connect(inputspec, 'realigned_file',
                             compcor, 'realigned_file')
    nuisance.connect(inputnode_nc, 'nc',
                             compcor, 'nc')
    nuisance.connect(inputspec, 'wm_mask',
                             compcor, 'wm_mask')
    nuisance.connect(inputspec, 'csf_mask',
                             compcor, 'csf_mask')
    nuisance.connect(inputspec, 'realigned_file',
                             glb_sig, 'realigned_file')
    nuisance.connect(inputspec, 'realigned_file',
                             gm_sig, 'realigned_file')
    nuisance.connect(inputspec, 'gm_mask',
                             gm_sig, 'mask_file')
    nuisance.connect(inputspec, 'realigned_file',
                             wm_sig, 'realigned_file')
    nuisance.connect(inputspec, 'realigned_file',
                             csf_sig, 'realigned_file')
    nuisance.connect(inputspec, 'wm_mask',
                             wm_sig, 'mask_file')
    nuisance.connect(inputspec, 'csf_mask',
                             csf_sig, 'mask_file')
    nuisance.connect(inputspec, 'realigned_file',
                             fp1_sig, 'realigned_file')
    nuisance.connect(inputspec, 'realigned_file',
                             lt_sig, 'realigned_file')
    nuisance.connect(glb_sig, 'global_component',
                             addoutliers, 'global_component')
    nuisance.connect(gm_sig, 'mask_mean_component',
                             addoutliers, 'gm_component')
    nuisance.connect(compcor, 'noise_components',
                             addoutliers, 'compcor_components')
    nuisance.connect(wm_sig, 'mask_mean_component',
                             addoutliers, 'wm_component')
    nuisance.connect(csf_sig, 'mask_mean_component',
                             addoutliers, 'csf_component')
    nuisance.connect(fp1_sig, 'firstprinc_component',
                             addoutliers, 'firstprinc_component')
    nuisance.connect(lt_sig, 'linear_component',
                             addoutliers, 'linear_component')
    nuisance.connect(inputspec, 'motion_components',
                             addoutliers, 'motion_components')
    nuisance.connect(inputnode_selector, 'selector',
                             addoutliers, 'selector')
    nuisance.connect(addoutliers, 'filter_file',
                             remove_noise, 'regressors_file')
    nuisance.connect(inputspec, 'realigned_file',
                             remove_noise, 'realigned_file')

    #Median angle correction on residual file
    nuisance.connect(remove_noise, 'residual_file',
                             median_angle, 'realigned_file')
    nuisance.connect(inputnode_target_angle_deg, 'target_angle_deg',
                             median_angle, 'target_angle_deg')

    nuisance.connect(median_angle, 'corrected_file',
                             outputspec, 'median_angle_corrected_file')
    nuisance.connect(median_angle, 'angles',
                             outputspec, 'angles')
    nuisance.connect(remove_noise, 'residual_file',
                             outputspec, 'residual_file')

    return nuisance
