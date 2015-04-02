import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants

#from nipype import logging
#logger = logging.getLogger('workflow')

def bandpass_voxels(realigned_file, bandpass_freqs, sample_period = None):
    """
    Performs ideal bandpass filtering on each voxel time-series.
    
    Parameters
    ----------
    realigned_file : string
        Path of a realigned nifti file.
    bandpass_freqs : tuple
        Tuple containing the bandpass frequencies. (LowCutoff, HighCutoff)
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

    print 'Frequency filtering using sample period: ', sample_period, 'sec'

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
                   wm_sig_file = None,
                   csf_sig_file = None,
                   gm_sig_file = None,
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
        Path of csv file of regressors used.  Filename corresponds to the name of each
        regressor in each column.
        
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
    import scipy
    from CPAC.nuisance import calc_compcor_components
    
    
    nii = nb.load(subject)
    data = nii.get_data().astype(np.float64)
    global_mask = (data != 0).sum(-1) != 0
    
    
    #Check and define regressors which are provided from files
    if wm_sig_file is not None:
        wm_sigs = np.load(wm_sig_file)
        if wm_sigs.shape[1] != data.shape[3]:
            raise ValueError('White matter signals length %d do not match data timepoints %d' % (wm_sigs.shape[1], data.shape[3]))
    if csf_sig_file is not None:
        csf_sigs = np.load(csf_sig_file)
        if csf_sigs.shape[1] != data.shape[3]:
            raise ValueError('CSF signals length %d do not match data timepoints %d' % (csf_sigs.shape[1], data.shape[3]))
    if gm_sig_file is not None:
        gm_sigs = np.load(gm_sig_file)
        if gm_sigs.shape[1] != data.shape[3]:
            raise ValueError('Grey matter signals length %d do not match data timepoints %d' % (gm_sigs.shape[1], data.shape[3]))
        
    if motion_file is not None:
        motion = np.genfromtxt(motion_file)
        if motion.shape[0] != data.shape[3]:
            raise ValueError('Motion parameters %d do not match data timepoints %d' % (motion.shape[0], data.shape[3]) )

    #Calculate regressors
    regressor_map = {'constant' : np.ones((data.shape[3],1))}
    if(selector['compcor']):
        print 'compcor_ncomponents ', compcor_ncomponents
        regressor_map['compcor'] = calc_compcor_components(data, compcor_ncomponents, wm_sigs, csf_sigs)
    
    if(selector['wm']):
        regressor_map['wm'] = wm_sigs.mean(0)
        
    if(selector['csf']):
        regressor_map['csf'] = csf_sigs.mean(0)
        
    if(selector['gm']):
        regressor_map['gm'] = gm_sigs.mean(0)
        
    if(selector['global']):
        regressor_map['global'] = data[global_mask].mean(0)
        
    if(selector['pc1']):
        bdata = data[global_mask].T
        bdatac = bdata - np.tile(bdata.mean(0), (bdata.shape[0], 1))
        U, S, Vh = np.linalg.svd(bdatac, full_matrices=False)
        regressor_map['pc1'] = U[:,0]
        
    if(selector['motion']):
        regressor_map['motion'] = motion
        
    if(selector['linear']):
        regressor_map['linear'] = np.arange(0, data.shape[3])
    
    if(selector['quadratic']):
        regressor_map['quadratic'] = np.arange(0, data.shape[3])**2
    
    print 'Regressors include: ', regressor_map.keys()
    
    X = np.zeros((data.shape[3], 1))
    csv_filename = ''
    for rname, rval in regressor_map.items():
        X = np.hstack((X, rval.reshape(rval.shape[0],-1)))
        csv_filename += '_' + rname
    X = X[:,1:]
    
    csv_filename = csv_filename[1:]
    csv_filename += '.csv'
    csv_filename = os.path.join(os.getcwd(), csv_filename)
    np.savetxt(csv_filename, X, delimiter='\t')
    
    print 'Regressors dim: ', X.shape, ' starting regression'
    
    Y = data[global_mask].T
    B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    Y_res = Y - X.dot(B)
    
    data[global_mask] = Y_res.T
    
    print 'Writing residual and regressors'
    img = nb.Nifti1Image(data, header=nii.get_header(), affine=nii.get_affine())
    residual_file = os.path.join(os.getcwd(), 'residual.nii.gz')
    img.to_filename(residual_file)
    
    #Easier to read for debugging purposes
    regressors_file = os.path.join(os.getcwd(), 'nuisance_regressors.mat')

    if scipy.__version__ == '0.7.0':
        scipy.io.savemat(regressors_file, regressor_map)                        ### for scipy v0.7.0
    else:
        scipy.io.savemat(regressors_file, regressor_map, oned_as='column')   ### for scipy v0.12: OK

    
    
    return residual_file, csv_filename



def extract_tissue_data(data_file,
                        ventricles_mask_file,
                        wm_seg_file, csf_seg_file, gm_seg_file,
                        wm_threshold=0.0, csf_threshold=0.0, gm_threshold=0.0):
    import numpy as np
    import nibabel as nb
    import os    
    from CPAC.nuisance import erode_mask
    from CPAC.utils import safe_shape
    
    print 'Tissues extraction thresholds wm %d, csf %d, gm %d' % (wm_threshold,
                                                                  csf_threshold,
                                                                  gm_threshold)

    try:
        data = nb.load(data_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % data_file)


    try:
        lat_ventricles_mask = nb.load(ventricles_mask_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % lat_ventricles_mask)


    if not safe_shape(data, lat_ventricles_mask):
        raise ValueError('Spatial dimensions for data and the lateral ventricles mask do not match')

    try:
        wm_seg = nb.load(wm_seg_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % wm_seg)


    if not safe_shape(data, wm_seg):
        raise ValueError('Spatial dimensions for data, white matter segment do not match')

    wm_mask = erode_mask(wm_seg > wm_threshold)
    wm_sigs = data[wm_mask]
    file_wm = os.path.join(os.getcwd(), 'wm_signals.npy')
    np.save(file_wm, wm_sigs)
    del wm_sigs

    try:
        csf_seg = nb.load(csf_seg_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % csf_seg)


    if not safe_shape(data, csf_seg):
        raise ValueError('Spatial dimensions for data, cerebral spinal fluid segment do not match')

    # Only take the CSF at the lateral ventricles as labeled in the Harvard
    # Oxford parcellation regions 4 and 43
    csf_mask = (csf_seg > csf_threshold)*(lat_ventricles_mask==1)
    csf_sigs = data[csf_mask]
    file_csf = os.path.join(os.getcwd(), 'csf_signals.npy')
    np.save(file_csf, csf_sigs)
    del csf_sigs


    try:
        gm_seg = nb.load(gm_seg_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % gm_seg)


    if not safe_shape(data, gm_seg):
        raise ValueError('Spatial dimensions for data, gray matter segment do not match')


    gm_mask = erode_mask(gm_seg > gm_threshold)
    gm_sigs = data[gm_mask]
    file_gm = os.path.join(os.getcwd(), 'gm_signals.npy')
    np.save(file_gm, gm_sigs)
    del gm_sigs



    nii = nb.load(wm_seg_file)
    wm_mask_file = os.path.join(os.getcwd(), 'wm_mask.nii.gz')
    csf_mask_file = os.path.join(os.getcwd(), 'csf_mask.nii.gz')
    gm_mask_file = os.path.join(os.getcwd(), 'gm_mask.nii.gz')
    nb.Nifti1Image(wm_mask, header=nii.get_header(), affine=nii.get_affine()).to_filename(wm_mask_file)
    nb.Nifti1Image(csf_mask, header=nii.get_header(), affine=nii.get_affine()).to_filename(csf_mask_file)
    nb.Nifti1Image(gm_mask, header=nii.get_header(), affine=nii.get_affine()).to_filename(gm_mask_file)

    return file_wm, file_csf, file_gm



def create_nuisance(use_ants, name='nuisance'):
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
        inputspec.mni_to_anat_linear_xfm : string (nifti file)
            Corresponding MNI to anatomical linear transformation 
        inputspec.func_to_anat_linear_xfm : string (nifti file)
            Corresponding EPI to anatomical linear transformation
        inputspec.harvard_oxford_mask : string (nifti file)
            Harvard Oxford parcellation for ventrical locations
        inputspec.motion_components : string (text file)
            Corresponding rigid-body motion parameters.  Matrix in the file should be of shape 
            (`T`, `R`), `T` timepoints and `R` motion parameters.
        inputspec.selector : dictionary
        inputspec.compcor_ncomponents : integer
        
    Workflow Outputs::

        outputspec.subject : string (nifti file)
            Path of residual file in nifti format
        outputspec.regressors : string (mat file)
            Path of csv file of regressors used.  Filename corresponds to the name of each
            regressor in each column.
            
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
                                                       'mni_to_anat_linear_xfm',
                                                       'anat_to_mni_initial_xfm',
                                                       'anat_to_mni_rigid_xfm',
                                                       'anat_to_mni_affine_xfm',
                                                       'func_to_anat_linear_xfm',
                                                       'lat_ventricles_mask',
                                                       'motion_components',
                                                       'selector',
                                                       'compcor_ncomponents',
                                                       'template_brain']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                        'regressors']),
                         name='outputspec')


    # Resampling the masks from 1mm to 2mm, but remaining in subject space
    wm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='wm_anat_to_2mm_flirt_applyxfm')
    wm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    wm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'wm_mask', wm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'wm_mask', wm_anat_to_2mm, 'reference')
 

    # Resampling the masks from 1mm to 2mm, but remaining in subject space
    csf_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='csf_anat_to_2mm_flirt_applyxfm')
    csf_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    csf_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'csf_mask', csf_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'csf_mask', csf_anat_to_2mm, 'reference')

    
    # Resampling the masks from 1mm to 2mm, but remaining in subject space
    gm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='gm_anat_to_2mm_flirt_applyxfm')
    gm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    gm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'gm_mask', gm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'gm_mask', gm_anat_to_2mm, 'reference')


    func_to_2mm = pe.Node(interface=fsl.FLIRT(), name='func_to_2mm_flirt_applyxfm')
    func_to_2mm.inputs.args = '-applyisoxfm 2'

    nuisance.connect(inputspec, 'subject', func_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'csf_mask', func_to_2mm, 'reference')
    nuisance.connect(inputspec, 'func_to_anat_linear_xfm', func_to_2mm, 'in_matrix_file')


    if use_ants == True:

        collect_linear_transforms = pe.Node(util.Merge(3), name='ho_mni_to_2mm_ants_collect_linear_transforms')

        ho_mni_to_2mm = pe.Node(interface=ants.ApplyTransforms(), name='ho_mni_to_2mm_ants_applyxfm')

        ho_mni_to_2mm.inputs.invert_transform_flags = [True, True, True]
        ho_mni_to_2mm.inputs.interpolation = 'NearestNeighbor'
        ho_mni_to_2mm.inputs.dimension = 3

        nuisance.connect(inputspec, 'anat_to_mni_initial_xfm', collect_linear_transforms, 'in1')
        nuisance.connect(inputspec, 'anat_to_mni_rigid_xfm', collect_linear_transforms, 'in2')
        nuisance.connect(inputspec, 'anat_to_mni_affine_xfm', collect_linear_transforms, 'in3')

        nuisance.connect(collect_linear_transforms, 'out', ho_mni_to_2mm, 'transforms')

        nuisance.connect(inputspec, 'lat_ventricles_mask', ho_mni_to_2mm, 'input_image')
        nuisance.connect(csf_anat_to_2mm, 'out_file', ho_mni_to_2mm, 'reference_image')

        #resample_to_2mm = pe.Node(interface=afni.Resample(), name='resample_to_2mm_ants_output'
        


    else:

        ho_mni_to_2mm = pe.Node(interface=fsl.FLIRT(), name='ho_mni_to_2mm_flirt_applyxfm')
        ho_mni_to_2mm.inputs.args = '-applyisoxfm 2'
        ho_mni_to_2mm.inputs.interp = 'nearestneighbour'

        nuisance.connect(inputspec, 'mni_to_anat_linear_xfm', ho_mni_to_2mm, 'in_matrix_file')
        nuisance.connect(inputspec, 'lat_ventricles_mask', ho_mni_to_2mm, 'in_file')
        nuisance.connect(inputspec, 'csf_mask', ho_mni_to_2mm, 'reference')


    tissue_masks = pe.Node(util.Function(input_names=['data_file',
                                                      'ventricles_mask_file',
                                                      'wm_seg_file', 'csf_seg_file', 'gm_seg_file',
                                                      'wm_threshold', 'csf_threshold', 'gm_threshold'],
                                         output_names=['file_wm', 'file_csf', 'file_gm'],
                                         function=extract_tissue_data),
                           name='tissue_masks')
    


    nuisance.connect(func_to_2mm, 'out_file', tissue_masks, 'data_file')
    nuisance.connect(wm_anat_to_2mm, 'out_file', tissue_masks, 'wm_seg_file')
    nuisance.connect(csf_anat_to_2mm, 'out_file', tissue_masks, 'csf_seg_file')
    nuisance.connect(gm_anat_to_2mm, 'out_file', tissue_masks, 'gm_seg_file')

    if use_ants == True:
        nuisance.connect(ho_mni_to_2mm, 'output_image', tissue_masks, 'ventricles_mask_file')

    else:
        nuisance.connect(ho_mni_to_2mm, 'out_file', tissue_masks, 'ventricles_mask_file')



    calc_r = pe.Node(util.Function(input_names=['subject',
                                                'selector',
                                                'wm_sig_file',
                                                'csf_sig_file',
                                                'gm_sig_file',
                                                'motion_file',
                                                'compcor_ncomponents'],
                                   output_names=['residual_file',
                                                'regressors_file'],
                                   function=calc_residuals),
                     name='residuals')
    nuisance.connect(inputspec, 'subject',
                     calc_r, 'subject')
    nuisance.connect(tissue_masks, 'file_wm',
                     calc_r, 'wm_sig_file')
    nuisance.connect(tissue_masks, 'file_csf',
                     calc_r, 'csf_sig_file')
    nuisance.connect(tissue_masks, 'file_gm',
                     calc_r, 'gm_sig_file')
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
