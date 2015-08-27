import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def joint_mask(subjects_file_list, mask_file):
    """
    Creates a joint mask (intersection) common to all the subjects in a provided list
    and a provided mask
    
    Parameters
    ----------
    subjects_file_list : list of strings
        A length `N` list of file paths of the nifti files of subjects
    mask_file : string
        Path to a mask file in nifti format
    
    Returns
    -------
    joint_mask : string
        Path to joint mask file in nifti format
    
    """
    import nibabel as nb
    import numpy as np
    import os
    
    from CPAC.utils import safe_shape
    
    nii = nb.load(mask_file)
    
    mask = nii.get_data().astype('bool')
    for subject_file in subjects_file_list:
        sdata = nb.load(subject_file).get_data().astype(np.float64).sum(-1)
        if not safe_shape(sdata, mask): raise ValueError('Subject %s with volume shape %s conflicts \
                                                          with mask shape %s' % ( subject_file,
                                                                                  str(sdata.shape),
                                                                                  str(mask.shape) ) )
        mask *= sdata.astype('bool')

    print '... joint subject/roi mask loaded'

    img = nb.Nifti1Image(mask, header=nii.get_header(), affine=nii.get_affine())
    img_file = os.path.join(os.getcwd(), 'joint_mask.nii.gz')
    img.to_filename(img_file)
    
    return img_file

def nifti_cwas(subjects_file_list, mask_file, regressor, cols, f_samples, 
               voxel_range, strata=None):
    """
    Performs CWAS for a group of subjects
    
    Parameters
    ----------
    subjects_file_list : list of strings
        A length `N` list of file paths of the nifti files of subjects
    mask_file : string
        Path to a mask file in nifti format
    regressor : ndarray
        Vector of shape (`S`) or (`S`, `1`), `S` subjects
    cols : list
        todo
    f_samples : integer
        Number of pseudo f values to sample using a random permutation test
    voxel_range : tuple
        (start, end) tuple specify the range of voxels (inside the mask) to perform cwas on.
        Index ordering is based on the np.where(mask) command
    strata : ndarray (optional)
        todo
    
    Returns
    -------
    F_file : string
        .npy file of pseudo-F statistic calculated for every voxel
    p_file : string
        .npy file of significance probabilities of pseudo-F values
    voxel_range : tuple
        Passed on by the voxel_range provided in parameters, used to make parallelization
        easier
        
    """
    import nibabel as nb
    import numpy as np
    import os
    from CPAC.cwas import calc_cwas
    
    #Check regressor is a column vector
    if(len(regressor.shape) == 1):
        regressor = regressor[:, np.newaxis]
    elif(len(regressor.shape) != 2):
        raise ValueError('Bad regressor shape: %s' % str(regressor.shape))
    
    if(len(subjects_file_list) != regressor.shape[0]):
        raise ValueError('Number of subjects does not match regressor size')
    
    #Load the data to produce the joint mask
    mask = nb.load(mask_file).get_data().astype('bool')
    mask_indices = np.where(mask)
    #batch_indices = tuple([mask_index[voxel_range[0]:voxel_range[1]] for mask_index in mask_indices])
    
    #Reload the data again to actually get the values, sacrificing CPU for smaller memory footprint
    subjects_data = [ nb.load(subject_file).get_data().astype('float64')[mask_indices].T 
                        for subject_file in subjects_file_list ]
    #subjects_data = np.array(subjects_data)
    print '... subject data loaded', len(subjects_data), 'batch voxel range', voxel_range
    
    F_set, p_set = calc_cwas(subjects_data, regressor, cols, f_samples, voxel_range, strata)
    
    print '... writing cwas data to disk'
    cwd = os.getcwd()
    F_file = os.path.join(cwd, 'pseudo_F.npy')
    p_file = os.path.join(cwd, 'significance_p.npy')

    np.save(F_file, F_set)
    np.save(p_file, p_set)
    
    return F_file, p_file, voxel_range

def merge_cwas_batches(cwas_batches, mask_file):
    import numpy as np
    import nibabel as nb
    import os
    
    def volumize(mask, data):
        volume = np.zeros_like(mask, dtype=data.dtype)
        volume[np.where(mask==True)] = data
        return volume
    
    F_files, p_files, voxel_range = zip(*cwas_batches)
    end_voxel = np.array(voxel_range).max()
    
    nii = nb.load(mask_file)
    mask = nii.get_data().astype('bool')
    
    F_set = np.zeros(end_voxel)
    p_set = np.zeros(end_voxel)
    for F_file, p_file, voxel_range in cwas_batches:
        F_batch = np.load(F_file)
        p_batch = np.load(p_file)
        F_set[voxel_range[0]:voxel_range[1]] = F_batch
        p_set[voxel_range[0]:voxel_range[1]] = p_batch
    
    F_vol = volumize(mask, F_set)
    p_vol = volumize(mask, p_set)
    
    cwd = os.getcwd()
    F_file = os.path.join(cwd, 'pseudo_F_volume.nii.gz')
    p_file = os.path.join(cwd, 'p_significance_volume.nii.gz')
    
    img = nb.Nifti1Image(F_vol, header=nii.get_header(), affine=nii.get_affine())
    img.to_filename(F_file)
    img = nb.Nifti1Image(p_vol, header=nii.get_header(), affine=nii.get_affine())
    img.to_filename(p_file)
    
    return F_file, p_file

def create_cwas_batches(mask_file, batches):
    import nibabel as nb
    import numpy as np
    mask = nb.load(mask_file).get_data().astype('bool')
    nVoxels = mask.sum()
    
    print 'voxels: ', nVoxels
    batch_size = nVoxels/batches
    #batch_indices = np.arange(0, nVoxels, batch_size)
    
    batch_list = []
    for i in range(batches):
        batch_list.append((i*batch_size, (i+1)*batch_size))
    # Add remainder voxels to last batch
    batch_list[-1] = (batch_list[-1][0], nVoxels)
    
    return batch_list

def compile_theano_functions():
    """
    Returns compiled theano functions.  
    
    Notes
    -----
    Originally used to speedup multiplication of large matrices and vectors.  Caused strange 
    issue in nipype where nipype unecessarily reran nodes that use these compiled functions.
    Not used in current implementation.
    """
    import theano.tensor as T
    import theano
    
    def TnormCols(X):
        """
        Theano expression which centers and normalizes columns of X `||x_i|| = 1`
        """
        Xc = X - X.mean(0)
        return Xc/T.sqrt( (Xc**2.).sum(0) )
    
    def TzscorrCols(Xn):
        """
        Theano expression which returns Fisher transformed correlation values between columns of a
        normalized input, `X_n`.  Diagonal is set to zero.
        """
        C_X = T.dot(Xn.T, Xn)-T.eye(Xn.shape[1])
        return 0.5*T.log((1+C_X)/(1-C_X))
    
    X,Y = T.dmatrices('X','Y')
    tdot = theano.function([X,Y], T.dot(X,Y))
    tnormcols = theano.function([X], TnormCols(X))

    return tdot, tnormcols

def create_cwas(name='cwas'):
    """
    Connectome Wide Association Studies
    
    This workflow performs CWAS on a group of subjects.
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
        
    Returns
    -------
    cwas : nipype.pipeline.engine.Workflow
        CWAS workflow.
        
    Notes
    -----
    
    Workflow Inputs::
    
        inputspec.roi : string (nifti file)
            Mask of region(s) of interest
        inputpsec.subjects : list (nifti files)
            4-D timeseries of a group of subjects normalized to MNI space
        inputspec.regressor : list (float)
            Corresponding list of the regressor variable of shape (`N`) or (`N`,`1`), `N` subjects
        inputspec.cols : list (int)
            todo
        inputspec.f_samples : int
            Number of permutation samples to draw from the pseudo F distribution
        inputspec.strata : None or ndarray
            todo
        inputspec.parallel_nodes : integer
            Number of nodes to create and potentially parallelize over
        
    Workflow Outputs::

        outputspec.F_map : string (nifti file)
            Pseudo F values of CWAS
        outputspec.p_map : string (nifti file)
            Significance p values calculated from permutation tests
            
    CWAS Procedure:
    
    1. Calculate spatial correlation of a voxel
    2. Correlate spatial z-score maps for every subject pair
    3. Convert matrix to distance matrix, `1-r`
    4. Calculate MDMR statistics for the voxel
    5. Determine significance of MDMR statistics with permutation tests
    
    Workflow Graph:
    
    .. image:: ../images/cwas.dot.png
        :width: 500
        
    Detailed Workflow Graph:
    
    .. image:: ../images/cwas_detailed.dot.png
        :width: 500
    
    References
    ----------
    .. [1] Shehzad Z, Kelly C, Reiss PT, Emerson JW, McMahon K, Copland DA, Castellanos FX, Milham MP. An Analytic Framework for Connectome-Wide Association Studies. Under Review.
    
    """
    
    inputspec = pe.Node(util.IdentityInterface(fields=['roi',
                                                       'subjects',
                                                       'regressor', 
                                                       'cols', 
                                                       'f_samples', 
                                                       'strata', 
                                                       'parallel_nodes']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['F_map',
                                                        'p_map']),
                         name='outputspec')
    
    cwas = pe.Workflow(name=name)
    
    ccb = pe.Node(util.Function(input_names=['mask_file',
                                             'batches'],
                                output_names=['batch_list'],
                                function=create_cwas_batches),
                  name='cwas_batches')
    
    ncwas = pe.MapNode(util.Function(input_names=['subjects_file_list',
                                                  'mask_file',
                                                  'regressor', 
                                                  'cols', 
                                                  'f_samples',
#                                                  'compiled_func',
                                                  'voxel_range', 
                                                  'strata'],
                                     output_names=['result_batch'],
                                     function=nifti_cwas),
                       name='cwas_batch',
                       iterfield=['voxel_range'])
    
    mcwasb = pe.Node(util.Function(input_names=['cwas_batches',
                                                'mask_file'],
                                   output_names=['F_file',
                                                 'p_file'],
                                   function=merge_cwas_batches),
                     name='cwas_volumes')

#    ctf = pe.Node(util.Function(input_names=[],
#                                output_names=['compiled_dot_norm'],
#                                function=compile_theano_functions),
#                  name='theano_functions')
    
    jmask = pe.Node(util.Function(input_names=['subjects_file_list', 
                                               'mask_file'],
                                  output_names=['joint_mask'],
                                  function=joint_mask),
                    name='joint_mask')
    
    #Compute the joint mask
    cwas.connect(inputspec, 'subjects',
                 jmask, 'subjects_file_list')
    cwas.connect(inputspec, 'roi',
                 jmask, 'mask_file')

    #Create batches based on the joint mask
    cwas.connect(jmask, 'joint_mask',
                 ccb, 'mask_file')
    cwas.connect(inputspec, 'parallel_nodes',
                 ccb, 'batches')
    
    #Compute CWAS over batches of voxels
    cwas.connect(jmask, 'joint_mask',
                 ncwas, 'mask_file')
    cwas.connect(inputspec, 'subjects',
                 ncwas, 'subjects_file_list')
    cwas.connect(inputspec, 'regressor',
                 ncwas, 'regressor')
    cwas.connect(inputspec, 'f_samples',
                 ncwas, 'f_samples')
    cwas.connect(inputspec, 'cols',
                 ncwas, 'cols')
#    cwas.connect(ctf, 'compiled_dot_norm',
#                 ncwas, 'compiled_func')
    cwas.connect(ccb, 'batch_list',
                 ncwas, 'voxel_range')
    cwas.connect(inputspec, 'strata',
                 ncwas, 'strata')
    
    #Merge the computed CWAS data
    cwas.connect(ncwas, 'result_batch',
                 mcwasb, 'cwas_batches')
    cwas.connect(jmask, 'joint_mask',
                 mcwasb, 'mask_file')
    
    cwas.connect(mcwasb, 'F_file',
                 outputspec, 'F_map')
    cwas.connect(mcwasb, 'p_file',
                 outputspec, 'p_map')
    
    return cwas