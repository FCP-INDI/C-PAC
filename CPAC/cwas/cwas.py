import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def nifti_cwas(subjects_file_list, mask_file, regressor, voxel_range):
    """
    Converts a dataset of subjects into a single `S` by `T` by `V` numpy array
    """
    import nibabel as nb
    import numpy as np
    import os
    from CPAC.cwas import calc_cwas

    mask = nb.load(mask_file).get_data().astype('bool')
    mask_indices = np.where(mask==True)
    batch_indices = tuple([mask_index[voxel_range[0]:voxel_range[1]] for mask_index in mask_indices])
    subjects_data = [nb.load(subject_file).get_data().astype('float64')[batch_indices].T for subject_file in subjects_file_list]
    subjects_data = np.array(subjects_data)
    print '... subject data loaded', subjects_data.shape, 'batch voxel range', voxel_range
    
    F_set, p_set = calc_cwas(subjects_data, regressor, 15000)

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
            Corresponding list of the regressor variable for subjects
        inputspec.parallel_nodes : integer
            Number of nodes to create and potentially parallelize over
        
    Workflow Outputs::

        outputspec.diff_map : string (nifti file)
    
    CWAS Procedure:
    
    1. Calculate spatial correlation of a voxel
    2. Correlate spatial z-score maps for every subject pair
    3. Convert matrix to distance matrix, `1-r`
    4. Calculate MDMR statistics for the voxel
    5. Determine significance of MDMR statistics with permutation tests
    
    References
    ----------
    .. [1] Shehzad Z, Kelly C, Reiss PT, Emerson JW, McMahon K, Copland DA, Castellanos FX, Milham MP. An Analytic Framework for Connectome-Wide Association Studies. Under Review.
    
    """
    
    inputspec = pe.Node(util.IdentityInterface(fields=['roi',
                                                       'subjects',
                                                       'regressor',
                                                       'parallel_nodes']),
                        name='inputspec')
    
    cwas = pe.Workflow(name=name)
    
    ccb = pe.Node(util.Function(input_names=['mask_file',
                                             'batches'],
                                output_names=['batch_list'],
                                function=create_cwas_batches),
                  name='cwas_batches')
    
    ncwas = pe.MapNode(util.Function(input_names=['subjects_file_list',
                                                  'mask_file',
                                                  'regressor',
                                                  'voxel_range'],
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
    
    
    cwas.connect(inputspec, 'roi',
                 ccb, 'mask_file')
    cwas.connect(inputspec, 'parallel_nodes',
                 ccb, 'batches')
    
    cwas.connect(inputspec, 'subjects',
                 ncwas, 'subjects_file_list')
    cwas.connect(inputspec, 'roi',
                 ncwas, 'mask_file')
    cwas.connect(inputspec, 'regressor',
                 ncwas, 'regressor')
    cwas.connect(ccb, 'batch_list',
                 ncwas, 'voxel_range')
    
    cwas.connect(ncwas, 'result_batch',
                 mcwasb, 'cwas_batches')
    cwas.connect(inputspec, 'roi',
                 mcwasb, 'mask_file')
    
    
    
    
    return cwas