import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from .cwas import (
    joint_mask,
    create_cwas_batches,
    merge_cwas_batches,
    nifti_cwas,
)


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
        inputspec.subjects : list (nifti files)
            4-D timeseries of a group of subjects normalized to MNI space
        inputspec.regressor : list (float)
            Corresponding list of the regressor variable of shape (`N`) or (`N`,`1`), `N` subjects
        inputspec.cols : list (int)
            todo
        inputspec.f_samples : int
            Number of permutation samples to draw from the pseudo F distribution
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
                                                       'parallel_nodes']),
                        name='inputspec')
                        
    outputspec = pe.Node(util.IdentityInterface(fields=['F_map',
                                                        'p_map']),
                         name='outputspec')
    
    cwas = pe.Workflow(name=name)
    
    ccb = pe.Node(util.Function(input_names=['mask_file',
                                             'batches'],
                                output_names='batch_list',
                                function=create_cwas_batches),
                  name='cwas_batches')
    
    ncwas = pe.MapNode(util.Function(input_names=['subjects_file_list',
                                                  'mask_file',
                                                  'regressor', 
                                                  'cols', 
                                                  'f_samples',
                                                  'voxel_range'],
                                     output_names=['result_batch'],
                                     function=nifti_cwas),
                       name='cwas_batch',
                       iterfield='voxel_range')
    
    jmask = pe.Node(util.Function(input_names=['subjects_file_list', 
                                               'mask_file'],
                                  output_names=['joint_mask'],
                                  function=joint_mask),
                    name='joint_mask')
    
    mcwasb = pe.Node(util.Function(input_names=['cwas_batches',
                                                'mask_file'],
                                   output_names=['F_file',
                                                 'p_file'],
                                   function=merge_cwas_batches),
                     name='cwas_volumes')
    
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
                 
    cwas.connect(ccb, 'batch_list',
                 ncwas, 'voxel_range')
    
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