
import os
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.utils.interfaces.function import Function


from .cwas import (
    joint_mask,
    create_cwas_batches,
    merge_cwas_batches,
    nifti_cwas,
)


def create_cwas(name='cwas', working_dir=None, crash_dir=None):
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
    
        inputspec.subjects : dict (subject id: nifti files)
            4-D timeseries of a group of subjects normalized to MNI space
        inputspec.roi : string (nifti file)
            Mask of region(s) of interest
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
    
    .. exec::
        from CPAC.cwas import create_cwas
        wf = create_cwas()
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/create_cwas.dot'
        )

    Workflow Graph:
    
    .. image:: ../../images/generated/cwas.png
        :width: 500
        
    Detailed Workflow Graph:
    
    .. image:: ../../images/generated/cwas_detailed.png
        :width: 500
    
    References
    ----------
    .. [1] Shehzad Z, Kelly C, Reiss PT, Emerson JW, McMahon K, Copland DA, Castellanos FX, Milham MP. An Analytic Framework for Connectome-Wide Association Studies. Under Review.
    
    """

    if not working_dir:
        working_dir = os.path.join(os.getcwd(), 'MDMR_work_dir')
    if not crash_dir:
        crash_dir = os.path.join(os.getcwd(), 'MDMR_crash_dir')

    workflow = pe.Workflow(name=name)
    workflow.base_dir = working_dir
    workflow.config['execution'] = {'hash_method': 'timestamp',
                                    'crashdump_dir': os.path.abspath(crash_dir)}

    inputspec = pe.Node(util.IdentityInterface(fields=['roi',
                                                       'subjects',
                                                       'regressor',
                                                       'participant_column',
                                                       'columns',
                                                       'permutations',
                                                       'parallel_nodes']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['F_map',
                                                        'p_map',
                                                        'neglog_p_map']),
                         name='outputspec')

    ccb = pe.Node(Function(input_names=['mask_file',
                                        'batches'],
                           output_names='batch_list',
                           function=create_cwas_batches,
                           as_module=True),
                  name='cwas_batches')

    ncwas = pe.MapNode(Function(input_names=['subjects',
                                             'mask_file',
                                             'regressor_file',
                                             'participant_column',
                                             'columns_string',
                                             'permutations',
                                             'voxel_range'],
                                output_names=['result_batch'],
                                function=nifti_cwas,
                                as_module=True),
                       name='cwas_batch',
                       iterfield='voxel_range')

    jmask = pe.Node(Function(input_names=['subjects',
                                          'mask_file'],
                             output_names=['joint_mask'],
                             function=joint_mask,
                             as_module=True),
                    name='joint_mask')

    mcwasb = pe.Node(Function(input_names=['cwas_batches',
                                           'mask_file'],
                              output_names=['F_file',
                                            'p_file',
                                            'neglog_p_file'],
                              function=merge_cwas_batches,
                              as_module=True),
                     name='cwas_volumes')

    #Compute the joint mask
    workflow.connect(inputspec, 'subjects',
                     jmask, 'subjects')
    workflow.connect(inputspec, 'roi',
                     jmask, 'mask_file')

    #Create batches based on the joint mask
    workflow.connect(jmask, 'joint_mask',
                     ccb, 'mask_file')
    workflow.connect(inputspec, 'parallel_nodes',
                     ccb, 'batches')

    #Compute CWAS over batches of voxels
    workflow.connect(jmask, 'joint_mask',
                     ncwas, 'mask_file')
    workflow.connect(inputspec, 'subjects',
                     ncwas, 'subjects')
    workflow.connect(inputspec, 'regressor',
                     ncwas, 'regressor_file')
    workflow.connect(inputspec, 'permutations',
                     ncwas, 'permutations')
    workflow.connect(inputspec, 'participant_column',
                     ncwas, 'participant_column')
    workflow.connect(inputspec, 'columns',
                     ncwas, 'columns_string')

    workflow.connect(ccb, 'batch_list',
                     ncwas, 'voxel_range')

    #Merge the computed CWAS data
    workflow.connect(ncwas, 'result_batch',
                     mcwasb, 'cwas_batches')
    workflow.connect(jmask, 'joint_mask',
                     mcwasb, 'mask_file')

    workflow.connect(mcwasb, 'F_file', outputspec, 'F_map')
    workflow.connect(mcwasb, 'p_file', outputspec, 'p_map')
    workflow.connect(mcwasb, 'neglog_p_file', outputspec, 'neglog_p_map')

    return workflow
