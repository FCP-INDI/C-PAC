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
        inputspec.phenotypes: list (float)
            Corresponding list of the phenotypic variable for subjects
            
    Workflow Outputs::

        outputspec.diff_map : string (nifti file)
    
    References
    ----------
    .. [1] Shehzad Z, Kelly C, Reiss PT, Emerson JW, McMahon K, Copland DA, Castellanos FX, Milham MP. An Analytic Framework for Connectome-Wide Association Studies. Under Review.
    
    """