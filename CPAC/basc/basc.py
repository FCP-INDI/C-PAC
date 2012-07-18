import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def group_stability_matrix(indiv_stability_list, n_bootstraps, k_clusters):
    """
    Calculate the group stability matrix of the entire dataset by bootstrapping the dataset
    
    Parameters
    ----------
    indiv_stability_list : array_like
        A length `N` list of matrices of shape (`T`, `V`), `N` subjects, `T` timepoints, `V` voxels
    n_bootstraps : integer
        Number of bootstrap datasets
    k_clusters : integer
        Number of clusters
    
    Returns
    -------
    G : array_like
        Group stability matrix of shape (`V`, `V`), `V` voxels
    clusters_G : array_like
        Length `V` array of cluster assignments for each voxel
    voxel_scores : array_like
        Length `V` array of within-cluster average values for each voxel
    """
    print 'Calculating group stability matrix for', len(indiv_stability_list), 'subjects.'
    
    from CPAC.basc import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average
    import numpy as np
    
    indiv_stability_set = np.asarray(indiv_stability_list)
    V = indiv_stability_set.shape[2]
    
    G = np.zeros((V,V))
    for bootstrap_i in range(n_bootstraps):
        J = standard_bootstrap(indiv_stability_set).mean(0)
        G += adjacency_matrix(cluster_timeseries(J, k_clusters, similarity_metric = 'data')[:,np.newaxis])
    G /= n_bootstraps

    clusters_G = cluster_timeseries(G, k_clusters, similarity_metric = 'data')
    voxel_scores = cluster_matrix_average(G, clusters_G)

    return G, clusters_G, voxel_scores

def nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, k_clusters, cbb_block_size = None):
    """
    Calculate the individual stability matrix for a single subject by using Circular Block Bootstrapping method
    for time-series data.
    
    Parameters
    ----------
    subject_file : string
        Nifti file of a subject
    roi_mask_file : string
        Region of interest (this method is too computationally intensive to perform on a whole-brain volume)
    n_bootstraps : integer
        Number of bootstraps
    k_clusters : integer
        Number of clusters
    cbb_block_size : integer, optional
        Size of the time-series block when performing circular block bootstrap
    
    Returns
    -------
    ism : array_like
        Individual stability matrix of shape (`V`, `V`), `V` voxels
    """
    print 'Calculating individual stability matrix of:', subject_file

    from CPAC.basc import individual_stability_matrix
    import nibabel as nb

    data = nb.load(subject_file).get_data().astype('float64')
    roi_mask_file = nb.load(roi_mask_file).get_data().astype('float64').astype('bool')
    Y = data[roi_mask_file].T
    print '(timepoints,voxels):', Y.shape
    print 'Circular bootstrap block size:', cbb_block_size
    
    ism = individual_stability_matrix(Y, n_bootstraps, k_clusters, cbb_block_size=cbb_block_size)
    
    return ism

def ndarray_to_vol(data_array, roi_mask_file, sample_file, filename):
    """
    Converts a numpy array to a nifti file given an roi mask
    
    Parameters
    ----------
    data_array : array_like
        A data array with the same column length and index alignment as the given roi_mask_file
    roi_mask_file : string
        Path of the roi_mask_file
    sample_file : string or list of strings
        Path of sample nifti file(s) to use for header of the output.  If list, the first file is chosen.
    filename : string
        Name of output file
        
    Returns
    -------
    img_file : string
        Path of the nifti file output
    
    """
    import nibabel as nb
    import numpy as np
    import os
    
    roi_mask_file = nb.load(roi_mask_file).get_data().astype('float64').astype('bool')
    out_vol = np.zeros_like(roi_mask_file, dtype=data_array.dtype)
    out_vol[roi_mask_file] = data_array.T

    nii = None
    if type(sample_file) is list:
        nii = nb.load(sample_file[0])
    else:
        nii = nb.load(sample_file)
        
    img = nb.Nifti1Image(out_vol, header=nii.get_header(), affine=nii.get_affine())
    img_file = os.path.join(os.getcwd(), filename)
    img.to_filename(img_file)
    
    return img_file

def create_basc(name='basc'):
    """
    Bootstrap Analysis of Stable Clusters (BASC)
    
    This workflow performs group-level BASC.  
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
    
    Returns
    -------
    basc : nipype.pipeline.engine.Workflow
        BASC workflow.
    
    Notes
    -----
    
    Workflow Inputs::
        
        inputspec.roi : string (nifti file)
            Mask of region(s) of interest
        inputpsec.subjects : list (nifti files)
            4-D timeseries of a group of subjects normalized to MNI space
    
    Workflow Outputs::
    
        outputspec.gsm : ndarray
            Group stability matrix
        outputspec.gsclusters: ndarray
            Matrix partitioning each cluster of the group stability matrix
        outputspec.gsmap: ndarray
            Group stability map using gsm and gscluster to calculate average within-cluster stability
        outputspec.gsclusters_img : string (nifti file)
            3-D volume of brain regions partitioned with gsclusters
        outputspec.gsmap_img : string (nifti file)
            3-D volume of brain regions associated with gs_map
    
    BASC Procedure:
    
    1. Generate individual stability matrices based on multiple clusterings of each bootstrap sample for a single subject
    2. Use stratified bootstrap to sample new datasets of subjects
    3. Calculate average stability matrix of each new dataset using individual stability matrices generated at step 1
    4. Cluster each average stabiilty matrix
    5. Average to create a group stability matrix
    6. Cluster the group stability matrix
    7. Calculate average within-cluster stability based on the clustering of step 6
    
    References
    ----------
    .. [1] P. Bellec, P. Rosa-Neto, O. C. Lyttelton, H. Benali, and A. C. Evans, "Multi-level bootstrap analysis of stable clusters in resting-state fMRI.," NeuroImage, vol. 51, no. 3, pp. 1126-39, Jul. 2010.
    
    Examples
    --------
    >>> from CPAC import basc

    """
    inputspec = pe.Node(util.IdentityInterface(fields=['subjects',
                                                       'roi',
                                                       'dataset_bootstraps',
                                                       'timeseries_bootstraps',
                                                       'k_clusters']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['gsm',
                                                        'gsclusters',
                                                        'gsmap',
                                                        'gsclusters_img',
                                                        'gsmap_img']),
                        name='outputspec')
    
    basc = pe.Workflow(name=name)
    
    nis = pe.MapNode(util.Function(input_names=['subject_file',
                                                'roi_mask_file',
                                                'n_bootstraps',
                                                'k_clusters',
                                                'cbb_block_size'],
                                   output_names=['individual_stability_matrices'],
                                   function=nifti_individual_stability),
                     name='nis',
                     iterfield=['subject_file'])
    
    gsm = pe.Node(util.Function(input_names=['indiv_stability_list',
                                             'n_bootstraps',
                                             'k_clusters'],
                                output_names=['group_stability_matrix',
                                              'group_stability_clusters',
                                              'group_stability_scores'],
                                function=group_stability_matrix),
                  name='gsm')

    gs_cluster_vol = pe.Node(util.Function(input_names=['data_array',
                                                        'roi_mask_file',
                                                        'sample_file',
                                                        'filename'],
                                           output_names=['img_file'],
                                           function=ndarray_to_vol),
                             name='gs_cluster_vol')

    gs_score_vol = pe.Node(util.Function(input_names=['data_array',
                                                      'roi_mask_file',
                                                      'sample_file',
                                                      'filename'],
                                         output_names=['img_file'],
                                         function=ndarray_to_vol),
                           name='gs_score_vol')

    # Gather outside workflow inputs
    basc.connect(inputspec, 'subjects',
                 nis, 'subject_file')
    basc.connect(inputspec, 'roi',
                 nis, 'roi_mask_file')
    basc.connect(inputspec, 'timeseries_bootstraps',
                 nis, 'n_bootstraps')
    basc.connect(inputspec, 'k_clusters',
                 nis, 'k_clusters')
    basc.connect(inputspec, 'dataset_bootstraps',
                 gsm, 'n_bootstraps')
    basc.connect(inputspec, 'k_clusters',
                 gsm, 'k_clusters')
    
    basc.connect(nis, 'individual_stability_matrices',
                 gsm, 'indiv_stability_list')

    basc.connect(inputspec, 'subjects',
                 gs_cluster_vol, 'sample_file')
    basc.connect(inputspec, 'roi',
                 gs_cluster_vol, 'roi_mask_file')
    gs_cluster_vol.inputs.filename = 'group_stability_clusters.nii.gz'
    
    basc.connect(gsm, 'group_stability_clusters',
                 gs_cluster_vol, 'data_array')

    basc.connect(inputspec, 'subjects',
                 gs_score_vol, 'sample_file')
    basc.connect(inputspec, 'roi',
                 gs_score_vol, 'roi_mask_file')
    gs_score_vol.inputs.filename = 'group_stability_scores.nii.gz'        
    
    basc.connect(gsm, 'group_stability_scores',
                 gs_score_vol, 'data_array')
    
    basc.connect(gsm, 'group_stability_matrix',
                 outputspec, 'gsm')
    basc.connect(gsm, 'group_stability_clusters',
                 outputspec, 'gsclusters')
    basc.connect(gsm, 'group_stability_scores',
                 outputspec, 'gsmap')

    basc.connect(gs_cluster_vol, 'img_file',
                 outputspec, 'gsclusters_img')
    basc.connect(gs_score_vol, 'img_file',
                 outputspec, 'gsmap_img')
    
    return basc