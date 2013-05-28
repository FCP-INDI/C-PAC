import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

def group_stability_matrix(indiv_stability_list, n_bootstraps, k_clusters, stratification=None):
    """
    Calculate the group stability matrix of the entire dataset by bootstrapping the dataset
    
    Parameters
    ----------
    indiv_stability_list : list of strings
        A length `N` list of file paths to numpy matrices of shape (`V`, `V`), `N` subjects, `V` voxels
    n_bootstraps : integer
        Number of bootstrap datasets
    k_clusters : integer
        Number of clusters
    stratification : array_like, optional
        List of integer entries denoting stratums for indiv_stability_list
    
    
    Returns
    -------
    G : array_like
        Group stability matrix of shape (`V`, `V`), `V` voxels
    clusters_G : array_like
        Length `V` array of cluster assignments for each voxel
    cluster_voxel_scores : array_like
        `K` by `V` matrix of within-cluster average values for each cluster of each voxel
    """
    print 'Calculating group stability matrix for', len(indiv_stability_list), 'subjects.'

    if stratification is not None:
        print 'Applying stratification to group dataset'
                
    from CPAC.basc import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average
    import numpy as np

    indiv_stability_set = np.asarray([np.load(ism_file) for ism_file in indiv_stability_list])
#    indiv_stability_set = indiv_stability_list
    print 'Individual stability list dimensions:', indiv_stability_set.shape
    
    V = indiv_stability_set.shape[2]
    
    G = np.zeros((V,V))
    for bootstrap_i in range(n_bootstraps):
        if stratification is not None:
            strata = np.unique(stratification)
            J = np.zeros((V,V))
            for stratum in strata:
                J += standard_bootstrap(indiv_stability_set[np.where(stratification == stratum)]).sum(0)
            J /= indiv_stability_set.shape[0]
        else:
            J = standard_bootstrap(indiv_stability_set).mean(0)
        G += adjacency_matrix(cluster_timeseries(J, k_clusters, similarity_metric = 'data')[:,np.newaxis])
    G /= n_bootstraps


    clusters_G = cluster_timeseries(G, k_clusters, similarity_metric = 'data')
    
    cluster_voxel_scores = cluster_matrix_average(G, clusters_G)
    
    # Cluster labels normally start from 0, start from 1 to provide contrast when viewing between 0 voxels
    clusters_G += 1

    return G, clusters_G, cluster_voxel_scores

def individual_group_clustered_maps(indiv_stability_list, clusters_G, roi_mask_file):
    """
    Calculate the individual stability maps of each subject based on the group stability clustering solution.
    
    Parameters
    ----------
    indiv_stability_list : list of strings
        A length `N` list of file paths to numpy matrices of shape (`V`, `V`), `N` subjects, `V` voxels
    clusters_G : array_like
        Length `V` array of cluster assignments for each voxel
        
    Returns
    -------
    individual_cluster_voxel_scores : list of strings
        A length `N` list of nifti files of the individual group clustered stability maps for each cluster.  Temporal
        dimension of each file corresponds to each subject.
    
    """
    import os
    import numpy as np
    from CPAC.basc import cluster_matrix_average, ndarray_to_vol
    
    indiv_stability_set = np.asarray([np.load(ism_file) for ism_file in indiv_stability_list])

    nSubjects = indiv_stability_set.shape[0]
    nVoxels = indiv_stability_set.shape[1]

    cluster_ids = np.unique(clusters_G)
    nClusters = cluster_ids.shape[0]
    
    cluster_voxel_scores = np.zeros((nClusters, nSubjects, nVoxels))
    for i in range(nSubjects):
        cluster_voxel_scores[:,i] = cluster_matrix_average(indiv_stability_set[i], clusters_G)
    
    icvs = []
    icvs_idx = 0
    for k in cluster_ids:
        icvs.append(ndarray_to_vol(cluster_voxel_scores[icvs_idx], roi_mask_file, roi_mask_file, 'individual_group_cluster%i_stability.nii.gz' % k))
        icvs_idx += 1
 
    return icvs

def nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, k_clusters, cbb_block_size = None, affinity_threshold = 0.5):
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
    affinity_threshold : float, optional
        Minimum threshold for similarity matrix based on correlation to create an edge
        
    Returns
    -------
    ism : array_like
        Individual stability matrix of shape (`V`, `V`), `V` voxels
    """
    print 'Calculating individual stability matrix of:', subject_file

    from CPAC.basc import individual_stability_matrix
    from CPAC.utils import safe_shape
    import nibabel as nb
    import numpy as np
    import os

    data = nb.load(subject_file).get_data().astype('float64')
    roi_mask_file = nb.load(roi_mask_file).get_data().astype('float64').astype('bool')
    
    if not safe_shape(roi_mask_file, data):
        raise ValueError('Subject %s with volume shape %s conflicts with mask shape %s' % (subject_file,
                                                                                           str(data.shape[:3]),
                                                                                           str(roi_mask_file.shape)) )
    
    Y = data[roi_mask_file].T
    print '(%i timepoints, %i voxels) and %i bootstraps' % (Y.shape[0], Y.shape[1], n_bootstraps)
    
    ism = individual_stability_matrix(Y, n_bootstraps, k_clusters, cbb_block_size=cbb_block_size, affinity_threshold=affinity_threshold)
    ism_file = os.path.join(os.getcwd(), 'individual_stability_matrix.npy')
    np.save(ism_file, ism)
    
    print 'Saving individual stability matrix %s for %s' % (ism_file, subject_file)
    
    return ism_file

def ndarray_to_vol(data_array, roi_mask_file, sample_file, filename):
    """
    Converts a numpy array to a nifti file given an roi mask
    
    Parameters
    ----------
    data_array : array_like
        A data array with the same column length and index alignment as the given roi_mask_file.  If data_array is two dimensional,
        first dimension is considered temporal dimension
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
    if(len(data_array.shape) == 1):
        out_vol = np.zeros_like(roi_mask_file, dtype=data_array.dtype)
        out_vol[roi_mask_file] = data_array
    elif(len(data_array.shape) == 2):
        out_vol = np.zeros((roi_mask_file.shape[0], roi_mask_file.shape[1], roi_mask_file.shape[2], data_array.shape[0]), dtype=data_array.dtype)
        out_vol[roi_mask_file] = data_array.T
    else:
        raise ValueError('data_array is %i dimensional, must be either 1 or 2 dimensional' % len(data_array.shape) )

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
        inputspec.dataset_bootstraps : integer
            Number of bootstrap samples of the dataset
        inputspec.timeseries_bootstraps : integer
            Number of bootstraps of each subject's timeseries
        inputspec.k_clusters : integer
            Number of clusters at both the individiual and group level
        inputspec.affinity_threshold : list (floats)
            Minimum threshold for similarity matrix based on correlation to create an edge
            
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
        outputspec.ismap_imgs : list of strings (nifti files)
            3-D volumes of stability scores of each cluster based on group clustering
            
    BASC Procedure:
    
    1. Generate individual stability matrices based on multiple clusterings of each bootstrap sample for a single subject
    2. Use stratified bootstrap to sample new datasets of subjects
    3. Calculate average stability matrix of each new dataset using individual stability matrices generated at step 1
    4. Cluster each average stabiilty matrix
    5. Average to create a group stability matrix
    6. Cluster the group stability matrix
    7. Calculate average within-cluster stability based on the clustering of step 6
    
    Workflow Graph:
    
    .. image:: ../images/basc.dot.png
        :width: 500
        
    Detailed Workflow Graph:
    
    .. image:: ../images/basc_detailed.dot.png
        :width: 500
        
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
                                                       'k_clusters',
                                                       'affinity_threshold']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['gsm',
                                                        'gsclusters',
                                                        'gsmap',
                                                        'gsclusters_img',
                                                        'gsmap_img',
                                                        'ismap_imgs']),
                        name='outputspec')
    
    basc = pe.Workflow(name=name)
    
    nis = pe.MapNode(util.Function(input_names=['subject_file',
                                                'roi_mask_file',
                                                'n_bootstraps',
                                                'k_clusters',
                                                'cbb_block_size',
                                                'affinity_threshold'],
                                   output_names=['individual_stability_matrices'],
                                   function=nifti_individual_stability),
                     name='individual_stability_matrices',
                     iterfield=['subject_file',
                                'affinity_threshold'])
    
    gsm = pe.Node(util.Function(input_names=['indiv_stability_list',
                                             'n_bootstraps',
                                             'k_clusters',
                                             'stratification'],
                                output_names=['group_stability_matrix',
                                              'group_stability_clusters',
                                              'group_stability_scores'],
                                function=group_stability_matrix),
                  name='group_stability_matrix')

    igcm = pe.Node(util.Function(input_names=['indiv_stability_list',
                                              'clusters_G',
                                              'roi_mask_file'],
                                 output_names=['individual_cluster_voxel_scores'],
                                 function=individual_group_clustered_maps),
                   name='individual_group_clustered_maps')

    gs_cluster_vol = pe.Node(util.Function(input_names=['data_array',
                                                        'roi_mask_file',
                                                        'sample_file',
                                                        'filename'],
                                           output_names=['img_file'],
                                           function=ndarray_to_vol),
                             name='group_stability_cluster_vol')

    gs_score_vol = pe.Node(util.Function(input_names=['data_array',
                                                      'roi_mask_file',
                                                      'sample_file',
                                                      'filename'],
                                         output_names=['img_file'],
                                         function=ndarray_to_vol),
                           name='group_stability_score_vol')

    # Gather outside workflow inputs
    basc.connect(inputspec, 'subjects',
                 nis, 'subject_file')
    basc.connect(inputspec, 'roi',
                 nis, 'roi_mask_file')
    basc.connect(inputspec, 'timeseries_bootstraps',
                 nis, 'n_bootstraps')
    basc.connect(inputspec, 'k_clusters',
                 nis, 'k_clusters')
    basc.connect(inputspec, 'affinity_threshold',
                 nis, 'affinity_threshold')
    
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
    
    basc.connect(nis, 'individual_stability_matrices',
                 igcm, 'indiv_stability_list')
    basc.connect(gsm, 'group_stability_clusters',
                 igcm, 'clusters_G')
    basc.connect(inputspec, 'roi',
                 igcm, 'roi_mask_file')
    
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
    
    basc.connect(igcm, 'individual_cluster_voxel_scores',
                 outputspec, 'ismap_imgs')
    return basc
