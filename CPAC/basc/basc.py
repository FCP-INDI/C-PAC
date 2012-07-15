import numpy as np
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
        
    """
    indiv_stability_set = np.asarray(indiv_stability_list)
    V = indiv_stability_set.shape[2]
    
    G = np.zeros((V,V))
    for bootstrap_i in range(n_bootstraps):
        J = standard_bootstrap(indiv_stability_set).mean(0)
        G += adjacency_matrix(cluster_timeseries(J, k_clusters, similarity_metric = 'correlation')[:,np.newaxis])
    G /= n_bootstraps

    return G

def nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, k_clusters, cbb_block_size = None):
    print 'Calculating individual stability matrix of:', subject_file

    from utils import individual_stability_matrix
    import nibabel as nb

    data = nb.load(subject_file).get_data().astype('float64')
    roi_mask_file = nb.load(roi_mask_file).get_data().astype('float64').astype('bool')
    Y = data[roi_mask_file].T
    print '(timepoints,voxels):', Y.shape
    print 'Circular bootstrap block size:', cbb_block_size
    
    ism = individual_stability_matrix(Y, n_bootstraps, k_clusters, cbb_block_size=cbb_block_size)
    
    return ism

def open_dataset(subject_list, roi_mask):
    print 'Inside dataset method'
    print subject_list
    nifti_files = [1, 2, 3, 4, 5, 6]
    return nifti_files

def dummy_method(x_variable):
    print 'Inside dummy method'
    print x_variable
    return x_variable*10

def create_basc(name='basc'):
    """
    Bootstrap Analysis of Stable Clusters (BASC)
    
    This workflow performs group-level BASC.  
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow
    
    Returns
    -------
    
    
    Notes
    -----
    
    Workflow Inputs::
        
        inputspec.roi : nifti file
            Mask of region(s) of interest
        inputpsec.subjects : list of nifti files
            4-D timeseries of a group of subjects normalized to MNI space
    
    Workflow Outputs::
    
        outputspec.gsm : ndarray
            Group stability matrix
        outputspec.gsclusters: ndarray
            Matrix partitioning each cluster of the group stability matrix
        outputspec.gsmap: ndarray
            Group stability map using gsm and gscluster to calculate average within-cluster stability
        outputspec.gsclusters_img : nifti file
            3-D volume of brain regions partitioned with gsclusters
        outputspec.gsmap_img : nifti file
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
    .. [1] P. Bellec, P. Rosa-Neto, O. C. Lyttelton, H. Benali, and A. C. Evans,
    "Multi-level bootstrap analysis of stable clusters in resting-state fMRI.," 
    NeuroImage, vol. 51, no. 3, pp. 1126-39, Jul. 2010.
    
    Examples
    --------
    
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
                                                        'gsmao_img']),
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
                                output_names=['group_stability_matrix'],
                                function=group_stability_matrix),
                     name='gsm')

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
    
    
    
    return basc