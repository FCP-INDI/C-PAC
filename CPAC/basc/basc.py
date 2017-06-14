import os
import numpy as np
import nibabel as nb
import sys
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import imp
from os.path import expanduser
#
#home2 = expanduser("~")
#
#
#CPAC = imp.load_source('CPAC', home2 + '/C-PAC')
#
#bascutils = imp.load_source('bascutils', home2 + '/C-PAC/CPAC/basc/utils.py')
#
##basc = imp.load_source('basc', home2 + '/C-PAC/CPAC/basc/basc.py')
#
#utilsutils = imp.load_source('utilsutils', home2 + '/C-PAC/CPAC/utils/utils.py')
#
#import CPAC
#from bascutils import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average, individual_stability_matrix
##from basc import group_stability_matrix, individual_group_clustered_maps, individual_stability_matrix, nifti_individual_stability, ndarray_to_vol, create_basc
#from utilsutils import safe_shape

#sys.path.insert(0, '~/C-PAC/CPAC')




#import CPAC

#you just do "import sys" (if not already imported), then "sys.path.insert(0, /path/to/install/dir)", then import CPAC

def group_stability_matrix(indiv_stability_list, n_bootstraps, n_clusters, stratification=None):
    """
    Calculate the group stability matrix of the entire dataset by bootstrapping the dataset

    Parameters
    ----------
    indiv_stability_list : list of strings
        A length `N` list of file paths to numpy matrices of shape (`V`, `V`), `N` subjects, `V` voxels
    n_bootstraps : integer
        Number of bootstrap datasets
    n_clusters : integer
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
#
    import os
    import numpy as np
    import nibabel as nb
    import utils
#    import sys
    #import CPAC
    #import nipype.pipeline.engine as pe
    #import nipype.interfaces.utility as util

#    from CPAC.basc.utils import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average, individual_stability_matrix
#    from CPAC.basc import group_stability_matrix, individual_group_clustered_maps, individual_stability_matrix, nifti_individual_stability, ndarray_to_vol, create_basc
#    from CPAC.utils.utils import safe_shape


    print 'Calculating group stability matrix for', len(indiv_stability_list), 'subjects.'

    if stratification is not None:
        print 'Applying stratification to group dataset'




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
                J += utils.standard_bootstrap(indiv_stability_set[np.where(stratification == stratum)]).sum(0)
            J /= indiv_stability_set.shape[0]
        else:
            J = utils.standard_bootstrap(indiv_stability_set).mean(0)
        G += utils.adjacency_matrix(utils.cluster_timeseries(J, n_clusters, similarity_metric = 'correlation')[:,np.newaxis])
    G /= n_bootstraps


    clusters_G = utils.cluster_timeseries(G, n_clusters, similarity_metric = 'correlation')

    cluster_voxel_scores = utils.cluster_matrix_average(G, clusters_G)

    # Cluster labels normally start from 0, start from 1 to provide contrast when viewing between 0 voxels
    clusters_G += 1

    gsm_file = os.path.join(os.getcwd(), 'group_stability_matrix.npy')
    np.save(gsm_file, G)
    
    clusters_G_file = os.path.join(os.getcwd(), 'clusters_G.npy')
    np.save(clusters_G_file, clusters_G)
    
    cluster_voxel_scores_file = os.path.join(os.getcwd(), 'cluster_voxel_scores.npy')
    np.save(cluster_voxel_scores_file, cluster_voxel_scores)
    #return G, clusters_G, cluster_voxel_scores

    return G, clusters_G, cluster_voxel_scores, gsm_file, clusters_G_file, cluster_voxel_scores_file



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
    import numpy as np
    import utils
    import basc

    indiv_stability_set = np.asarray([np.load(ism_file) for ism_file in indiv_stability_list])

    nSubjects = indiv_stability_set.shape[0]
    nVoxels = indiv_stability_set.shape[1]

    cluster_ids = np.unique(clusters_G)
    nClusters = cluster_ids.shape[0]

    cluster_voxel_scores = np.zeros((nClusters, nSubjects, nVoxels))
    for i in range(nSubjects):
        cluster_voxel_scores[:,i] = utils.cluster_matrix_average(indiv_stability_set[i], clusters_G)

    icvs = []
    icvs_idx = 0
    for k in cluster_ids:
        icvs.append(basc.ndarray_to_vol(cluster_voxel_scores[icvs_idx], roi_mask_file, roi_mask_file, 'individual_group_cluster%i_stability.nii.gz' % k))
        icvs_idx += 1

    return icvs

#def individual_stability_matrix(Y1, n_bootstraps, k_clusters, Y2=None, cross_cluster=False, cbb_block_size = None, affinity_threshold = 0.5):
#    """
#    Calculate the individual stability matrix of a single subject by bootstrapping their time-series
#
#    Parameters
#    ----------
#    Y1 : array_like
#        A matrix of shape (`V`, `N`) with `V` voxels `N` timepoints
#    Y2 : array_like
#        A matrix of shape (`V`, `N`) with `V` voxels `N` timepoints
#        For Cross-cluster solutions- this will be the matrix by which Y1 is clustered
#    n_bootstraps : integer
#        Number of bootstrap samples
#    k_clusters : integer
#        Number of clusters
#    cbb_block_size : integer, optional
#        Block size to use for the Circular Block Bootstrap algorithm
#    affinity_threshold : float, optional
#        Minimum threshold for similarity matrix based on correlation to create an edge
#
#    Returns
#    -------
#    S : array_like
#        A matrix of shape (`V1`, `V1`), each element v1_{ij} representing the stability of the adjacency of voxel i with voxel j
#    """
#    if affinity_threshold < 0.0:
#        raise ValueError('affinity_threshold %d must be non-negative value' % affinity_threshold)
#
#    #flipped the N and V values bc originally data was being put in transposed
#    N1 = Y1.shape[1]
#    V1 = Y1.shape[0]
#
#    if(cbb_block_size is None):
#        cbb_block_size = int(np.sqrt(N1))
#
#    S = np.zeros((V1, V1))
#
#    if (cross_cluster is True):
#        for bootstrap_i in range(n_bootstraps):
#            N2 = Y2.shape[1]
#            cbb_block_size2 = int(np.sqrt(N2))
#            Y_b1 = timeseries_bootstrap(Y1, cbb_block_size)
#            Y_b2 = timeseries_bootstrap(Y2, cbb_block_size2)
#            S += adjacency_matrix(cross_cluster_timeseries(Y_b1, Y_b2, k_clusters, similarity_metric = 'correlation'))
#        S /= n_bootstraps
#    else:
#        for bootstrap_i in range(n_bootstraps):
#            Y_b1 = timeseries_bootstrap(Y1, cbb_block_size)
#            S += adjacency_matrix(cluster_timeseries(Y_b1.T, k_clusters, similarity_metric = 'correlation', affinity_threshold = affinity_threshold)[:,np.newaxis])
#        S /= n_bootstraps
#
#    return S




def nifti_individual_stability(subject_file, roi_mask_file, n_bootstraps, n_clusters, output_size, cross_cluster=False, roi2_mask_file=None, cbb_block_size=None, affinity_threshold=0.5):
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
    n_clusters : integer
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


    import os
    import numpy as np
    import nibabel as nb
    import utils
    #from utils import individual_stability_matrix
    #from CPAC.basc.utils import standard_bootstrap, adjacency_matrix, cluster_timeseries, cluster_matrix_average, individual_stability_matrix


    print 'Calculating individual stability matrix of:', subject_file


    data = nb.load(subject_file).get_data().astype('float64')


    if (roi2_mask_file != None):
        roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float64').astype('bool')
        roi2_mask_nparray = nb.load(roi2_mask_file).get_data().astype('float64').astype('bool')

        #safe shape may be broken?
#        if not safe_shape(roi_mask_file, data):
#            raise ValueError('Subject %s with volume shape %s conflicts with mask shape %s' % (subject_file,
#                                                                                               str(data.shape[:3]),
#                                     s                                                          str(roi_mask_file.shape)) )
#        if not safe_shape(roi2_mask_file, data):
#            raise ValueError('Subject %s with volume shape %s conflicts with mask shape %s' % (subject_file,
#                                                                                           str(data.shape[:3]),
#                                                                                           str(roi2_mask_file.shape)) )
        roi_mask_file_nb=nb.load(roi_mask_file)
        roi2_mask_file_nb=nb.load(roi2_mask_file)
        
        
        roi1data = data[roi_mask_nparray]
        print '(%i voxels, %i timepoints and %i bootstraps)' % (roi1data.shape[0], roi1data.shape[1], n_bootstraps)
        roi2data = data[roi2_mask_nparray]
        print '(%i voxels, %i timepoints and %i bootstraps)' % (roi2data.shape[0], roi2data.shape[1], n_bootstraps)

        Y1_compressed = utils.data_compression(roi1data, roi_mask_file_nb, output_size).T
        Y2_compressed = utils.data_compression(roi2data, roi2_mask_file_nb, output_size).T

        ism = utils.individual_stability_matrix(Y1_compressed, n_bootstraps, n_clusters, Y2_compressed, cross_cluster, cbb_block_size, affinity_threshold)


        #def individual_stability_matrix(Y1, n_bootstraps, k_clusters, Y2=None, cross_cluster=False, cbb_block_size = None, affinity_threshold = 0.5):

        ism_file = os.path.join(os.getcwd(), 'individual_stability_matrix.npy')
        np.save(ism_file, ism)
        print 'Saving individual stability matrix %s for %s' % (ism_file, subject_file)
        return ism_file

    else:

        roi_mask_nparray = nb.load(roi_mask_file).get_data().astype('float64').astype('bool')

#        if not safe_shape(roi_mask_file, data):
#            raise ValueError('Subject %s with volume shape %s conflicts with mask shape %s' % (subject_file,
#                                                                                               str(data.shape[:3]),
#                                                                                               str(roi_mask_file.shape)) )

        roi1data = data[roi_mask_nparray]
        print '(%i voxels, %i timepoints and %i bootstraps' % (roi1data.shape[0], roi1data.shape[1], n_bootstraps)
        
        roi_mask_file_nb=nb.load(roi_mask_file)
        Y1_compressed = data_compression(roi1data, roi_mask_file_nb, output_size).T


        ism = utils.individual_stability_matrix(Y1_compressed, n_bootstraps, n_clusters, cbb_block_size, affinity_threshold)
        ism_file = os.path.join(os.getcwd(), 'individual_stability_matrix.npy')
        np.save(ism_file, ism)
        print 'Saving individual stability matrix %s for %s' % (ism_file, subject_file)
        return ism_file
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
        inputspec.n_clusters : integer
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
    import nipype.pipeline.engine as pe


    inputspec = pe.Node(util.IdentityInterface(fields=['subject_file_list',
                                                       'roi_mask_file',
                                                       'dataset_bootstraps',
                                                       'timeseries_bootstraps',
                                                       'n_clusters',
                                                       'output_size',
                                                       'cross_cluster',
                                                       'roi2_mask_file',
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
                                                'n_clusters',
                                                'output_size',
                                                'cross_cluster',
                                                'roi2_mask_file',
                                                'cbb_block_size',
                                                'affinity_threshold'],
                                   output_names=['ism_file'],
                                   function=nifti_individual_stability),
                     name='individual_stability_matrices',
                     iterfield=['subject_file',
                                'affinity_threshold'])

    nis.inputs.cbb_block_size=None
    #return G, clusters_G, clusters_voxel_scores, gsm_file, clusters_G_file, clusters_voxel_scores_file


    gsm = pe.Node(util.Function(input_names=['indiv_stability_list',
                                             'n_bootstraps',
                                             'n_clusters',
                                             'stratification'],
                                output_names=['G',
                                              'clusters_G',
                                              'cluster_voxel_scores',
                                              'gsm_file',
                                              'clusters_G_file',
                                              'cluster_voxel_scores_file'],
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



    #run_basc_workflow(subject_file_list, roi_mask_file, dataset_bootstraps, timeseries_bootstraps, k_clusters, cross_cluster=cross_cluster, roi2_mask_file=roi2_mask_file, affinity_threshold=affinity_threshold, out_dir=out_dir, run=run)


    # Gather outside workflow inputs
    basc.connect(inputspec, 'subject_file_list',
                 nis, 'subject_file')
    basc.connect(inputspec, 'roi_mask_file',
                 nis, 'roi_mask_file')
    basc.connect(inputspec, 'timeseries_bootstraps',
                 nis, 'n_bootstraps')
    basc.connect(inputspec, 'output_size',
                 nis, 'output_size')
    basc.connect(inputspec, 'roi2_mask_file',
                 nis, 'roi2_mask_file')
    basc.connect(inputspec, 'n_clusters',
                 nis, 'n_clusters')
    basc.connect(inputspec, 'affinity_threshold',
                 nis, 'affinity_threshold')
    basc.connect(inputspec, 'cross_cluster',
                 nis, 'cross_cluster')
   
    

    basc.connect(inputspec, 'dataset_bootstraps',
                 gsm, 'n_bootstraps')
    basc.connect(inputspec, 'n_clusters',
                 gsm, 'n_clusters')

    basc.connect(nis, 'ism_file',
                 gsm, 'indiv_stability_list')

    basc.connect(inputspec, 'subject_file_list',
                 gs_cluster_vol, 'sample_file')
    basc.connect(inputspec, 'roi_mask_file',
                 gs_cluster_vol, 'roi_mask_file')
    gs_cluster_vol.inputs.filename = 'group_stability_clusters.nii.gz'

    basc.connect(nis, 'ism_file',
                 igcm, 'indiv_stability_list')
    basc.connect(gsm, 'clusters_G',
                 igcm, 'clusters_G')
    basc.connect(inputspec, 'roi_mask_file',
                 igcm, 'roi_mask_file')

    basc.connect(gsm, 'clusters_G',
                 gs_cluster_vol, 'data_array')

    basc.connect(inputspec, 'subject_file_list',
                 gs_score_vol, 'sample_file')
    basc.connect(inputspec, 'roi_mask_file',
                 gs_score_vol, 'roi_mask_file')
    gs_score_vol.inputs.filename = 'group_stability_scores.nii.gz'

        #return G, clusters_G, clusters_voxel_scores, gsm_file, clusters_G_file, clusters_voxel_scores_file


    basc.connect(gsm, 'cluster_voxel_scores',
                 gs_score_vol, 'data_array')

    basc.connect(gsm, 'gsm_file',
                 outputspec, 'gsm')
    basc.connect(gsm, 'clusters_G_file',
                 outputspec, 'gsclusters')
    basc.connect(gsm, 'cluster_voxel_scores_file',
                 outputspec, 'gsmap')

    basc.connect(gs_cluster_vol, 'img_file',
                 outputspec, 'gsclusters_img')
    basc.connect(gs_score_vol, 'img_file',
                 outputspec, 'gsmap_img')

    basc.connect(igcm, 'individual_cluster_voxel_scores',
                 outputspec, 'ismap_imgs')
    return basc
