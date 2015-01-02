# Import packages
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
# Import CPAC functions
from CPAC.network_centrality import *
from CPAC.network_centrality.core import *


# Function to create the network centrality workflow
def create_resting_state_graphs(allocated_memory = None,
                                wf_name = 'resting_state_graph'):
    '''
    Workflow to calculate degree and eigenvector centrality as well as 
    local functional connectivity density (lfcd) measures for the 
    resting state data.
    
    Parameters
    ----------
    generate_graph : boolean
        when true the workflow plots the adjacency matrix graph 
        and converts the adjacency matrix into compress sparse 
        matrix and stores it in a .mat file. By default its False
    wf_name : string
        name of the workflow
        
    Returns 
    -------
    wf : workflow object
        resting state graph workflow object
          
    Notes
    -----
    
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/network_centrality/resting_state_centrality.py>`_
    
    Workflow Inputs::
    
        inputspec.subject: string (nifti file)
            path to resting state input data for which centrality measure is to be calculated
            
        inputspec.template : string (existing nifti file)
            path to mask/parcellation unit 
        
        inputspec.method_option: string (int)
            0 for degree centrality, 1 for eigenvector centrality, 2 for lFCD
        
        inputspec.threshold: string (float)
            pvalue/sparsity_threshold/threshold value
        
        inputspec.threshold_option: string (int)
            threshold options:  0 for probability p_value, 1 for sparsity threshold, any other for threshold value
           
        centrality_options.weight_options : string (list of boolean)
            list of two booleans for binarize and weighted options respectively
        
        centrality_options.method_options : string (list of boolean)
            list of two booleans for Degree and Eigenvector centrality method options respectively
        
    Workflow Outputs::
    
        outputspec.centrality_outputs : string (list of nifti files)
            path to list of centrality outputs for binarized or/and weighted and
            degree or/and eigen_vector 
        
        outputspec.threshold_matrix : string (numpy file)
            path to file containing thresholded correlation matrix
        
        outputspec.correlation_matrix : string (numpy file)
            path to file containing correlation matrix
        
        outputspec.graph_outputs : string (mat and png files)
            path to matlab compatible sparse adjacency matrix files 
            and adjacency graph images 
    
    Order of commands:
    
    - load the data and template, based on template type (parcellation unit ar mask)
      extract timeseries
    
    - Calculate the correlation matrix for the image data for each voxel in the mask or node
      in the parcellation unit
    
    - Based on threshold option (p_value or sparsity_threshold), calculate the threshold value
    
    - Threshold the correlation matrix
     
    - Based on weight options for edges in the network (binarize or weighted), calculate Degree 
      or Vector Based centrality measures
     
    
    High Level Workflow Graph:
    
    .. image:: ../images/resting_state_centrality.dot.png
       :width: 1000
    
    
    Detailed Workflow Graph:
    
    .. image:: ../images/resting_state_centrality_detailed.dot.png
       :width: 1000
    
    Examples
    --------
    
    >>> import resting_state_centrality as graph
    >>> wflow = graph.create_resting_state_graphs()
    >>> wflow.inputs.centrality_options.method_options=[True, True]
    >>> wflow.inputs.centrality_options.weight_options=[True, True]
    >>> wflow.inputs.inputspec.subject = '/home/work/data/rest_mc_MNI_TR_3mm.nii.gz'
    >>> wflow.inputs.inputspec.template = '/home/work/data/mask_3mm.nii.gz'
    >>> wflow.inputs.inputspec.threshold_option = 1
    >>> wflow.inputs.inputspec.threshold = 0.0744
    >>> wflow.base_dir = 'graph_working_directory'
    >>> wflow.run()
    
    '''
    
    # Instantiate workflow with input name
    wf = pe.Workflow(name = wf_name)
    
    # Instantiate inputspec node
    inputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                       'template',
                                                       'method_option',
                                                       'threshold_option',
                                                       'threshold',
                                                       'weight_options']),
                        name='inputspec')
    
    # Instantiate calculate_centrality main function node
    calculate_centrality = pe.Node(util.Function(input_names = ['datafile',
                                                                'template',
                                                                'method_option',
                                                                'threshold_option',
                                                                'threshold',
                                                                'weight_options',
                                                                'allocated_memory'],
                                                 output_names = ['out_list'],
                                                 function = calc_centrality),
                                   name = 'calculate_centrality')
    
    # Connect inputspec node to main function node
    wf.connect(inputspec, 'subject', 
               calculate_centrality, 'datafile')
    wf.connect(inputspec, 'template', 
               calculate_centrality, 'template')
    wf.connect(inputspec, 'method_option',
               calculate_centrality, 'method_option')
    wf.connect(inputspec, 'threshold_option', 
               calculate_centrality, 'threshold_option')
    wf.connect(inputspec, 'threshold', 
               calculate_centrality, 'threshold')
    wf.connect(inputspec,'weight_options',
               calculate_centrality,'weight_options')
    
    # Specify allocated memory from workflow input to function node
    calculate_centrality.inputs.allocated_memory = allocated_memory
    
    # Instantiate outputspec node
    outputspec = pe.Node(util.IdentityInterface(fields=['centrality_outputs',
                                                        'threshold_matrix',
                                                        'correlation_matrix',
                                                        'graph_outputs']),
                         name = 'outputspec')
    
    # Connect function node output list to outputspec node
    wf.connect(calculate_centrality, 'out_list',
               outputspec, 'centrality_outputs')
    
    # Return the connected workflow
    return wf


# Function to load in nifti files and extract info for centrality calculation
def load(datafile, template=None):
    '''
    Method to read data from datafile and mask/parcellation unit
    and store the mask data, timeseries, affine matrix, mask type
    and scans. The output of this method is used by all other nodes.
    
    Note that this function also will internally compute it's own 
    brain mask by getting all voxels with non-zero variance in the
    timeseries.
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to subject data file
    template : string (nifti file) or None (default: None)
        path to mask/parcellation unit
        if none, then will be mask with all 1s
        
    Returns
    -------
    timeseries_data : ndarray
        Masked timeseries of the input data. 
    affine : ndarray
        Affine matrix of the input data
    final_mask : ndarray
        Mask/parcellation unit matrix
    template_type : string 
        0 for mask, 1 for parcellation unit 
    scans : string (int)
        total no of scans in the input data
        
    Raises
    ------
    Exception
    '''

    import os
    import nibabel as nib
    import numpy as np

    
    try:
        if isinstance(datafile, list):
            img = nib.load(datafile[0])
        else:
            img = nib.load(datafile) 
        
        data    = img.get_data().astype(np.float32)
        aff     = img.get_affine()    
        scans   = data.shape[3]
        
        datmask     = data.var(axis=3).astype('bool')
        if template is None:
            mask    = np.ones((data.shape[:3]))
        else:
            mask    = nib.load(template).get_data().astype(np.float32)
        
    except:
        print "Error in loading images for graphs"
        raise
    
    
    if mask.shape != data.shape[:3]:
        raise Exception("Invalid Shape Error. mask and data file have"\
                        "different shape please check the voxel size of the two files")
    
    #check for parcellation
    nodes = np.unique(mask).tolist()
        
    #extract timeseries
    if len(nodes)>2:
        nodes.sort()
        print "sorted nodes", nodes
        
        flag=1
        for n in nodes:
            if n > 0:
                node_array = data[(mask == n) & datmask]
                avg = np.mean(node_array, axis =0)
                if flag:
                    timeseries = avg
                    flag=0
                else:
                    timeseries = np.vstack((timeseries, avg))
            final_mask  = datmask
        #template_type is 1 for parcellation
        template_type = 1
    else:
        #template_type is 0 for mask
        template_type = 0
        mask = mask.astype('bool')
        final_mask = mask & datmask
        timeseries = data[final_mask]
    
    return timeseries, aff, final_mask, template_type, scans


# Function to calculate centrality using a correlation threshold 
def get_centrality_by_rvalue(ts_normd, 
                             template, 
                             method_option, 
                             weight_options, 
                             r_value, 
                             block_size):
    '''
    Method to calculate degree/eigenvector centrality and lFCD
    via correlation (r-value) threshold
    
    Parameters
    ----------
    ts_normd : ndarray (float)
        timeseries of shape (ntpts x nvoxs) that is normalized; i.e. the data 
        is demeaned and divided by its L2-norm
    template : ndarray
        three dimensional array with non-zero elements corresponding to the
        indices at which the lFCD metric is analyzed
    method_option : integer
        0 - degree centrality calculation, 
        1 - eigenvector centrality calculation, 
        2 - lFCD calculation
    weight_options : list (boolean)
        weight_options[0] - True or False to perform binary counting
        weight_options[1] - True or False to perform weighted counting
    threshold : a float
        threshold (as correlation r) value
    block_size : an integer
        the number of rows (voxels) to compute timeseries correlation over
        at any one time
    
    Returns
    -------
    out_list : list (string, ndarray)
        list of (string,ndarray) elements corresponding to:
        string - the name of the metric
        ndarray - the array of values to be mapped for that metric
    '''
    
    # Import packages
    from CPAC.network_centrality.utils import cluster_data
    
    # Init variables
    out_list = []
    nvoxs = ts_normd.shape[1]
    # ntpts = timeseries.shape[0]
    calc_degree = False
    calc_eigen = False
    calc_lfcd = False
    # Select which method we're going to perform
    if method_option == 0:
        calc_degree = True
    elif method_option == 1:
        calc_eigen = True
    elif method_option == 2:
        calc_lfcd = True
    # Weighting
    out_binarize = weight_options[0]
    out_weighted = weight_options[1]
    
    # Init degree centrality outputs
    if calc_degree:
        # If binary weighting, init output map
        if out_binarize:
            degree_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('degree_centrality_binarize', degree_binarize))
        # If connection weighting, init output map
        if out_weighted:
            degree_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('degree_centrality_weighted', degree_weighted))
    # Init eigenvector centrality outputs
    if calc_eigen:
        r_matrix = np.zeros((nvoxs,nvoxs), dtype=ts_normd.dtype)
        # If binary weighting, init output map
        if out_binarize:
            eigen_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_binarize', eigen_binarize))
        # If connection weighting, init output map
        if out_weighted:
            eigen_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
    # Init lFCD outputs
    if calc_lfcd:
        # If binary weighting, init output map
        if out_binarize:
            lfcd_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('lfcd_binarize', lfcd_binarize))
        # If connection weighting, init output map
        if out_weighted:
            lfcd_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('lfcd_weighted', lfcd_weighted))
    
    # Prepare to loop through and calculate correlation matrix
    n = 0
    m = block_size
    block_no = 1
    
    # Run as long as our last row index is <= nvoxs
    while m <= nvoxs:
        # First, compute block of correlation matrix
        print 'running block %d: rows %d thru %d' % (block_no, n, m)
        rmat_block = np.dot(ts_normd[:,n:m].T, ts_normd)
        
        # Degree centrality calculation
        if calc_degree:
            if weight_options[0]:
                degree_centrality(rmat_block, r_value, method='binarize', 
                                  out=degree_binarize[n:m])
            if weight_options[1]:
                degree_centrality(rmat_block, r_value, method='weighted', 
                                  out=degree_weighted[n:m])
        
        # Eigenvector centrality - append global corr. matrix
        if calc_eigen:
            r_matrix[n:m] = rmat_block
        
        # lFCD - perform lFCD algorithm
        if calc_lfcd:
            xyz_a = np.argwhere(template)
            krange = rmat_block.shape[0]
            print '...iterating through seeds in block - lfcd'
            for k in range (0,krange):
                corr_seed = rmat_block[k,:]
                labels = cluster_data(corr_seed,r_value,xyz_a)
                seed_label = labels[n+k]
                if out_binarize:
                    if seed_label > 0:
                        lfcd = np.sum(labels==seed_label)
                    else:
                        lfcd = 1
                    lfcd_binarize[n+k] = lfcd
                if out_weighted:
                    if seed_label > 0:
                        lfcd = np.sum(corr_seed*(labels==seed_label))
                    else:
                        lfcd = 1
                    lfcd_weighted[n+k] = lfcd
        
        # Delete block of corr matrix and increment indices
        del rmat_block
        
        # Move next block start point up to last block finish point
        n = m
        # If we finished at nvoxs last time, break the loop
        if n == nvoxs:
            break
        # Else, if our next block runs over nvoxs, limit it to nvoxs
        elif (m+block_size) > nvoxs:
            m = nvoxs
        # Else, just increment end of next block by block_size
        else:
            m += block_size
        # Increment block number
        block_no += 1
    
    # Correct for self-correlation in degree centrality
    if calc_degree:
        if out_binarize:
            idx = np.where(degree_binarize)
            degree_binarize[idx] = degree_binarize[idx]-1
        if out_weighted:
            idx = np.where(degree_weighted)
            degree_weighted[idx] = degree_weighted[idx]-1
    
    # Perform eigenvector measures
    try:
        if calc_eigen:
            if out_binarize:
                print '...calculating binarize eigenvector'
                eigen_binarize[:] = eigenvector_centrality(r_matrix, 
                                                           r_value, 
                                                           method='binarize').squeeze()
            if out_weighted:
                print '...calculating weighted eigenvector'
                eigen_weighted[:] = eigenvector_centrality(r_matrix, 
                                                           r_value, 
                                                           method='weighted').squeeze()
            del r_matrix
        
    except Exception:
        print 'Error in calcuating eigen vector centrality'
        raise
    
    # Return list of outputs
    return out_list


# Function to calculate centrality with a sparsity threhold
def get_centrality_by_sparsity(ts_normd, 
                               method_option, 
                               weight_options, 
                               threshold, 
                               block_size):
    '''
    Method to calculate degree/eigenvector centrality via sparsity threshold
    
    Parameters
    ----------
    ts_normd : ndarray
        timeseries of shape (ntpts x nvoxs) that is normalized; i.e. the data 
        is demeaned and divided by its L2-norm
    method_option : integer
        0 - degree centrality calculation, 
        1 - eigenvector centrality calculation, 
        2 - lFCD calculation
    weight_options : list (boolean)
        weight_options[0] - True or False to perform binary counting
        weight_options[1] - True or False to perform weighted counting
    threshold : a float
        sparsity_threshold value
    block_size : an integer
        the number of rows (voxels) to compute timeseries correlation over
        at any one time
    
    Returns
    -------
    out_list : list (string, ndarray)
        list of (string,ndarray) elements corresponding to:
        string - the name of the metric
        ndarray - the array of values to be mapped for that metric
    '''
    
    # Import packages
    import scipy as sp

    # Init variables
    out_list = []
    nvoxs = ts_normd.shape[1]
    # ntpts = timeseries.shape[0]
    calc_degree = False
    calc_eigen = False
    
    # Select which method we're going to perform
    if method_option == 0:
        calc_degree = True
    elif method_option == 1:
        calc_eigen = True
    
    # Weighting
    out_binarize = weight_options[0]
    out_weighted = weight_options[1]
    
    # Init degree centrality outputs
    if calc_degree:
        # If binary weighting, init output map
        if out_binarize:
            degree_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('degree_centrality_binarize', degree_binarize))
        # If connection weighting, init output map
        if out_weighted:
            degree_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('degree_centrality_weighted', degree_weighted))
    
    # Init eigenvector centrality outputs
    if calc_eigen:
        r_matrix = np.zeros((nvoxs,nvoxs), dtype = ts_normd.dtype)
        # If binary weighting, init output map
        if out_binarize:
            eigen_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_binarize', eigen_binarize))
        # If connection weighting, init output map
        if out_weighted:
            eigen_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
    
    # Get the number of connections to keep
    sparse_num = np.round((nvoxs**2-nvoxs)*threshold/2.0)
    
    # Prepare to loop through and calculate correlation matrix
    n = 0
    m = block_size
    block_no = 1
    r_value = -1
    # Init wij list
    z = np.array([])
    wij_global = np.rec.fromarrays([z.astype(ts_normd.dtype), 
                                    z.astype('int32'), 
                                    z.astype('int32')])
    # Form the initial blockwise mask (to only grab upper triangle of data)
    block_triu = np.triu(np.ones((block_size,block_size)), k=1).astype('bool')
    block_rect = np.ones((block_size,nvoxs-block_size), dtype='bool')
    block_mask = np.concatenate((block_triu,block_rect), axis=1)
    # Delete matrices to save memory
    del block_triu, block_rect
    # Calculate correlations step - prune connections for degree
    while n <= nvoxs:
        # First, compute block of correlation matrix
        print 'running block %d: rows %d thru %d' % (block_no, n, m-1)
        # Calculate wij over entire matrix by block
        # Do this for both deg and eig, more efficient way to compute r_value
        rmat_block = np.dot(ts_normd[:,n:m].T, 
                            ts_normd[:,n:])
        # Shrink block_mask down on every iteration after the first block
        if n > 0:
            block_mask = block_mask[:,:-block_size]
        # Get elements as an array
        rmat_block = rmat_block[block_mask]
        thr_idx = np.where(rmat_block >= r_value)   
        rmat_block = rmat_block[thr_idx]
        print 'number of passing correlations is %d' % len(rmat_block)
        # Add global offset
        idx = np.where(block_mask)
        i = idx[0][thr_idx].astype('int32') + n
        j = idx[1][thr_idx].astype('int32') + n
        # Free some memory
        del idx, thr_idx
        w_global = np.concatenate([wij_global.f0,rmat_block])
        i_global = np.concatenate([wij_global.f1,i])
        j_global = np.concatenate([wij_global.f2,j])
        # Free some memory
        del i,j,rmat_block
        # Grab indices and weights that pass and combine into list
        wij_global = np.rec.fromarrays([w_global,i_global,j_global])
        # Free some memory
        del w_global, i_global, j_global
        # Pass these into the global set and sort (ascending) by correlation
        print 'sorting list...'
        wij_global.sort()
        # And trim list if it's greater than the number of connections we want
        if len(wij_global) > sparse_num:
            wij_global = wij_global[-sparse_num:]
        r_value = wij_global[0][0]
        
        # If we're doing eigen, store block into full matrix
        if calc_eigen:
            r_matrix[n:m] = np.dot(ts_normd[:,n:m].T, ts_normd)
        
        # Move next block start point up to last block finish point
        n = m
        # If we finished at nvoxs last time, break the loop
        if n == nvoxs:
            break
        # Else, if our next block runs over nvoxs, limit it to nvoxs
        elif (m+block_size) > nvoxs:
            m = nvoxs
        # Else, just increment end of next block by block_size
        else:
            m += block_size
        # Increment block number
        block_no += 1
    
    # Calculate centrality step
    # Degree - use ijw list to create a sparse matrix
    if calc_degree:
        # Create sparse (symmetric) matrix of all correlations that survived
        print 'creating sparse matrix'
        # Extract the weights and indices from the global list
        w = wij_global.f0
        i = wij_global.f1
        j = wij_global.f2
        del wij_global
        # And compute degree centrality on sparse matrix
        if weight_options[0]:
            # Create the sparse correlation matrix (upper triangle) from wij's
            Rsp = sp.sparse.coo_matrix((np.ones(len(w)),(i,j)),
                                        shape=(nvoxs,nvoxs))
            # Make it symmetric
            Rsp = Rsp + Rsp.T
            Rcsr = Rsp.tocsr()
            degree_binarize[:] = np.array(Rcsr.sum(axis=0))
        if weight_options[1]:
            # Create the sparse correlation matrix (upper triangle) from wij's
            Rsp = sp.sparse.coo_matrix((w,(i,j)), shape=(nvoxs,nvoxs))
            del w, i, j
            # Make it symmetric
            Rsp = Rsp + Rsp.T
            Rcsr = Rsp.tocsr()
            degree_weighted[:] = np.array(Rcsr.sum(axis=0))
        del Rsp
    
    # Eigenvector - compute the r value from entire matrix
    if calc_eigen:
        del wij_global
        # Finally compute centrality using full matrix and r_value
        if out_binarize:
            print '...calculating binarize eigenvector'
            eigen_binarize[:] = eigenvector_centrality(r_matrix, r_value, 
                                                       method='binarize').squeeze()
        if out_weighted:
            print '...calculating weighted eigenvector'
            eigen_weighted[:] = eigenvector_centrality(r_matrix, r_value, 
                                                       method='weighted').squeeze()
        del r_matrix
    
    # Return list of outputs
    return out_list


# Function to calculated a quick centrality measure
def get_centrality_fast(timeseries,
                        method_options):
    '''
    Method to calculate degree and eigen vector centrality. 
    Relative to `get_centrality_opt`, it runs fast by not directly computing 
    the correlation matrix. As a consequence, there are several differences/
    limitations from the standard approach:
    
    1. Cannot specify a correlation threshold
    2. As a consequence, the weighted dense matrix centrality is computed
    3. No memory limit is specified since it is assumed to be ok
    
    Note that the approach doesn't directly calculate the complete correlation
    matrix.
    
    Parameters
    ----------
    timeseries_data : numpy array
        timeseries of the input subject
    method_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    
    Returns
    -------
    out_list : string (list of tuples)
        list of tuple containing output name to be used to store nifti image
        for centrality and centrality matrix. this will only be weighted since
        the fast approaches are limited to this type of output.
    
    Raises
    ------
    Exception
    '''
    
    # Import packages
    from CPAC.network_centrality import fast_degree_centrality,\
                                        fast_eigenvector_centrality
    from CPAC.cwas.subdist import norm_cols
    
    try:
        out_list    = []
        nvoxs       = timeseries.shape[0]
        ntpts       = timeseries.shape[1]
        
        # It's assumed that there is enough memory
        # So a block size isn't set here
        
        calc_degree  = method_options[0]
        calc_eigen   = method_options[1]
        
        print "Normalize Time-series"
        timeseries = norm_cols(timeseries.T)
        
        print "Computing centrality across %i voxels" % nvoxs
        
        if calc_degree:
            print "...calculating degree"
            degree_weighted = fast_degree_centrality(timeseries)
            out_list.append(('degree_centrality_weighted', degree_weighted))
        
        if calc_eigen:
            print "...calculating eigen"
            eigen_weighted = fast_eigenvector_centrality(timeseries)
            out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
        
        return out_list   
    
    except Exception: 
        print "Error in calcuating centrality"
        raise


# Main centrality function utilized by the centrality workflow
def calc_centrality(datafile,
                    template,
                    method_option,
                    threshold_option,
                    threshold,
                    weight_options,
                    allocated_memory):
    '''
    Method to calculate centrality and map them to a nifti file
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to subject data file
    template : string (nifti file)
        path to mask/parcellation unit
    method_option : integer
        0 - degree centrality calculation, 1 - eigenvector centrality calculation, 2 - lFCD calculation
    threshold_option : an integer
        0 for probability p_value, 1 for sparsity threshold, 
        2 for actual threshold value, and 3 for no threshold and fast approach
    threshold : a float
        pvalue/sparsity_threshold/threshold value
    weight_options : list (boolean)
        list of booleans, where, weight_options[0] corresponds to binary counting 
        and weight_options[1] corresponds to weighted counting (e.g. [True,False]) 
    allocated_memory : string
        amount of memory allocated to degree centrality
    
    Returns
    -------
    out_list : list
        list containing out mapped centrality images
    '''
    
    # Import packages
    from CPAC.network_centrality import load,\
                                        get_centrality_by_rvalue,\
                                        get_centrality_by_sparsity,\
                                        get_centrality_fast,\
                                        map_centrality_matrix,\
                                        calc_blocksize,\
                                        convert_pvalue_to_r
    from CPAC.cwas.subdist import norm_cols
    
    # Check for input errors
    if weight_options.count(True) == 0:
        raise Exception("Invalid values in weight options" \
                        "At least one True value is required")
    # If it's sparsity thresholding, check for (0,1]
    if threshold_option == 1:
        if threshold <= 0 or threshold > 1:
            raise Exception('Threshold value must be a positive number'\
                            'greater than 0 and less than or equal to 1.'\
                            '\nCurrently it is set at %d' % threshold)
    if method_option == 2 and threshold_option != 2:
        raise Exception('lFCD must use correlation-type thresholding.'\
                         'Check the pipline configuration has this setting')
    import time
    start = time.clock()
    
    # Init variables
    out_list = []
    ts, aff, mask, t_type, scans = load(datafile, template)
    
    # If we're doing eigenvectory centrality, need entire correlation matrix
    if method_option == 0 and threshold_option == 1:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    sparsity_thresh=threshold)
    elif method_option == 1:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    include_full_matrix=True)
    # Otherwise, compute blocksize with regards to available memory
    else:
        block_size = calc_blocksize(ts, memory_allocated=allocated_memory,
                                    include_full_matrix=False)
    # Normalize the timeseries for easy dot-product correlation calc.
    ts_normd = norm_cols(ts.T)
    
    # P-value threshold centrality
    if threshold_option == 0:
        r_value = convert_pvalue_to_r(scans, threshold)
        centrality_matrix = get_centrality_by_rvalue(ts_normd, 
                                                     mask, 
                                                     method_option, 
                                                     weight_options, 
                                                     r_value, 
                                                     block_size)
    # Sparsity threshold
    elif threshold_option == 1:
        centrality_matrix = get_centrality_by_sparsity(ts_normd, 
                                                       method_option, 
                                                       weight_options, 
                                                       threshold, 
                                                       block_size)
    # R-value threshold centrality
    elif threshold_option == 2:
        centrality_matrix = get_centrality_by_rvalue(ts_normd, 
                                                     mask, 
                                                     method_option, 
                                                     weight_options, 
                                                     threshold, 
                                                     block_size)
    # For fast approach (no thresholding)
    elif threshold_option == 3:
        centrality_matrix = get_centrality_fast(ts, method_option)
    # Otherwise, incorrect input for threshold_option
    else:
        raise Exception('Option must be between 0-3 and not %s, check your '\
                        'pipeline config file' % str(threshold_option))
    
    # Print timing info
    print 'Timing:', time.clock() - start
 
    # Map the arrays back to images
    for mat in centrality_matrix:
        centrality_image = map_centrality_matrix(mat, aff, mask, t_type)
        out_list.append(centrality_image)
    
    # Finally return
    return out_list

