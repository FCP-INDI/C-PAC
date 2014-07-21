import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.network_centrality import *
from CPAC.network_centrality.core import *

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

    return wf

# # Function to combine outputs for different centrality methods
# def create_merging_graph(wf_name = 'merging_graph'):
#     # Instantiate workflow
#     wf = pe.Node(name = wf_name)
#     inputspec = pe.Node(util.IdentityInterface(fields=['merge_list']),
#                         name = 'inputspec')
#     merge_node = pe.Node(util.Function(input_names=['deg_list',
#                                                     'eig_list',
#                                                     'lfcd_list'],
#                                           output_names = ['out_list'],
#                                           function = merge_lists),
#                             name = 'merge_node')
#     wf.connect(inputspec, 'merge_list',
#                merge_node, 'in_list')
#     outputspec = pe.Node(util.IdentityInterface(fields=['merged_list']),
#                          name = 'outputspec')
#     wf.connect(merge_node,'out_list',
#                outputspec,'merged_list')
#     return wf

    
def load(datafile, template=None):
    
    """
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
    timeseries_data: ndarray
        Masked timeseries of the input data. 
    affine: ndarray
        Affine matrix of the input data
    mask_data: ndarray
        Mask/parcellation unit matrix
    template_type: string 
        0 for mask, 1 for parcellation unit 
    scans: string (int)
        total no of scans in the input data
        
    Raises
    ------
    Exception
    """

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


def get_centrality_by_sparsity(timeseries, 
                   method_option,
                   weight_options,
                   threshold,
                   memory_allocated):
    
    """
    Method to calculate degree and eigen vector centrality
    
    Parameters
    ----------
    timeseries : numpy array
        timeseries of the input subject
    method_options : string (list of boolean)
        list of two booleans for degree and eigen options respectively
    weight_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    threshold : float
        sparsity threshold for the correlation values
    memory_allocated : a string
        amount of memory allocated to degree centrality
    
    Returns
    -------
    out_list : string (list of tuples)
        list of tuple containing output name to be used to store nifti image
        for centrality and centrality matrix
    
    Raises
    ------
    Exception
    """
    
    import os
    import numpy as np
    from CPAC.network_centrality import calc_blocksize,\
                                        convert_sparsity_to_r,\
                                        degree_centrality,\
                                        eigenvector_centrality
    from CPAC.cwas.subdist import norm_cols
    
    out_list=[]
    
    try:
        # Calculate the block size (i.e., number of voxels) to compute part of the
        # connectivity matrix at once.
        # 
        # We still use a block size to calculate the whole correlation matrix
        # because of issues in numpy that lead to extra memory usage when
        # computing the dot product.
        # See https://cmi.hackpad.com/Numpy-Memory-Issues-BlV9Pg5nRDM.
        block_size  = calc_blocksize(timeseries, memory_allocated, include_full_matrix=True)
        
        nvoxs = timeseries.shape[0]
        ntpts = timeseries.shape[1]
        
        calc_degree = False         # init degree measure flag to false
        calc_eigen = False          # init eigen measure flag to false
        calc_lfcd= False            # init lFCD measure flag to false
        
        # Select which method we're going to perform
        if method_option == 0:
            calc_degree = True
        elif method_option == 1:
            calc_eigen = True
        elif method_option == 2:
            calc_lfcd = True
        
        # Set weighting parameters
        out_binarize = weight_options[0]
        out_weighted = weight_options[1]
        
        corr_matrix = np.zeros((nvoxs, nvoxs), dtype = timeseries.dtype)
        
        
        print "Normalize TimeSeries"
        timeseries  = norm_cols(timeseries.T)
        
        
        print "Computing centrality across %i voxels" % nvoxs
        j = 0
        i = block_size
        while i <= timeseries.shape[1]:
            print "running block ->", i,j
            
            print "...correlating"
            np.dot(timeseries[:,j:i].T, timeseries, out=corr_matrix[j:i])
            
            j = i
            if i == nvoxs:
                break
            elif (i+block_size) > nvoxs:
                i = nvoxs
            else:
                i += block_size
        
        
        print "Calculating threshold"
        r_value = convert_sparsity_to_r(corr_matrix, threshold, full_matrix = True)
        print "r_value ->", r_value
        
        
        if calc_degree:
            if out_binarize:
                print "...calculating binarize degree"
                degree_binarize = degree_centrality(corr_matrix, r_value, method="binarize")
                out_list.append(('degree_centrality_binarize', degree_binarize))
            if out_weighted:
                print "...calculating weighted degree"
                degree_weighted = degree_centrality(corr_matrix, r_value, method="weighted")
                out_list.append(('degree_centrality_weighted', degree_weighted))
        
        
        if calc_eigen:
            if out_binarize:
                print "...calculating binarize eigenvector"
                eigen_binarize = eigenvector_centrality(corr_matrix, r_value, method="binarize")
                out_list.append(('eigenvector_centrality_binarize', eigen_binarize))
            if out_weighted:
                print "...calculating weighted eigenvector"
                eigen_weighted = eigenvector_centrality(corr_matrix, r_value, method="weighted")
                out_list.append(('eigenvector_centrality_weighted', eigen_weighted))            
            
    except Exception:
        print "Error while calculating centrality"
        raise
    
    return out_list


def calc_via_rvalue(ts_normd, 
                    template, 
                    method_option, 
                    weight_options, 
                    r_value, 
                    block_size):
    
    # Import packages
    #from CPAC.network_centrality import convert_pvalue_to_r,\
    #                                    degree_centrality
    
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
            out_list.append(('degree_centrality_binarize', degree_binarize)
        # If connection weighting, init output map
        if out_weighted:
            degree_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('degree_centrality_weighted', degree_weighted))
    # Init eigenvector centrality outputs
    if calc_eigen:
        r_matrix = np.zeros((nvoxs,nvoxs), dtyp=ts_normd.dtype))
        # If binary weighting, init output map
        if out_binarize:
            eigen_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_binarize', eigen_binarize)
        # If connection weighting, init output map
        if out_weighted:
            eigen_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
    # Init lFCD outputs
    if calc_lfcd:
        # If binary weighting, init output map
        if out_binarize:
            lfcd_binarize = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('lfcd_binarize', lfcd_binarize)
        # If connection weighting, init output map
        if out_weighted:
            lfcd_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('lfcd_weighted', lfcd_weighted))
    
    # Prepare to loop through and calculate correlation matrix
    n = 0
    m = (n+1)*block_size
    exit_flg = 0
    
    while not exit_flg:
        # If we've gone past the number of voxels, limit m and set exit flag
        if m > nvoxs:
            m = nvoxs
            exit_flg = 1
        # First, compute block of correlation matrix
        print 'running block %d: rows %d thru %d' % (n+1, n*block_size, m)
        rmat_block = np.dot(ts_normd[:,n*block_size:m].T, ts_normd)
        
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
        n += 1
        m = (n+1)*block_size
        
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
        except Exception:
            print 'Error in calcuating eigen vector centrality'
            raise
        
    return out_list


def calc_via_sparsity(ts_normd, 
                      template, 
                      method_option, 
                      weight_options, 
                      threshold, 
                      block_size):
    
    # Import packages
    #from CPAC.network_centrality import convert_pvalue_to_r,\
    #                                    degree_centrality
    
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
            out_list.append(('degree_centrality_binarize', degree_binarize)
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
            out_list.append(('eigenvector_centrality_binarize', eigen_binarize)
        # If connection weighting, init output map
        if out_weighted:
            eigen_weighted = np.zeros(nvoxs, dtype=ts_normd.dtype)
            out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
    
    # Get the number of connections to keep
    no_conn = int((nvoxs**2-nvoxs)*threshold/2)
    
    # Prepare to loop through and calculate correlation matrix
    n = 0
    m = (n+1)*block_size
    exit_flg = 0
    r_value = -1
    ijw_global = []
    
    while not exit_flg:
        # If we've gone past the number of voxels, limit m and set exit flag
        if m > nvoxs:
            m = nvoxs
            exit_flg = 1
        # First, compute block of correlation matrix
        print 'running block %d: rows %d thru %d' % (n+1, n*block_size, m)
        rmat_block = np.dot(ts_normd[:,n*block_size:m].T, 
                            ts_normd[:,n*block_size:])
        # Get upper-right chunk of block shifted by 1 (off main diagonal)
        rmat_triu = np.triu(rmat_block,1)
        # Get elements as an array
        idx = np.where(rmat_triu)
        rarr = rmat_triu[idx]
        global_row_idx = idx[0] + n*block_size
        global_col_idx = idx[1] + n*block_size
        
        # Threshold elements
        thr_idx = np.where(rarr >= threshold)
        print 'number of passing correlations id %d' % len(thr[0])
        
        # Grab indices and weights that pass and combine into list
        ijw = (global_row_idx[thr_idx], global_col_idx[thr_idx], rarr[thr_idx])
        ijw_block = [(ijw[0][i], ijw[1][i], ijw[2][i]) \
                     for i in range(len(thr[0]))]
        
        # Pass these into the global set and sort (ascending) by correlation
        ijw_global.extend(ijw_block)
        print 'sorting list...'
        ijw_global = sorted(ijw_global, key=lambda x: x[2])
        # And trim list if it's greater than the number of connections we want
        if len(ijw_global) > no_conns:
            ijw_global = ijw_global[-no_conns:]
        r_value = ijw_global[0][2]
        
        # Increment indices
        n += 1
        m = (n+1)*block_size
        if calc_eigen:
            r_matrix[n:m] = np.dot(ts_normd[:,n*block_size:m].T, ts_normd)
    
    if calc_degree:
        # Create sparse (symmetric) matrix of all correlations that survived
        print 'creating sparse matrix'
        i = []
        j = []
        w = []
        for ii,jj,ww in ijw_global
            i.append(ii)
            j.append(jj)
            w.append(ww)
        Rsp = sp.sparse.coo_matrix((w,(i,j)), shape=(nvoxs,nvoxs))
        Rsp = Rsp + Rsp.T
        
        # And compute degree centrality on sparse matrix
        if weight_options[0]:
            degree_centrality(Rsp, r_value, method='binarize', 
                              out=degree_binarize[n:m])
        if weight_options[1]:
            degree_centrality(Rsp, r_value, method='weighted', 
                              out=degree_weighted[n:m])
    if calc_eigen:
        if out_binarize:
            print '...calculating binarize eigenvector'
            eigen_binarize = eigenvector_centrality(r_matrix, r_value, 
                                                    method='binarize')
        if out_weighted:
            print '...calculating weighted eigenvector'
            eigen_weighted = eigenvector_centrality(r_matrix, r_value, 
                                                    method='weighted')
        
        
    return out_list
        

def get_centrality_by_thresh(timeseries,
                             template,
                             method_option,
                             weight_options,
                             threshold,
                             r_value,
                             memory_allocated):
    """
    Method to calculate degree and eigen vector centrality. 
    This method takes into consideration the amount of memory
    allocated by the user to calculate degree centrality.
    
    Parameters
    ----------
    timeseries_data : numpy array
        timeseries of the input subject
    template : numpy array
        Mask/ROI template for timeseries of subject
    method_option : integer
        0 - degree centrality calculation, 1 - eigenvector centrality calculation, 2 - lFCD calculation
    weight_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    threshold : float
        p-value threshold for the correlation values (ignored if the r_value option is specified)
    r_value : float
        threshold value in terms of the correlation (this will override the threshold option)
    memory_allocated : a string
        amount of memory allocated to degree centrality
        
    Returns
    -------
    out_list : string (list of tuples)
        list of tuple containing output name to be used to store nifti image
        for centrality and centrality matrix 
    
    Raises
    ------
    Exception
    """
    
    
    import numpy as np
    import os
    from CPAC.network_centrality import calc_blocksize,\
                                        cluster_data,\
                                        convert_pvalue_to_r,\
                                        degree_centrality,\
                                        eigenvector_centrality
    from CPAC.cwas.subdist import norm_cols
    
    try:                         
        # Init variables for use
        out_list = []
        nvoxs = timeseries.shape[0]
        ntpts = timeseries.shape[1]
        
        r_matrix = None             # init correlation matrix
        calc_degree = False         # init degree measure flag to false
        calc_eigen = False          # init eigen measure flag to false
        calc_lfcd= False            # init lFCD measure flag to false
        
        # Select which method we're going to perform
        if method_option == 0:
            calc_degree = True
        elif method_option == 1:
            calc_eigen = True
        elif method_option == 2:
            calc_lfcd = True
        
        # Set weighting parameters
        out_binarize = weight_options[0]
        out_weighted = weight_options[1]
        
        # Calculate the block size (i.e., number of voxels) to compute part of the
        # connectivity matrix at once.
        if calc_eigen:
            # We still use a block size to calculate the whole correlation matrix
            # because of issues in numpy that lead to extra memory usage when
            # computing the dot product.
            # See https://cmi.hackpad.com/Numpy-Memory-Issues-BlV9Pg5nRDM.
            block_size = calc_blocksize(timeseries, memory_allocated, include_full_matrix=True)
        else:
            block_size = calc_blocksize(timeseries, memory_allocated)
        
        if r_value == None:
            print "Calculating threshold"
            r_value = convert_pvalue_to_r(ntpts, threshold)
            print "...%s -> %s" % (threshold, r_value)
        
        print "Setup Intermediates/Outputs"
        # Degree matrix init
        if calc_degree:
            print "...degree"
            if out_binarize:
                degree_binarize = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('degree_centrality_binarize', degree_binarize))
            if out_weighted:
                degree_weighted = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('degree_centrality_weighted', degree_weighted))
        # Eigen matrix init
        if calc_eigen:
            print "...eigen"
            r_matrix = np.zeros((nvoxs, nvoxs), dtype=timeseries.dtype)
            if out_binarize:
                eigen_binarize = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('eigenvector_centrality_binarize', eigen_binarize))
            if out_weighted:
                eigen_weighted = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
        # lFCD matrix init
        if calc_lfcd:
            print "...degree"
            if out_binarize:
                lfcd_binarize = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('lFCD_binarize', lfcd_binarize))
            if out_weighted:
                lfcd_weighted = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('lFCD_weighted', lfcd_weighted))
        
        # Normalize the timeseries columns for simple correlation calc via dot product later..
        print "Normalize TimeSeries"
        timeseries = norm_cols(timeseries.T)
        
        # Init blocking indices for correlation matrix calculation
        print "Computing centrality across %i voxels" % nvoxs
        i = block_size
        j = 0
        # Calculate correlation matrix in blocks while loop
        while i <= nvoxs:
            print "running block ->", i, j
           
            try:
                print "...correlating"
                corr_matrix = np.dot(timeseries[:,j:i].T, timeseries)
            except:
                raise Exception("Error in calcuating block wise correlation for the block %i,%i"%(j,i))
                      
            if calc_eigen:
                print "...storing correlation matrix"
                r_matrix[j:i] = corr_matrix
            
            if calc_degree:
                if out_binarize:
                    print "...calculating binarize degree"
                    degree_centrality(corr_matrix, r_value, method="binarize", out=degree_binarize[j:i])
                if out_weighted:
                    print "...calculating weighted degree"
                    degree_centrality(corr_matrix, r_value, method="weighted", out=degree_weighted[j:i])
            
            if calc_lfcd:
                xyz_a = np.argwhere(template)
                krange = corr_matrix.shape[0]
                print "...iterating through seeds in block - lfcd"
                for k in range (0,krange):
                    corr_seed = corr_matrix[k,:]
                    labels = cluster_data(corr_seed,r_value,xyz_a)
                    seed_label = labels[j+k]
                    if out_binarize:
                        if seed_label > 0:
                            lfcd = np.sum(labels==seed_label)
                        else:
                            lfcd = 1
                        lfcd_binarize[j+k] = lfcd
                    if out_weighted:
                        if seed_label > 0:
                            lfcd = np.sum(corr_seed*(labels==seed_label))
                        else:
                            lfcd = 1
                        lfcd_weighted[j+k] = lfcd
                            
            print "...removing temporary correlation matrix"
            del corr_matrix
           
            j = i
            if i == nvoxs:
                break
            elif (i+block_size) > nvoxs:
                i = nvoxs
            else:
                i += block_size
        
        # In case there are any zeros in lfcd matrix, set them to 1
        if calc_lfcd:
            if out_binarize:
                lfcd_binarize[np.argwhere(lfcd_binarize == 0)] = 1
            if out_weighted:
                lfcd_weighted[np.argwhere(lfcd_weighted == 0)] = 1
        
        # Perform eigenvector measures if necessary
        try:
            if calc_eigen:
                if out_binarize:
                    print "...calculating binarize eigenvector"
                    eigen_binarize[:] = eigenvector_centrality(r_matrix, r_value, method="binarize").squeeze()
                if out_weighted:
                    print "...calculating weighted eigenvector"
                    eigen_weighted[:] = eigenvector_centrality(r_matrix, r_value, method="weighted").squeeze()
        except Exception:
            print "Error in calcuating eigen vector centrality"
            raise
        
        if calc_degree:
            print "...removing effect of auto-correlation on degree"
            if out_binarize:
                degree_binarize[degree_binarize!=0] = degree_binarize[degree_binarize!=0] - 1
            if out_weighted:
                degree_weighted[degree_weighted!=0] = degree_weighted[degree_weighted!=0] - 1
        
        return out_list
    
    except Exception: 
        print "Error in calcuating Centrality"
        raise


def get_centrality_fast(timeseries,
                        method_options):
    """
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
    """
    
    
    import numpy as np
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
                                        get_centrality_by_thresh,\
                                        get_centrality_by_sparsity,\
                                        get_centrality_fast,\
                                        map_centrality_matrix
    from CPAC.cwas.subdist import norm_cols
    
    # Raise exceptions if no weights/methods are selected
     if method_options.count(True) == 0:  
         raise Exception("Invalid values in method_options " \
                         "At least one True value is required")
   
    if weight_options.count(True) == 0:
        raise Exception("Invalid values in weight options" \
                        "At least one True value is required")
    
    import time
    start = time.clock()
    
    # Init variables
    out_list = []
    ts, aff, mask, t_type, scans = load(datafile, template)
    
    # If we're doing eigenvectory centrality, need entire correlation matrix
    if method_option == 1:
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
        r_value = convert_pvalue_to_r(threshold)
        centrality_matrix = calc_via_rvalue(ts_normd, 
                                            mask, 
                                            method_option, 
                                            weight_options, 
                                            r_value, 
                                            block_size)
    # Sparsity threshold
    elif threshold_option == 1:
        centrality_matrix = calc_via_sparsity(ts_normd, 
                                              mask, 
                                              method_option, 
                                              weight_options, 
                                              threshold, 
                                              block_size)
    # R-value threshold centrality
    elif threshold_option == 2:
        centrality_matrix = calc_via_rvalue(ts_normd, 
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
    print 'Timing:', (time.clock() - start)
 
    # Map the arrays back to images
    for mat in centrality_matrix:
        centrality_image = map_centrality_matrix(mat, aff, mask, t_type)
        out_list.append(centrality_image)
    
    return out_list

