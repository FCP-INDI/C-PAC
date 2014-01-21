import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util

from CPAC.network_centrality import *
from CPAC.network_centrality.core import *

def create_resting_state_graphs(allocated_memory = None,
                                wf_name = 'resting_state_graph'):
    """
    Workflow to calculate degree and eigenvector centrality
    measures for the resting state data.
    
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
        
        inputspec.threshold_option: string (int)
            threshold options:  0 for probability p_value, 1 for sparsity threshold, any other for threshold value
        
        inputspec.threshold: string (float)
            pvalue/sparsity_threshold/threshold value
           
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
    
    """
    wf = pe.Workflow(name = wf_name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                       'template',
                                                       'threshold',
                                                       'threshold_option']),
                        name='inputspec')
    
    outputspec = pe.Node(util.IdentityInterface(fields=['centrality_outputs',
                                                        'threshold_matrix',
                                                        'correlation_matrix',
                                                        'graph_outputs']),
                         name = 'outputspec')
    
    centrality_options = pe.Node(util.IdentityInterface(fields = ['weight_options',
                                                                  'method_options']),
                                 name = 'centrality_options')
            
    
    calculate_centrality = pe.Node(util.Function(input_names = ['datafile', 
                                                                'template', 
                                                                'method_options', 
                                                                'weight_options',
                                                                'option',
                                                                'threshold',
                                                                'allocated_memory'],
                                                 output_names = ['out_list'],
                                                 function = calc_centrality),
                                   name = 'calculate_centrality')
    
    wf.connect(inputspec, 'subject', 
               calculate_centrality, 'datafile')
    wf.connect(inputspec, 'template', 
               calculate_centrality, 'template')
    
    wf.connect(centrality_options, 'method_options',
               calculate_centrality, 'method_options')
    wf.connect(centrality_options, 'weight_options',
               calculate_centrality, 'weight_options')
    
    wf.connect(inputspec, 'threshold', 
               calculate_centrality, 'threshold')
    wf.connect(inputspec, 'threshold_option', 
               calculate_centrality, 'option')
    
    calculate_centrality.inputs.allocated_memory = allocated_memory
    
    wf.connect(calculate_centrality, 'out_list',
               outputspec, 'centrality_outputs')
    
    return wf



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
                   method_options,
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
        block_size  = calc_blocksize(timeseries, memory_allocated, include_full_matrix=True)
        
        nvoxs = timeseries.shape[0]
        ntpts = timeseries.shape[1]
        
        calc_degree  = method_options[0]
        calc_eigen   = method_options[1]
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
                out_list.append(('degree_centrality_binarize', degree_matrix))
            if out_weighted:
                print "...calculating weighted degree"
                degree_weighted = degree_centrality(corr_matrix, r_value, method="weighted")
                out_list.append(('degree_centrality_weighted', degree_matrix))
        
        
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



def get_centrality_by_thresh(timeseries,
                          method_options,
                          weight_options,
                          memory_allocated,
                          threshold, 
                          r_value = None):
    """
    Method to calculate degree and eigen vector centrality. 
    This method takes into consideration the amount of memory
    allocated by the user to calculate degree centrality.
    
    Parameters
    ----------
    timeseries_data : numpy array
        timeseries of the input subject
    method_options : string (list of boolean)
        list of two booleans for degree and eigen options respectively
    weight_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    memory_allocated : a string
        amount of memory allocated to degree centrality
    threshold : float
        p-value threshold for the correlation values (ignored if the r_value option is specified)
    r_value : float
        threshold value in terms of the correlation (this will override the threshold option)
    
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
                                        convert_pvalue_to_r,\
                                        degree_centrality,\
                                        eigenvector_centrality
    from CPAC.cwas.subdist import norm_cols
    
    try:                         
        out_list = []
        nvoxs = timeseries.shape[0]
        ntpts = timeseries.shape[1]
                
        r_matrix = None
        
        calc_degree  = method_options[0]
        calc_eigen   = method_options[1]
        out_binarize = weight_options[0]
        out_weighted = weight_options[1]
        
        if calc_eigen:
            block_size = calc_blocksize(timeseries, memory_allocated, include_full_matrix=True)
        else:
            block_size = calc_blocksize(timeseries, memory_allocated)
        
        if r_value == None:
            print "Calculating threshold"
            r_value = convert_pvalue_to_r(ntpts, threshold)
            print "...%s -> %s" % (threshold, r_value)
        
        
        print "Setup Intermediates/Outputs"
        
        if calc_degree:
            print "...degree"
            if out_binarize:
                degree_binarize = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('degree_centrality_binarize', degree_binarize))
            if out_weighted:
                degree_weighted = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('degree_centrality_weighted', degree_weighted))
        
        if calc_eigen:
            print "...eigen"
            r_matrix = np.zeros((nvoxs, nvoxs), dtype=timeseries.dtype)
            if out_binarize:
                eigen_binarize = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('eigenvector_centrality_binarize', eigen_binarize))
            if out_weighted:
                eigen_weighted = np.zeros(nvoxs, dtype=timeseries.dtype)
                out_list.append(('eigenvector_centrality_weighted', eigen_weighted))
                
        
        print "Normalize TimeSeries"
        timeseries = norm_cols(timeseries.T)
        
                
        print "Computing centrality across %i voxels" % nvoxs
        i = block_size
        j = 0
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
            
           print "...removing temporary correlation matrix"
           del corr_matrix
           
           j = i
           if i == nvoxs:
               break
           elif (i+block_size) > nvoxs:
               i = nvoxs
           else:
               i += block_size
        
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
            degree_binarize[degree_binarize!=0] = degree_binarize[degree_binarize!=0] - 1
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
                    method_options, 
                    weight_options,
                    option,
                    threshold,
                    allocated_memory):
    
    """
    Method to calculate centrality and map them to a nifti file
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to subject data file
    template : string (nifti file)
        path to mask/parcellation unit
    method_options : list (boolean)
        list of two booleans for degree and eigen options respectively
    weight_options : list (boolean)
        list of two booleans for binarize and weighted options respectively
    option : an integer
        0 for probability p_value, 1 for sparsity threshold, 
        2 for actual threshold value, and 3 for no threshold and fast approach
    threshold : a float
        pvalue/sparsity_threshold/threshold value
    allocated_memory : string
        amount of memory allocated to degree centrality
    
    
    Returns
    -------
    out_list : list
        list containing out mapped centrality images
        
    """
    
    from CPAC.network_centrality import load,\
                                        get_centrality_by_thresh,\
                                        get_centrality_by_sparsity,\
                                        get_centrality_fast,\
                                        map_centrality_matrix
    
    out_list = []
    
    if method_options.count(True) == 0:  
        raise Exception("Invalid values in method_options " \
                        "At least one True value is required")
   
    if weight_options.count(True) == 0:
        raise Exception("Invalid values in weight options" \
                        "At least one True value is required")
   
    
    ts, aff, mask, t_type, scans = load(datafile, template)
    
    
    import time
    start = time.clock()
    
    #p-value threshold centrality
    if option == 0:
        centrality_matrix = get_centrality_by_thresh(ts,
                                               method_options, 
                                               weight_options,
                                               allocated_memory,
                                               threshold)
                                               
    #for sparsity threshold
    elif option == 1: 
        centrality_matrix = get_centrality_by_sparsity(ts, 
                                           method_options,
                                           weight_options,
                                           threshold,
                                           allocated_memory)
                                               
    #r-value threshold centrality
    elif option == 2:
        centrality_matrix = get_centrality_by_thresh(ts, 
                                           method_options,
                                           weight_options,
                                           allocated_memory, 
                                           None,
                                           threshold)
    
    #for fast approach (no thresholding)
    elif option == 3:
        centrality_matrix = get_centrality_fast(ts, method_options)
    
    else:
        raise Exception("option must be between 0-3 and not %s" % str(option))
    
    print "timing:", (time.clock() - start)
    
    def get_image(matrix, template_type):
        centrality_image = map_centrality_matrix(matrix, 
                                                 aff, 
                                                 mask,
                                                 template_type)
        out_list.append(centrality_image) 
         
    for mat in centrality_matrix:
        get_image(mat, t_type)
               
    return out_list
