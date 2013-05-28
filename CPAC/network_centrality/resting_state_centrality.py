import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
from CPAC.network_centrality import *

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
    
   
    read_data = pe.Node(util.Function( input_names = ['datafile', 
                                                      'template'],
                                       output_names= ['timeseries_data', 
                                                      'affine', 
                                                      'mask_data',
                                                      'mask_type',
                                                      'scans'],
                                      function = load),
                        name = 'read_data')
    
    wf.connect(inputspec, 'subject', 
               read_data, 'datafile')
    wf.connect(inputspec, 'template', 
               read_data, 'template')
    
    
    calculate_centrality = pe.Node(util.Function(input_names = ['method_options', 
                                                                'weight_options',
                                                                'option',
                                                                'threshold',
                                                                'timeseries_data',
                                                                'scans',
                                                                'template_type',
                                                                'template_data', 
                                                                'affine',
                                                                'allocated_memory'],
                                                 output_names = ['out_list'],
                                                 function = calc_centrality),
                                   name = 'calculate_centrality')
    
    
    wf.connect(centrality_options, 'method_options',
               calculate_centrality, 'method_options')
    wf.connect(centrality_options, 'weight_options',
               calculate_centrality, 'weight_options')
    
    wf.connect(inputspec, 'threshold', 
               calculate_centrality, 'threshold')
    wf.connect(inputspec, 'threshold_option', 
               calculate_centrality, 'option')
    
    wf.connect(read_data, 'timeseries_data', 
               calculate_centrality, 'timeseries_data')
    wf.connect(read_data, 'scans', 
               calculate_centrality, 'scans')
    wf.connect(read_data, 'mask_data',
                calculate_centrality, 'template_data')
    wf.connect(read_data, 'affine',
               calculate_centrality, 'affine')
    wf.connect(read_data, 'mask_type',
               calculate_centrality, 'template_type')
    
    calculate_centrality.inputs.allocated_memory = allocated_memory
    
    wf.connect(calculate_centrality, 'out_list',
               outputspec, 'centrality_outputs')
    
    return wf



def load(datafile, template):
    
    """
    Method to read data from datafile and mask/parcellation unit
    and store the mask data, timeseries, affine matrix, mask type
    and scans. The output of this method is used by all other nodes.
    
    Parameters
    ----------
    datafile : string (nifti file)
        path to subject data file
    template : string (nifti file)
        path to mask/parcellation unit
        
    Returns
    -------
    timeseries_data: string (numpy npy file)
        path to file containing timeseries of the input data
    affine: string (numpy npy file)
        path to file containing affine matrix of teh input data
    mask_data: string (numpy npy file)
        path to file containing mask/parcellation unit matrix
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
        
        data = img.get_data().astype(np.float32)
        aff = img.get_affine()    
        mask = nib.load(template).get_data().astype(np.float32)  
        scans = data.shape[3]
        
    except:
        print "Error in loading images for graphs"
        raise
    
    
    
    if mask.shape != data.shape[:3]:
        raise Exception("Invalid Shape Error. mask and data file have"\
                        "different shape please check the voxel size of the two files")

    
    #check for parcellation
    nodes = np.unique(mask).tolist()
    nodes.sort()
    print "sorted nodes", nodes
    
    #extract timeseries
    if len(nodes)>2:
        flag=1
        for n in nodes:
            if n > 0:
                node_array = data[mask == n]
                avg = np.mean(node_array, axis =0)
                if flag:
                    timeseries = avg
                    flag=0
                else:
                    timeseries = np.vstack((timeseries, avg))
        #template_type is 1 for parcellation
        template_type = 1
    else:
        #template_type is 0 for mask
        template_type = 0
        mask = mask.astype('bool')
        timeseries = data[mask]
        
    #saving files
    img_name = os.path.splitext(os.path.splitext(os.path.basename(datafile))[0])[0] + '.npy'
    mask_name = os.path.splitext(os.path.splitext(os.path.basename(template))[0])[0]  + '.npy'
    
    mask_data = os.path.join(os.getcwd(), mask_name)
    timeseries_data = os.path.join(os.getcwd(),img_name)
    affine = os.path.join(os.getcwd(),'affine.npy')
        
    np.save(mask_data, mask)
    np.save(timeseries_data, timeseries)
    np.save(affine, aff)
    
    return timeseries_data, affine, mask_data, template_type, scans


def get_centrality(timeseries_data, 
                   method_options,
                   weight_options,
                   threshold,
                   option,
                   scans,
                   memory_allocated):
    
    """
    Method to calculate degree and eigen vector centrality
    
    Parameters
    ----------
    weight_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    method_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    threshold_matrix : string (numpy npy file)
        path to file containing thresholded correlation matrix 
    correlation_matrix : string (numpy npy file)
        path to file containing correlation matrix
    template_data : string (numpy npy file)
        path to file containing mask or parcellation unit data    
    
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
    from CPAC.network_centrality import load_mat,\
                                        calc_corrcoef,\
                                        calc_blocksize,\
                                        calc_threshold,\
                                        calc_eigenV
    
    from scipy.sparse import csc_matrix
    
    out_list=[]
    
    try:
        
        timeseries = load_mat(timeseries_data)
        shape = timeseries.shape
        block_size = calc_blocksize(shape, memory_allocated)
        corr_matrix = np.zeros((shape[0], shape[0]), dtype = np.float16)
        j=0
        i = block_size
        
        while i <= timeseries.shape[0]:
            print "block ->", i,j 
            temp_matrix = np.nan_to_num(calc_corrcoef(timeseries[j:i].T, timeseries.T))
            corr_matrix[j:i] = temp_matrix
            j = i   
            if i == timeseries.shape[0]:
                break
            elif (i+block_size) > timeseries.shape[0]: 
                i = timeseries.shape[0] 
            else:
                i += block_size
        
        r_value = calc_threshold(option, 
                                 threshold, 
                                 scans, 
                                 corr_matrix,
                                 full_matrix = True)
        
        print "r_value ->", r_value
                
        if method_options[0]:
            
            print "calculating binarize degree centrality matrix..."
            degree_matrix = np.sum( corr_matrix > r_value , axis = 1)  -1
            out_list.append(('degree_centrality_binarize', degree_matrix))
            
            print "calculating weighted degree centrality matrix..."
            degree_matrix = np.sum( corr_matrix*(corr_matrix > r_value), axis= 1) -1
            out_list.append(('degree_centrality_weighted', degree_matrix))
            
        
        if method_options[1]:
            out_list.extend(calc_eigenV(corr_matrix, 
                                           r_value, 
                                           weight_options))
    
    except Exception:
        print "Error while calculating centrality"
        raise
    
    return out_list



def get_centrality_opt(timeseries_data,
                       method_options,
                       weight_options,
                       memory_allocated,
                       threshold,
                       scans,
                       r_value = None):
    
    """
    Method to calculate degree and eigen vector centrality. 
    This method takes into consideration the amount of memory
    allocated by the user to calculate degree centrality.
    
    Parameters
    ----------
    timeseries_data : numpy array
        timeseries of the input subject
    weight_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    method_options : string (list of boolean)
        list of two booleans for binarize and weighted options respectively
    memory_allocated : a string
        amount of memory allocated to degree centrality
    scans : an integer
        number of scans in the subject
    r_value :a float
        threshold value
    
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
    from CPAC.network_centrality import load_mat,\
                                        calc_corrcoef,\
                                        calc_blocksize,\
                                        calc_eigenV,\
                                        calc_threshold
    #from scipy.sparse import dok_matrix
    
    try:                                    
        out_list =[]
        timeseries = load_mat(timeseries_data)
        shape = timeseries.shape
        try:
            block_size = calc_blocksize(shape, memory_allocated)
        except:
           raise Exception("Error in calculating block size")
        
        r_matrix = None
        
        if method_options[0]:
            if weight_options[0]:
                degree_mat_binarize = np.zeros(shape[0], dtype= np.float32)
                out_list.append(('degree_centrality_binarize', degree_mat_binarize))
    
            if weight_options[1]:
                degree_mat_weighted = np.zeros(shape[0], dtype = np.float32)
                out_list.append(('degree_centrality_weighted', degree_mat_weighted))
            
        if method_options[1]:
            r_matrix = np.zeros((shape[0], shape[0]), dtype = np.float32)
    
        j=0
        i = block_size
        
        while i <= timeseries.shape[0]:
           
           print "running block ->", i, j 
           try:
               corr_matrix = np.nan_to_num(calc_corrcoef(timeseries[j:i].T, timeseries.T))
           except:
               raise Exception("Error in calcuating block wise correlation for the block %,%"%(j,i))
           
           if r_value == None:
                r_value = calc_threshold(1, threshold, scans, corr_matrix, full_matrix = False)
    
           if method_options[1]:
               r_matrix[j:i] = corr_matrix 
    
           if method_options[0]:
               if weight_options[0]:
                   degree_mat_binarize[j:i] = np.sum((corr_matrix > r_value).astype(np.float32), axis = 1) -1
               if weight_options[1]:
                   degree_mat_weighted[j:i] = np.sum(corr_matrix*(corr_matrix > r_value).astype(np.float32), axis = 1) -1
        
           j = i   
           if i == timeseries.shape[0]:
               break
           elif (i+block_size) > timeseries.shape[0]: 
               i = timeseries.shape[0] 
           else:
               i += block_size    
        
        try:
            if method_options[1]:
                out_list.extend(calc_eigenV(r_matrix, r_value, weight_options))
        except Exception:
            print "Error in calcuating eigen vector centrality"
            raise
        
        return out_list   
    
    except Exception: 
        print "Error in calcuating Centrality"
        raise
 
 
def calc_eigenV(r_matrix, 
                r_value, 
                weight_options):
    
    """
    Method to calculate Eigen Vector Centrality
    
    Parameters
    ----------
    
    r_matrix : numpy array
        correlation matrix
    r_value : a float
        threshold value
    weight_options : list (boolean)
        list of two booleans for binarize and weighted options respectively
        
        
    Returns
    -------
    out_list : list
        list containing eigen vector centrality maps
        
    
    """
    
    
    import scipy.sparse.linalg as LA
    import numpy as np
    from scipy.sparse import csc_matrix    

    out_list =[]
    
    try:
        def getEigenVectorCentrality(matrix):
            """
            from numpy import linalg as LA
            w, v = LA.eig(a)
            index = np.argmax(w)
            eigenValue = w.max()
            eigenvector= v[index]
            """
            #using scipy method, which is a wrapper to the ARPACK functions
            #http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
            eigenValue, eigenVector= LA.eigsh(matrix, k=1, which='LM', maxiter=1000)
            print "eigenValues : ", eigenValue
            eigen_matrix=(matrix.dot(np.abs(eigenVector)))/eigenValue[0]
            return eigen_matrix
    except:
        raise Exception("Exception in calculating eigenvector centrality")
    
    if weight_options[0]:
        print "calculating eigen vector centrality matrix..."
        eigen_matrix_binarize = getEigenVectorCentrality(csc_matrix((r_matrix> r_value).astype(np.float32)))        
        out_list.append(('eigenvector_centrality_binarize', eigen_matrix_binarize))
    
    if weight_options[1]:
        print "calculating weighted eigen vector centrality matrix..."
        eigen_matrix_weighted = getEigenVectorCentrality(csc_matrix(r_matrix*(r_matrix > r_value).astype(np.float32)))        
        out_list.append(('eigenvector_centrality_weighted', eigen_matrix_weighted))
        
    return out_list


def calc_centrality(method_options, 
                    weight_options,
                    option,
                    threshold,
                    timeseries_data,
                    scans,
                    template_type,
                    template_data, 
                    affine,
                    allocated_memory):
    
    """
    Method to calculate centrality and map them to a nifti file
    
    Parameters
    ----------
        
    method_options : list (boolean)
        list of two booleans for binarize and weighted options respectively
    weight_options : list (boolean)
        list of two booleans for binarize and weighted options respectively
    option : an integer
        0 for probability p_value, 1 for sparsity threshold, 
        any other for threshold value
    threshold : a float
        pvalue/sparsity_threshold/threshold value
    timeseries_data : string (numpy filepath)
        timeseries of the input subject
    scans : an integer
        number of scans in the subject
    template_type : an integer
        0 for mask, 1 for roi
    template_data : string (numpy filepath)
        path to file containing mask/parcellation unit matrix
    affine : string (filepath)
        path to file containing affine matrix of the input data
    allocated_memory : string
        amount of memory allocated to degree centrality
    
    
    Returns
    -------
    out_list : list
        list containing out mapped centrality images
        
    """
    
    from CPAC.network_centrality import map_centrality_matrix,\
                                        get_centrality, \
                                        get_centrality_opt,\
                                        calc_threshold
    
    out_list = []
    
    if method_options.count(True) == 0:  
        raise Exception("Invalid values in method_options " \
                        "Atleast one True value is required")
   
    if weight_options.count(True) == 0:
        raise Exception("Invalid values in weight options" \
                        "Atleast one True value is required")
   
    #for sparsity threshold
    if option == 1 and allocated_memory == None:
        
        centrality_matrix = get_centrality(timeseries_data, 
                                           method_options,
                                           weight_options,
                                           threshold,
                                           option,
                                           scans,
                                           allocated_memory)
    #optimized centrality
    else:
        
        if option ==1 :
            r_value = None
        else:
            r_value = calc_threshold(option, 
                                     threshold,
                                     scans)
        
        print "inside optimized_centraltity, r_value ->", r_value
        
        centrality_matrix = get_centrality_opt(timeseries_data,
                                               method_options, 
                                               weight_options,
                                               allocated_memory,
                                               threshold,
                                               scans,
                                               r_value)     
        
    def get_image(matrix, template_type):
        
        centrality_image = map_centrality_matrix(matrix, 
                                                 affine, 
                                                 template_data,
                                                 template_type)
        out_list.append(centrality_image) 
         
    for mat in centrality_matrix:
        get_image(mat, template_type)
               
    return out_list