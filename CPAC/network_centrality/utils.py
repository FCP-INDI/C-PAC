
def convert_pvalue_to_r(scans, threshold):
        
    import scipy.stats as s
    import math
    
    """
    Method to calculate correlation threshold from p_value
    
    Parameters
    ----------
    scans : int
        Total number of scans in the data
    threshold : float
        input p_value
    
    Returns
    -------
    rvalue : float
        correlation threshold value 
    """
    
    #p_value =0.05
    print "p_value ->", threshold
    x = 1-threshold/2
    dof = scans-2
    #Inverse Survival Function (Inverse of SF)
    tvalue = s.t.isf(x, dof)
    rvalue = math.sqrt(math.pow(tvalue, 2)/(dof+ math.pow(tvalue,2)))
    
    return rvalue
    

def convert_sparsity_to_r(rmatrix, threshold, full_matrix):
    import numpy as np
    """
    Method to calculate correlation threshold from sparsity threshold
    
    Parameters
    ----------
    rmatrix : array_like
        correlation matrix
    threshold : float
        input sparsity threshold
    full_matrix : boolean
        True, if sparsity threshold is calculated on the entire matrix
        False, if sparsity threshold is calculated only for the block 
    
    Returns
    -------
    rvalue : float
        correlation threshold value 
     
    """

    #SparsityThreshold=0.0744
    print "Sparsity threshold ->", threshold
    
    def get_upper_triangle(matrix):
        s = matrix.shape[0]
        upperT = np.triu(np.ones([s,s]) - np.eye(s)).astype('bool')
        #getting only the upper triangle, since it is an symmetric matrix
        val = matrix[upperT]
        return val
    
    
    if full_matrix:
        val = get_upper_triangle(rmatrix)
        val.sort()
        size = np.round(val.size*threshold)
        rvalue= val[-size]
    else:
        #split the rmatrix into a square matrix and a rectangle block
        val1= get_upper_triangle(rmatrix[:,:rmatrix.shape[0]])
        val2 = rmatrix[:, rmatrix.shape[0]:].flatten()
        #concatenate two arrays
        val3 = np.concatenate([val1,val2])
        #sort the array
        val3.sort()
        #calculating sparsity threshold for a block
        size = np.round(val3.size*threshold)
        rvalue = val3[-size]
    
    return rvalue

    
def load_mat(mat_file):
    """
    Simple method to load a npy file
    
    Parameters
    ----------
    mat_file : string (numpy file or list of numpy file)
        any numpy image
    
    Returns
    -------
    matrix : numpy matrix
        
    """
    
    import numpy as np
        
    if isinstance(mat_file, list):
        matrix = np.load(mat_file[0])
    else:
        matrix = np.load(mat_file)
    return matrix 
    


def calc_threshold(option, 
                   threshold,
                   scans,
                   corr_matrix= None,
                   full_matrix = True):  
    
    """
    Method to calculate threshold based
    on the threshold method chosen
    
    Parameters
    ----------
    option : an integer
        threshold option
    threshold : a float
        thrshold value
    scans : an integer
        no of timepoints
    corr_matrix : numpy array
        correlation matrix
    full_matrix : boolean
        True, if full matrix is considered.
        False, if only upper triangle is considered.
    
    Return 
    ------
    r_value : a float
        threshold value
    
    
    """
        
    print "threshold_option -->", option
     
    try:
         if option == 0:
             r_value = convert_pvalue_to_r(scans, threshold)
         elif option == 1:
             r_value = convert_sparsity_to_r(corr_matrix, threshold, full_matrix)
         else:
             r_value = threshold
    except:
         print "Exception in calculating threshold value"
         raise
     
    print "r_value --> ", r_value
     
    return r_value
 

def map_centrality_matrix(centrality_matrix, affine, template_data, template_type):
    """
    Method to map centrality matrix to a nifti image
    
    Parameters
    ----------
    centrality_matrix : tuple (string, array_like)
        tuple containing matrix name and degree/eigenvector centrality matrix
    affine : string (numpy mat file)
        path to file containing image affine matrix
    template_data : string (numpy mat file)
        path to file containing mask or roi data matrix
    template_type : int
        type of template: 0 for mask, 1 for roi
    
    Returns
    -------
    out_file : string (nifti image)
        nifti image mapped from the centrality matrix
    
    Raises
    ------
    Exception
    
    """
    
    import nibabel as nib
    import os
    import numpy as np
    
    try:
        
        mask = load_mat(template_data)   
        aff = load_mat(affine)
        out_file, matrix = centrality_matrix
       
        out_file = os.path.join(os.getcwd(), out_file + ".nii.gz")
        sparse_m = np.zeros((mask.shape), dtype=float)
     
        print "mapping centrality matrix to nifti image...", out_file
            
        if int(template_type) == 0:
            cords = np.argwhere(mask)        
            index=0
            for val in cords:
                x,y,z=val
                sparse_m[x,y,z]= matrix[index]
                index+=1
        
        elif int(template_type) == 1:
            nodes = np.unique(mask).tolist()
            nodes.sort()
            index = 0
            for n in nodes:
                if n> 0:
                    cords = np.argwhere(mask==n)
                    for val in cords:
                        x,y,z = val
                        if isinstance(matrix[index], list):
                            sparse_m[x,y,z] = matrix[index][0]
                        else:
                            sparse_m[x,y,z]=matrix[index]
                    index+=1
                        
    
        nifti_img = nib.Nifti1Image(sparse_m, aff)
        nifti_img.to_filename(out_file)
        
        return out_file
    except:
        print "Error in mapping centrality matrix to nifti image"
        raise
   
   
def calc_corrcoef(X, Y=None):
    """
    Method to calculate correlation 
    Each of the columns in X will be correlated 
    with each of the columns in Y. Each column 
    represents a variable, with the rows containing the observations.
    
    Parameters
    ----------
    X : numpy array
       array of shape x1, x2
    Y : numpy array
      array of shape y1, y2
    
    Returns
    -------
    r : numpy array
      array containing correlation values of shape x2, y2
    
    """
    import numpy as np
    
    if Y is None:
        Y = X
    
    if X.shape[0] != Y.shape[0]:
        raise Exception("X and Y must have the same number of rows.")
    
    X = X.astype(float)
    Y = Y.astype(float)
    
    X -= X.mean(axis=0)[np.newaxis,...]
    Y -= Y.mean(axis=0)
    
    xx = np.sum(X**2, axis=0)
    yy = np.sum(Y**2, axis=0)
    
    r = np.dot(X.T, Y)/np.sqrt(np.multiply.outer(xx,yy))
    
    return r


def calc_blocksize (shape, memory_allocated = None):
    """
    Method to calculate blocksize to calculate correlation matrix
    as per the memory allocated by the user. By default, the block
    size is 1000. 

    Parameters
    ----------
    shape : tuple
       shape of array
    memory_allocated : float
       memory allocated in GB for degree centrality
    
    Returns
    -------
    block_size : an integer
      size of block for matrix calculation
    """
    
    import warnings
    
    block_size = 1000
    
    def get_size(num, unit):
        
        for x in range(3):
            if unit == 'GB':
                num /= 1024.0
            elif unit == 'bytes':
                num *= 1024.0
        return float(num)
    
    if memory_allocated:
        block_size =  int(0.8*(get_size(memory_allocated, 'bytes') - shape[0]*shape[1]*8 - shape[0]*8*2)/(shape[0]*8*4 + shape[1]*8))
        
    if block_size > shape[0]:
        block_size = shape[0]
    elif block_size < 1:
        raise MemoryError(" Not enough memory available to perform degree centrality")
            
    print "block_size -> ", block_size
    
    
    return block_size
    
    
def check_timeseries(data):
    """
    Method to check if the array contains
    any zeros values. If it contains zeros
    then return the indices of those points. 
    
    Parameters
    ----------
    data : numpy array
    
    Returns
    -------
    index : list
        indices of all where a
    data : numpy array
    """
    index= np.where(np.all(data==0, axis=1))[0].tolist()
    print "index where timeseries is zero ", index
    
    if index:
        data = data[~np.all(data == 0, axis=1)]
        print "new shape", data.shape
        
    return index, data 

    
    
