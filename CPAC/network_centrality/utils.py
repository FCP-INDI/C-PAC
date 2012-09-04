
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
    

def convert_sparsity_to_r(rmatrix, threshold):
    import numpy as np
    """
    Method to calculate correlation threshold from sparsity threshold
    
    Parameters
    ----------
    rmatrix : array_like
        correaltion matrix
    threshold : float
        input sparsity threshold
    
    Returns
    -------
    rvalue : float
        correlation threshold value 
     
    """

    #SparsityThreshold=0.0744
    print "Sparsity threshold ->", threshold
    s = rmatrix.shape[0]
    upperT = np.triu(np.ones([s,s]) - np.eye(s)).astype('bool')
    #getting only the upper triangle, since it is an symmetric matrix
    val = rmatrix[upperT]
    del upperT
    val.sort()
    size = np.round(val.size*threshold)
    rvalue= val[-size:][0]
    
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
    
    
def get_centrality_matrix(threshold_matrix, correlation_matrix, 
                          weight_options, method_options):
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
    import scipy.sparse.linalg as LA
    import numpy as np
    
    corr_matrix = load_mat(correlation_matrix)
    adjacency_matrix = load_mat(threshold_matrix)
    out_list=[]
    
    try:
        
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
                print "eigenValue : ", eigenValue
                eigen_matrix=(matrix.dot(np.abs(eigenVector)))/eigenValue[0]
                return eigen_matrix
        except:
            raise Exception("Exception in calculating eigenvector centrality")
        
        if weight_options[0]:
            
            if method_options[0]:
                print "calculating binarize degree centrality matrix..."
                degree_matrix = adjacency_matrix.sum(axis=1) -1
                out_list.append(('degree_centrality_binarize', degree_matrix))
            
            if method_options[1]:
                print "calculating eigen vector centrality matrix..."
                eigen_matrix_binarize = getEigenVectorCentrality(adjacency_matrix.astype(np.float))
                out_list.append(('eigenvector_centrality_binarize', eigen_matrix_binarize))
        
        if weight_options[1]:
            
            weighted_adjacency_matrix = adjacency_matrix*corr_matrix
            
            if method_options[0]:
                print "calculating weighted degree centrality matrix..."
                degree_matrix = (weighted_adjacency_matrix).sum(axis=1) -1
                out_list.append(('degree_centrality_weighted', degree_matrix))
            
            if method_options[1]:
                print "calculating weighted eigen vector centrality matrix..."
                eigen_matrix_weighted = getEigenVectorCentrality(weighted_adjacency_matrix)
                out_list.append(('eigenvector_centrality_weighted', eigen_matrix_weighted))
    
    
    except Exception:
        print "Error while calculating centrality"
        raise
    
    return out_list


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
        print "Error in mapping centrality matrix to nifti image", out_file, matrix.shape
        raise
        



        

    
    
