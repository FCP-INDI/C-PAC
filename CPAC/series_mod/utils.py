#in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
#mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')


def compute_ROI_corr(in_file, mask_file):

    from CPAC.series_mod import gen_roi_timeseries
    from CPAC.series_mod import corr

    ROI_data = gen_roi_timeseries(in_file, mask_file)
    corr_mat = corr(ROI_data)

    return corr_mat
    
def compute_MI(in_file, mask_file, bins):

    from CPAC.series_mod import gen_roi_timeseries
    from CPAC.series_mod import transform
    from CPAC.series_mod import mutual_information
    
    import numpy as np

    ROI_data = gen_roi_timeseries(in_file, mask_file)    
    ROI_data = transform(ROI_data,bins).astype(int)
    
    n_var = ROI_data.shape[0];
    
    MI_mat = np.zeros((n_var,n_var))    
    
    for i_ in range(n_var):
        for j_ in range(n_var):
            MI_mat[i_,j_] = mutual_information(ROI_data[i_,:],ROI_data[j_,:])
        
    ## CHECK THE MATRICES SHAPE AND RESULTS

    return MI_mat

def gen_voxel_timeseries(in_file):

    """
    Extracts voxelwise timeseries and return a np array with them.
    
    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    Returns
    -------

    data_array =  voxel(x,y,z) * timepoints

    """
    import numpy as np    
    import nibabel as nb
    
    img_data = nb.load(in_file).get_data()
    #TR = datafile.header['pixdim'][4]

    #print img_data.shape
    (n_x, n_y, n_z, n_t) = img_data.shape
    voxel_data_array = np.reshape(img_data, (n_x*n_y*n_z, n_t), order='F')
    
    
    return  voxel_data_array

def gen_roi_timeseries(in_file, mask_file):

    """
    Extracts ROI timeseries and return a np array with them.
    

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 
        
    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)
        Must be 3D    

    Returns
    -------

    data_array = ROI_number * timepoints

    """
    import numpy as np    
    import nibabel as nb
    
    img_data = nb.load(in_file).get_data()
    #TR = datafile.header['pixdim'][4]
    n_samples = img_data.shape[3]

    mask_data = nb.load(mask_file).get_data()
    # Cast as rounded-up integer
    mask_data = np.int64(np.ceil(mask_data)) #ROI numbers, int64 is enough

    if mask_data.shape != img_data.shape[:3]:
        raise Exception('Invalid Shape Error.'\
                        'Please check the voxel dimensions.'\
                        'Data and roi should have'\
                        'same shape')

    nodes = np.unique(mask_data).tolist()
    nodes.remove(0) #quits the ROI number '0'
    roi_data_array = np.zeros((len(nodes),n_samples)) #,dtype=float change?

    # Be carefull with number of ROIs and np-arrays
    nodes.sort()
    for n in nodes:
        node_array = img_data[mask_data == n]
        avg = np.mean(node_array, axis=0)
        roi_data_array[n-1] = np.round(avg, 6)
    
    return  roi_data_array

def corr(timeseries):

    """
    Computes the Network Correlation Matrix for a timeseries * timepoints
    nparray

    Parameters
    ----------

    in_file : nparray: timeseries * timepoints


    Returns
    -------

    out_file : correlation matrix

    """

    import numpy as np    
 
    corr_matrix=np.corrcoef(timeseries)
    
    
    ## IN CASE WE WOULD LIKE TO ADD A THRESHOLD FEATURE (set inside the function 
    # or from input)
    # C.corrcoef[C.corrcoef<0.7] = 0
    return  corr_matrix
  

  
#def autocorr(timeseries):
#    
#    from scipy.signal import correlation   
#    
#    autocorr_matrix = numpy.corrcoef(timeseries,'full')
#    
#    return result[autocorr_matrix.size/2:]

 

#be careful, it gets an array and takes a variable as a column!! 
def partial_corr(C):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in C, controlling 
    for the remaining variables in C.


    Parameters
    ----------
    C : array-like, shape (n, p)
        Array with the different variables. Each column of C is taken as a variable


    Returns
    -------
    P : array-like, shape (p, p)
        P[i, j] contains the partial correlation of C[:, i] and C[:, j] controlling
        for the remaining variables in C.
        
    %%%
    Partial Correlation in Python (clone of Matlab's partialcorr)
    
    This uses the linear regression approach to compute the partial 
    correlation (might be slow for a huge number of variables). The 
    algorithm is detailed here:
    
        http://en.wikipedia.org/wiki/Partial_correlation#Using_linear_regression
    
    Taking X and Y two variables of interest and Z the matrix with all the variable minus {X, Y},
    the algorithm can be summarized as
    
        1) perform a normal linear least-squares regression with X as the target and Z as the predictor
        2) calculate the residuals in Step #1
        3) perform a normal linear least-squares regression with Y as the target and Z as the predictor
        4) calculate the residuals in Step #3
        5) calculate the correlation coefficient between the residuals from Steps #2 and #4; 
    
        The result is the partial correlation between X and Y while controlling for the effect of Z
    
    
    Date: Nov 2014
    Author: Fabian Pedregosa-Izquierdo, f@bianp.net
    Testing: Valentina Borghesani, valentinaborghesani@gmail.com    
            
    """
    import numpy as np
    from scipy import stats, linalg
    
    C = np.asarray(C)
    C = C.T ## HERE PROBLEM SOLVED, but take this in mind !!
    p = C.shape[1]
    P_corr = np.zeros((p, p), dtype=np.float)
    for i in range(p):
        P_corr[i, i] = 1
        for j in range(i+1, p):
            idx = np.ones(p, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = linalg.lstsq(C[:, idx], C[:, j])[0]
            beta_j = linalg.lstsq(C[:, idx], C[:, i])[0]
 
            res_j = C[:, j] - C[:, idx].dot( beta_i)
            res_i = C[:, i] - C[:, idx].dot(beta_j)
            
            corr = stats.pearsonr(res_i, res_j)[0]
            P_corr[i, j] = corr
            P_corr[j, i] = corr
        
    return P_corr


def compute_te(in_file, mask_file):

    """
    Computes the Pairwise Transfer Entropy Matrix for ROIs in the mask
    For now, we are computing GC, since GC and TE are equivalent for Gaussian
    signals.
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.103.238701

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    mask_file : nifti file
        Mask of the EPI File(Only Compute Correlation of voxels in the mask)
        Must be 3D

    Returns
    -------

    out_file : Transfer Entropy matrix

    import numpy as np
    import nibabel as nb
    
    #NEED TO  WORK ON THIS
    
    return  g1
    
    """
    
def transform(x_old, Nbins):
    
    '''
    TRANSFORM This funcion computes transforms in a matrix
    to obtain a matrix scaled between the especified number of bins
    INPUT:
      x_old: nobservations * nvariables matrix
      Nbins: Number of bins of the transformed matrix (NBins=2 values betwen
     -1, NBins=3 values between 0-2...)

    OUTPUT:  
    x_new: New scaled matrix
    '''
     
    import numpy as np

    [npoints, num_vals] = x_old.shape
    xmax = Nbins-1
    xmin = 0
    ymax = x_old.max()
    ymin = x_old.min()
    
    #x_new = ((xmax - xmin)/(ymax - ymin) )* x_old - ( (xmax - xmin) / (ymax - ymin) ) * ymin + xmin;
    
    x = (xmax-xmin)/(ymax-ymin)
    x_new = x * x_old
    x_new = x_new - (x*ymin+xmin)
    x_new = np.round(x_new)
    
    return x_new

def entropy(*X):
    
    """
    Same function for entropy as for joint entropy
    """

    import numpy as np
    import itertools 
    
    n_insctances = len(X[0])
    H = 0
    for classes in itertools.product(*[set(x) for x in X]):
        v = np.array([True] * n_insctances)
        for predictions, c in zip(X, classes):
            v = np.logical_and(v, predictions == c)
        p = np.mean(v)
        H += -p * np.log2(p) if p > 0 else 0
    return H  
    
def mutual_information(X,Y):

    from CPAC.series_mod import entropy

    Hx = entropy(X)
    Hy = entropy(Y)
    Hxy = entropy(X,Y)
    MI = Hx + Hy - Hxy
    
    return MI

def cond_entropy(X,Y):
    
    """
    Conditional entropy H(X|Y) = H(Y,X) - H(Y). X conditioned on Y
    """
    from CPAC.series_mod import entropy

    Hy = entropy(Y)
    Hyx = entropy(Y,X)
    CE = Hyx - Hy    
    
    return CE
    
#def entropy(*X):
#    n_insctances = len(X[0])
#    H = 0
#    for classes in itertools.product(*[set(x) for x in X]):
#        v = reduce(np.logical_and, (predictions, c for predictions, c in zip(X, classes)))
#        p = np.mean(v)
#        H += -p * np.log2(p) if p > 0 else 0
#    return H   
#   
#def entropy(*X):
#    return = np.sum(-p * np.log2(p) if p > 0 else 0 for p in \
#    (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product (*[set(x) for x in X])))

    
