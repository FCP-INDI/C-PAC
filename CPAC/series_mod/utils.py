#in_file = ('/home/asier/git/C-PAC/CPAC/series_mod/Standard-clean_func_preproc.nii.gz')
#mask_file = ('/home/asier/git/C-PAC/CPAC/series_mod/AAL_Contract_90_3MM.nii.gz')

# THIS SCRIPT USES Nvar * Ntimepoints LIKE MATRIX STRUCTURES

#from matplotlib.pylab import *
#matshow(corr_mat);colorbar();clim(-1,1);show()
#matshow(pcorr_mat);colorbar();clim(-1,1);show()
#matshow(MI_mat);colorbar();clim(0,np.max(MI_mat));show()


def compute_ROI_corr(in_file, mask_file):

    from CPAC.series_mod import gen_roi_timeseries
    from CPAC.series_mod import corr  
    import numpy as np
    
    ROI_data = gen_roi_timeseries(in_file, mask_file)
    corr_mat = corr(ROI_data)
    
    np.savetxt(in_file[:-7]+'_corr.txt', corr_mat)
    make_image_from_bin( in_file, corr_mat, mask_file) 
    
#    img = nb.Nifti1Image(K, header=res_img.get_header(), affine=res_img.get_affine())
#
#    reho_file = os.path.join(os.getcwd(), 'ReHo.nii.gz')
#
#    img.to_filename(reho_file)
#
#    out_file = reho_file

    return corr_mat
    
def compute_ROI_pcorr(in_file, mask_file):

    from CPAC.series_mod import gen_roi_timeseries
    from CPAC.series_mod import partial_corr  
    import numpy as np
    
    ROI_data = gen_roi_timeseries(in_file, mask_file)
    pcorr_mat = partial_corr(ROI_data)
    
    np.savetxt(in_file[:-7]+'partial_corr.txt', pcorr_mat)
    make_image_from_bin( in_file, pcorr_mat, mask_file)    

    return pcorr_mat    
    
def compute_MI(in_file, mask_file):

    from CPAC.series_mod import gen_roi_timeseries
    from CPAC.series_mod import transform
    from CPAC.series_mod import mutual_information
    import numpy as np
    import math
    


    ROI_data = gen_roi_timeseries(in_file, mask_file)    
    
    n_var = ROI_data.shape[0]
    points = ROI_data.shape[1]
    bins = math.pow(points, 1/3.) #to the 3rd due to Equiquantization formula
    # Proposed by Milan Palus. n+1 where n is the number of vars in the computation
    # as it is pairwise, n+1 is 3
    bins = np.round(bins)
    
    ROI_data = transform(ROI_data,bins).astype(int)
    
 
    
    MI_mat = np.zeros((n_var,n_var))    
    
    for i_ in range(n_var):
        for j_ in range(n_var):
            MI_mat[i_,j_] = mutual_information(ROI_data[i_,:],ROI_data[j_,:])
        
    np.savetxt(in_file[:-7]+'_MI.txt', MI_mat)
    make_image_from_bin( in_file, MI_mat, mask_file)    
        
    ## CHECK THE MATRICES SHAPE AND RESULTS

    return MI_mat
    
def compute_TE(in_file, mask_file):

    from CPAC.series_mod import gen_roi_timeseries
    from CPAC.series_mod import transform
    from CPAC.series_mod import transfer_entropy
    import numpy as np
    import math
    


    ROI_data = gen_roi_timeseries(in_file, mask_file)    
    
    n_var = ROI_data.shape[0]
    points = ROI_data.shape[1]
    bins = math.pow(points, 1/3.) #to the 3rd due to Equiquantization formula
    # Proposed by Milan Palus. n+1 where n is the number of vars in the computation
    # as it is pairwise, n+1 is 3
    bins = np.round(bins)
    
    ROI_data = transform(ROI_data,bins).astype(int)
    
    TE_mat = np.zeros((n_var,n_var))    
    
    for i_ in range(n_var):
        for j_ in range(n_var):
            TE_mat[i_,j_] = transfer_entropy(ROI_data[i_,:],ROI_data[j_,:],1)
    
    
    np.savetxt(in_file[:-7]+'_TE.txt', TE_mat)  
    make_image_from_bin( in_file, TE_mat, mask_file)
    
    return TE_mat

    
def compute_ApEn(in_file, m_param, r_param):

    from CPAC.series_mod import gen_voxel_timeseries
    from CPAC.series_mod import ap_entropy 
    import numpy as np
    
    data = gen_voxel_timeseries(in_file)
    ApEn_vector = ap_entropy(data,m_param,r_param)
    
    np.savetxt(in_file[:-7]+'_ApEn.txt', ApEn_vector)
    #make_image_from_bin( in_file, ApEn_vector, mask_file )    


    return ApEn_vector




#
#
#
#
# _          _                    __                  _   _                 
#| |        | |                  / _|                | | (_)                
#| |__   ___| |_ __   ___ _ __  | |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
#| '_ \ / _ \ | '_ \ / _ \ '__| |  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#| | | |  __/ | |_) |  __/ |    | | | |_| | | | | (__| |_| | (_) | | | \__ \
#|_| |_|\___|_| .__/ \___|_|    |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#             | |                                                           
#             |_|      
#
#
#
#




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
    for index, n in enumerate(nodes):
        node_array = img_data[mask_data == n]
        avg = np.mean(node_array, axis=0)
        roi_data_array[index] = np.round(avg, 6)
        
#    for n in nodes:
#        node_array = img_data[mask_data == n]
#        avg = np.mean(node_array, axis=0)
#        roi_data_array[n-1] = np.round(avg, 6)        
        
    
    return  roi_data_array

def make_image_from_bin( image, binfile, mask ):

    import numpy as np 
    import nibabel as nb

    # read in the mask
    nim=nb.load(mask)

    # read in the binary data    
    if( binfile.endswith(".npy") ):
        print "Reading",binfile,"as a npy filetype"
        a = np.load(binfile)
    else:
        print "Reading",binfile,"as a binary file of doubles"
        a = np.fromfile(binfile)

    imdat=nim.get_data()
    print "shape",np.shape(a)
    print "sum",sum(imdat)

    # map the binary data to mask
    mask_voxels=(imdat.flatten()>0).sum()
    print "shape2",np.shape(a[0:mask_voxels])
    imdat[imdat>0]=np.short(a[0:mask_voxels].flatten())

    # write out the image as nifti
    thdr=nim.get_header()
    thdr['scl_slope']=1
    
    nim_aff = nim.get_affine()

    nim_out = nb.Nifti1Image(imdat, nim_aff, thdr)
    #nim_out.set_data_dtype('int16')
    nim_out.to_filename(image)

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
        Array with the different variables. Each row of C is taken as a variable !!


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
    C = C.T ## HERE PROBLEM SOLVED, but take this in mind !! Now, row = variable
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
    
def entropy_cc(X,Y): #ENTROPY CORRELATION COEFFICIENT
    
    """
    Entropy correlation coefficient p(H) = sqrt(I(X,Y)/0.5*(H(X)+H(Y)))
    """
    import numpy as np
    from CPAC.series_mod import entropy
    from CPAC.series_mod import mutual_information
    
    Hx = entropy(X)
    Hy = entropy(Y)
    Ixy = mutual_information(Y,X)
    ECC = np.sqrt(Ixy/(0.5*(Hx+Hy)))
    
    return ECC   
    

def transfer_entropy(X, Y, lag):
    #========================================================
    # TRANSFER_ENTROPY This funcion computes the transfer entropy for two given 
    # signals
    #  te =transfer_entropy(data,lag)
    #  INPUT:
    #   X: target (points of the signal)
    #   Y: source (points of the signal)
    #   lag: The number of samples to lag to obtain future series
    #  OUTPUT:  
    #   te: Raw Transfer entropy 

    import numpy as np
    from CPAC.series_mod import cond_entropy
    from CPAC.series_mod import entropy
    
    # future of i
    Fi = np.roll(X, -lag)
    # past of i
    Pi = X
    # past of j
    Pj = Y
    
    #Transfer entropy
    Inf_from_Pi_to_Fi = cond_entropy(Fi, Pi)

    # same as cond_entropy(Fi, Pi_Pj)
    Hy = entropy(Pi,Pj)
    Hyx = entropy(Fi,Pj,Pi)  
    Inf_from_Pi_Pj_to_Fi = Hyx - Hy     

    TE_from_j_to_i = Inf_from_Pi_to_Fi-Inf_from_Pi_Pj_to_Fi         
        
    return TE_from_j_to_i    
    
def gtetruncated(Y,X,Z):

    # X=N x 1  Y=N x m  Z=N x m'
    # te = gaussian tranfer entropy  Y -> X | Z
    # Y=series[:,:,0].T
    # X=Y[:,0].copy()
    # Z=series[:,:,1].T

    import numpy as np

    #Z = Z.T
    #Y = Y.T
    
    [Nobs, nvar] = np.shape(Z)
    
    
    
    xzy = np.dot(X.T, np.array(np.hstack((Y, Z)))) / (Nobs-1)
    #xzy = np.expand_dims(xzy,0)
    
    covxy = np.cov(np.array(np.hstack((Y, Z))).T)    
    xzy_lstsq_covxy = np.linalg.lstsq(covxy.T,xzy.T)[0].T
    xzy_lstsq_covxy_dot_zxy = np.dot(xzy_lstsq_covxy,xzy.T)
        
    sigmazy = np.var(X) - xzy_lstsq_covxy_dot_zxy
    #sigmazy=var(X)-xzy/cov([Y Z])*xzy';

    xz = np.dot(X.T, Z) / (Nobs-1)
    #xz = np.expand_dims(xz,0)

    covz = np.cov(Z.T)   
    xz_lstsq_covz = np.linalg.lstsq(covz.T,xz.T)[0].T
    xz_lstsq_covz_dot_xz = np.dot(xz_lstsq_covz,xz.T)
    
    sigmaz = (np.var(X) - xz_lstsq_covz_dot_xz)
    #sigmaz=var(X)-xz/cov(Z)*xz'; 
    
    te = 0.5 * np.log(sigmaz / sigmazy)
    
    return te, sigmaz, sigmazy   
    
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




def ap_entropy(X, M, R):
        
    import numpy as np
	
    N = len(X)
    
    Em = embed_seq(X, 1, M)	
    Emp = embed_seq(X, 1, M + 1)
    
    Cm, Cmp = np.zeros(N - M + 1), np.zeros(N - M)

    for i in xrange(0, N - M):
        for j in xrange(i, N - M):
            if in_range(Em[i].any(), Em[j].any(), R):
                Cm[i] += 1												
                Cm[j] += 1
                if abs(Emp[i][-1] - Emp[j][-1]) <= R:
                    Cmp[i] += 1
                    Cmp[j] += 1
        if in_range(Em[i].any(), Em[N-M].any(), R):
            Cm[i] += 1
            Cm[N-M] += 1
	
    Cm[N - M] += 1 
    Cm /= (N - M +1 )
    Cmp /= ( N - M )
    Phi_m, Phi_mp = sum(np.log(Cm)),  sum(np.log(Cmp))

    ApEn = (Phi_m - Phi_mp) / (N - M)

    return ApEn

def embed_seq(X,Tau,D):
    
    import numpy as np
    
    N = len(X)

    Y=np.zeros((N - (D - 1) * Tau, D))
    for i in xrange(0, N - (D - 1) * Tau):
        for j in xrange(0, D):
            Y[i][j] = X[i + j * Tau]
    return Y


def in_range(i,j,k):
    """
    Returns True if i in [j, k[
    * 0 <= i, j, k < MAX
    * no order is assumed between j and k: we can have k < j
    """
    if j <= k:
        return j <= i < k
    # j > k :
    return j <= i or i < k



def butter_bandpass(lowcut, highcut, fs, order=5):
    
    from scipy.signal import butter    
    
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    
    from scipy.signal import lfilter 
    
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    
    return y
    
def phase_sync(X):
    
    # Computes the PS_value as explained in Ponce et al. 2015
    import numpy as np
    import scipy.signal.signaltools as sigtool

    htx = sigtool.hilbert(X)
    # This converts (real, imag) to an angle
    px = np.angle(htx)
        
    return np.sqrt(np.mean(np.cos(px))**2 + np.mean(np.sin(px))**2) #PS_value      
    
    
def PLV(X, Y):
    """
    FUNC: phasediff - phase locking value
    DESCR: Compute the phase difference using the hilbert transform of s0
    and s1 and then the phase lock

    return
        PLV    : the phase difference
    
    """
    import numpy as np
    import scipy.signal.signaltools as sigtool


    htx = sigtool.hilbert(X)
    hty = sigtool.hilbert(Y)
    px = np.angle(htx)
    py = np.angle(hty)
    
    psi = px - py
    
    return np.sqrt(np.mean(np.cos(psi))**2 + np.mean(np.sin(psi))**2) #PLV_value

