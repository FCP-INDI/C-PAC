import scipy.io  
import numpy as np
import nibabel as nib
import os
from scipy import signal
from scipy.sparse import lil_matrix
from scipy.ndimage.filters import gaussian_filter
from numpy import ndarray
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
#Loading the data file, which is now a matfile. this returns a matlab dictionary with variables names as keys and loaded matrices as values.

def qpp_wf(img,nd,window_length,number_randomPermutations,cth,n_itr_threshold,max_itr,path_for_saving):
    # TODOconsolidate this to one function as Behnaz suggested
    """This code is adapted from the paper
       "Quasi-periodic patterns(QP):Large-scale dynamics in resting state fMRI that correlate"\
       with local infraslow electrical activity" Shella Keilholz,D et al.NeuroImage,Volume 84, 1 January 2014."\
       The paper implemnts the algorithms for finding QPP in resting state fMRI using matlab"\
       This project is an attempt to adopt the algorithm in python, and to integrate into C-PAC.
       Input:
       ------
       B: 2D nifti image 
       msk: mask of the 2D nifti image
       nd: number of subjects*number of runs per subject
       wl: window length
       nrp: number of repetitions 
       cth: threshold
       n_itr_th: number of iterations
       mx_itr: maximum number of repetitions 
       pfs: path to save the template, FTP, ITP and iter files
       
       
       Returns:
       -------
       time_course_file: 2D array of time points where QPP is detected in .npy format
       ftp_file: 1D array of Final Time Points in .npy format
       itp_file: 1D array of Final Time points in .npy format
       iter_file: 1D array of iterations in .npy format 
       
       Notes:
       -----
       i) If using a .mat file as an input, save only the image with flag 'v7.0' to make it scipy.io loadmat compatible
       (This functionality will soon be replaced by importing with NifTi format only)
       
       ii) To show the peaks found in the signal, add a show=True boolean values in the "find peaks" command.
       A "True" value plots the peaks that are found in the signal.
       
       Examples:
       --------
       >> python detectqppv.py '/path/to/Data/file.mat'
       'path/to/mask/file/' 30 6 0.2 0.3 1 15 'path/to/save/results/' 6 1
    """
    #define these variables so they make more sense --> didn't change it in the detectqppv function because laziness


    n_timePoints = img.shape[1]
    #shape of time dimension
    n_xaxis = img.shape[0] #shape of x dimensions
    print(img.shape)
    n_tempDim = int(n_timePoints/nd) #number of temporal dimensions
    #use int to prevent floating point errors during initializations
    n_inspect_segment = n_tempDim-window_length+1 #no.of inspectable segments
    nTf = (n_xaxis*window_length)
    #make it a boolean mask - all values with entries greater than zeros will become 1 the rest will be zero
    #no real use of mask anywhere else?
    #msk = np.zeros((n_xaxis,1))
    #msk[(np.sum(abs(img)) > 0)] = 1
    #a = np.where(msk[:,0]==1)
    #img=img[a[0],:]
    #defining 3D arrays here. Each array within the 2D array will finally be a nX*wl shape column vector, which will store the flattened template segment values

    flattened_segment_array = lil_matrix((n_timePoints,n_xaxis*window_length))
    flattened_segment_array_2 = lil_matrix((n_timePoints,n_xaxis*window_length))

    #for each subject*run store the template into the flattened_segment_array array. Instead of using transpose and multiplication, just us dot product of the template square to be stored in bchfn,
    #This step. Presumably is done to maximize the peaks that are found within the arrays(eplained below)
    for i in range(nd):
        for ich in range(n_inspect_segment):
            template=img[:,(i)*n_tempDim+ich:(i)*n_tempDim+window_length+ich]
            #change template from a row vector to a column vector
            template = ndarray.flatten(template)
            # insert the template into the bchfn array (this template will be a 1D array)
            flattened_segment_array[i*n_tempDim+ich] = template
            #normalize
            template=template-np.sum(template)/nTf
            #get dot product
            #template_trans = np.transpose(template)
            temp_dot = np.dot(template,template)
            template_sqrt = np.sqrt(temp_dot)
            template=template/template_sqrt
            #add said template into bchfn
            flattened_segment_array[(i)*n_tempDim+ich] = template
            #removing nan values and making them 0 to prevent further issues in calculations
            A = np.isnan(flattened_segment_array_2)
            flattened_segment_array_2[A] = 0

    #array initialized to later be deleted from the random ITP array
    random_selection_array=np.zeros((nd,window_length-1))
    #filling the sequence with range of numbers from wl+2 to nt
    for i in range(1,nd+1):
        random_selection_array[i-1,:] = range(i*n_tempDim-window_length+2,i*n_tempDim+1)
    #delete instances of ITP from i2x
    initial_timePoints=np.arange(1,n_timePoints+1)
    random_selection_array = ndarray.flatten(random_selection_array)
    initial_timePoints = np.delete(initial_timePoints,i2x-1,0)

    #permute the numbers within ITP

    initial_timePoints = np.random.permutation(initial_timePoints)
    initial_timePoints = initial_timePoints[0:n_randomPermutations]

    #Initialize the time course that will later on be saved
    time_course=np.zeros((n_randomPermutations,n_timePoints))
    #this is the list of final timepoints
    final_timePoints = [[None]]*n_randomPermutations
    #iteration is an array
    iteration = np.zeros(n_randomPermutations)

    for irp in range(n_randomPermutations):
        #initialize a matrix template_holder which will hold the templates
        template_holder=np.zeros(n_timePoints)

        for i in range(nd):
            for ich in range(n_inspect_segment):
                #this confusing flattened_segment_array_2_1 is only so we can finally do dot product
                flattened_segment_array_2_1 =flattened_segment_array_2[initial_timePoints[irp]]
                flattened_segment_array_2_2 =flattened_segment_array_2[i*n_tempDim+ich]
                template_holder[(i)*n_tempDim+ich]= np.dot(flattened_segment_array_2_1,flattened_segment_array_2_2)

        #using MARCUS DEUTRE'S awesome detect_peaks.py function which is a replica of the matlab find peaks function
        #switching off show true until it is necessary, in order to test code.
        peaks= detect_peaks(template_holder,mph=cth[0],mpd=window_length)
                            #show=True)
        #indexes = pu.indexes(c, thresh=c[0])
        #You're deleting the first and last instances of the peaks that are now in the 'peaks' array
        for i in range(nd):
            if i*n_tempDim in peaks:
                peaks = np.delete(peaks,np.where(peaks==(i)*n_tempDim))
            if i*n_tempDim+n_inspect_segment in peaks:
                peaks = np.delete(peaks,np.where(peaks==i*n_tempDim+n_inspect_segment))

        #house three copies of templates (inefficient) which is then used to decide between the correlation coefficient in the next loop
        template_holder_0 = template_holder
        template_holder_00 = template_holder
        template_holder_000 = template_holder
        itr = 1
        #loops over specified number of dimensions
        while itr<=max_itr:
            #applying a gaussian filter to smooth the template file. Default threshold of 0.5
            template_holder = gaussian_filter(template_holder,0.5)
            #setting the threshold value for the local maxima
            if itr<=n_iter_threshold:
                initial_threshold=0
            else:
                initial_threshold=1
            #you can now set the threshold for the local maxima
            threshold=cth[initial_threshold]

            n_signals=np.size(peaks)
            if n_signals<=1:
                break
            template = [peaks[0]]
            for i in range(1,n_signals):
                template=template+flattened_segment_array[peaks[i]]
            template=template/n_signals
            #perform a repeate of the operations in order to find peaks in the template
            #template_trans2=np.transpose(template)
            template=template-np.sum(template)/nTf
            template=template/np.sqrt(np.dot(template,template))
            for i in range(nd):
                for ich in range(n_inspect_segment):
                    template_holder[i*n_tempDim+ich]=np.dot(template,flattened_segment_array[(i)*n_tempDim+ich])
            peaks=detect_peaks(template_holder,mph=cth[1],mpd=window_length)
            for i in range(nd):
                if i * n_tempDim in peaks:
                    peaks = np.delete(peaks, np.where(peaks == (i) * n_tempDim))
                if i * n_tempDim + n_inspect_segment in peaks:
                    peaks = np.delete(peaks, np.where(peaks == i * n_tempDim + n_inspect_segment))
            template_holder_0_norm = (template_holder_0 - np.mean(template_holder_0))/(np.std(template_holder_0))
            #use the correlation coefficient. It returns a matrix and therefore, the first entry of that matrix will be the correlation coefficient value
            if (np.corrcoef(template_holder_0,template_holder)[0,1]>0.9999) or (np.corrcoef(template_holder_00,template_holder)[0,1]>0.9999) or (np.corrcoef(template_holder_000,template_holder)[0,1]>0.9999):
                break

            template_holder_000=template_holder_00
            template_holder_00=template_holder_0
            template_holder_0=template_holder
            itr=itr+1
        if n_signals>1:
            time_course[irp,:]=template_holder
            final_timePoints[irp] = signals.tolist()
            iteration[irp]=itr
    #save everything!!
    plt.plot(template,'b')
    plt.title('Template of QPP(nd=6,wl=30,subjects=7)')
    plt.xlabel('avg of func.data of length WL(30)')
    plt.savefig("{0}/Temple_QPP.png".format(pfs))
    mdict = {}
    mdict["C"] = time_course
    mdict["FTP"] = final_timePoints
    mdict["ITER"] = iteration
    mdict["ITP"] = initial_timePoints
    np.save('{0}/template_file'.format(pfs),template)
    np.save('{0}/time_course_file'.format(pfs),time_course)
    np.save('{0}/ftp_file'.format(pfs),final_timePoints)
    np.save('{0}/iter_file'.format(pfs),iteration)
    np.save('{0}/itp_file'.format(pfs),initial_timePoints)
    #finding the best template, T1 or the  QPP, out of nRP templates the template with the maximum sum of correlation at the supra-threshold local
    #maxima is selected as T1, hence T1 would have higher correlation and more  occurance compared to  other templates

    n_randomPermutations = time_course.shape[0]
    sum_correlation = np.zeros(n_randomPermutations)
    if len(final_timePoints) < n_randomPermutations:
        print("Error: There were no peaks found in your data.Please try to run with a different dataset")
    for i in range(n_randomPermutations):
        if np.any(final_timePoints[i]):
            temp_ftp = final_timePoints[i]
            temp_ftp = [int(x) for x in temp_ftp]
            #check out the list index out of range error that pops up
            sum_correlation[i] =np.sum(time_course[i,temp_ftp])

    is_sum_correlation = np.argsort(sum_correlation)[::-1]
    T1 = is_sum_correlation[0]
    time_course_sum_correlation = time_course[T1,:] #time course sum correlation is an array that contains
    if ftp:
        final_timePoints_1 = final_timePoints[T1]
    else:
        raise Exception("The program will end now, because we could not find any signal correlation and the Final Time Point "
                        " was not found with QPP. Please try running with more data.")
    #basic metrics of T1
    metrics_T1 = np.empty(3)
    final_timePoints_1 = [int(x) for x in final_timePoints_1]
    metrics_T1[0] = np.median(time_course_sum_correlation[final_timePoints_1])
    metrics_T1[1] = np.median(np.diff(final_timePoints_1))
    metrics_T1[2] = len(final_timePoints_1)
    # plots
    # QPP correlation timecourse and metrics

    plt.plot(time_course_sum_correlation,'b')
    plt.plot(final_timePoints_1,time_course_sum_correlation[final_timePoints_1],'g^')
    plt.axis([0,nd*n_tempDim,-1,1])
    plt.xticks(np.arange(n_tempDim,n_temPoints,step=n_tempDim))
    plt.yticks(np.arange(-1,1,step=0.2))
    plt.xlabel('Time points of functional data,TR(s)')
    plt.title('QPP 2D array')
    plt.savefig("{0}/QPP_2D_array.png".format(pfs))

    #% building T1, by averaging segments of B with length 2*WL starting at FTP-WL/2; extra WL/2 at each end is primarily to have die-off effect; it
    #is also used in fine-phase-matching two QPPs when comparing them

    n_final_timePoints = len(final_timePoints_1)

    window_length_starting=round(window_length/2)

    window_length_ending= window_length_starting-np.remainder(window_length,2)

    best_template = np.zeros((n_xaxis,2*window_length))

    for i in range(n_final_timePoints):
        start_time=final_timePoints_1[i]-window_length_starting
        start_time = int(start_time)
        end_time=final_timePoints_1[i]+window_length-1+window_length_ending
        end_time = int(end_time)
        z_starting=None
        zs_flag=False
        if start_time <= 0:
            z_starting=np.zeros((n_xaxis,abs(start_time)+1))
            start_time=1
            zs_flag = True
        z_ending=None
        ze_flag=False
        if end_time>n_timePoints:
            z_ending = np.zeros((n_xaxis,end_time-n_timePoints))
            end_time=n_timePoints
            ze_flag=True
        if zs_flag:
            conct_array = np.concatenate((z_starting,img[:,start_time-1:end_time]),axis=1)
        else:
            conct_array = img[:,start_time-1:end_time]
        if ze_flag:
            conct_array2 = np.concatenate((conct_array,z_ending),axis=1)
        else:
            conct_array2 = conct_array

        best_template = best_template+conct_array2
    best_template=best_template/nFTP

    return img,nd,best_template,time_course_sum_correlation,path_for_saving

def regressqpp(img,nd,best_template,time_course_sum_correlation,path_for_saving):
    #to do: check shape of c in loop
    window_length = np.round(best_template.shape[1]/2)
    window_length_starting = np.round(window_length/2)
    window_length_ending=np.round(window_length/2)+window_length
    best_template_c=best_template_1[0,window_length_starting:window_length_ending]
    best_template_c_new =best_template[:,window_length_starting:window_length_ending]
    n_xaxis = img.shape[0]
    n_timePoints = img.shape[1]
    n_tempDim = n_timePoints/nd
    nTf = (n_xaxis * window_length)
    regressed_img=np.zeros((n_xaxis,n_timePoints))
    for i in range(nd):
        start_time=(i)*n_tempDim
        template_holder = time_course_sum_correlation[start_time:start_time+n_tempDim]
        for ix in range(nd):
            x = np.convolve(template_holder,best_template_c,mode='valid')
            y = img[ix,start_time+window_length-1:start_time+n_tempDims]
            x_dot = np.dot(x,x)
            y_dot = np.dot(x,y)
            #adding a small value to prevent zero division error
            if x_dot == 0:
                x_dot = x_dot + 1e-25
            if y_dot == 0:
                y_dot = y_dot + 1e-25
            beta=y_dot/x_dot
            regressed_img[ix,start_time+window_length-1:start_time+n_tempDim]=y-x*beta

    regressed_timecourse=np.zeros((1,n_timePoints))
    ntf=n_xaxis*window_length
    template=np.array(best_template_c_new.reshape(best_template_c_new.shape[0]*best_template_c_new.shape[1]))
    template=template-np.sum(template)/nTf
    t_dot = np.dot(template,template)
    template=template/np.sqrt(t_dot)
    best_template_n = template.reshape(1,-1)

    for i in range(nd):
        start_time=(i)*n_timePoints
        for ich in range(n_tempDim-window_length):
            template = regressed_img[:,start_time+ich:start_time+ich+window_length]
            template=template.flatten()
            template=template-np.sum(template)/ntf
            if np.any(template==0):
                template=template+1e-25
            template_for_regression = template/np.sqrt(np.dot(template,template))
            regressed_timecourse[:,start_time+ich]=np.dot(best_template_n,template_for_regression)

    C1r_plt = regressed_timecourse.flatten()
    np.save('{0}/regressor file'.format(path_for_saving),regressed_timecourse)
    plt.plot(time_course_sum_correlation,'b')
    plt.plot(C1r_plt,'r')
    plt.axis([0,nd*n_tempDim,-1,1])
    plt.xticks(np.arange(n_tempDim,n_timePoints,step=n_tempDim))
    plt.yticks(np.arange(-1,1,step=0.2))
    plt.title('Time course overlapped with regressor')
    plt.savefig("{0}/Time_course_and_regressor.png".format(path_for_saving))
#    if glassr_360:
#        indu=np.nonzero(np.triu(np.ones(nX),1))

#        indl=np.nonzero(np.tril(np.ones(nX),-1))

#        FCF=np.zeros((nX,nX))

#        bo,nind,p4lb,ylb = RDRG2Y7(B)
#        FC=np.corrcoef(bo)
#        FCF[indu] = FC[indu]
#        plt.imshow(X=FCF, cmap='jet', vmin=-0.6, vmax=0.8)
#        plt.gca().set_aspect('equal',adjustable='box')
#        plt.colorbar
#        plt.grid()
#        plt.show()

#        bo_r, nind_r, p4lb,ylb = RDRG2Y7(Br)
#        FC=np.corrcoef(bo_r)
#        FCF[indl]=FC[indl]

#        plt.imshow(X=FCF,cmap='jet',vmin=-0.6,vmax=0.8)
#        plt.gca().set_aspect('equal',adjustable='box')
#        plt.colorbar
#        plt.grid()
#        plt.show()

#        for i in range(2,7):
#            arr_1=np.array([0,nX])
#            arr_2=np.array([nind[i],nind[i]])
#            plt.plot(arr_1,arr_2,'k')
#            plt.plot(arr_2,arr_1,'k')
#        plt.gca()
#        plt.xticks(p4lb)
#        plt.yticks(p4lb)
#        plt.show()
    return Br, C1r

def RDRG2Y7(bi):
    glassr_file = scipy.io.loadmat('myGlssr360.mat')
    for keys in glassr_file:
        g2y7=glassr_file['G2Y7']
        ylb1=glassr_file['YLB1']
        ylb=glassr_file['YLB']

    g2y7=g2y7.tolist()
    g2y7=[item for sublist in g2y7 for item in sublist]
    n=7
    ind=[]
    #print(g2y7.shape)
    nind=np.zeros(n+1)
    #There are 7 regions in the brain, so we're just trying to mark those regions in the input image B
    #and copy it into another image bo
    #ind has the ordered list of each occurrence of the 1-7 indices
    #nind is the number of times each n (1 to 7) occur in the
    for i in range(n):
        ind.append([p for p,x in enumerate(g2y7) if x==i+1])
        nind[i+1]=len(ind[i])


    #type of ind, change ind to a list of lists
    # each element of each list, will now
    # be int values
#    for x in range(len(ind)):
#        print(ind[x])
    nind=np.cumsum(nind)

    nind=nind.astype(int)

    bo=np.zeros(bi.shape)
    for i in range(n):
        bo[nind[i]:nind[i+1],:]=bi[ind[i],:]
    p4lb=np.zeros(n)
    for i in range(n):
        p4lb[i]=(nind[i+1]+nind[i])/2
    return bo,nind,p4lb,ylb

#if __name__ == '__main__':

#    qpp_wf(d,msk,nd,wl,nrp,cth,n_itr_th,mx_itr,pfs)




