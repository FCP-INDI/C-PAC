import scipy.io  
import numpy as np
import nibabel as nib
import os
from scipy import signal
from scipy.signal import find_peaks
from scipy.ndimage.filters import gaussian_filter
from numpy import ndarray
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
#Loading the data file, which is now a matfile. this returns a matlab dictionary with variables names as keys and loaded matrices as values.

def qpp_wf(B,msk,nd,wl,nrp,cth,n_itr_th,mx_itr,pfs):
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
    #get parameters of the image shape to decide the
    #shape of the cells,arrays,etc
    nT = B.shape[1] #smaller value
    nX = B.shape[0] #larger value
    #print(nT,nX)
    nt = int(nT/nd) #to prevent floating point errors during initializations
    nch = nt-wl+1
    nTf = (nX*wl)
    #make it a boolean mask - all valies with entries greater than zeros will become 1 the rest will be zero
    #no real use of mask anywhere else?
    msk = np.zeros((nX,1))
    msk[(np.sum(abs(B)) > 0)] = 1
    a = np.where(msk[:,0]==1)
    B=B[a[0],:]
    #defining 3D arrayshere. Each array within the 2D array will finally be a nX*wl shape column vector, which will store the template values
    bchf = np.zeros((nT,nX*wl))
    bchfn = np.zeros((nT,nX*wl))
    #for each subject*run store the template into the bchf array. Instead of using transpose and multiplication, just us dot product of the template square to be stored in bchfn,
    #This step. Presumably is done to maximize the peaks that are found within the arrays(eplained below)
    for i in range(nd):
        for ich in range(nch):
            template=B[:,(i)*nt+ich:(i)*nt+wl+ich]
            #change template from a row vector to a column vector
            template = ndarray.flatten(template)
            # insert the template into the bchfn array (this template will be a 1D array)
            bchf[i*nt+ich] = template
            #normalize
            template=template-np.sum(template)/nTf
            #get dot product
            #template_trans = np.transpose(template)
            temp_dot = np.dot(template,template)
            template_sqrt = np.sqrt(temp_dot)
            template=template/template_sqrt
            #add said template into bchfn
            bchfn[(i)*nt+ich] = template
            #removing nan values and making them 0 to prevent further issues in calculations
            A = np.isnan(bchfn)
            bchfn[A] = 0
        #todo: have to make args.nd and other args.something as just the variable name
    #array initialized to later be deleted from the random ITP array
    i2x=np.zeros((nd,wl-1))
    #filling the sequence with range of numbers from wl+2 to nt
    for i in range(1,nd+1):
        i2x[i-1,:] = range(i*nt-wl+2,i*nt+1)
    #delete instances of ITP from i2x
    itp=np.arange(1,nT+1)
    i2x = ndarray.flatten(i2x)
    itp = np.delete(itp,i2x-1,0)
    #permute the numbers within ITP
    #itp = np.random.permutation(itp)
    #itp = np.random.permutation(itp)
    itp = itp[0:nrp]
    #Initialize the time course that will later on be saved
    time_course=np.zeros((nrp,nT))
    ftp = [[None]]*nrp
    iter = np.zeros(nrp)
    for irp in range(nrp):
        #initialize a matrix c which will hold the templates
        c=np.zeros(nT)
        for i in range(nd):
            for ich in range(nch):
                #bchfn_transpose = np.transpose(bchfn[itp[irp]])
                bchfn_1 =bchfn[itp[irp]]
                bchfn_2 =bchfn[i*nt+ich]
                c[(i)*nt+ich]= np.dot(bchfn_1,bchfn_2)
                #print(c.shape)
        #using MARCUS DEUTRE'S awesome detect_peaks.py function which is a replica of the matlab find peaks function
        #switching off show true until it is necessary, in order to test code.
        peaks= detect_peaks(c,mph=cth[0],mpd=wl)
                            #show=True)
        #indexes = pu.indexes(c, thresh=c[0])
        #You're deleting the first and last instances of the peaks that are now in the 'peaks' array
        for i in range(nd):
            if i*nt in peaks:
                peaks = np.delete(peaks,np.where(peaks==(i)*nt))
            if i*nt+nch in peaks:
                peaks = np.delete(peaks,np.where(peaks==i*nt+nch))
        #house three copies of templates (inefficient) which is then used to decide between the correlation coefficient in the next loop
        c_0 = c
        c_00 = c
        c_000 = c
        itr = 1
        #peaks_size = peaks.size
        #print(peaks)
        #print(peaks.shape)
        #print(peaks.size)
        while itr<=mx_itr:
            c = gaussian_filter(c,0.5)
            if itr<=n_itr_th:
                ith=0
            else:
                ith=1
            th=cth[ith]
            tpsgth=peaks
            n_tpsgth=np.size(tpsgth)
            if n_tpsgth<=1:
                break
            template = bchf[tpsgth[0]]
            for i in range(1,n_tpsgth):
                template=template+bchf[tpsgth[i]]
            template=template/n_tpsgth
            #perform a repeate of the operations in order to find peaks in the template
            #template_trans2=np.transpose(template)
            template=template-np.sum(template)/nTf
            template=template/np.sqrt(np.dot(template,template))
            for i in range(nd):
                for ich in range(nch):
                    c[i*nt+ich]=np.dot(template,bchfn[(i)*nt+ich])
            peaks=detect_peaks(c,mph=cth[1],mpd=wl)
            for i in range(nd):
                if i * nt in peaks:
                    peaks = np.delete(peaks, np.where(peaks == (i) * nt))
                if i * nt + nch in peaks:
                    peaks = np.delete(peaks, np.where(peaks == i * nt + nch))
            c_0_norm = (c_0 - np.mean(c_0))/(np.std(c_0))
            #use the correlation coefficient. It returns a matrix and therefore, the first entry of that matrix will be the correlation coefficient value
            if (np.corrcoef(c_0,c)[0,1]>0.9999) or (np.corrcoef(c_00,c)[0,1]>0.9999) or (np.corrcoef(c_000,c)[0,1]>0.9999):
                break

            c_000=c_00
            c_00=c_0
            c_0=c
            itr=itr+1
        if n_tpsgth>1:
            time_course[irp,:]=c
            ftp[irp] = tpsgth.tolist()
            iter[irp]=itr
    #save everything!!
    plt.plot(template,'b')
    plt.title('Template of QPP(nd=6,wl=30,subjects=7)')
    plt.xlabel('avg of func.data of length WL(30)')
    plt.savefig("{0}/Temple_QPP.png".format(pfs))
    mdict = {}
    mdict["C"] = time_course
    mdict["FTP"] = ftp
    mdict["ITER"] = iter
    mdict["ITP"] = itp
    np.save('{0}/template_file'.format(pfs),template)
    np.save('{0}/time_course_file'.format(pfs),time_course)
    np.save('{0}/ftp_file'.format(pfs),ftp)
    np.save('{0}/iter_file'.format(pfs),iter)
    np.save('{0}/itp_file'.format(pfs),itp)

    return time_course,ftp,itp,iter


def z(x,y):
    x =x[:]-np.mean(x[:])
    y =y[:]-np.mean(y[:])
    if np.norm(x)==0 or np.norm(y)==0:
        z = nan
    else:
        x_trans= np.transpose(X)
        z = (x_trans*y)/norm(x)/norm(y)
    return z

def BSTT(time_course,ftp,nd,B,pfs):
    #load the important files into this

    nT = B.shape[1]  # smaller value
    nX = B.shape[0]  # larger value
    nt = int(nT / nd)  # to prevent floating point errors during initializations
    nRp = time_course.shape[0]
    scmx = np.zeros(nRp)
    if len(ftp) < nRp:
        print("Error: There were no peaks found in your data.Please try to run with a different dataset")
    for i in range(nRp):
        if np.any(ftp[i]):
            p = ftp[i]
            p = [int(x) for x in p]
            #check out the list index out of range error that pops up
            scmx[i] =np.sum(time_course[i,p])

    isscmx = np.argsort(scmx)[::-1]
    T1 = isscmx[0]
    C_1 = time_course[T1,:]
    if ftp:
        FTP1 = ftp[T1]
    else:
        raise Exception("The program will end now, because we could not find any signal correlation and the Final Time Point "
                        " was not found with QPP. Please try running with more data.")
    Met1 = np.empty(3)

    FTP1 = [int(x) for x in FTP1]
    Met1[0] = np.median(C_1[FTP1])
    Met1[1] = np.median(np.diff(FTP1))
    Met1[2] = len(FTP1)
    # plots
    # QPP correlation timecourse and metrics

    plt.plot(C_1,'b')
    plt.plot(FTP1,C_1[FTP1],'g^')
    plt.axis([0,nd*nt,-1,1])
    plt.xticks(np.arange(nt,nT,step=nt))
    plt.yticks(np.arange(-1,1,step=0.2))
    plt.xlabel('Time points of functional data,TR(s)')
    plt.title('QPP 2D array')
    plt.savefig("{0}/QPP_2D_array.png".format(pfs))

    return C_1,FTP1,Met1

def TBLD2WL(B,wl,FTP1,pfs):
    nT = B.shape[1]  # smaller value
    nX = B.shape[0]
    nFTP = len(FTP1)

    WLhs0=round(wl/2)

    WLhe0= WLhs0-np.remainder(wl,2)

    T = np.zeros((nX,2*wl)) #shape is 360,60

    for i in range(nFTP):
        ts=FTP1[i]-WLhs0
        ts = int(ts)
        te=FTP1[i]+wl-1+WLhe0
        te = int(te)
        zs=None
        zs_flag=False
        if ts <= 0:
            zs=np.zeros((nX,abs(ts)+1))
            ts=1
            zs_flag = True
        ze=None
        ze_flag=False
        if te>nT:
            ze = np.zeros((nX,te-nT))
            te=nT
            ze_flag=True
        if zs_flag:
            conct_array = np.concatenate((zs,B[:,ts-1:te]),axis=1)
        else:
            conct_array = B[:,ts-1:te]
        if ze_flag:
            conct_array2 = np.concatenate((conct_array,ze),axis=1)
        else:
            conct_array2 = conct_array

        T = T+conct_array2
    T=T/nFTP

    return T

def regressqpp(B,nd,T1,C_1,pfs):
    #to do: check shape of c in loop
    wl = np.round(T1.shape[1]/2)
    wlhs = np.round(wl/2)
    wlhe=np.round(wl/2)+wl
    T1c=T1[0,wlhs:wlhe]
    T1c_new =T1[:,wlhs:wlhe]
    nX = B.shape[0]
    nT = B.shape[1]
    nt = nT/nd
    nTf = (nX * wl)
    Br=np.zeros((nX,nT))
    for i in range(nd):
        ts=(i)*nt
        c = C_1[ts:ts+nt]
        for ix in range(nd):
            x = np.convolve(c,T1c,mode='valid')
            y = B[ix,ts+wl-1:ts+nt]
            x_dot = np.dot(x,x)
            y_dot = np.dot(x,y)
            #adding a small value to prevent zero division error
            if x_dot == 0:
                x_dot = x_dot + 1e-25
            if y_dot == 0:
                y_dot = y_dot + 1e-25
            beta=y_dot/x_dot
            Br[ix,ts+wl-1:ts+nt]=y-x*beta

    C1r=np.zeros((1,nT))
    ntf=nX*wl
    T=np.array(T1c_new.reshape(T1c_new.shape[0]*T1c_new.shape[1]))
    T=T-np.sum(T)/nTf
    t_dot = np.dot(T,T)
    T=T/np.sqrt(t_dot)
    T1n = T.reshape(1,-1)

    for i in range(nd):
        ts=(i)*nt
        for ich in range(nt-wl):
            T = Br[:,ts+ich:ts+ich+wl]
            T=T.flatten()
            T=T-np.sum(T)/ntf
            if np.any(T==0):
                T=T+1e-25
            bch = T/np.sqrt(np.dot(T,T))
            C1r[:,ts+ich]=np.dot(T1n,bch)

    C1r_plt = C1r.flatten()
    np.save('{0}/regressor file'.format(pfs),C1r)
    plt.plot(C_1,'b')
    plt.plot(C1r_plt,'r')
    plt.axis([0,nd*nt,-1,1])
    plt.xticks(np.arange(nt,nT,step=nt))
    plt.yticks(np.arange(-1,1,step=0.2))
    plt.title('Time course overlapped with regressor')
    plt.savefig("{0}/Time_course_and_regressor.png".format(pfs))
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




