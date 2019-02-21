import scipy.io
import nibabel as nib
import numpy as np





def qpp_wf(img,mask,wl,nrp,cth,n_itr_th,mx_itr,pfs,nsubj,nrn):

    if img.endswith('.mat'):
        D_file = scipy.io.loadmat(img)
        for keys in D_file:
            D = D_file['D']
        D = np.array(D)

    else:
        ##This is the function to import the img into an array object
        ##D_file is now an array object
        ##nib is nibabel package
        D_file = nib.load(img)
        D_img = D_file.dataobj
        D_img=np.array(D_img)

        ##shape of D_img is (61,73,61,7)
        D_img = D_img.reshape(D_img.shape[0]*D_img.shape[1],D_img.shape[2],D_img.shape[3])

        ##D_img is now (4453,61,7)
        D = [[None]*nrn]*nsubj
        #each element of D[i] should be of size (D_img.shape[0],D_img.shape[1]*D_img.shape[2])
        ##initializing D, which is a list of lists
        for i in range(nsubj):
            for j in range(nrn):
                D[i][j] = D_img[:,:,i+j*nsubj]



    if  mask.endswith('.nii'):
        data1 = nib.load(mask)
        msk_img = np.array(data1.dataobj)
        #import msk

        #reshape for masks
        msk_shape = msk_img.shape[:-1]
        m_voxels = np.prod(msk_img.shape[:-1])
        msk = msk_img.reshape(m_voxels,msk_img.shape[-1])

    else:
        msk_file = h5py.File(mask)
        msk_img = msk_file['M']
        msk_img = np.array(msk_img)
        msk_shape = msk_img.shape[:-1]
        m_voxels = np.prod(msk_img.shape[:-1])
        msk = msk_img.reshape(m_voxels,msk_img.shape[-1])

    nx = D[0][0].shape[0]
    nt = D[0][0].shape[1]
    nd = nsubj*nrn
    nt_new = nt * nd
    B = np.zeros((nx,nt_new))
    id =1
    for isbj in range(nsubj):
        for irn in range(nrn):
                B[:,(id-1)*nt:id*nt] = (stats.zscore(D[isbj][irn],axis=1))
                id += 1
    B=np.around(B, decimals=4)
    msk = np.zeros((nx,1))
    msk[(np.sum(abs(B)) > 0)] = 1
    A = np.isnan(B)
    B[A] = 0

    with open('b.txt','w') as f:
        for line in B:
            f.write("%s\n" %line)


    start_time = time.time()
    #generate qpp
    time_course, ftp, itp, iter = qppv(B, msk, nd, wl, nrp, cth, n_itr_th, mx_itr, pfs)
    #choose best template
    C_1,FTP1,Met1 = BSTT(time_course,ftp,nd,B)
    #regress QPP
    T =TBLD2WL(B,wl,FTP1)
    Br, C1r=regressqpp(B, nd, T, C_1,args.glassr_360)
    print("-----%s seconds ----"%(time.time() - start_time))

if __name__ == "__main__":

    #import argparse
    #parser = argparse.ArgumentParser()

    #parser.add_argument("img", type=str,help='Provide the path to the 2D nifti file')

    #parser.add_argument("mask", type=str, help='provide the path to the mask of 2D nifti file')


    #parser.add_argument("wl", type=int,help='provide the length of window you would like to search for the template in')

    #parser.add_argument("nrp", type=int, help='provide the number of random permutations you would like to perform')

    #parser.add_argument("cth", nargs= '+',type=float,help='provide the threshold value, as a list')

    #parser.add_argument("n_itr_th", type=int, help='provide the number of scans contatenated')

    #parser.add_argument("mx_itr", type=int, help='provide the maximum number of iterations')

    #parser.add_argument("pfs", type=str, help='provide the path to the directory you would like to save the files in')

    #parser.add_argument("nsubj", type=int, help='provide the number of subjects')

    #parser.add_argument("nrn",type=int,help='provide the number of runs per subject')

    #parser.add_argument("glassr_360",type=bool,help='If you have data organized with glasser 360 parcels, set this option to true,else set it to False.')


    #args = parser.parse_args()
    qpp_wf(img,mask,wl,nrp,cth,n_itr_th,mx_itr,pfs,nsubj,nrn)



