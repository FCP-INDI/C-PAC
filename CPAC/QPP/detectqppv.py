import nibabel as nib
import numpy as np
import scipy.io
import h5py
import numpy as np
import os
import nibabel as nib
from scipy import stats
import scipy.io
from CPAC.QPP.QPPv0418 import qpp_wf,regressqpp
import time
import sys
from nilearn.masking import compute_epi_mask

def check_merge_list(merge_list):

    char_1=merge_list[0]
    equality_flag=True
    for chars in merge_list[1:]:
        if chars.shape != char_1.shape:
            equality_flag = False
            break
    if equality_flag == True:
        return True
    else:
        return False
    #merged_empty = np.empty((nsubj, nrn, r_subject.shape[0], r_subject.shape[1]))


def qppv(img,mask,flag_3d_4d,wl,cth,n_itr_th,mx_itr,pfs,nsubj,nrn):

    nrn = int(nrn)
    nsubj=int(nsubj)

    if flag_3d_4d == False:
        mask=nib.load(mask)

        ##This is the function to import the img into an array object
        ##D_file is now an array object
        ##nib is nibabel package
        ##This is the function to import the img into an array object
        ##D_file is now an array object
        ##nib is nibabel package
        D_file = nib.load(img)
        D_img = D_file.dataobj
        D_img = np.array(D_img)

        ##shape of D_img is (61,73,61,7)
        D_img = D_img.reshape(D_img.shape[0] * D_img.shape[1], D_img.shape[2], D_img.shape[3])

        ##D_img is now (4453,61,7)
        D = [[None] * nrn] * nsubj
        # each element of D[i] should be of size (D_img.shape[0],D_img.shape[1]*D_img.shape[2])
        ##initializing D, which is a list of lists
        for i in range(nsubj):
            for j in range(nrn):
                D[i][j] = D_img[:, :, i + j * nsubj]
        ##initializing D, which is a list of lists
        nx = D[0][0].shape[0]
        nt = D[0][0].shape[1]
        print(nx,nt)
        nd = nsubj * nrn
        nrp = nd
        nt_new = nt * nd
        B = np.zeros((nx, nt_new))
        id = 1
        for isbj in range(nsubj):
            for irn in range(nrn):
                B[:, (id - 1) * nt:id * nt] = stats.zscore(D[isbj][irn], axis=1)
                id += 1
        B = np.around(B, decimals=4)
        A = np.isnan(B)
        B[A] = 0
    else:
        if img.endswith('.mat'):
            D_file = scipy.io.loadmat(img)
            for keys in D_file:
                D = D_file['D']
            nx = D[0][0].shape[0]
            nt = D[0][0].shape[1]
            nd = nsubj * nrn
            nt_new = nt * nd
            nrp=nd
            B = np.zeros((nx, nt_new))
            id = 1
            for isbj in range(nsubj):
                for irn in range(nrn):
                    B[:, (id - 1) * nt:id * nt] = (stats.zscore(D[isbj][irn], axis=1))
                    id += 1
            B = np.around(B, decimals=4)
            A = np.isnan(B)
            mask = np.zeros((nx, nt))
            mask[(np.sum(abs(B)) > 0)] = 1

        else:
            sub_img=nib.load(img)
            sub_img = sub_img.dataobj
            sub_img=np.array(sub_img)
            sub_img = sub_img.reshape(sub_img.shape[0] * sub_img.shape[1] * sub_img.shape[2], sub_img.shape[3])

            sub_img=stats.zscore(sub_img,axis=1)
            print(type(sub_img))
            nx = sub_img.shape[0]
            nt = sub_img.shape[1]
            nd = nsubj*nrn
            nrp=nd
            mask = nib.load(mask)
            mask_array = mask.dataobj
            mask = np.array(mask_array)
            print(mask.shape)

    start_time = time.time()
    #generate qpp
    img,nd,best_template,time_course_sum_correlation,path_for_saving=qpp_wf(img, mask, nd, wl, nrp, cth, n_itr_th, mx_itr, pfs)
 
    print("-----%s seconds ----"%(time.time() - start_time))
if __name__ == "__main__":


    qppv(img,mask,flag_3d_4d, wl, cth, n_itr_th, mx_itr, pfs, nsubj, nrn)


