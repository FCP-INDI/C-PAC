import numpy as np
import nibabel as nb
import os
import sys


dtseries = sys.argv[0] # dtseries file
mask = sys.argv[1] # surface mask file
cortex_file = sys.argv[2] # cortex file 
surf = sys.argv[3] # surface file
scalar_dt =  sys.argv[4] # mean_timeseries
structure_name = sys.argv[5]
surf_reho = sys.argv[6] # Path to reho file

def get_neighbours(faces, vertex_id, depth=1):
    fois = np.where(faces == vertex_id)[0]
    nbrs = list(np.unique(faces[fois]))
    if depth > 1:
        depth -= 1
        new_nbrs = []
        for n in nbrs:
            new_nbrs += get_neighbours(faces, n, depth=depth)
        nbrs = list(np.unique(new_nbrs))
    return nbrs


def ccs_ReHo(dt_file, surf_file):

    #Finding neighbors
    surf = np.array(surf_file.agg_data()) # 2 separate arrays
    faces = np.array(surf_file.agg_data('triangle'))
    vertex = np.array(surf_file.agg_data('pointset'))
    voi = range(0,(vertex.shape)[0])
    nbrs = []
    for voi in voi:
        output = get_neighbours(faces, voi, depth = 1)
        nbrs.append(output)

    #Ting ReHo Calculation
    tsmat = np.array(dt_file.agg_data())
    nsp = (tsmat.shape)[1] #vertex
    ntp = (tsmat.shape)[0] #timepoints

    tmp_ts = []
    cReHo = []
    for i in range(0,nsp):
        tmp_ts = np.squeeze(tsmat[:,i])
        tmp_ts = tmp_ts.reshape(-1,1)
        if np.std(tmp_ts) > 0:
            tmp_nbrs = nbrs[i]
            nbrs_ts = np.squeeze(tsmat[:,tmp_nbrs]) 
            std = np.std(nbrs_ts, axis = 0)
            where = (np.where(std > 0))[0]
            nbrs_ts_arr = nbrs_ts[:, where]
            ts = np.concatenate((tmp_ts, nbrs_ts_arr), axis = 1)
            m = (ts.shape)[1]
            I = np.argsort(ts, axis = 0, kind = 'mergesort')
            R = np.argsort(I, axis = 0, kind = 'mergesort')
            S = (np.sum((np.sum(R, axis = 1)**2))) - (ntp*np.mean(np.sum(R,axis=1))**2)
            F = m*m*(ntp*ntp*ntp-ntp)
            cReHo = np.append(cReHo, (12 * (S/F)))
            cReHo = cReHo.reshape(-1,1)

    return(cReHo)    

cReHo = ccs_ReHo(cortex_file, surf)

## Axes and Header stuff ##
axes = [dtseries.header.get_axis(i) for i in range(dtseries.ndim)]
axes1 = [scalar_dt.header.get_axis(i) for i in range(scalar_dt.ndim)]

time_axis, brain_model_axis = axes #dtseries
time_axis1, brain_axis1 = axes1 #dscalar

# Select the structures you want
structure_names = [structure_name]  # List of structure names

brain_models = [bm for bm in brain_model_axis.iter_structures()
                if bm[0] in structure_names]

new_dataobj = np.concatenate([dtseries.dataobj[0, bm[1]] for bm in brain_models], axis=0)
new_dataobj = np.transpose(new_dataobj.reshape(-1,1))

new_brain_model_axis = sum(
    (bm[2] for bm in brain_models[1:]), brain_models[0][2])

new_cifti = nb.Cifti2Image(new_dataobj,
                       header=(time_axis1, new_brain_model_axis),
                       nifti_header=dtseries.nifti_header)

## Saving image ##                       
img = nb.Cifti2Image(np.transpose(cReHo), header = new_cifti.header, nifti_header=new_cifti.nifti_header)
reho_file = surf_reho
img.to_filename(reho_file)