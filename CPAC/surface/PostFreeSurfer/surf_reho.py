
def get_neighbours(faces, vertex_id, depth=1):
    import numpy as np
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

    import numpy as np

    #Finding neighbors
    #surf = np.array(surf_file.agg_data()) # 2 separate arrays
    faces = np.array(surf_file.agg_data('triangle'))
    vertex = np.array(surf_file.agg_data('pointset'))
    voi = range(0,(vertex.shape)[0])
    nbrs = []
    for voi in voi:
        output = get_neighbours(faces, voi, depth = 1)
        nbrs.append(output)

    #ReHo Calculation
    tsmat = np.array(dt_file.agg_data())
    nsp = (tsmat.shape)[1] #vertex
    ntp = (tsmat.shape)[0] #timepoints

    tmp_ts = []
    cReHo = np.zeros(nsp)
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
            cReHo[i] = 12 * S / F
            cReHo = cReHo.reshape(-1,1)

    return(cReHo)    

def run_surf_reho(subject, dtseries, mask, cortex_file, \
                surface_file, mean_timeseries, reho_filename, structure_name):

    import os
    import nibabel as nb
    import numpy as np
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    from CPAC.surface.PostFreeSurfer.surf_reho import ccs_ReHo
    from CPAC.surface.PostFreeSurfer.surf_reho import get_neighbours

    dtseries = nb.load(dtseries) # dtseries file
    mask = nb.load(mask) # surface mask file
    cortex_file = nb.load(cortex_file) # cortex file 
    surf_file = nb.load(surface_file) # surface file
    scalar_dt =  nb.load(mean_timeseries) # mean_timeseries
    surf_reho = os.path.join(os.getcwd(), f'{subject}_{reho_filename}')

    cReHo = ccs_ReHo(cortex_file, surf_file)
    structure_names = [structure_name]

    ## Get header information from dtseries and mean timeseries ##
    axes = [dtseries.header.get_axis(i) for i in range(dtseries.ndim)]
    axes1 = [scalar_dt.header.get_axis(i) for i in range(scalar_dt.ndim)]

    time_axis1, brain_axis1 = axes1 #dscalar
    brain_model_axis = nb.cifti2.cifti2_axes.BrainModelAxis.from_surface(
                    vertices=cReHo[:, 0], nvertex=cReHo.shape[0], name=structure_name)

    # Only use data from specific region
    brain_models = [bm for bm in brain_model_axis.iter_structures()
                    if bm[0] in structure_names]

    # Extract data from dtseries for every element in brain models and make into dataobj
    new_dataobj = np.concatenate([dtseries.dataobj[0, bm[1]] for bm in brain_models], axis=0)
    new_dataobj = np.transpose(new_dataobj.reshape(-1,1))

    # Get axis information for new dataobj
    new_brain_model_axis = sum(
        (bm[2] for bm in brain_models[1:]), brain_models[0][2])

    new_cifti = nb.Cifti2Image(new_dataobj,
                        header=(time_axis1, new_brain_model_axis),
                        nifti_header=dtseries.nifti_header)

    ## Saving image ##                       
    img = nb.Cifti2Image(np.transpose(cReHo), header=new_cifti.header, nifti_header=new_cifti.nifti_header)                     
    reho_file = surf_reho
    img.to_filename(reho_file)

    return surf_reho