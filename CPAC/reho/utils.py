

def getOpString(mean, std_dev):

    """
    Generate the Operand String to be used in workflow nodes to supply 
    mean and std deviation to alff workflow nodes

    Parameters
    ----------

    mean: string
        mean value in string format
    
    std_dev : string
        std deviation value in string format


    Returns
    -------

    op_string : string


    """

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))

    op_string = str1 + " -mas %s"

    return op_string


def f_kendall(timeseries_matrix):

    """
    Calculates the Kendall's coefficient of concordance for a number of
    time-series in the input matrix

    Parameters
    ----------

    timeseries_matrix : ndarray
        A matrix of ranks of a subset subject's brain voxels

    Returns
    ------- 

    kcc : float
        Kendall's coefficient of concordance on the given input matrix 
  
    """

    import numpy as np
    nk = timeseries_matrix.shape

    n = nk[0]
    k = nk[1]

    sr = np.sum(timeseries_matrix, 1)
    sr_bar = np.mean(sr)

    s = np.sum(np.power(sr, 2)) - n*np.power(sr_bar, 2)

    kcc = 12 *s/np.power(k, 2)/(np.power(n, 3) - n)

    return kcc

def compute_reho(in_file, mask_file, cluster_size):

    """
    Computes the ReHo Map, by computing tied ranks of the timepoints,
    followed by computing Kendall's coefficient concordance(KCC) of a timeseries with its neighbours

    Parameters
    ----------

    in_file : nifti file
        4D EPI File 

    mask_file : nifti file
        Mask of the EPI File(Only Compute ReHo of voxels in the mask)

    cluster_size : integer
        for a brain voxel the number of neighbouring brain voxels to use for KCC.


    Returns
    -------

    out_file : nifti file
        ReHo map of the input EPI image

    """

    import nibabel as nb
    import numpy as np
    import os
    import sys
    from CPAC.reho.utils import f_kendall

    out_file = None

    res_fname = (in_file)
    res_mask_fname = (mask_file)
    CUTNUMBER = 10

#    nvoxel = 27

    if not (cluster_size == 27 or cluster_size == 19 or cluster_size == 7):
        cluster_size = 27

    nvoxel = cluster_size

    res_img = nb.load(res_fname)
    res_mask_img = nb.load(res_mask_fname)


    res_data = res_img.get_data()
    res_mask_data = res_mask_img.get_data()

    print res_data.shape
    (n_x, n_y, n_z, n_t) = res_data.shape

    res_data = np.reshape(res_data, (n_x*n_y*n_z, n_t), order='F').T

    Ranks_res_data = np.tile((np.zeros((1, (res_data.shape)[1]))), [(res_data.shape)[0], 1])
    #Ranks_res_data = np.tile(np.uint16(np.zeros((1, (res_data.shape)[1]))), [(res_data.shape)[0], 1])

    segment_length = np.ceil(float((res_data.shape)[1])/float(CUTNUMBER))

    for icut in range(0, CUTNUMBER):

        segment = None

        if not (icut == (CUTNUMBER - 1)):

            segment = np.array(np.arange(icut * segment_length, (icut+1) * segment_length))

        else:

            segment = np.array(np.arange(icut * segment_length, ((res_data.shape)[1])))

        segment = np.int64(segment[np.newaxis])

        res_data_piece = res_data[:, segment[0]]
        nvoxels_piece = (res_data_piece.shape)[1]

        res_data_sorted = np.sort(res_data_piece, 0, kind='mergesort')
        sort_index = np.argsort(res_data_piece, axis=0, kind='mergesort')

        db = np.diff(res_data_sorted, 1, 0)

        db = db == 0

        sumdb = np.sum(db, 0)

        temp_array = np.array(np.arange(0, n_t))
        temp_array = temp_array[:, np.newaxis]

        sorted_ranks = np.tile(temp_array, [1, nvoxels_piece])

        if np.any(sumdb[:]):

            tie_adjust_index = np.flatnonzero(sumdb)

            for i in range(0, len(tie_adjust_index)):

                ranks = sorted_ranks[:, tie_adjust_index[i]]

                ties = db[:, tie_adjust_index[i]]

                tieloc = np.append(np.flatnonzero(ties), n_t + 2)
                maxties = len(tieloc)
                tiecount = 0

                while(tiecount < maxties -1):
                    tiestart = tieloc[tiecount]
                    ntied = 2
                    while(tieloc[tiecount + 1] == (tieloc[tiecount] + 1)):
                        tiecount += 1
                        ntied += 1

                    ranks[tiestart:tiestart + ntied ] = np.ceil(np.float32(np.sum(ranks[tiestart:tiestart + ntied ]))/np.float32(ntied))
                    tiecount += 1

                sorted_ranks[:, tie_adjust_index[i]] = ranks

        del db, sumdb
        sort_index_base = np.tile(np.multiply(np.arange(0, nvoxels_piece), n_t), [n_t, 1])
        sort_index += sort_index_base
        del sort_index_base

        ranks_piece = np.zeros((n_t, nvoxels_piece))

        ranks_piece = ranks_piece.flatten(order='F')
        sort_index = sort_index.flatten(order='F')
        sorted_ranks = sorted_ranks.flatten(order='F')

        ranks_piece[sort_index] = np.array(sorted_ranks)

        ranks_piece = np.reshape(ranks_piece, (n_t, nvoxels_piece), order='F')

#        ranks_piece = np.uint16(ranks_piece)

        del sort_index, sorted_ranks

        Ranks_res_data[:, segment[0]] = ranks_piece

        sys.stdout.write('.')


    Ranks_res_data = np.reshape(Ranks_res_data, (n_t, n_x, n_y, n_z), order='F')


    K = np.zeros((n_x, n_y, n_z))

    mask_cluster = np.ones((3, 3, 3))

    if nvoxel == 19:
        mask_cluster[0, 0, 0] = 0
        mask_cluster[0, 2, 0] = 0
        mask_cluster[2, 0, 0] = 0
        mask_cluster[2, 2, 0] = 0
        mask_cluster[0, 0, 2] = 0
        mask_cluster[0, 2, 2] = 0
        mask_cluster[2, 0, 2] = 0
        mask_cluster[2, 2, 2] = 0

    elif nvoxel == 7:

        mask_cluster[0, 0, 0] = 0
        mask_cluster[0, 1, 0] = 0
        mask_cluster[0, 2, 0] = 0
        mask_cluster[0, 0, 1] = 0
        mask_cluster[0, 2, 1] = 0
        mask_cluster[0, 0, 2] = 0
        mask_cluster[0, 1, 2] = 0
        mask_cluster[0, 2, 2] = 0
        mask_cluster[1, 0, 0] = 0
        mask_cluster[1, 2, 0] = 0
        mask_cluster[1, 0, 2] = 0
        mask_cluster[1, 2, 2] = 0
        mask_cluster[2, 0, 0] = 0
        mask_cluster[2, 1, 0] = 0
        mask_cluster[2, 2, 0] = 0
        mask_cluster[2, 0, 1] = 0
        mask_cluster[2, 2, 1] = 0
        mask_cluster[2, 0, 2] = 0
        mask_cluster[2, 1, 2] = 0
        mask_cluster[2, 2, 2] = 0


    for i in range(1, n_x - 1):

        for j in range(1, n_y -1):

            for k in range(1, n_z -1):



                block = Ranks_res_data[:, i-1:i+2, j-1:j+2, k-1:k+2]
                mask_block = res_mask_data[i-1:i+2, j-1:j+2, k-1:k+2]


                if not(int(mask_block[1, 1, 1]) == 0):


                    if nvoxel == 19 or nvoxel == 7:

                        mask_block = np.multiply(mask_block, mask_cluster)

                    R_block = np.reshape(block, (block.shape[0], 27), order='F')
                    mask_R_block = R_block[:, np.argwhere(np.reshape(mask_block, (1, 27), order='F') > 0)[:, 1]]

                    K[i, j, k] = f_kendall(mask_R_block)


    img = nb.Nifti1Image(K, header=res_img.get_header(), affine=res_img.get_affine())

    reho_file = os.path.join(os.getcwd(), 'ReHo.nii.gz')

    img.to_filename(reho_file)

    out_file = reho_file

    return out_file


