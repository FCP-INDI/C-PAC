

def calc_compcor_components(data, nComponents, wm_sigs, csf_sigs):
    import scipy.signal as signal
    import numpy as np
    
    wmcsf_sigs = np.vstack((wm_sigs, csf_sigs))

    # filter out any voxels whose variance equals 0
    print 'Removing zero variance components'
    wmcsf_sigs = wmcsf_sigs[wmcsf_sigs.std(1)!=0,:]

    if wmcsf_sigs.shape.count(0):
        print 'No wm or csf signals left after removing those with zero variance'
        raise IndexError
    
    print 'Detrending and centering data'
    Y = signal.detrend(wmcsf_sigs, axis=1, type='linear').T
    Yc = Y - np.tile(Y.mean(0), (Y.shape[0], 1))
    Yc = Yc / np.tile(np.array(Y.std(0)).reshape(1, Y.shape[1]), (Y.shape[0], 1))
    
    print 'Calculating SVD decomposition of Y*Y\''
    U, S, Vh = np.linalg.svd(Yc)
    
    return U[:, :nComponents]


def erode_mask(data):
    import numpy as np

    mask = data != 0
    eroded_mask = np.zeros_like(data, dtype='bool')
    max_x, max_y, max_z = data.shape
    x, y, z = np.where(data != 0)
    for i in range(x.shape[0]):
        if (max_x-1) == x[i] or \
           (max_y-1) == y[i] or \
           (max_z-1) == z[i] or \
           x[i] == 0 or \
           y[i] == 0 or \
           z[i] == 0:
            eroded_mask[x[i], y[i], z[i]] = False
        else:
            eroded_mask[x[i],y[i],z[i]] = mask[x[i], y[i], z[i]] * \
                                          mask[x[i] + 1, y[i], z[i]] * \
                                          mask[x[i], y[i] + 1, z[i]] * \
                                          mask[x[i], y[i], z[i] + 1] * \
                                          mask[x[i] - 1, y[i], z[i]] * \
                                          mask[x[i], y[i] - 1, z[i]] * \
                                          mask[x[i], y[i], z[i] - 1]

    eroded_data = np.zeros_like(data)
    eroded_data[eroded_mask] = data[eroded_mask]
    
    return eroded_data


def find_offending_time_points(thresh_metric, out_file_path="censors.tsv", fd_file_path=None, dvars_file_path=None,
                               fd_threshold=None, dvars_threshold=None, number_of_previous_trs_to_remove=0,
                               number_of_subsequent_trs_to_remove=0):
    """

    Applies criterion in method to find time points whose FD or DVARS (or both) are above threshold

    :param thresh_metric: metric for determining offending time points, either 'FD', 'DVARS', or 'FD+DVARS'. In the last
        case time points will be chosen from the intersection of the FD and DVARS criteria
    :param out_file_path: name of output TSV, which will contain the indices of offending time points list in a
        single column
    :param fd_file_path: path to TSV containing framewise displacement as a single column
    :param dvars_file_path: path to TSV containing DVARS as a single column
    :param fd_threshold: threshold to apply to framewise displacement, can be a value such as 0.2 or a floating
      point multiple of the standard deviation specified as, e.g. '1.5 SD'
    :param dvars_threshold: threshold to apply to DVARS, can be a value such as 0.5 or a floating
      point multiple of the standard deviation specified as, e.g. '1.5 SD'
    :param number_of_previous_trs_to_remove: extent of censorship window before the censor
    :param number_of_subsequent_trs_to_remove: extent of censorship window after the censor

    :return: file path to output file
    """
    import numpy as np
    import os
    import re

    if not thresh_metric or thresh_metric not in ['FD', 'DVARS', 'FD+DVARS']:
        raise ValueError("Improper value for method ({0}), should be one of ['FD', 'DVARS', "
                         "'FD+DVARS']".format(thresh_metric))

    offending_time_points = []
    time_course_len = 0

    if thresh_metric in ['FD', 'FD+DVARS']:

        if not fd_file_path:
            raise ValueError("Method {0} requires the specification of a framewise displacement file path, "
                             "none received".format(thresh_metric))

        if not os.path.isfile(fd_file_path):
            raise ValueError("Framewise displacement file {0} could not be found.".format(fd_file_path))

        if not fd_threshold:
            raise ValueError("Method {0} requires the specification of a framewise displacement threshold, "
                             "none received".format(thresh_metric))

        framewise_displacement = np.loadtxt(fd_file_path)

        time_course_len = framewise_displacement.shape[0]

        try:
            fd_threshold_sd = re.match("([0-9]*\.*[0-9]*)\s*SD", fd_threshold)
            if fd_threshold_sd:
                fd_threshold_sd = float(fd_threshold_sd.groups()[0])
                fd_threshold = framewise_displacement.mean() + fd_threshold_sd * framewise_displacement.std()
            else:
                fd_threshold = float(fd_threshold)
        except:
            print("Could not translate fd_threshold {0} into a meaningful value".format(fd_threshold))
            raise

        offending_time_points = np.where(framewise_displacement > fd_threshold)[0].tolist()

    if thresh_metric in ['DVARS', 'FD+DVARS']:

        if not dvars_file_path:
            raise ValueError("Method {0} requires the specification of a DVARS file path, "
                             "none received".format(thresh_metric))

        if not os.path.isfile(dvars_file_path):
            raise ValueError("DVARS file {0} could not be found.".format(dvars_file_path))

        if not dvars_threshold:
            raise ValueError("Method {0} requires the specification of a DVARS threshold, "
                             "none received".format(thresh_metric))

        dvars = np.loadtxt(dvars_file_path)

        time_course_len = dvars.shape[0]

        try:
            dvars_threshold_sd = re.match("([0-9]*\.*[0-9]*)\s*SD", dvars_threshold)
            if dvars_threshold_sd:
                dvars_threshold_sd = float(dvars_threshold_sd.groups()[0])
                dvars_threshold = dvars.mean() + dvars_threshold_sd * dvars.std()
            else:
                dvars_threshold = float(dvars_threshold)
        except:
            print("Could not translate dvars_threshold {0} into a meaningful value".format(dvars_threshold))
            raise

        if not offending_time_points:
            offending_time_points = np.where(dvars > dvars_threshold)[0].tolist()
        else:
            offending_time_points = list(set(offending_time_points).intersection(
                np.where(dvars > dvars_threshold)[0].tolist()))

    extended_censors = []
    for censor in offending_time_points:
        extended_censors += range((censor - number_of_previous_trs_to_remove),
                                  (censor + number_of_subsequent_trs_to_remove + 1))

    extended_censors = [censor for censor in np.unique(extended_censors) if 0 <= censor < time_course_len]

    censor_vector = np.ones((time_course_len, 1))
    censor_vector[extended_censors] = 0

    out_file_path = os.path.join(os.getcwd(), out_file_path)
    np.savetxt(out_file_path, censor_vector)

    return out_file_path
