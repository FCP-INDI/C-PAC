

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


def create_temporal_variance_mask(functional_data_file_path, mask_file_path, output_file_name,
                                  threshold, by_slice=False):
    """
    Create a mask by applying threshold to the temporal variance of 4D nifti file in functional_data_file_path. Only
    non-zero voxels in mask will be considered for inclusion.

    :param functional_data_file_path: 4D nifti file containing functional data
    :param output_file_name: name of 3D nifti file containing created mask, the current directory will be prepended
        to this value, and the mask will be written to the resulting location.
    :param mask_file_path: name of 3D nifti file containing mask to use to restrict the voxels considered by the masking
        operation
    :param threshold: only voxels whose temporal variance meet the threshold criterion will be included in the created
        mask. Appropriate values are:
         - a floating point value, values whose temporal variance is greater than this value will be included in mask
         - a floating point value followed by SD, (1.5SD), values whose temporal variance is greater that 1.5 standard
              deviations of voxels will be included in the mask
         - a floating point value followed by PCT, (2PCT), values whose temporal variance is in the specified percentile
              _from the top_ will be included in the mask. For example 2pct results in the top 2% voxels being included
    :param by_slice: indicates whether threshold criterion should be applied by slice, or to all data, only changes
        result for thresholds expressed in terms of SD or pct
    :return: the full path of the 3D nifti file containing the mask created by this operation.
    """

    import os
    import re
    import nibabel as nb
    import numpy as np

    # begin by verifying the input parameters
    if not (functional_data_file_path and (functional_data_file_path.endswith(".nii") or
                                           functional_data_file_path.endswith(".nii.gz"))):
        raise ValueError("Improper functional file specified ({0}), should be a 4D nifti file.")

    if mask_file_path:
        if not (mask_file_path.endswith(".nii") or mask_file_path.endswith(".nii.gz")):
            raise ValueError("Improper mask file specified ({0}), should be a 3D nifti file.")

    if not output_file_name:
        raise ValueError("Output file name must be specified. Received None")

    if not threshold:
        raise ValueError("Threshold must be specified. Received None")

    # first just assume that we are using a variance cutoff, and then change it if we learn differently, this avoids
    # the else clauses for some of the ifs below
    threshold_method = "VAR"
    threshold_value = threshold

    if isinstance(threshold, str):
        regex_match = re.match("([0-9]*\.*[0-9]*)\s*SD", threshold)

        if regex_match:
            threshold_method = "SD"
            threshold_value = regex_match.groups()[0]
        else:
            regex_match = re.match("([0-9]*\.*[0-9]*)\s*PCT", threshold)

            if regex_match:
                threshold_method = "PCT"
                threshold_value = regex_match.groups()[0]

    try:
        threshold_value = float(threshold_value)
    except:
        print("Error converting threshold value {0} from {1} to a floating point number. The threshold value can "
              "contain SD or PCT for selecting a threshold based on the variance distribution, otherwise it should "
              "be a floating point number.".format(threshold_value, threshold))
        raise

    if threshold_value < 0:
        raise ValueError("Threshold value should be positive, instead of {0}.".format(threshold_value))

    if threshold_method is "PCT" and threshold_value >= 100.0:
        raise ValueError("Percentile should be less than 100, received {0}.".format(threshold_value))

    if not isinstance(by_slice, bool):
        raise ValueError("Parameter by_slice should be a boolean.")

    print("Calculating variance mask from {0}, using {1}, {2} threshold method, and a cutoff of {3}".format(
        functional_data_file_path, mask_file_path, threshold_method, threshold_value))

    # load in functional data and create the variance map
    functional_data_img = nb.load(functional_data_file_path)

    if len(functional_data_img.shape) == 4 and functional_data_img.shape[3] < 3:
        raise ValueError("Functional data used to create mask ({0}) should be 4D and should contain 3 or more "
                         "time points.".format(functional_data_file_path))

    functional_data_variance = functional_data_img.get_data().var(3)

    # if a mask was provided read it in and make sure it is appropriate
    mask_image = []
    mask_data = []

    if mask_file_path:
        mask_image = nb.load(mask_file_path)

        if not(np.all(mask_image.shape == functional_data_img.shape[0:3]) and
               np.all(mask_image.affine == functional_data_img.affine)):
            raise ValueError("Shape and affine of mask image {0} ({1} {2}) should match those of the functional data"
                             " {3} ({4} {5})".format(mask_file_path, mask_image.shape, mask_image.affine,
                                                     functional_data_file_path, functional_data_img.shape,
                                                     functional_data_img.affine))

        # get the data and make sure that it is binary
        mask_data = mask_image.get_data()
        mask_data[mask_data > 0] = 1
        mask_data[mask_data != 1] = 0
    else:
        mask_data = np.uint8(functional_data_variance > 0)

    print("data variance shape {0}".format(functional_data_variance.shape))

    if by_slice is True:
        functional_data_variance = functional_data_variance.reshape((np.prod(functional_data_variance.shape[0:2]),
                                                                     functional_data_variance.shape[2]))
    else:
        functional_data_variance = functional_data_variance.reshape((np.prod(functional_data_variance.shape[0:3]), 1))

    # conform output file and mask to functional data shape
    output_variance_mask = np.zeros(functional_data_variance.shape, dtype=np.uint8)
    mask_data = mask_data.reshape(functional_data_variance.shape)

    for slice_number in range(0, functional_data_variance.shape[1]):

        # make sure that there are some voxels at this slice, if not, move on
        if np.sum(mask_data[:, slice_number]) == 0:
            continue

        if threshold_method is "PCT":
            slice_threshold_value = np.percentile(functional_data_variance[mask_data[:, slice_number] == 1,
                                                                           slice_number], 100.0 - threshold_value)

        elif threshold_method is "SD":
            slice_threshold_value = functional_data_variance[mask_data[:, slice_number] == 1, slice_number].mean() + \
                                    threshold_value * functional_data_variance[mask_data[:, slice_number] == 1,
                                                                               slice_number].std()

        else:
            slice_threshold_value = threshold_value

        output_variance_mask[:, slice_number] = mask_data[:, slice_number] & \
                                                (functional_data_variance[:, slice_number] > slice_threshold_value)

    # make sure that the output mask is the correct shape and format
    output_variance_mask = np.uint8(output_variance_mask)
    output_variance_mask = output_variance_mask.reshape(mask_image.shape)

    # now write it out!
    output_file_path = os.path.join(os.getcwd(), output_file_name)
    output_img = nb.Nifti1Image(output_variance_mask, mask_image.affine)
    output_img.to_filename(output_file_path)

    return output_file_path


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
