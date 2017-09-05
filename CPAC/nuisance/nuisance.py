import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
from nipype.interfaces import afni
import nuisance_afni_interfaces
from CPAC.nuisance import find_offending_time_points, create_temporal_variance_mask


def mask_summarize_time_course(functional_file_path, mask_file_path, output_file_path, method="DetrendNormMean",
                               mask_vol_index=None, mask_label=None, num_pcs=1):
    """
    Calculates summary time course for voxels specified by a mask. Methods for
    summarizing the the time courses include mean and principle component
    analysis.
                               
    :param functional_file_path: path to nifti file corresponding to the functional data that will be used in the mean 
            calculation.
    :param mask_file_path: a nifti file indicating the voxels that should be included in the summary. The mask may 
            include multiple regions and multiple volumes, which can be specified using the 'mask_vol_index' and 
            'mask_label' parameters. Default is to calculate a separate summary for each volume of the mask file 
            and to combine all non-zero regions of each volume into the same summary. The orientation, voxel sizes, and
            'space' of the mask should match those of the functional data.
    :param output_file_path: full path and name of the output TSV
    :param method: the method used for calculating the summary. suitable values are: 
            'Mean' - the average of voxel indices
            'NormMean' - z-score normalize voxels prior to calculating summary
            'DetrendNormMean' - detrend (degree=2) and normalize voxels prior to calculating summary
            'PCA' - summarize data as top num_pcs from a principle components analysis (i.e. compcor)
            default is 'DetrendNormMean'
    :param mask_vol_index: (int) if the mask contains multiple volumes, index of the volume to use. default behavior is 
            to calculate a separate summary for each volume, which are output to separate columns of the output TSV
    :param mask_label: (list of tuples) if a mask contains multiple labelled regions, the desired regions can be 
            selected by providing a list containing a tuple for each volume to be considered, if more than one label is
            provided, the union of the specified regions will be used 
    :param num_pcs: If PCA is chosen, the number of the largest PCs that should be returned. Default is the first.
    :return: name of the TSV file containing the output.
    """

    import os
    import nibabel as nb
    import numpy as np

    # first check inputs to make sure that they are OK
    if not functional_file_path or \
            (not functional_file_path.endswith(".nii.gz") and not functional_file_path.endswith(".nii")) \
            or not os.path.isfile(functional_file_path):
        raise ValueError("Invalid value for functional_data ({}). Should be the filename to an input nifti "
                         "file".format(functional_file_path))

    if not mask_file_path or (not mask_file_path.endswith(".nii.gz") and not mask_file_path.endswith(".nii")) or \
            not os.path.isfile(mask_file_path):
        raise ValueError("Invalid value for mask ({}). Should be the path to an existing nifti "
                         "file".format(mask_file_path))

    if not output_file_path:
        raise ValueError("Invalid value for output_file_path ({}).".format(output_file_path))

    if method and method not in ['Mean', 'NormMean', 'DetrendNormMean', 'PCA']:
        raise ValueError("Invalid value for method ({}). Should be one of ['Mean', 'NormMean', 'DetrendNormMean',"
                         " 'PCA']".format(method))

    if mask_vol_index and not (isinstance(mask_vol_index, (int, long)) and mask_vol_index > 0):
        raise ValueError("Invalid value for mask_vol_index ({}). Should be an integer > 0".format(mask_vol_index))

    if mask_label:
        if not (isinstance(mask_label, list)):
            raise ValueError("Invalid type for mask_label ({}). Should be a list.".format(type(mask_label), mask_label))

        for mask_label_index, mask_label_val in enumerate(mask_label):
            if not isinstance(mask_label_val, (int, long)) and not isinstance(mask_label_val, tuple) \
                    and not isinstance(mask_label_val, list):
                raise ValueError("Invalid type ({0}) for mask label #{1}: {2}, should be integer, "
                                 "list, or tuple".format(type(mask_label_val), mask_label_index, mask_label_val))

    if num_pcs and not (isinstance(num_pcs, (int, long)) and num_pcs > 0):
        raise ValueError("Invalid value for num_pcs ({0}). Should be an integer > 0".format(num_pcs))

    functional_image = nb.load(functional_file_path)

    if len(functional_image.shape) < 4 or functional_image.shape[3] == 1:
        raise ValueError("Expecting 4D file nifti file, but {0} has shape {1}.".format(functional_file_path,
                                                                                       functional_image.shape))

    functional_image_data = (functional_image.get_data()).reshape((np.prod(functional_image.shape[0:3]),
                                                                   functional_image.shape[3]))

    time_course_length = functional_image.shape[3]

    print("Mask image data comes from {0}".format(mask_file_path))

    # read in the mask data
    mask_image = nb.load(mask_file_path)

    if mask_image.shape[0:3] != functional_image.shape[0:3] or \
            not np.allclose(mask_image.affine, functional_image.affine):
        raise ValueError("Mask ({0}) and functional image ({1}) must be in the same space. Please check the header "
                         "and verify that the shape and transform are the same.".format(mask_file_path,
                                                                                        functional_file_path))

    if len(mask_image.shape) > 3:

        mask_image_data = mask_image.get_data().reshape((np.prod(mask_image.shape[0:3]), mask_image.shape[3]))

        if mask_vol_index and mask_vol_index >= mask_image.shape[3]:
            raise ValueError("Requested mask volume {0}, from file that only contains "
                             "volumes 0 ... {1}".format(mask_vol_index, mask_vol_index-1))

        mask_image_data = mask_image_data[:, [mask_vol_index]]

    else:
        mask_image_data = mask_image.get_data().reshape((np.prod(mask_image.shape[0:3]), 1))

        if mask_vol_index and mask_vol_index != 0:
            print("Requested mask volume {0}, from file that only contains a single volume (0), "
                  "ignoring.".format(mask_vol_index))

    if mask_label and len(mask_label) > 1 and len(mask_label) != mask_image_data.shape[1]:
        raise ValueError("Length of mask labels {0} should either match the number of mask volumes {1}, "
                         "or should be 1".format(len(mask_label), mask_image_data.shape[1]))

    number_of_mask_volumes = mask_image_data.shape[1]

    # the calculated roi summaries will be held in a list of lists
    time_course_summaries = []

    for mask_index in range(0, mask_image_data.shape[1]):

        # extract the time series for the intended voxels
        if not mask_label:
            voxel_indices = np.where(mask_image_data[:, mask_index] > 0)[0].tolist()

            print("{0} nonzero positive voxels in mask {1}".format(len(voxel_indices), mask_file_path))

        else:
            # if we only have one set of mask labels, but many mask volumes,
            # we apply the same labels to each volume, otherwise we use a
            # different set of mask labels for each volume
            mask_label_this_mask_vol = []
            if len(mask_label) == 1:
                mask_label_this_mask_vol = mask_label[0]
            else:
                mask_label_this_mask_vol = mask_label[mask_index]

            if isinstance(mask_label_this_mask_vol, int):
                mask_label_this_mask_vol = [mask_label_this_mask_vol]

            # calculate the union between the different masks
            voxel_indices = []
            for mask_label_val in mask_label_this_mask_vol:
                voxel_indices += np.where(mask_image_data[:, mask_index] == mask_label_val)[0].tolist()

            print("{0} voxels in mask {1} with labels {2}".format(len(voxel_indices), mask_file_path,
                                                                  mask_label_this_mask_vol))

        # make sure that the voxel indices are unique
        voxel_indices = np.unique(voxel_indices)

        if len(voxel_indices) == 0:
            raise ValueError("Time series extraction failed, no voxels in mask {0} match label(s) {1}".format(
                mask_file_path, mask_label))

        time_courses = functional_image_data[voxel_indices, :]

        # exclude time courses with zero variance to avoid wasting our time on
        # meaningless data, and to avoid NaNs
        time_courses = time_courses[np.isclose(time_courses.var(1), 0.0) == False, :]

        if time_courses.shape[0] == 0:
            raise ValueError("None of the {0} in-mask voxels have non-zero variance time"
                             " courses.".format(len(voxel_indices)))

        print("{0} voxels survived variance filter".format(time_courses.shape[0]))

        # Make sure we have a 2D array, even if it only contains a single time
        # course
        if len(time_courses.shape) == 1:
            time_courses = time_courses.reshape((time_courses.shape[0], 1))

        if method in ["PCA"]:

            # compcor begins with linearly detrending the columns of the
            # matrix
            def linear_detrend_columns(image_array_2d):
                """
                perform quadratic detrending on each row of a 2D numpy array representing a functional image

                :param image_array_2d: 2D numpy array containing functional data to be processed
                :return: 2D numpy array of residuals
                """
                column_len = image_array_2d.shape[0]
                polynomial_design_matrix = np.array([range(0, column_len), [1] * column_len])
                polynomial_coefficients = np.polyfit(range(0, column_len), image_array_2d, deg=1)
                return image_array_2d - polynomial_coefficients.transpose().dot(polynomial_design_matrix).transpose()

            time_courses = linear_detrend_columns(time_courses)

        # normalize data as requested
        if method in ["PCA", "DetrendNormMean"]:

            def quadratic_detrend_rows(image_array_2d):
                """
                perform quadratic detrending on each row of a 2D numpy array representing a functional image
                
                :param image_array_2d: 2D numpy array containing functional data to be processed
                :return: 2D numpy array of residuals
                """

                print("2D image shape {0} {1}".format(image_array_2d.shape, len(image_array_2d.shape)))

                row_len = image_array_2d.shape[1]
                polynomial_design_matrix = np.array([[x * x for x in range(0, row_len)], range(0, row_len),
                                                     [1]*row_len])
                polynomial_coefficients = np.polyfit(range(0, row_len), image_array_2d.transpose(), deg=2)
                return image_array_2d - polynomial_coefficients.transpose().dot(polynomial_design_matrix)

            time_courses = quadratic_detrend_rows(time_courses)

        if method in ["PCA", "DetrendNormMean", "NormMean"]:
            time_courses = time_courses - np.tile(time_courses.mean(1).reshape(time_courses.shape[0], 1),
                                                  (1, time_courses.shape[1]))
            time_courses = time_courses / np.tile(time_courses.std(1).reshape(time_courses.shape[0], 1),
                                                  (1, time_courses.shape[1]))

        # now summarise
        if method in ["DetrendNormMean", "NormMean", "Mean"]:
            time_course_summaries.append(time_courses.mean(0))

        elif method in ["PCA"]:
            [u, s, v] = np.linalg.svd(time_courses, full_matrices=False)
            time_course_summaries.append(v[0:num_pcs, :])

    if len(time_course_summaries) != number_of_mask_volumes:
        raise ValueError("Expected {0} summaries, one for mask volume, "
                         "but received {1}".format(number_of_mask_volumes,
                                                   len(time_course_summaries)))

    output_file_path = os.path.join(os.getcwd(), output_file_path)
    with open(output_file_path, "w") as ofd:

        for time_point in range(0, time_course_length):

            row_values = []

            if time_point == 0:

                # write out the header information
                ofd.write("# CPAC version {0}\n".format("1.010"))
                ofd.write("# Time courses extracted from {0} using mask {1} and "
                          "method {2}\n".format(functional_file_path, mask_file_path, method))

                for time_course_summary_index, time_course_summary in enumerate(time_course_summaries):
                    if method in ["PCA"]:
                        if time_course_summary.shape[0] != num_pcs:
                            raise ValueError("Time course summary {0} expected {1} pcs, but received {2}".format(
                                time_course_summary_index, num_pcs, time_course_summary.shape[0]))
                        row_values += ["mask#{0}_{1}#{2}".format(time_course_summary_index, method, pc)
                                       for pc in range(0, num_pcs)]
                    else:
                        row_values += ["mask#{0}_{1}".format(time_course_summary_index, method)]

                ofd.write("#"+"\t".join(row_values)+"\n")

                row_values = []

            for time_course_summary_index, time_course_summary in enumerate(time_course_summaries):
                if method in ["PCA"]:
                    row_values += ["{0}".format(time_course_summary[pc, time_point]) for pc in range(0, num_pcs)]
                else:
                    row_values += ["{0}".format(time_course_summary[time_point])]

            ofd.write("\t".join(row_values) + "\n")

    return output_file_path


def gather_nuisance(functional_file_path, selector, output_file_path, grey_matter_summary_file_path=None,
                    white_matter_summary_file_path=None, csf_summary_file_path=None, acompcorr_file_path=None,
                    tcompcorr_file_path=None, global_summary_file_path=None, motion_parameters_file_path=None,
                    dvars_file_path=None, framewise_displacement_file_path=None, censor_file_path=None):
    """
    Gathers the various nuisance regressors together into a single tab
    separated values file that is an appropriate input into 3dTproject
    
    :param selector: Dictionary that indicates which nuisance regressors should be included in the model, along with
        parameters that indicate how the regressors should be constructed from the inputs. Example selector:
            selector = {'aCompCorr' : None | {num_pcs = <number of components to retain>, 
                                            tissues = 'WM' | 'CSF' | 'WM+CSF',
                                            include_delayed = True | False, 
                                            include_squared = True | False,
                                            include_delayed_squared = True | False},
                        'WhiteMatter' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       include_delayed = True | False, 
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'Ventricles' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       include_delayed = True | False, 
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'GreyMatter' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       include_delayed = True | False, 
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'GlobalSignal' : None | {summary_method = 'PC', 'Mean', 'NormMean' or 'DetrendNormMean',
                                           num_pcs = <number of components to retain>,
                                           include_delayed = True | False, 
                                           include_squared = True | False,
                                           include_delayed_squared = True | False},
                        'Motion' : None | {include_delayed = True | False, 
                                           include_squared = True | False,
                                           include_delayed_squared = True | False},
                        'DVARS' : None | {include_delayed = True | False,
                                          include_squared = True | False,
                                          include_delayed_squared = True | False},
                        'FD' : None | {include_delayed = True | False,
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'SpikeRegression' : None | {number_of_previous_trs_to_remove = True | False,
                                                     number_of_subsequent_trs_to_remove = True | False}
                        }
                        
    :param functional_file_path: path to file that the regressors are being calculated for, is used to calculate
           the length of the regressors for error checking and in particular for calculating spike regressors
    :param output_file_path: path to output TSV that will contain the various nuisance regressors as columns
    :param grey_matter_summary_file_path: path to TSV that includes summary of grey matter time courses, e.g. output of 
        mask_summarize_time_course
    :param white_matter_summary_file_path: path to TSV that includes summary of white matter time courses, e.g. output 
        of mask_summarize_time_course
    :param csf_summary_file_path: path to TSV that includes summary of csf time courses, e.g. output 
        of mask_summarize_time_course
    :param acompcorr_file_path: path to TSV that includes acompcorr time courses, e.g. output
        of mask_summarize_time_course
    :param tcompcorr_file_path: path to TSV that includes tcompcorr time courses, e.g. output
        of mask_summarize_time_course
    :param global_summary_file_path: path to TSV that includes summary of global time courses, e.g. output 
        of mask_summarize_time_course 
    :param motion_parameters_file_path: path to TSV that includes motion parameters
    :param dvars_file_path: path to TSV that includes DVARS as a single column
    :param framewise_displacement_file_path: path to TSV that includes framewise displacement as a single column
    :param censor_file_path: path to TSV with a single column with 1's for indices that should be retained and 0's
              for indices that should be censored
    :return: 
    """

    import os
    import numpy as np
    import nibabel as nb

    # make a mapping for the input files to make things easier loop over
    input_files_dict = {'aCompCorr': acompcorr_file_path,
                        'tCompCorr': tcompcorr_file_path,
                        'Ventricles': csf_summary_file_path,
                        'GlobalSignal': global_summary_file_path,
                        'GreyMatter': grey_matter_summary_file_path,
                        'WhiteMatter': white_matter_summary_file_path,
                        'Motion': motion_parameters_file_path,
                        'DVARS': dvars_file_path,
                        'FD': framewise_displacement_file_path,
                        'Censor': censor_file_path}

    # some very basic error checking on the primary arguments

    if not functional_file_path or not os.path.isfile(functional_file_path) or \
            (not functional_file_path.endswith(".nii") and not functional_file_path.endswith(".nii.gz")):
        raise ValueError("Invalid value for input_file ({}). Should be a nifti file "
                         "and should exist".format(functional_file_path))

    if not isinstance(selector, dict):
        raise ValueError("Invalid type for selector {0}, expecting dict".format(type(selector)))

    # now go through and put'er together
    column_names = []
    nuisance_regressors = []

    functional_image = nb.load(functional_file_path)

    if len(functional_image.shape) < 4 or functional_image.shape[3] < 2:
        raise ValueError("Invalid input_file ({}). Expected 4D file.".format(functional_file_path))

    regressor_length = functional_image.shape[3]

    for regressor_type in ['aCompCorr', 'tCompCorr', 'Ventricles', 'GlobalSignal', 'GreyMatter', 'WhiteMatter',
                           'Motion', 'DVARS', 'FD']:
        if regressor_type in selector and selector[regressor_type]:
            if not input_files_dict[regressor_type]:
                raise ValueError('Regressor type {0} specified in selector but the corresponding '
                                 'file was not found!'.format(regressor_type))

            try:
                regressors = np.loadtxt(input_files_dict[regressor_type])
            except:
                print("Could not read regressor {0} from {1}.".format(regressor_type, input_files_dict[regressor_type]))
                raise

            if regressors.shape[0] != regressor_length:
                raise ValueError("Number of time points in {0} ({1}) is inconsistent with "
                                 "length of functional file {2} ({3})".format(input_files_dict[regressor_type],
                                                                              regressors.shape[0], functional_file_path,
                                                                              regressor_length))

            if regressor_type is "Motion":
                num_regressors = 6
            elif not selector[regressor_type]['num_pcs']:
                num_regressors = 1
            else:
                num_regressors = selector[regressor_type]['num_pcs']

            if len(regressors.shape) == 1:
                regressors = regressors.reshape((regressors.shape[0], 1))

            if regressors.shape[1] != num_regressors:
                raise ValueError('Expecting {0} regressors for {1}, but only found {2} in '
                                 'file {3}.'.format(selector[regressor_type]['num_pcs'],
                                                    regressor_type,
                                                    regressors.shape[1],
                                                    input_files_dict[regressor_type]))

            summary_method = ''
            if regressor_type in ['aCompCorr', 'tCompCorr', 'Ventricles', 'GlobalSignal', 'GreyMatter', 'WhiteMatter']:
                # sanity check on method
                if 'summary_method' not in selector[regressor_type] or \
                        not selector[regressor_type]['summary_method']:
                    raise ValueError('Summarization method for {0} not specified. Should be one of: "PCA", "Mean", '
                                     '"NormMean" or "DetrendNormMean".'.format(regressor_type))

                if selector[regressor_type]['summary_method'] not in ["PCA", "Mean", "NormMean", "DetrendNormMean"]:
                    raise ValueError('Invalid summarization method for {0}. {1} should be one of: "PCA", "Mean", '
                                     '"NormMean" or "DetrendNormMean".'.format(regressor_type,
                                                                               selector[regressor_type][
                                                                                   'summary_method']))
                summary_method = selector[regressor_type]['summary_method']

            # if other parameters are missing, assume they are False by
            # default
            for regressor_derivative in ['include_delayed', 'include_squared', 'include_delayed_squared']:
                if regressor_derivative not in selector[regressor_type] or \
                        not selector[regressor_type][regressor_derivative]:

                    selector[regressor_type][regressor_derivative] = False

            motion_labels = ["RotY", "RotX", "RotZ", "Y", "X", "Z"]

            # add in the regressors, making sure to also add in the column
            # name
            for regressor_index in range(0, regressors.shape[1]):
                if regressor_type is "Motion":
                    regressor_name = motion_labels[regressor_index]
                else:
                    regressor_name = "{0}{1}{2}".format(regressor_type, summary_method, regressor_index)

                column_names.append(regressor_name)
                nuisance_regressors.append(regressors[:, regressor_index])

                if selector[regressor_type]['include_delayed']:
                    column_names.append("{0}Back".format(regressor_name))
                    nuisance_regressors.append(np.append(regressors[1:, regressor_index], [0.0]))

                if selector[regressor_type]['include_squared']:
                    column_names.append("{0}Sq".format(regressor_name))
                    nuisance_regressors.append(np.square(regressors[:, regressor_index]))

                if selector[regressor_type]['include_delayed_squared']:
                    column_names.append("{0}BackSq".format(regressor_name))
                    nuisance_regressors.append(np.square(np.append(regressors[1:, regressor_index], [0.0])))

    print("Finished tissue and motion regressors, moving on to Censor")

    # now we add in the spike regressors
    regressor_type = 'Censor'

    if regressor_type in selector and 'censor_method' in selector[regressor_type] and \
            selector[regressor_type]['censor_method'] is 'SpikeRegression':

        if not input_files_dict[regressor_type]:
            raise ValueError('Regressor type {0} specified in selector but the corresponding '
                             'file was not found!'.format(regressor_type))

        try:
            censor_vector = np.loadtxt(input_files_dict[regressor_type])

        except:
            print("Could not read regressor {0} from {1}.".format(regressor_type, input_files_dict[regressor_type]))
            raise

        if len(censor_vector.shape) > 1 and censor_vector.shape[1] > 1:
            raise ValueError("Invalid format for censor file {0}, should be a single column containing "
                             "1s for volumes to keep and 0s for volumes to censor. {0} has {1} columns.")

        # get the indices
        censor_indices = np.where(censor_vector == 0)[0]

        # first check to make sure that we have some indices
        if len(censor_indices) > 0:

            # if number_of_previous_trs_to_remove and number_of_subsequent_trs_to_remove are not set, assume they should
            # be zero
            if not ('number_of_previous_trs_to_remove' in selector[regressor_type] and
                        selector[regressor_type]['number_of_previous_trs_to_remove']):
                selector[regressor_type]['number_of_previous_trs_to_remove'] = 0

            if not ('number_of_subsequent_trs_to_remove' in selector[regressor_type] and
                    selector[regressor_type]['number_of_subsequent_trs_to_remove']):
                selector[regressor_type]['number_of_subsequent_trs_to_remove'] = 0

            for censor_list_index, censor_index in enumerate(censor_indices):
                if censor_index >= regressor_length:
                    raise ValueError("Censor #{0}: {1} from {2} is out of range, calculated"
                                     " regressor length is {3}".format(censor_list_index, censor_index,
                                                                       input_files_dict[regressor_type],
                                                                       regressor_length))
                censor_begin_index = censor_index - selector[regressor_type]['number_of_previous_trs_to_remove']
                censor_end_index = censor_index + selector[regressor_type]['number_of_subsequent_trs_to_remove']

                # enforce boundary conditions
                if censor_begin_index < 0:
                    censor_begin_index = 0

                if censor_end_index >= regressor_length:
                    censor_end_index = regressor_length - 1

                spike_regressor = np.zeros((regressor_length, 1))
                spike_regressor[censor_begin_index:censor_end_index+1, 0] = 1

                column_names.append("SpikeRegression{0}".format(censor_index))
                nuisance_regressors.append(spike_regressor.flatten())

        else:
            print("Censor file {0} was empty, spike regressors will not be included in"
                  " nuisance model.".format(input_files_dict[regressor_type]))

    # finally write out regressor file
    output_file_path = os.path.join(os.getcwd(), output_file_path)
    with open(output_file_path, "w") as ofd:

        for time_point in range(0, regressor_length):

            row_values = []

            if time_point == 0:
                # write out the header information
                ofd.write("# CPAC version {0}\n".format("1.010"))
                ofd.write("# Nuisance regressors\n")
                ofd.write("#"+"\t".join(column_names)+"\n")

            for regressor_index, regressor in enumerate(nuisance_regressors):
                row_values += ["{0}".format(regressor[time_point])]

            ofd.write("\t".join(row_values) + "\n")

    return output_file_path


def create_nuisance_workflow(use_ants, selector, name='nuisance'):
    """
    Workflow for the removal of various signals considered to be noise from
    resting state fMRI data. The residual signals for linear regression
    denoising is performed in a single model.  Therefore the residual
    time-series will be orthogonal to all signals.

    Parameters
    ---------- 
    :param use_ants: flag indicating whether FNIRT or ANTS is used
    :param selector: Dictionary containing configuration parameters for nuisance regression
            selector = {'Anaticor' : None | {radius = <radius in mm>},
                        'aCompCor' : None | {num_pcs = <number of components to retain>,
                                            tissues = 'WM' | 'CSF' | 'WM+CSF',
                                            include_delayed = True | False,
                                            include_squared = True | False,
                                            include_delayed_squared = True | False},
                        'tCompCor' : None | {num_pcs = <number of components to retain>,
                                            threshold = <floating point number = cutoff as raw variance value,
                                                         floating point number followed by SD (ex. 1.5SD) = mean + a
                                                             multiple of the SD,
                                                         floating point number followed by PCT (ex. 2PCT) = percentile
                                                             from the top (ex is top 2%),
                                            by_slice = boolean, whether or not the threshold criterion should be applied
                                                       by slice or across the entire volume, makes most sense for SD or
                                                       PCT
                                            tissues = 'WM' | 'CSF' | 'WM+CSF',
                                            include_delayed = True | False,
                                            include_squared = True | False,
                                            include_delayed_squared = True | False},
                        'WhiteMatter' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       include_delayed = True | False, 
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'Ventricles' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       include_delayed = True | False, 
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'GreyMatter' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       include_delayed = True | False, 
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'GlobalSignal' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                           num_pcs = <number of components to retain>,
                                           include_delayed = True | False, 
                                           include_squared = True | False,
                                           include_delayed_squared = True | False},
                        'Motion' : None | {include_delayed = True | False, 
                                           include_squared = True | False,
                                           include_delayed_squared = True | False},
                        'Censor' : None | { thresh_metric = 'RMSD','DVARS', or 'RMSD+DVARS',
                                            threshold = <threshold to be applied to metric, if using 
                                              RMSD+DVARS, this should be a tuple (RMSD thresh, DVARS thresh)>,
                                            number_of_previous_trs_to_remove = True | False,
                                            number_of_subsequent_trs_to_remove = True | False,
                                            method = 'Kill', 'Zero', 'Interpolate', 'SpikeRegression'},
                        'PolyOrt' : None | { degree = <polynomial degree up to which will be removed, e.g. 2 means 
                                                       constant + linear + quadratic, practically that is probably,
                                                       the most that will be need esp. if band pass filtering>},
                        'Bandpass' : None | { bottom_frequency = <frequency in hertz of the highpass part of the pass 
                                                                  band, frequencies below this will be removed>,
                                              top_frequency = <frequency in hertz of the lowpass part of the pass 
                                                               band, frequencies above this will be removed>},
                        }
    :param name: Name of the workflow, defaults to 'nuisance'
    :return: nuisance : nipype.pipeline.engine.Workflow
        Nuisance workflow.
        
    Notes
    -----

    Workflow Inputs::

        inputspec.functional_file_path : string (nifti file)
            Path to realigned and motion corrected functional image (nifti) file.
        inputspec.wm_mask_file_path : string (nifti file)
            Corresponding white matter mask.
        inputspec.csf_mask_file_path : string (nifti file)
            Corresponding cerebral spinal fluid mask.
        inputspec.gm_mask_file_path : string (nifti file)
            Corresponding grey matter mask.
        inputspec.functional_brain_mask_file_path : string (nifti file)
            Whole brain mask corresponding to the functional data.
        inputspec.mni_to_anat_linear_xfm_file_path: string (nifti file)
            FLIRT Linear MNI to Anat transform
        inputspec.anat_to_mni_initial_xfm_file_path: string (nifti file)
            ANTS initial transform from anat to MNI
        inputspec.anat_to_mni_rigid_xfm_file_path: string (nifti file)
            ANTS rigid (6 parameter, no scaling) transform from anat to MNI
        inputspec.anat_to_mni_affine_xfm_file_path: string (nifti file)
            ANTS affine (13 parameter, scales and shears) transform from anat to MNI
        inputspec.func_to_anat_linear_xfm_file_path: string (nifti file)
            FLIRT Linear Transform between functional and anatomical spaces 
        inputspec.brain_template_file_path: string (nifti file)
            Template used to define MNI space
        inputspec.mni_to_anat_linear_xfm_file_path : string (nifti file)
            Corresponding MNI to anatomical linear transformation 
        inputspec.func_to_anat_linear_xfm_file_path : string (nifti file)
            Corresponding EPI to anatomical linear transformation
        inputspec.lat_ventricles_mask_file_path : string (nifti file)
            Mask of lateral ventricles calculated from the Harvard Oxford Atlas.
        inputspec.motion_parameter_file_path : string (text file)
            Corresponding rigid-body motion parameters.  Matrix in the file should be of shape 
            (`T`, `R`), `T` time points and `R` motion parameters.
        inputspec.fd_file_path : string (text file)
            Framewise displacement calculated from the motion parameters.
        inputspec.dvars_file_path : string (text file)
            DVARS calculated from the functional data.
        inputspec.selector : dictionary

    Workflow Outputs::

        outputspec.residual_file_path : string (nifti file)
            Path of residual file in nifti format
        outputspec.regressors_file_path : string (TSV file)
            Path of TSV file of regressors used. Column name indicates the regressors included .

    Nuisance Procedure:

    1. Compute nuisance regressors based on input selections.
    2. Calculate residuals with respect to these nuisance regressors in a
       single model for every voxel.

    Workflow Graph:

    .. image:: ../images/nuisance.dot.png
        :width: 500

    Detailed Workflow Graph:

    .. image:: ../images/nuisance_detailed.dot.png
        :width: 500    

    """

    nuisance = pe.Workflow(name=name)

    # should these be created with traits for stricter typing?
    inputspec = pe.Node(util.IdentityInterface(fields=['functional_file_path',
                                                       'wm_mask_file_path',
                                                       'csf_mask_file_path',
                                                       'gm_mask_file_path',
                                                       'functional_brain_mask_file_path',
                                                       'mni_to_anat_linear_xfm_file_path',
                                                       'anat_to_mni_initial_xfm_file_path',
                                                       'anat_to_mni_rigid_xfm_file_path',
                                                       'anat_to_mni_affine_xfm_file_path',
                                                       'func_to_anat_linear_xfm_file_path',
                                                       'lat_ventricles_mask_file_path',
                                                       'motion_parameters_file_path',
                                                       'fd_file_path',
                                                       'dvars_file_path',
                                                       'brain_template_file_path',
                                                       'selector']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['residual_file_path',
                                                        'regressors_file_path']),
                         name='outputspec')

    # Resample the white matter mask to 2mm space in subject space
    wm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='wm_anat_to_2mm_flirt_applyxfm')
    wm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    wm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'wm_mask_file_path', wm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'wm_mask_file_path', wm_anat_to_2mm, 'reference')

    # erode the white matter mask by 1 voxel to avoid overlap with grey matter
    wm_anat_2mm_erode = pe.Node(interface=afni.Calc(), name='wm_anat_2mm_erode')
    wm_anat_2mm_erode.inputs.args = '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    wm_anat_2mm_erode.inputs.expr = 'a*(1-amongst(0,b,c,d,e,f,g))'
    wm_anat_2mm_erode.inputs.outputtype = 'NIFTI_GZ'
    wm_anat_2mm_erode.inputs.out_file = 'wm_mask_2mm_eroded.nii.gz'

    nuisance.connect(wm_anat_to_2mm, 'out_file', wm_anat_2mm_erode, 'in_file_a')

    # Resample grey matter masks from 1mm to 2mm, but remaining in subject
    # space
    gm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='gm_anat_to_2mm_flirt_applyxfm')
    gm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    gm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'gm_mask_file_path', gm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'gm_mask_file_path', gm_anat_to_2mm, 'reference')

    # erode the grey matter mask by 1 voxel to avoid overlap with other
    # tissues
    gm_anat_2mm_erode = pe.Node(interface=afni.Calc(), name='gm_anat_2mm_erode')
    gm_anat_2mm_erode.inputs.args = '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    gm_anat_2mm_erode.inputs.expr = 'a*(1-amongst(0,b,c,d,e,f,g))'
    gm_anat_2mm_erode.inputs.outputtype = 'NIFTI_GZ'
    gm_anat_2mm_erode.inputs.out_file = 'gm_mask_2mm_eroded.nii.gz'

    nuisance.connect(gm_anat_to_2mm, 'out_file', gm_anat_2mm_erode, 'in_file_a')

    # multiply the CSF mask by a lateral ventricles to reduce it to the CSF in
    # the lateral ventricles

    # first transform the lat ventricles to match the CSF
    if use_ants is True:
        # perform the transform using ANTS
        collect_linear_transforms = pe.Node(util.Merge(3), name='ho_mni_to_2mm_ants_collect_linear_transforms')

        nuisance.connect(inputspec, 'anat_to_mni_initial_xfm_file_path', collect_linear_transforms, 'in1')
        nuisance.connect(inputspec, 'anat_to_mni_rigid_xfm_file_path', collect_linear_transforms, 'in2')
        nuisance.connect(inputspec, 'anat_to_mni_affine_xfm_file_path', collect_linear_transforms, 'in3')

        lat_ven_mni_to_anat = pe.Node(interface=ants.ApplyTransforms(), name='lat_ven_mni_to_anat_ants')
        lat_ven_mni_to_anat.inputs.invert_transform_flags = [True, True, True]
        lat_ven_mni_to_anat.inputs.interpolation = 'NearestNeighbor'
        lat_ven_mni_to_anat.inputs.dimension = 3

        nuisance.connect(collect_linear_transforms, 'out', lat_ven_mni_to_anat, 'transforms')
        nuisance.connect(inputspec, 'lat_ventricles_mask_file_path', lat_ven_mni_to_anat, 'input_image')
        nuisance.connect(inputspec, 'csf_mask_file_path', lat_ven_mni_to_anat, 'reference_image')

    else:
        # perform the transform using FLIRT
        lat_ven_mni_to_anat = pe.Node(interface=fsl.FLIRT(), name='lat_ven_mni_to_anat_flirt')
        lat_ven_mni_to_anat.inputs.args = '-applyisoxfm 2'
        lat_ven_mni_to_anat.inputs.interp = 'nearestneighbour'

        nuisance.connect(inputspec, 'mni_to_anat_linear_xfm', lat_ven_mni_to_anat, 'in_matrix_file')
        nuisance.connect(inputspec, 'lat_ventricles_mask_file_path', lat_ven_mni_to_anat, 'in_file')
        nuisance.connect(inputspec, 'csf_mask_file_path', lat_ven_mni_to_anat, 'reference')

    # reduce CSF mask to the lateral ventricles
    mask_cfs_with_lat_ven = pe.Node(interface=afni.Calc(), name='mask_cfs_with_lat_ven')
    mask_cfs_with_lat_ven.inputs.expr = 'step(a)*step(b)'
    mask_cfs_with_lat_ven.inputs.outputtype = 'NIFTI_GZ'
    mask_cfs_with_lat_ven.inputs.out_file = 'cfs_lat_ven_mask.nii.gz'

    nuisance.connect(inputspec, 'csf_mask_file_path', mask_cfs_with_lat_ven, 'in_file_a')
    nuisance.connect(lat_ven_mni_to_anat, 'output_image', mask_cfs_with_lat_ven, 'in_file_b')

    # Resample the csf restricted lateral ventricle masks from 1mm to 2mm, but
    # remaining in subject space
    lat_ven_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='lat_ven_anat_to_2mm_flirt_applyxfm')
    lat_ven_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    lat_ven_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(mask_cfs_with_lat_ven, 'out_file', lat_ven_anat_to_2mm, 'in_file')
    nuisance.connect(mask_cfs_with_lat_ven, 'out_file', lat_ven_anat_to_2mm, 'reference')

    # erode the lateral ventricle mask by 1 voxel to avoid overlap with other
    # tissues
    lat_ven_2mm_erode = pe.Node(interface=afni.Calc(), name='lat_ven_2mm_erode')
    lat_ven_2mm_erode.inputs.args = '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    lat_ven_2mm_erode.inputs.expr = 'a*(1-amongst(0,b,c,d,e,f,g))'
    lat_ven_2mm_erode.inputs.outputtype = 'NIFTI_GZ'
    lat_ven_2mm_erode.inputs.out_file = 'lat_ven_2mm_eroded.nii.gz'

    nuisance.connect(lat_ven_anat_to_2mm, 'out_file', lat_ven_2mm_erode, 'in_file_a')

    # combine the three tissue files, GM, WM, and lat ven into a single file
    # to make them easier to use in the future
    # the BIDS spec says that GM = 1, WM = 2 and LV = 3
    combine_masks = pe.Node(interface=afni.Calc(), name='combine_masks')
    combine_masks.inputs.expr = '1*step(a) + 2*step(b) + 3*step(c)'
    combine_masks.inputs.outputtype = 'NIFTI_GZ'
    combine_masks.inputs.out_file = 'combined_masks_2mm.nii.gz'

    nuisance.connect(gm_anat_2mm_erode, 'out_file', combine_masks, 'in_file_a')
    nuisance.connect(wm_anat_2mm_erode, 'out_file', combine_masks, 'in_file_b')
    nuisance.connect(lat_ven_2mm_erode, 'out_file', combine_masks, 'in_file_c')

    # now upsample the functional data to 2mm, this is a compromise to try and
    # avoid partial voluming between tissues but not going to 1mm space which
    # would use more memory than we have available.
    func_to_2mm = pe.Node(interface=fsl.FLIRT(), name='func_to_2mm_flirt_applyxfm')
    func_to_2mm.inputs.args = '-applyisoxfm 2'
    func_to_2mm.interface.estimated_memory_gb = 2.0

    nuisance.connect(inputspec, 'functional_file_path', func_to_2mm, 'in_file')
    nuisance.connect(combine_masks, 'out_file', func_to_2mm, 'reference')
    nuisance.connect(inputspec, 'func_to_anat_linear_xfm_file_path', func_to_2mm, 'in_matrix_file')

    if 'tCompCor' in selector and selector['tCompCor']:

        # for tCompCor calculate a temporal variance mask using the input
        # parameters
        create_tcompcor_mask = pe.Node(util.Function(input_names=['functional_data_file_path',
                                                                  'mask_file_path',
                                                                  'threshold',
                                                                  'output_file_name',
                                                                  'by_slice'],
                                                     output_names=['tCompCor_mask_file_path'],
                                                     function=create_temporal_variance_mask),
                                       name='create_tCompCor_mask')

        if not ('threshold' in selector['tCompCor'] and selector['tCompCor']['threshold']):
            raise ValueError("tCompCor requires a threshold value, but none received.")

        create_tcompcor_mask.inputs.threshold = selector['tCompCor']['threshold']

        if not ('by_slice' in selector['tCompCor'] and selector['tCompCor']['by_slice'] is not None):
            raise ValueError("tCompCor requires a value for by_slice, but none received.")

        create_tcompcor_mask.inputs.by_slice = selector['tCompCor']['by_slice']

        nuisance.connect(inputspec, 'functional_file_path', create_tcompcor_mask, 'functional_data_file_path')
        nuisance.connect(inputspec, 'functional_brain_mask_file_path', create_tcompcor_mask, 'mask_file_path')

        create_tcompcor_mask.inputs.output_file_name = 'variance_mask.nii.gz'

    # lets build up a resource pool to make it easier to find what we need
    # when we are putting everything together
    tissue_mask_resource_pool = {
        'tCompCor': {
            'functional_file_path': (inputspec, 'functional_file_path'),
            'mask_file_path': (create_temporal_variance_mask, 'out_file'),
            'mask_label': 1
        },
        'aCompCor_WM': {
            'functional_file_path': (func_to_2mm, 'out_file'),
            'mask_file_path': (combine_masks, 'out_file'),
            'mask_label': 2
        },
        'aCompCor_CSF': {
            'functional_file_path': (func_to_2mm, 'out_file'),
            'mask_file_path': (combine_masks, 'out_file'),
            'mask_label': 3
        },
        'aCompCor_WM+CSF': {
            'functional_file_path': (func_to_2mm, 'out_file'),
            'mask_file_path': (combine_masks, 'out_file'),
            'mask_label': (2, 3)
        },
        'GlobalSignal': {
            'functional_file_path': (inputspec, 'functional_file_path'),
            'mask_file_path': (inputspec, 'functional_brain_mask_file_path'),
            'mask_label': 1
        },
        'GreyMatter': {
            'functional_file_path': (func_to_2mm, 'out_file'),
            'mask_file_path': (combine_masks, 'out_file'),
            'mask_label': 1
        },
        'WhiteMatter': {
            'functional_file_path': (func_to_2mm, 'out_file'),
            'mask_file_path': (combine_masks, 'out_file'),
            'mask_label': 2
        },
        'Ventricles': {
            'functional_file_path': (func_to_2mm, 'out_file'),
            'mask_file_path': (combine_masks, 'out_file'),
            'mask_label': 3
        }
    }

    tissue_regressor_nodes = {}

    # now that we have the tissue masks, we can extract the various tissue
    # regressors
    for tissue_regressor in ['aCompCor', 'GreyMatter', 'WhiteMatter', 'Ventricles', 'GlobalSignal']:

        if tissue_regressor in selector and selector[tissue_regressor]:

            if tissue_regressor is 'aCompCor':
                if 'tissues' not in selector[tissue_regressor]:
                    selector[tissue_regressor]['tissues'] = "WM+CSF"
                tissue_mask_resource_pool_key = "_".join([tissue_regressor, selector[tissue_regressor]['tissues']])
            else:
                tissue_mask_resource_pool_key = tissue_regressor

            tissue_regressor_nodes[tissue_regressor] = pe.Node(util.Function(input_names=['functional_file_path',
                                                                                          'mask_file_path',
                                                                                          'output_file_path',
                                                                                          'summary_method',
                                                                                          'mask_vol_index',
                                                                                          'mask_label',
                                                                                          'num_pcs'],
                                                                             output_names=['regressor_file_path'],
                                                                             function=mask_summarize_time_course),
                                                               name='summarise_regressor_{}'.format(tissue_regressor))

            tissue_regressor_nodes[tissue_regressor].interface.estimated_memory_gb = 1.0

            # probably seems redundant to the user to specify CompCor and PCA
            # summarization, since it is redundant
            # lets just set here to keep the logic the same as for other
            # tissue without a bunch of branching
            if tissue_regressor in ['aCompCor', 'tCompCor']:
                selector[tissue_regressor]['summary_method'] = "PCA"

            if 'summary_method' not in selector[tissue_regressor]:
                raise ValueError("Missing method for summarizing {0} tissue voxels into a nuisance"
                                 " regressor".format(tissue_regressor))
            tissue_regressor_nodes[tissue_regressor].inputs.method = selector[tissue_regressor]['summary_method']

            if selector[tissue_regressor]['summary_method'] is 'PCA':
                if 'num_pcs' not in selector[tissue_regressor]:
                    raise ValueError("Summarization method for {0} is PCA, but num_pcs is not specified in "
                                     "selector".format(tissue_regressor))
                tissue_regressor_nodes[tissue_regressor].inputs.num_pcs = selector[tissue_regressor]['num_pcs']

            tissue_regressor_nodes[tissue_regressor].inputs.mask_vol_index = 0
            tissue_regressor_nodes[tissue_regressor].inputs.mask_label = \
                [tissue_mask_resource_pool[tissue_mask_resource_pool_key]['mask_label']]

            tissue_regressor_nodes[tissue_regressor].inputs.output_file_path = \
                '{}_regressor.tsv'.format(tissue_regressor)

            nuisance.connect(tissue_mask_resource_pool[tissue_mask_resource_pool_key]['functional_file_path'][0],
                             tissue_mask_resource_pool[tissue_mask_resource_pool_key]['functional_file_path'][1],
                             tissue_regressor_nodes[tissue_regressor], 'functional_file_path')

            nuisance.connect(tissue_mask_resource_pool[tissue_mask_resource_pool_key]['mask_file_path'][0],
                             tissue_mask_resource_pool[tissue_mask_resource_pool_key]['mask_file_path'][1],
                             tissue_regressor_nodes[tissue_regressor], 'mask_file_path')

    # now add in anaticor if its requested
    anaticor_regressor_to_functional_space = None
    if 'Anaticor' in selector and selector['Anaticor']:

        # make sure we have a radius
        if 'radius' not in selector['Anaticor']:
            raise ValueError('Anaticor specified in selector, but not radius. Radius is a required parameter.')

        # construct the regressors for anaticor
        # '3dLocalstat -prefix __WMeLOCAL_r${r} -nbhd 'SPHERE('${r}')' \
        #    -stat mean -mask  __mask_WMe${view} \
        #    -use_nonmask ${fn_epi}'
        construct_anaticor_regressor = pe.Node(interface=nuisance_afni_interfaces.Localstat(),
                                               name="construct_anaticor_regressor")

        construct_anaticor_regressor.interface.num_threads = 4
        construct_anaticor_regressor.inputs.neighborhood = 'SPHERE({0})'.format(selector['Anaticor']['radius'])
        construct_anaticor_regressor.inputs.statistic = 'mean'
        construct_anaticor_regressor.inputs.use_nonmask = True
        construct_anaticor_regressor.inputs.output = 'anaticor_regressor_2mm.nii.gz'
        construct_anaticor_regressor.inputs.output_type = 'NIFTI_GZ'

        nuisance.connect(wm_anat_2mm_erode, 'out_file', construct_anaticor_regressor, 'mask')
        nuisance.connect(func_to_2mm, 'out_file', construct_anaticor_regressor, 'mask')

        # down sample the data to match functional image space
        anaticor_regressor_to_functional_space = pe.Node(interface=fsl.FLIRT(),
                                                         name='anaticor_regressor_to_functional_space')
        anaticor_regressor_to_functional_space.inputs.args = '-applyxfm'

        nuisance.connect(construct_anaticor_regressor, 'out_file', anaticor_regressor_to_functional_space, 'in_file')
        nuisance.connect(inputspec, 'functional_file_path', anaticor_regressor_to_functional_space, 'reference')

    # deal with censoring
    find_censors = None
    if 'Censor' in selector and selector['Censor']:

        if 'censor_method' not in selector['Censor'] or not selector['Censor']['censor_method']:
            raise ValueError('Censoring requested, but method not provided.')

        if selector['Censor']['censor_method'] not in ['Kill', 'Zero', 'Interpolate', 'SpikeRegression']:
            raise ValueError("Improper censoring method specified ({0}), should be one of ['Kill', 'Zero', "
                             "'Interpolate', 'SpikeRegression'].".format(selector['Censor']['censor_method']))

        find_censors = pe.Node(util.Function(input_names=['thresh_metric',
                                                          'out_file_path',
                                                          'fd_file_path',
                                                          'dvars_file_path',
                                                          'fd_threshold',
                                                          'dvars_threshold',
                                                          'number_of_previous_trs_to_remove',
                                                          'number_of_subsequent_trs_to_remove'],
                                             output_names=['out_file'],
                                             function=find_offending_time_points),
                               name="find_censors_fd_dvars_extend")

        if 'thresh_metric' not in selector['Censor'] or not selector['Censor']['thresh_metric']:
            raise ValueError('Censoring requested, but thresh_metric not provided.')
        find_censors.inputs.thresh_metric = selector['Censor']['thresh_metric']

        find_censors.inputs.out_file_path = "censors.1D"

        if selector['Censor']['thresh_metric'] in ['FD', 'FD+DVARS']:
            if 'fd_threshold' not in selector['Censor'] or not selector['Censor']['fd_threshold']:
                raise ValueError('Censoring thresh_metric {} requires fd_threshold but it was not '
                                 'provided.'.format(selector['Censor']['thresh_metric']))
            find_censors.inputs.fd_threshold = selector['Censor']['fd_threshold']

        if selector['Censor']['thresh_metric'] in ['DVARS', 'FD+DVARS']:
            if 'dvars_threshold' not in selector['Censor'] or not selector['Censor']['dvars_threshold']:
                raise ValueError('Censoring thresh_metric {} requires dvars_threshold but it was not '
                                 'provided.'.format(selector['Censor']['thresh_metric']))
            find_censors.inputs.dvars_threshold = selector['Censor']['dvars_threshold']

        if 'number_of_previous_trs_to_remove' in selector['Censor'] and \
                selector['Censor']['number_of_previous_trs_to_remove'] and \
                selector['Censor']['censor_method'] is not 'SpikeRegression':
            find_censors.inputs.number_of_previous_trs_to_remove = \
                selector['Censor']['number_of_previous_trs_to_remove']
        else:
            find_censors.inputs.number_of_previous_trs_to_remove = 0

        if 'number_of_subsequent_trs_to_remove' in selector['Censor'] and \
                selector['Censor']['number_of_subsequent_trs_to_remove'] and \
                selector['Censor']['censor_method'] is not 'SpikeRegression':
            find_censors.inputs.number_of_subsequent_trs_to_remove = \
                selector['Censor']['number_of_subsequent_trs_to_remove']
        else:
            find_censors.inputs.number_of_subsequent_trs_to_remove = 0

        nuisance.connect(inputspec, "fd_file_path", find_censors, "fd_file_path")
        nuisance.connect(inputspec, "dvars_file_path", find_censors, "dvars_file_path")

    gather_nuisance_input_map = {'DVARS': ('dvars_file_path', inputspec, 'dvars_file_path'),
                                 'FD': ('framewise_displacement_file_path', inputspec,
                                        'framewise_displacement_file_path'),
                                 'Motion': ('motion_parameters_file_path', inputspec, 'motion_parameters_file_path'),
                                 }

    for tissue_regressor in [('aCompCor', 'acompcorr_file_path'),
                             ('tCompCor', 'tcompcorr_file_path'),
                             ('GlobalSignal', 'global_summary_file_path'),
                             ('GreyMatter', 'grey_matter_summary_file_path'),
                             ('Ventricles', 'csf_summary_file_path'),
                             ('WhiteMatter', 'white_matter_summary_file_path')]:
        if tissue_regressor[0] in tissue_regressor_nodes:
            gather_nuisance_input_map[tissue_regressor[0]] = (tissue_regressor[1],
                                                              tissue_regressor_nodes[tissue_regressor[0]],
                                                              'regressor_file_path')

    if 'Censor' in selector and selector['Censor']:
        gather_nuisance_input_map['Censor'] = ('censor_file_path', find_censors, 'out_file')

    build_nuisance_regressors = None

    # first check to see if we have any regressors, if so we need to combine
    # them into a single file
    build_nuisance_regressors_flag = False
    for regressor_key in gather_nuisance_input_map:
        if regressor_key in selector and selector[regressor_key]:
            build_nuisance_regressors_flag = True
            break

    if build_nuisance_regressors_flag is True:

        # make spike regressors if requested and combine all of the regressors
        # into a single file
        build_nuisance_regressors = pe.Node(util.Function(input_names=['functional_file_path',
                                                                       'selector',
                                                                       'output_file_path',
                                                                       'grey_matter_summary_file_path',
                                                                       'white_matter_summary_file_path',
                                                                       'csf_summary_file_path',
                                                                       'acompcorr_file_path',
                                                                       'tcompcorr_file_path',
                                                                       'global_summary_file_path',
                                                                       'motion_parameters_file_path',
                                                                       'framewise_displacement_file_path',
                                                                       'dvars_file_path',
                                                                       'censor_file_path'],
                                                          output_names=['out_file'],
                                                          function=gather_nuisance),
                                            name="build_nuisance_regressors")

        build_nuisance_regressors.inputs.selector = selector
        build_nuisance_regressors.inputs.output_file_path = "nuisance_regressors.1D"

        nuisance.connect(inputspec, 'functional_file_path', build_nuisance_regressors, 'functional_file_path')

        for regressor_key, regressor_val in gather_nuisance_input_map.iteritems():

            if regressor_key in selector and selector[regressor_key]:

                if regressor_key is 'Censor' and selector[regressor_key]['censor_method'] is not 'SpikeRegression':
                    continue

                try:
                    nuisance.connect(regressor_val[1], regressor_val[2], build_nuisance_regressors, regressor_val[0])
                except:
                    print("Error trying to connect {0} to build_nuisance_regressors".format(regressor_key))
                    raise

    # the finale, invoke 3dTproject to perform nuisance variable regression
    nuisance_regression = pe.Node(interface=nuisance_afni_interfaces.Tproject(),
                                  name='nuisance_regression')

    nuisance_regression.inputs.out_file = 'residuals.nii.gz'
    nuisance_regression.inputs.outputtype = 'NIFTI_GZ'
    nuisance_regression.inputs.normalize = True

    nuisance.connect(inputspec, 'functional_file_path', nuisance_regression,
                     'in_file')

    nuisance.connect(inputspec, 'functional_brain_mask_file_path',
                     nuisance_regression, 'mask')

    if build_nuisance_regressors_flag is True:
        nuisance.connect(build_nuisance_regressors, 'out_file',
                         nuisance_regression, 'orthogonalize_file')

    if 'Censor' in selector and 'censor_method' in selector['Censor']:
        if 'censor_method' not in selector['Censor'] or not selector['Censor']['censor_method']:
            raise ValueError('Censoring requested, but method not provided.')

        if selector['Censor']['censor_method'] in ['Kill', 'Zero', 'Interpolate']:
            if selector['Censor']['censor_method'] is 'Interpolate':
                nuisance_regression.inputs.censor_mode = 'NTRP'
            else:
                nuisance_regression.inputs.censor_mode = selector['Censor']['censor_method'].upper()

            nuisance.connect(find_censors, 'out_file', nuisance_regression, 'censor_file')

    if 'PolyOrt' in selector or selector['PolyOrt']:
        if 'degree' not in selector['PolyOrt'] or not selector['PolyOrt']['degree']:
            raise ValueError('Polynomial orthogonalization requested, but degree not provided.')

        nuisance_regression.inputs.orthogonalize_polynomial = selector['PolyOrt']['degree']

    if 'Anaticor' in selector and selector['Anaticor']:
        nuisance.connect(anaticor_regressor_to_functional_space, 'out_file', nuisance_regression,
                         'orthogonalize_dataset')

    if 'Bandpass' in selector and selector['Bandpass']:
        if 'bottom_frequency' in selector['Bandpass'] and selector['Bandpass']['bottom_frequency']:
            raise ValueError('Bandpass filtering requested, but bottom_frequency not provided. Set to 0 if you would '
                             'like a lowpass only filter')
        if 'top_frequency' in selector['Bandpass'] and selector['Bandpass']['top_frequency']:
            raise ValueError('Bandpass filtering requested, but top_frequency not provided. Set to 9999 if you would '
                             'like a lowpass only filter')

        nuisance_regression.inputs.bandpass = [float(selector['Bandpass']['bottom_frequency']),
                                               float(selector['Bandpass']['top_frequency'])]

    nuisance.connect(nuisance_regression, 'out_file', outputspec, 'residual_file_path')
    nuisance.connect(build_nuisance_regressors, 'out_file', outputspec, 'regressors_file_path')

    return nuisance
