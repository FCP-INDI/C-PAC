import os
import numpy as np
import nibabel as nb

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
from nipype.interfaces import afni

from CPAC.nuisance import find_offending_time_points, create_temporal_variance_mask
import CPAC.utils as utils

from scipy.fftpack import fft, ifft

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants


def gather_nuisance(functional_file_path, selector, output_file_path, grey_matter_summary_file_path=None,
                    white_matter_summary_file_path=None, csf_summary_file_path=None, acompcorr_file_path=None,
                    tcompcorr_file_path=None, global_summary_file_path=None, motion_parameters_file_path=None,
                    dvars_file_path=None, framewise_displacement_file_path=None, censor_file_path=None):
    """
    Gathers the various nuisance regressors together into a single tab separated values file that is an appropriate for
    input into 3dTproject
    
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

    motion_labels = ["RotY", "RotX", "RotZ", "Y", "X", "Z"]

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

            # if other parameters are missing, assume they are False by default
            for regressor_derivative in ['include_delayed', 'include_squared', 'include_delayed_squared']:
                if regressor_derivative not in selector[regressor_type] or \
                        not selector[regressor_type][regressor_derivative]:

                    selector[regressor_type][regressor_derivative] = False

            # add in the regressors, making sure to also add in the column name
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


def calc_residuals(subject,
                   selector,
                   despiking=False,
                   wm_sig_file = None,
                   csf_sig_file = None,
                   gm_sig_file = None,
                   motion_file = None,
                   compcor_ncomponents = 0,
                   frames_ex=None):
    """
    Calculates residuals of nuisance regressors for every voxel for a subject.
    
    Parameters
    ----------
    subject : string
        Path of a subject's realigned nifti file.
    selector : dictionary
        Dictionary of selected regressors.  Keys are  represented as a string of the regressor name and keys 
        are True/False.  See notes for an example.
    wm_sig_file : string, optional
        Path to subject's white matter mask (in the same space as the subject's functional file)
    csf_sig_file : string, optional
        Path to subject's cerebral spinal fluid mask (in the same space as the subject's functional file)
    gm_sig_file : string, optional
        Path to subject's grey matter mask (in the same space as the subject's functional file)
    compcor_ncomponents : integer, optional
        The first `n` principal of CompCor components to use as regressors.  Default is 0.
    frames_ex : string, optional
        Filepath to the 1D file describing the volumes to be excluded (for
        de-spiking), selected via the threshold set for excessive motion.
        
    Returns
    -------
    residual_file : string
        Path of residual file in nifti format
    regressors_file : string
        Path of csv file of regressors used.  Filename corresponds to the name of each
        regressor in each column.
        
    Notes
    -----
    
    Example of selector parameter:
    
    >>> selector = {'compcor' : True,
    >>> 'wm' : True,
    >>> 'csf' : True,
    >>> 'gm' : True,
    >>> 'global' : True,
    >>> 'pc1' : True,
    >>> 'motion' : True,
    >>> 'linear' : True,
    >>> 'quadratic' : True}
    """
    
    nii = nb.load(subject)
    data = nii.get_data().astype(np.float64)
    global_mask = (data != 0).sum(-1) != 0
    
    # Check and define regressors which are provided from files
    if wm_sig_file is not None:
        wm_sigs = np.load(wm_sig_file)
        if wm_sigs.shape[1] != data.shape[3]:
            raise ValueError('White matter signals length {0} do not match '
                             'data timepoints {1}'.format(wm_sigs.shape[1], 
                                                          data.shape[3]))
        if wm_sigs.size == 0:
            raise ValueError('White matter signal file {0} is '
                             'empty'.format(wm_sig_file))
        
    if csf_sig_file is not None:
        csf_sigs = np.load(csf_sig_file)
        if csf_sigs.shape[1] != data.shape[3]:
            raise ValueError('CSF signals length {0} do not match data '
                             'timepoints {1}'.format(csf_sigs.shape[1], 
                                                     data.shape[3]))
        if csf_sigs.size == 0:
            raise ValueError('CSF signal file {0} is '
                             'empty'.format(csf_sig_file))
        
    if gm_sig_file is not None:
        gm_sigs = np.load(gm_sig_file)
        if gm_sigs.shape[1] != data.shape[3]:
            raise ValueError('Grey matter signals length {0} do not match '
                             'data timepoints {1}'.format(gm_sigs.shape[1], 
                                                          data.shape[3]))
        if gm_sigs.size == 0:
            raise ValueError('Grey matter signal file {0} is '
                             'empty'.format(gm_sig_file))
        
    if motion_file is not None:
        motion = np.genfromtxt(motion_file)
        if motion.shape[0] != data.shape[3]:
            raise ValueError('Motion parameters {0} do not match data '
                             'timepoints {1}'.format(motion.shape[0], 
                                                     data.shape[3]))
        if motion.size == 0:
            raise ValueError('Motion signal file {0} is '
                             'empty'.format(motion_file))

    # Calculate regressors
    regressor_map = {'constant': np.ones((data.shape[3], 1))}

    if selector['compcor']:
        if not wm_sig_file:
            err = "\n\n[!] CompCor cannot be run because the white matter " \
                  "mask was not generated.\n\n"
            raise Exception(err)
        if not csf_sig_file:
            err = "\n\n[!] CompCor cannot be run because the CSF mask " \
                  "was not generated.\n\n"
            raise Exception(err)

        regressor_map['compcor'] = \
            calc_compcor_components(data, compcor_ncomponents,
                                    wm_sigs, csf_sigs)
    
    if selector['wm']:
        regressor_map['wm'] = wm_sigs.mean(0)
        
    if selector['csf']:
        regressor_map['csf'] = csf_sigs.mean(0)
        
    if selector['gm']:
        regressor_map['gm'] = gm_sigs.mean(0)
        
    if selector['global']:
        regressor_map['global'] = data[global_mask].mean(0)
        
    if selector['pc1']:
        bdata = data[global_mask].T
        bdatac = bdata - np.tile(bdata.mean(0), (bdata.shape[0], 1))
        U, S, Vh = np.linalg.svd(bdatac, full_matrices=False)
        regressor_map['pc1'] = U[:, 0]
        
    if selector['motion']:
        regressor_map['motion'] = motion
        
    if selector['linear']:
        regressor_map['linear'] = np.arange(0, data.shape[3])
    
    if selector['quadratic']:
        regressor_map['quadratic'] = np.arange(0, data.shape[3])**2

    # insert the de-spiking regressor matrix here, if running de-spiking
    if frames_ex:
        despike_mat = create_despike_regressor_matrix(frames_ex, nii.shape[3])
        # this needs to be "is not None" instead of "if despike_mat:" because
        # despike_mat could be either a Numpy array or None
        if despike_mat is not None:
            regressor_map['despike'] = despike_mat

    X = np.zeros((data.shape[3], 1))
    csv_filename = ''
    for rname, rval in regressor_map.items():
        X = np.hstack((X, rval.reshape(rval.shape[0],-1)))
        csv_filename += '_' + rname
    X = X[:,1:]
    
    csv_filename = csv_filename[1:]
    csv_filename += '.csv'
    csv_filename = os.path.join(os.getcwd(), csv_filename)
    np.savetxt(csv_filename, X, delimiter='\t')
    
    if np.isnan(X).any() or np.isnan(X).any():
        raise ValueError('Regressor file contains NaN')

    Y = data[global_mask].T

    try:
        B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    except np.linalg.LinAlgError as e:
        if "Singular matrix" in e:
            raise Exception("Error details: {0}\n\nSingular matrix error: "
                            "The nuisance regression configuration you "
                            "selected may have been too stringent, and the "
                            "regression could not be completed. Ensure your "
                            "parameters (such as the motion threshold for "
                            "de-spiking or scrubbing) are not too "
                            "extreme.\n\n".format(e))
        else:
            raise Exception("Error details: {0}\n\nSomething went wrong with "
                            "nuisance regression.\n\n".format(e))

    Y_res = Y - X.dot(B)
    
    data[global_mask] = Y_res.T
    
    img = nb.Nifti1Image(data, header=nii.get_header(),
                         affine=nii.get_affine())
    residual_file = os.path.join(os.getcwd(), 'residual.nii.gz')
    img.to_filename(residual_file)
    
    # Easier to read for debugging purposes
    regressors_file = os.path.join(os.getcwd(), 'nuisance_regressors.mat')

    if scipy.__version__ == '0.7.0':
        # for scipy v0.7.0
        scipy.io.savemat(regressors_file, regressor_map)
    else:
        # for scipy v0.12: OK
        scipy.io.savemat(regressors_file, regressor_map, oned_as='column')
    
    return residual_file, regressors_file


def extract_tissue_data(data_file,
                        ventricles_mask_file,
                        wm_seg_file, csf_seg_file, gm_seg_file):
    import numpy as np
    import nibabel as nb
    import os    
    from CPAC.nuisance import erode_mask
    from CPAC.utils import safe_shape

    try:
        data = nb.load(data_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % data_file)


    try:
        lat_ventricles_mask = nb.load(ventricles_mask_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % lat_ventricles_mask)

    if not safe_shape(data, lat_ventricles_mask):
        raise ValueError('Spatial dimensions for data and the lateral '
                         'ventricles mask do not match')

    try:
        wm_seg = nb.load(wm_seg_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % wm_seg)

    if not safe_shape(data, wm_seg):
        raise ValueError('Spatial dimensions for data, white matter segment '
                         'do not match')

    wm_mask = erode_mask(wm_seg > 0)
    wm_sigs = data[wm_mask]
    file_wm = os.path.join(os.getcwd(), 'wm_signals.npy')
    np.save(file_wm, wm_sigs)
    del wm_sigs

    try:
        csf_seg = nb.load(csf_seg_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % csf_seg)

    if not safe_shape(data, csf_seg):
        raise ValueError('Spatial dimensions for data, cerebral spinal '
                         'fluid segment do not match')

    # Only take the CSF at the lateral ventricles as labeled in the Harvard
    # Oxford parcellation regions 4 and 43
    csf_mask = (csf_seg > 0)*(lat_ventricles_mask==1)
    csf_sigs = data[csf_mask]
    file_csf = os.path.join(os.getcwd(), 'csf_signals.npy')
    np.save(file_csf, csf_sigs)
    del csf_sigs

    try:
        gm_seg = nb.load(gm_seg_file).get_data().astype('float64')
    except:
        raise MemoryError('Unable to load %s' % gm_seg)

    if not safe_shape(data, gm_seg):
        raise ValueError('Spatial dimensions for data, gray matter '
                         'segment do not match')

    gm_mask = erode_mask(gm_seg > 0)
    gm_sigs = data[gm_mask]
    file_gm = os.path.join(os.getcwd(), 'gm_signals.npy')
    np.save(file_gm, gm_sigs)
    del gm_sigs

    nii = nb.load(wm_seg_file)
    wm_mask_file = os.path.join(os.getcwd(), 'wm_mask.nii.gz')
    csf_mask_file = os.path.join(os.getcwd(), 'csf_mask.nii.gz')
    gm_mask_file = os.path.join(os.getcwd(), 'gm_mask.nii.gz')
    nb.Nifti1Image(wm_mask, header=nii.get_header(), affine=nii.get_affine()).to_filename(wm_mask_file)
    nb.Nifti1Image(csf_mask, header=nii.get_header(), affine=nii.get_affine()).to_filename(csf_mask_file)
    nb.Nifti1Image(gm_mask, header=nii.get_header(), affine=nii.get_affine()).to_filename(gm_mask_file)

    return file_wm, file_csf, file_gm


def create_nuisance(use_ants, name='nuisance'):
    """
    Workflow for the removal of various signals considered to be noise in resting state
    fMRI data.  The residual signals for linear regression denoising is performed in a single
    model.  Therefore the residual time-series will be orthogonal to all signals.
    
    Parameters
    ----------
    name : string, optional
        Name of the workflow.
    
    Returns
    -------
    nuisance : nipype.pipeline.engine.Workflow
        Nuisance workflow.
        
    Notes
    -----
    
    Workflow Inputs::

        inputspec.subject : string (nifti file)
            Path of the subject's realigned nifti file.
        inputspec.wm_mask : string (nifti file)
            Corresponding white matter mask.
        inputspec.csf_mask : string (nifti file)
            Corresponding cerebral spinal fluid mask.
        inputspec.gm_mask : string (nifti file)
            Corresponding grey matter mask.
        inputspec.mni_to_anat_linear_xfm : string (nifti file)
            Corresponding MNI to anatomical linear transformation 
        inputspec.func_to_anat_linear_xfm : string (nifti file)
            Corresponding EPI to anatomical linear transformation
        inputspec.harvard_oxford_mask : string (nifti file)
            Harvard Oxford parcellation for ventrical locations
        inputspec.motion_components : string (text file)
            Corresponding rigid-body motion parameters.  Matrix in the file should be of shape 
            (`T`, `R`), `T` timepoints and `R` motion parameters.
        inputspec.selector : dictionary
        inputspec.compcor_ncomponents : integer
        
    Workflow Outputs::

        outputspec.subject : string (nifti file)
            Path of residual file in nifti format
        outputspec.regressors : string (mat file)
            Path of csv file of regressors used.  Filename corresponds to the name of each
            regressor in each column.
            
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
    
    inputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                       'wm_mask',
                                                       'csf_mask',
                                                       'gm_mask',
                                                       'mni_to_anat_linear_xfm',
                                                       'anat_to_mni_initial_xfm',
                                                       'anat_to_mni_rigid_xfm',
                                                       'anat_to_mni_affine_xfm',
                                                       'func_to_anat_linear_xfm',
                                                       'lat_ventricles_mask',
                                                       'motion_components',
                                                       'selector',
                                                       'compcor_ncomponents',
                                                       'template_brain',
                                                       'frames_ex']),
                        name='inputspec')
    outputspec = pe.Node(util.IdentityInterface(fields=['subject',
                                                        'regressors']),
                         name='outputspec')

    # Resampling the masks from 1mm to 2mm, but remaining in subject space
    wm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='wm_anat_to_2mm_flirt_applyxfm')
    wm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    wm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'wm_mask', wm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'wm_mask', wm_anat_to_2mm, 'reference')

    # Resampling the masks from 1mm to 2mm, but remaining in subject space
    csf_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='csf_anat_to_2mm_flirt_applyxfm')
    csf_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    csf_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'csf_mask', csf_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'csf_mask', csf_anat_to_2mm, 'reference')

    # Resampling the masks from 1mm to 2mm, but remaining in subject space
    gm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='gm_anat_to_2mm_flirt_applyxfm')
    gm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    gm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'gm_mask', gm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'gm_mask', gm_anat_to_2mm, 'reference')

    func_to_2mm = pe.Node(interface=fsl.FLIRT(), name='func_to_2mm_flirt_applyxfm')
    func_to_2mm.inputs.args = '-applyisoxfm 2'
    func_to_2mm.interface.estimated_memory_gb = 2.0

    nuisance.connect(inputspec, 'subject', func_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'csf_mask', func_to_2mm, 'reference')
    nuisance.connect(inputspec, 'func_to_anat_linear_xfm', func_to_2mm, 'in_matrix_file')

    if use_ants:
        collect_linear_transforms = pe.Node(util.Merge(3), name='ho_mni_to_2mm_ants_collect_linear_transforms')

        ho_mni_to_2mm = pe.Node(interface=ants.ApplyTransforms(), name='ho_mni_to_2mm_ants_applyxfm')

        ho_mni_to_2mm.inputs.invert_transform_flags = [True, True, True]
        ho_mni_to_2mm.inputs.interpolation = 'NearestNeighbor'
        ho_mni_to_2mm.inputs.dimension = 3

        nuisance.connect(inputspec, 'anat_to_mni_initial_xfm', collect_linear_transforms, 'in1')
        nuisance.connect(inputspec, 'anat_to_mni_rigid_xfm', collect_linear_transforms, 'in2')
        nuisance.connect(inputspec, 'anat_to_mni_affine_xfm', collect_linear_transforms, 'in3')

        nuisance.connect(collect_linear_transforms, 'out', ho_mni_to_2mm, 'transforms')

        nuisance.connect(inputspec, 'lat_ventricles_mask', ho_mni_to_2mm, 'input_image')
        nuisance.connect(csf_anat_to_2mm, 'out_file', ho_mni_to_2mm, 'reference_image')

        # resample_to_2mm = pe.Node(interface=afni.Resample(), name='resample_to_2mm_ants_output'
        
    else:
        ho_mni_to_2mm = pe.Node(interface=fsl.FLIRT(), name='ho_mni_to_2mm_flirt_applyxfm')
        ho_mni_to_2mm.inputs.args = '-applyisoxfm 2'
        ho_mni_to_2mm.inputs.interp = 'nearestneighbour'

        nuisance.connect(inputspec, 'mni_to_anat_linear_xfm', ho_mni_to_2mm, 'in_matrix_file')
        nuisance.connect(inputspec, 'lat_ventricles_mask', ho_mni_to_2mm, 'in_file')
        nuisance.connect(inputspec, 'csf_mask', ho_mni_to_2mm, 'reference')

    tissue_masks = pe.Node(util.Function(input_names=['data_file',
                                                      'ventricles_mask_file',
                                                      'wm_seg_file', 'csf_seg_file', 'gm_seg_file'],
                                         output_names=['file_wm', 'file_csf', 'file_gm'],
                                         function=extract_tissue_data),
                           name='tissue_masks')
    tissue_masks._interface.estimated_memory_gb = 3.0

    nuisance.connect(func_to_2mm, 'out_file', tissue_masks, 'data_file')
    nuisance.connect(wm_anat_to_2mm, 'out_file', tissue_masks, 'wm_seg_file')
    nuisance.connect(csf_anat_to_2mm, 'out_file', tissue_masks, 'csf_seg_file')
    nuisance.connect(gm_anat_to_2mm, 'out_file', tissue_masks, 'gm_seg_file')

    if use_ants:
        nuisance.connect(ho_mni_to_2mm, 'output_image', tissue_masks, 'ventricles_mask_file')
    else:
        nuisance.connect(ho_mni_to_2mm, 'out_file', tissue_masks, 'ventricles_mask_file')

    calc_imports = ['import os', 'import scipy', 'import numpy as np',
                    'import nibabel as nb', 
                    'from CPAC.nuisance import calc_compcor_components',
                    'from CPAC.nuisance.utils import create_despike_regressor_matrix']
    calc_r = pe.Node(util.Function(input_names=['subject',
                                                'selector',
                                                'wm_sig_file',
                                                'csf_sig_file',
                                                'gm_sig_file',
                                                'motion_file',
                                                'compcor_ncomponents',
                                                'frames_ex'],
                                   output_names=['residual_file',
                                                 'regressors_file'],
                                   function=calc_residuals,
                                   imports=calc_imports),
                     name='residuals')
    
    nuisance.connect(inputspec, 'subject', calc_r, 'subject')
    nuisance.connect(tissue_masks, 'file_wm', calc_r, 'wm_sig_file')
    nuisance.connect(tissue_masks, 'file_csf', calc_r, 'csf_sig_file')
    nuisance.connect(tissue_masks, 'file_gm', calc_r, 'gm_sig_file')
    nuisance.connect(inputspec, 'motion_components', calc_r, 'motion_file')
    nuisance.connect(inputspec, 'selector', calc_r, 'selector')
    nuisance.connect(inputspec, 'compcor_ncomponents', 
                     calc_r, 'compcor_ncomponents')
    nuisance.connect(inputspec, 'frames_ex', calc_r, 'frames_ex')

    nuisance.connect(calc_r, 'residual_file', outputspec, 'subject')
    nuisance.connect(calc_r, 'regressors_file', outputspec, 'regressors')
    
    return nuisance


def create_nuisance_workflow(pipeline_resource_pool, nuisance_configuration_selector, functional_specifier, use_ants,
                             name='nuisance'):
    """
    Workflow for the removal of various signals considered to be noise from resting state
    fMRI data.  The residual signals for linear regression denoising is performed in a single
    model.  Therefore the residual time-series will be orthogonal to all signals.

    Parameters
    ----------
    :param pipeline_resource_pool: dictionary of pipeline resources and their source nodes
    :param nuisance_configuration_selector: dictionary describing nuisance regression to be performed
    :param functional_specifier: key corresponding to functional data that should be used
    :param use_ants: flag indicating whether FNIRT or ANTS is used
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
        inputspec.selector : Dictionary containing configuration parameters for nuisance regression
            selector = {'aCompCor' : None | {num_pcs = <number of components to retain>,
                                            tissue = 'WM' | 'CSF' | 'WM+CSF',
                                            extraction_resolution = None | floating point value indicating isotropic
                                                resolution (ex. 2 for 2mm x 2mm x 2mm that data should be extracted at,
                                                the corresponding tissue mask will be resampled to this resolution. The
                                                functional data will also be resampled to this resolution, and the
                                                extraction will occur at this new resolution. The goal is to avoid
                                                contamination from undesired tissue components when extracting nuisance
                                                regressors,
                                             erode_mask = True | False, whether or not the mask should be eroded to
                                                further avoid a mask overlapping with a different tissue class,
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
                                                       PCT,
                                            include_delayed = True | False,
                                            include_squared = True | False,
                                            include_delayed_squared = True | False},
                        'WhiteMatter' : None | {summary_method = 'Anaticor', 'PCA', 'Mean', 'NormMean' or
                                           'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       anaticor_radius = <radius in mm>,
                                       extraction_resolution = None | floating point value (same as for aCompCor),
                                       erode_mask = True | False (same as for aCompCor),
                                       include_delayed = True | False,
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'Ventricles' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       extraction_resolution = None | floating point value (same as for aCompCor),
                                       erode_mask = True | False (same as for aCompCor),
                                       include_delayed = True | False,
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'GreyMatter' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                       num_pcs = <number of components to retain>,
                                       extraction_resolution = None | floating point value (same as for aCompCor),
                                       erode_mask = True | False (same as for aCompCor),
                                       include_delayed = True | False,
                                       include_squared = True | False,
                                       include_delayed_squared = True | False},
                        'GlobalSignal' : None | {summary_method = 'PCA', 'Mean', 'NormMean' or 'DetrendNormMean',
                                           num_pcs = <number of components to retain>,
                                           extraction_resolution = None | floating point value (same as for aCompCor),
                                           erode_mask = True | False (same as for aCompCor),
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
                                                               band, frequencies above this will be removed>}
                        }

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

    # build regressors and combine them into a single file
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

    build_nuisance_regressors.inputs.selector = nuisance_configuration_selector
    build_nuisance_regressors.inputs.output_file_path = "nuisance_regressors.1D"

    nuisance.connect(inputspec, 'functional_file_path', build_nuisance_regressors, 'functional_file_path')

    # regressor map to simplify construction of the needed regressors
    regressor_map = {'aCompCor': ('acompcorr_file_path', []),
                     'tCompCor': ('tcompcorr_file_path', []),
                     'GlobalSignal': ('global_summary_file_path', []),
                     'GreyMatter': ('grey_matter_summary_file_path', []),
                     'WhiteMatter': ('white_matter_summary_file_path', []),
                     'Ventricles': ('csf_summary_file_path', []),
                     'DVARS': ('dvars_file_path', (inputspec, 'dvars_file_path')),
                     'FD': ('framewise_displacement_file_path', (inputspec, 'framewise_displacement_file_path')),
                     'Motion': ('motion_parameters_file_path', (inputspec, 'motion_parameters_file_path'))}

    for regressor_type, regressor_specification in regressor_map.items():

        # go ahead and set summary_method for tCompCor and aCompCor so make the logic easier
        if regressor_type in ['tCompCor', 'aCompCor']:
            nuisance_configuration_selector[regressor_type]['summary_method'] = 'PCA'

        if regressor_type in nuisance_configuration_selector and nuisance_configuration_selector[regressor_type]:
            if not isinstance(regressor_specification, tuple):
                raise ValueError("Regressor specification should be a tuple")

            if len(regressor_specification) != 2 or not regressor_specification[1]:

                # we don't have the regressor, look for it in the resource pool, build a corresponding key,
                # this is seperated in to a mask key and an extraction key, which when concatenated provide
                # the resource key for the regressor
                regressor_descriptor = {'tissue': regressor_type}

                if regressor_type is 'aCompCor':
                    if 'tissue' in nuisance_configuration_selector[regressor_type] and \
                            nuisance_configuration_selector[regressor_type]['tissue']:
                        regressor_descriptor = {'tissue': nuisance_configuration_selector[regressor_type]['tissue']}
                    else:
                        raise ValueError("Tissue type required for aCompCor, but none specified")

                if regressor_type is 'tCompCor':
                    if 'threshold' in nuisance_configuration_selector[regressor_type] and \
                            nuisance_configuration_selector[regressor_type]['threshold']:
                        if 'by_slice' in nuisance_configuration_selector[regressor_type] and \
                                nuisance_configuration_selector[regressor_type]['by_slice']:
                            regressor_descriptor = {'tissue': 'functional-variance-{0}-bySlice'.format(
                                nuisance_configuration_selector[regressor_type]['threshold'])}
                        else:
                            regressor_descriptor = {'tissue': 'functional-variance-{0}'.format(
                                nuisance_configuration_selector[regressor_type]['threshold'])}
                    else:
                        raise ValueError("Threshold required for tCompCor, but none specified.")

                if 'extraction_resolution' in nuisance_configuration_selector[regressor_type] and \
                        nuisance_configuration_selector[regressor_type]['extraction_resolution']:
                    regressor_descriptor['resolution'] = \
                        nuisance_configuration_selector[regressor_type]['extraction_resolution']

                if 'erode_mask' in nuisance_configuration_selector[regressor_type]:
                    regressor_descriptor['erosion'] = nuisance_configuration_selector[regressor_type]['erode_mask']

                # now combine with extraction information to get the regressor file
                if 'summary_method' in nuisance_configuration_selector[regressor_type]and \
                        nuisance_configuration_selector[regressor_type]['summary_method']:
                    regressor_descriptor['extraction'] = \
                        nuisance_configuration_selector[regressor_type]['summary_method']

                    if nuisance_configuration_selector[regressor_type]['summary_method'] is 'PCA' and 'num_pcs' in \
                            nuisance_configuration_selector[regressor_type] and \
                            nuisance_configuration_selector[regressor_type]['num_pcs']:
                        regressor_descriptor['extraction'] += \
                            '_{0}pcs'.format(nuisance_configuration_selector[regressor_type]['num_pcs'])
                    else:
                        raise ValueError("Summary method PCA requires num_pcs, but received none.")
                else:
                    raise ValueError("Summary method required for {0}, but none specified".format(regressor_type))

                regressor_file_resource_key = "_".join([regressor_descriptor[key] for key in ['tissue', 'resolution',
                                                                                              'erosion', 'extraction']
                                                        if key in regressor_descriptor])

                if not (regressor_file_resource_key in pipeline_resource_pool and
                        pipeline_resource_pool[regressor_file_resource_key]):

                    mask_resample_erode_file_resource_key = "_".join([regressor_descriptor[key]
                                                                      for key in ['tissue', 'resolution', 'erosion']
                                                                      if key in regressor_descriptor])

                    if "erosion" in regressor_descriptor and regressor_descriptor["erosion"] and \
                            not (mask_resample_erode_file_resource_key in pipeline_resource_pool and
                                 pipeline_resource_pool[mask_resample_erode_file_resource_key]):

                        mask_resample_file_resource_key = "_".join([regressor_descriptor[key]
                                                                    for key in ['tissue', 'resolution']
                                                                    if key in regressor_descriptor])

                        if not (mask_resample_file_resource_key in pipeline_resource_pool and
                                pipeline_resource_pool[mask_resample_file_resource_key]):

                            mask_file_resource_key = regressor_descriptor['tissue']

                            if not (mask_file_resource_key in pipeline_resource_pool and
                                    pipeline_resource_pool[mask_file_resource_key]):

                                if 'functional-variance' in mask_file_resource_key:
                                    pipeline_resource_pool[mask_file_resource_key] = \
                                        insert_create_variance_mask_node(nuisance, pipeline_resource_pool['functional'],
                                                                         nuisance_configuration_selector)

                                else:
                                    raise ValueError("Nuisance regression requires a {0} mask, but not found in the "
                                                     "pipeline resource pool. Perhaps there is an error in the "
                                                     "configuration?".format(mask_file_resource_key))

                            anatomical_to_epi_transform_key = "{0}_to_anat_linear_xform.mat".format(
                                functional_specifier)

                            if not (anatomical_to_epi_transform_key in pipeline_resource_pool and
                                    pipeline_resource_pool[anatomical_to_epi_transform_key]):
                                raise ValueError("Resampling mask {0} requires transform {1}, but not found".format(
                                    mask_file_resource_key, anatomical_to_epi_transform_key))

                            pipeline_resource_pool[mask_resample_file_resource_key] = \
                                utils.mask_nodes.insert_node_apply_flirt_linear_transform(
                                    nuisance, pipeline_resource_pool[mask_file_resource_key],
                                    pipeline_resource_pool[functional_specifier],
                                    pipeline_resource_pool[anatomical_to_epi_transform_key],
                                    regressor_descriptor["resolution"], mask_resample_file_resource_key)

                        if "censor" in regressor_descriptor:
                            pipeline_resource_pool[mask_resample_erode_file_resource_key] = \
                                utils.mask_nodes.insert_erode_mask_node(nuisance,
                                                                        pipeline_resource_pool[
                                                                            mask_resample_file_resource_key],
                                                                        mask_resample_file_resource_key)

                    functional_file_at_resolution_key = functional_specifier + "_res-{0}".format(
                        regressor_descriptor["resolution"])

                    if not (functional_file_at_resolution_key in pipeline_resource_pool and
                            pipeline_resource_pool[functional_file_at_resolution_key]):

                        if not (functional_specifier in pipeline_resource_pool and
                                pipeline_resource_pool[functional_specifier]):
                            raise ValueError(
                                "Unable to find functional file corresponding to {0} in resource pool.".format(
                                    functional_specifier))

                        # also find the epi - anat xform that will be needed to copy functional into mask space
                        epi_anatomical_transform_key = "func_to_anat_linear_xfm"

                    if nuisance_configuration_selector[regressor_type] is not "anaticor":
                        pipeline_resource_pool[regressor_file_resource_key] = \
                            insert_mask_summarize_time_course_node(
                                nuisance,
                                pipeline_resource_pool[functional_file_key],
                                pipeline_resource_pool[mask_resample_erode_file_resource_key],
                                mask_vol_index=0,  # we always use a 3D mask
                                mask_label=1,  # we explicitly construct mask so that it is binary
                                summarization_method=nuisance_configuration_selector["method"],
                                num_pcs=nuisance_configuration_selector[regressor_type]["num_pcs"],
                                node_name=regressor_type)

                    sys.exit()


                    # couldn't find regressor in resource pool, so lets look for the mask
                    if mask_file_resource_key in pipeline_resource_pool and \
                            pipeline_resource_pool[mask_file_resource_key]:
                        mask_source = pipeline_resource_pool[mask_file_resource_key]
                    else:
                        # couldn't find exact mask, so lets check for a less processed mask that meets our needs
                        if 1:

                            regressor_specification[1] = pipeline_resource_pool[regressor_file_resource_key]

    if tissue_regressor[0] in tissue_regressor_nodes:
            gather_nuisance_input_map[tissue_regressor[0]] = (tissue_regressor[1],
                                                              tissue_regressor_nodes[tissue_regressor[0]],
                                                              'regressor_file_path')

    if 'Censor' in selector and selector['Censor']:
        gather_nuisance_input_map['Censor'] = ('censor_file_path', find_censors, 'out_file')


    build_nuisance_regressors = None


    # first check to see if we have any regressors, if so we need to combine them into a single file
    build_nuisance_regressors_flag = any(
        selector.get(regressor_key)
        for regressor_key in gather_nuisance_input_map
    )

    if build_nuisance_regressors_flag:

        for regressor_key, regressor_val in gather_nuisance_input_map.iteritems():

            if regressor_key in selector and selector[regressor_key]:

                if regressor_key is 'Censor' and selector[regressor_key]['censor_method'] is not 'SpikeRegression':
                    continue

                nuisance.connect(regressor_val[1], regressor_val[2], build_nuisance_regressors, regressor_val[0])


    # for regressor_type in [""]
    if ('WhiteMatter' in nuisance_configuration_selector and nuisance_configuration_selector['WhiteMatter']):

        if not ('gm_wm_ventricle_mask' in pipeline_resource_pool and pipeline_resource_pool['gm_wm_ventricle_mask']):
            if not ('ventricle_mask' in pipeline_resource_pool and pipeline_resource_pool['ventricle_mask']):
                # Resample grey matter masks from 1mm to 2mm, but remaining in subject space
                gm_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='gm_anat_to_2mm_flirt_applyxfm')
                gm_anat_to_2mm.inputs.args = '-applyisoxfm 2'
                gm_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(inputspec, 'gm_mask_file_path', gm_anat_to_2mm, 'in_file')
    nuisance.connect(inputspec, 'gm_mask_file_path', gm_anat_to_2mm, 'reference')

    # erode the grey matter mask by 1 voxel to avoid overlap with other tissues
    gm_anat_2mm_erode = pe.Node(interface=afni.Calc(), name='gm_anat_2mm_erode')
    gm_anat_2mm_erode.inputs.args = '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    gm_anat_2mm_erode.inputs.expr = 'a*(1-amongst(0,b,c,d,e,f,g))'
    gm_anat_2mm_erode.inputs.outputtype = 'NIFTI_GZ'
    gm_anat_2mm_erode.inputs.out_file = 'gm_mask_2mm_eroded.nii.gz'

    nuisance.connect(gm_anat_to_2mm, 'out_file', gm_anat_2mm_erode, 'in_file_a')

    # multiply the CSF mask by a lateral ventricles to reduce it to the CSF in the lateral ventricles

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

    # Resample the csf restricted lateral ventricle masks from 1mm to 2mm, but remaining in subject space
    lat_ven_anat_to_2mm = pe.Node(interface=fsl.FLIRT(), name='lat_ven_anat_to_2mm_flirt_applyxfm')
    lat_ven_anat_to_2mm.inputs.args = '-applyisoxfm 2'
    lat_ven_anat_to_2mm.inputs.interp = 'nearestneighbour'

    nuisance.connect(mask_cfs_with_lat_ven, 'out_file', lat_ven_anat_to_2mm, 'in_file')
    nuisance.connect(mask_cfs_with_lat_ven, 'out_file', lat_ven_anat_to_2mm, 'reference')

    # erode the lateral ventricle mask by 1 voxel to avoid overlap with other tissues
    lat_ven_2mm_erode = pe.Node(interface=afni.Calc(), name='lat_ven_2mm_erode')
    lat_ven_2mm_erode.inputs.args = '-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k'
    lat_ven_2mm_erode.inputs.expr = 'a*(1-amongst(0,b,c,d,e,f,g))'
    lat_ven_2mm_erode.inputs.outputtype = 'NIFTI_GZ'
    lat_ven_2mm_erode.inputs.out_file = 'lat_ven_2mm_eroded.nii.gz'

    nuisance.connect(lat_ven_anat_to_2mm, 'out_file', lat_ven_2mm_erode, 'in_file_a')

    # combine the three tissue files, GM, WM, and lat ven into a single file to make them easier to use in the future
    # the BIDS spec says that GM = 1, WM = 2 and LV = 3
    combine_masks = pe.Node(interface=afni.Calc(), name='combine_masks')
    combine_masks.inputs.expr = '1*step(a) + 2*step(b) + 3*step(c)'
    combine_masks.inputs.outputtype = 'NIFTI_GZ'
    combine_masks.inputs.out_file = 'combined_masks_2mm.nii.gz'

    nuisance.connect(gm_anat_2mm_erode, 'out_file', combine_masks, 'in_file_a')
    nuisance.connect(wm_anat_2mm_erode, 'out_file', combine_masks, 'in_file_b')
    nuisance.connect(lat_ven_2mm_erode, 'out_file', combine_masks, 'in_file_c')

    # now upsample the functional data to 2mm, this is a compromise to try and avoid partial voluming between tissues
    # but not going to 1mm space which would use more memory than we have available.
    func_to_2mm = pe.Node(interface=fsl.FLIRT(), name='func_to_2mm_flirt_applyxfm')
    func_to_2mm.inputs.args = '-applyisoxfm 2'
    func_to_2mm.interface.estimated_memory_gb = 2.0

    nuisance.connect(inputspec, 'functional_file_path', func_to_2mm, 'in_file')
    nuisance.connect(combine_masks, 'out_file', func_to_2mm, 'reference')
    nuisance.connect(inputspec, 'func_to_anat_linear_xfm_file_path', func_to_2mm, 'in_matrix_file')


    if selector.get('tCompCor'):

        # for tCompCor calculate a temporal variance mask using the input parameters
        create_tcompcor_mask = pe.Node(util.Function(input_names=['functional_data_file_path',
                                                                  'mask_file_path',
                                                                  'threshold',
                                                                  'output_file_name',
                                                                  'by_slice'],
                                                     output_names=['tCompCor_mask_file_path'],
                                                     function=create_temporal_variance_mask),
                                       name='create_tCompCor_mask')

        if not selector['tCompCor'].get('threshold'):
            raise ValueError("tCompCor requires a threshold value, but none received.")

        create_tcompcor_mask.inputs.threshold = selector['tCompCor']['threshold']

        if selector['tCompCor'].get('by_slice') is None:
            raise ValueError("tCompCor requires a value for by_slice, but none received.")

        create_tcompcor_mask.inputs.by_slice = selector['tCompCor']['by_slice']

        nuisance.connect(inputspec, 'functional_file_path', create_tcompcor_mask, 'functional_data_file_path')
        nuisance.connect(inputspec, 'functional_brain_mask_file_path', create_tcompcor_mask, 'mask_file_path')

        create_tcompcor_mask.inputs.output_file_name = 'variance_mask.nii.gz'


    # lets build up a resource pool to make it easier to find what we need when we are putting everything together
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


    # now that we have the tissue masks, we can extract the various tissue regressors
    for tissue_regressor in ['aCompCor', 'GreyMatter', 'WhiteMatter', 'Ventricles', 'GlobalSignal']:
        if tissue_regressor in selector and selector[tissue_regressor]:
            if tissue_regressor is 'aCompCor':
                if 'tissues' not in selector[tissue_regressor]:
                    selector[tissue_regressor]['tissues'] = "WM+CSF"
                tissue_mask_resource_pool_key = "_".join([tissue_regressor, selector[tissue_regressor]['tissues']])
            else:
                tissue_mask_resource_pool_key = tissue_regressor



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
        construct_anaticor_regressor = pe.Node(interface=afni.Localstat(),
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

    if selector.get('Censor'):

        if not selector['Censor'].get('censor_method'):
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


    # the finale, invoke 3dTproject to perform nuisance variable regression
    nuisance_regression = pe.Node(interface=afni.TProject(), name='nuisance_regression')

    nuisance_regression.inputs.out_file = 'residuals.nii.gz'
    nuisance_regression.inputs.outputtype = 'NIFTI_GZ'
    nuisance_regression.inputs.normalize = True

    nuisance.connect(inputspec, 'functional_file_path', nuisance_regression, 'in_file')
    nuisance.connect(inputspec, 'functional_brain_mask_file_path', nuisance_regression, 'mask')

    if build_nuisance_regressors_flag:
        nuisance.connect(build_nuisance_regressors, 'out_file', nuisance_regression, 'orthogonalize_file')



    if selector.get('Censor'):

        if selector['Censor']['censor_method'] in ['Kill', 'Zero', 'Interpolate']:
            if selector['Censor']['censor_method'] is 'Interpolate':
                nuisance_regression.inputs.censor_mode = 'NTRP'
            else:
                nuisance_regression.inputs.censor_mode = selector['Censor']['censor_method'].upper()

            nuisance.connect(find_censors, 'out_file', nuisance_regression, 'censor_file')



    if selector.get('PolyOrt'):

        if not selector['PolyOrt'].get('degree'):
            raise ValueError('Polynomial orthogonalization requested, but degree not provided.')

        nuisance_regression.inputs.orthogonalize_polynomial = selector['PolyOrt']['degree']



    if selector.get('Anaticor'):
        nuisance.connect(anaticor_regressor_to_functional_space, 'out_file',
                         nuisance_regression, 'orthogonalize_dataset')



    if selector.get('Bandpass'):
        if not selector['Bandpass'].get('bottom_frequency'):
            raise ValueError('Bandpass filtering requested, but bottom_frequency not provided. Set to 0 if you would '
                             'like a lowpass only filter')

        if not elector['Bandpass'].get('top_frequency'):
            raise ValueError('Bandpass filtering requested, but top_frequency not provided. Set to 9999 if you would '
                             'like a lowpass only filter')

        nuisance_regression.inputs.bandpass = [float(selector['Bandpass']['bottom_frequency']),
                                               float(selector['Bandpass']['top_frequency'])]



    nuisance.connect(nuisance_regression, 'out_file', outputspec, 'residual_file_path')
    nuisance.connect(build_nuisance_regressors, 'out_file', outputspec, 'regressors_file_path')

    return nuisance

