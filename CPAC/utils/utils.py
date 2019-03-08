import os
import fnmatch
import threading
from inspect import currentframe, getframeinfo , stack

global_lock = threading.Lock()

files_folders_wf = {
    'anatomical_brain': 'anat',
    'anatomical_brain_mask': 'anat',
    'qc': 'qc',
    'anatomical_reorient': 'anat',
    'anatomical_to_mni_linear_xfm': 'anat',
    'mni_to_anatomical_linear_xfm': 'anat',
    'mni_to_anatomical_nonlinear_xfm': 'anat',
    'anatomical_to_mni_nonlinear_xfm': 'anat',
    'anatomical_gm_mask': 'anat',
    'anatomical_csf_mask': 'anat',
    'anatomical_wm_mask': 'anat',
    'ants_initial_xfm': 'anat',
    'ants_rigid_xfm': 'anat',
    'ants_affine_xfm': 'anat',
    'mean_functional': 'func',
    'functional_preprocessed_mask': 'func',
    'functional_to_spatial_map': 'func',
    'functional_mask_to_spatial_map': 'func',
    'fmap_phase_diff': 'func',
    'fmap_magnitude': 'func',
    'functional_distortion_corrected': 'func',
    'despiked_fieldmap': 'func',
    'prepared_fieldmap_map': 'func',
    'fieldmap_mask': 'func',
    'slice_time_corrected': 'func',
    'slice_timing_corrected': 'func',
    'movement_parameters': 'parameters',
    'max_displacement': 'parameters',
    'xform_matrix': 'parameters',
    'output_means': 'parameters',
    'functional_preprocessed': 'func',
    'functional_brain_mask': 'func',
    'motion_correct': 'func',
    'motion_correct_smooth': 'func',
    'motion_correct_to_standard': 'func',
    'motion_correct_to_standard_smooth': 'func',
    'mean_functional_in_anat': 'func',
    'coordinate_transformation': 'func',
    'raw_functional': 'func',
    'selected_func_volume': 'func',
    'anatomical_wm_edge': 'registration',
    'anatomical_to_functional_xfm': 'registration',
    'inverse_anatomical_to_functional_xfm': 'registration',
    'functional_gm_mask': 'segmentation',
    'functional_wm_mask': 'segmentation',
    'functional_csf_mask': 'segmentation',
    'frame_wise_displacement_power': 'parameters',
    'frame_wise_displacement_jenkinson': 'parameters',
    'functional_nuisance_residuals': 'func',
    'functional_nuisance_regressors': 'func',
    'functional_median_angle_corrected': 'func',
    'power_spectrum_distribution': 'alff',
    'functional_freq_filtered': 'func',
    'scrubbing_movement_parameters': 'parameters',
    'despiking_frames_included': 'parameters',
    'despiking_frames_excluded': 'parameters',
    'scrubbing_frames_included': 'parameters',
    'scrubbing_frames_excluded': 'parameters',
    'motion_params': 'parameters',
    'power_params': 'parameters',
    'scrubbed_preprocessed': 'func',
    'functional_to_standard': 'func',
    'functional_brain_mask_to_standard': 'func',
    'mean_functional_to_standard': 'func',
    'functional_to_anat_linear_xfm': 'registration',
    'functional_to_mni_linear_xfm': 'registration',
    'mni_to_functional_linear_xfm': 'registration',
    'ants_symmetric_initial_xfm': 'registration',
    'ants_symmetric_rigid_xfm': 'registration',
    'ants_symmetric_affine_xfm': 'registration',
    'anatomical_to_symmetric_mni_nonlinear_xfm': 'registration',
    'symmetric_mni_to_anatomical_nonlinear_xfm': 'registration',
    'symmetric_mni_to_anatomical_linear_xfm': 'registration',
    'anat_to_symmetric_mni_ants_composite_xfm': 'registration',
    'symmetric_anatomical_to_standard': 'registration',
    'anatomical_to_symmetric_mni_linear_xfm': 'registration',
    'anatomical_to_standard': 'anat',
    'leaf_node_to_standard': 'func',
    'vmhc_raw_score': 'vmhc',
    'vmhc_fisher_zstd': 'vmhc',
    'vmhc_fisher_zstd_zstat_map': 'vmhc',
    'alff': 'alff',
    'alff_input_functional': 'alff',
    'falff': 'alff',
    'alff_smooth': 'alff',
    'falff_smooth': 'alff',
    'alff_to_standard': 'alff',
    'falff_to_standard': 'alff',
    'alff_to_standard_smooth': 'alff',
    'falff_to_standard_smooth': 'alff',
    'alff_to_standard_zstd': 'alff',
    'falff_to_standard_zstd': 'alff',
    'alff_to_standard_smooth_zstd': 'alff',
    'falff_to_standard_smooth_zstd': 'alff',
    'alff_to_standard_zstd_smooth': 'alff',
    'falff_to_standard_zstd_smooth': 'alff',
    'alff_input_functional': 'alff',
    'reho': 'reho',
    'reho_smooth': 'reho',
    'reho_to_standard': 'reho',
    'reho_to_standard_smooth': 'reho',
    'reho_to_standard_zstd': 'reho',
    'reho_to_standard_smooth_zstd': 'reho',
    'reho_to_standard_zstd_smooth': 'reho',
    'connectome_PearsonCorr': 'connectome',
    'connectome_PartialCorr': 'connectome',
    'voxel_timeseries': 'timeseries',
    'roi_timeseries': 'timeseries',
    'roi_timeseries_for_SCA': 'timeseries',
    'roi_timeseries_for_SCA_multreg': 'timeseries',
    'sca_roi_files': 'sca_roi',
    'sca_roi_files_smooth': 'sca_roi',
    'sca_roi_files_to_standard': 'sca_roi',
    'sca_roi_files_to_standard_smooth': 'sca_roi',
    'sca_roi_files_to_standard_fisher_zstd': 'sca_roi',
    'sca_roi_files_to_standard_smooth_fisher_zstd': 'sca_roi',
    'sca_roi_files_to_standard_fisher_zstd_smooth': 'sca_roi',
    'bbregister_registration': 'surface_registration',
    'left_hemisphere_surface': 'surface_registration',
    'right_hemisphere_surface': 'surface_registration',
    'vertices_timeseries': 'timeseries',
    'centrality': 'centrality',
    'centrality_smooth': 'centrality',
    'centrality_zstd': 'centrality',
    'centrality_smooth_zstd': 'centrality',
    'centrality_zstd_smooth': 'centrality',
    'centrality_graphs': 'centrality',
    'seg_probability_maps': 'anat',
    'seg_mixeltype': 'anat',
    'seg_partial_volume_map': 'anat',
    'seg_partial_volume_files': 'anat',
    'spatial_map_timeseries': 'timeseries',
    'spatial_map_timeseries_for_DR': 'timeseries',
    'dr_tempreg_maps_files': 'spatial_regression',
    'dr_tempreg_maps_files_smooth': 'spatial_regression',
    'dr_tempreg_maps_zstat_files': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_smooth': 'spatial_regression',
    'dr_tempreg_maps_files_to_standard': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_to_standard': 'spatial_regression',
    'dr_tempreg_maps_files_to_standard_smooth': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_to_standard_smooth': 'spatial_regression',
    'sca_tempreg_maps_files': 'sca_roi',
    'sca_tempreg_maps_files_smooth': 'sca_roi',
    'sca_tempreg_maps_zstat_files': 'sca_roi',
    'sca_tempreg_maps_zstat_files_smooth': 'sca_roi',
    'lesion_reorient': 'lesion',
    'lesion_mask': 'lesion'
}


def get_zscore(input_name, map_node=False, wf_name='z_score'):
    """
    Workflow to calculate z-scores

    Parameters
    ----------
    wf_name : string
        name of the workflow

    Returns
    -------
    wf : workflow object

    Notes
    -----
    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/network_centrality/z_score.py>`_


    Workflow Inputs::

        inputspec.input_file : string
            path to input functional derivative file for which z score has to be calculated
        inputspec.mask_file : string
            path to whole brain functional mask file required to calculate zscore

    Workflow Outputs::

        outputspec.z_score_img : string
             path to image containing Normalized Input Image Z scores across full brain.

    High Level Workflow Graph:

    .. image:: ../images/zscore.dot.png
       :width: 500


    Detailed Workflow Graph:

    .. image:: ../images/zscore_detailed.dot.png
       :width: 500

    Example
    -------
    >>> import get_zscore as z
    >>> wf = z.get_zscore()
    >>> wf.inputs.inputspec.input_file = '/home/data/graph_working_dir/calculate_centrality/degree_centrality_binarize.nii.gz'
    >>> wf.inputs.inputspec.mask_file = '/home/data/graphs/GraphGeneration/new_mask_3m.nii.gz'
    >>> wf.run()

    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['input_file',
                                                       'mask_file']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['z_score_img']),
                         name='outputspec')

    if map_node:
        mean = pe.MapNode(interface=fsl.ImageStats(),
                          name='mean',
                          iterfield=['in_file'])

        standard_deviation = pe.MapNode(interface=fsl.ImageStats(),
                                        name='standard_deviation',
                                        iterfield=['in_file'])

        op_string = pe.MapNode(util.Function(input_names=['mean', 'std_dev'],
                                             output_names=['op_string'],
                                             function=get_operand_string),
                               name='op_string',
                               iterfield=['mean', 'std_dev'])

        z_score = pe.MapNode(interface=fsl.MultiImageMaths(),
                             name='z_score',
                             iterfield=['in_file', 'op_string'])

    else:
        mean = pe.Node(interface=fsl.ImageStats(), name='mean')

        standard_deviation = pe.Node(interface=fsl.ImageStats(),
                                     name='standard_deviation')

        op_string = pe.Node(util.Function(input_names=['mean', 'std_dev'],
                                          output_names=['op_string'],
                                          function=get_operand_string),
                            name='op_string')

        z_score = pe.Node(interface=fsl.MultiImageMaths(), name='z_score')

    # calculate the mean
    mean.inputs.op_string = '-k %s -m'
    wflow.connect(inputNode, 'input_file', mean, 'in_file')
    wflow.connect(inputNode, 'mask_file', mean, 'mask_file')

    # calculate the standard deviation
    standard_deviation.inputs.op_string = '-k %s -s'
    wflow.connect(inputNode, 'input_file', standard_deviation, 'in_file')
    wflow.connect(inputNode, 'mask_file', standard_deviation, 'mask_file')

    # calculate the z-score
    wflow.connect(mean, 'out_stat', op_string, 'mean')
    wflow.connect(standard_deviation, 'out_stat', op_string, 'std_dev')

    #z_score.inputs.out_file = input_name + '_zstd.nii.gz'

    wflow.connect(op_string, 'op_string', z_score, 'op_string')
    wflow.connect(inputNode, 'input_file', z_score, 'in_file')
    wflow.connect(inputNode, 'mask_file', z_score, 'operand_files')

    wflow.connect(z_score, 'out_file', outputNode, 'z_score_img')

    return wflow


def get_fisher_zscore(input_name, map_node=False, wf_name='fisher_z_score'):
    """
    Runs the compute_fisher_z_score function as part of a one-node workflow.
    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl

    wflow = pe.Workflow(name=wf_name)

    inputNode = pe.Node(util.IdentityInterface(fields=['correlation_file',
                                                       'timeseries_one_d']),
                        name='inputspec')

    outputNode = pe.Node(
        util.IdentityInterface(fields=['fisher_z_score_img']),
        name='outputspec')

    if map_node:
        # node to separate out
        fisher_z_score = pe.MapNode(
            util.Function(input_names=['correlation_file',
                                       'timeseries_one_d',
                                       'input_name'],
                          output_names=['out_file'],
                          function=compute_fisher_z_score),
            name='fisher_z_score',
            iterfield=['correlation_file'])
    else:
        fisher_z_score = pe.Node(
            util.Function(input_names=['correlation_file',
                                       'timeseries_one_d',
                                       'input_name'],
                          output_names=['out_file'],
                          function=compute_fisher_z_score),
            name='fisher_z_score')

    fisher_z_score.inputs.input_name = input_name

    wflow.connect(inputNode, 'correlation_file',
                  fisher_z_score, 'correlation_file')
    wflow.connect(inputNode, 'timeseries_one_d',
                  fisher_z_score, 'timeseries_one_d')
    wflow.connect(fisher_z_score, 'out_file',
                  outputNode, 'fisher_z_score_img')

    return wflow


def compute_fisher_z_score(correlation_file, timeseries_one_d, input_name):
    """
    Computes the fisher z transform of the input correlation map
    If the correlation map contains data for multiple ROIs then
    the function returns z score for each ROI as a seperate nifti
    file


    Parameters
    ----------

    correlation_file: string
        Input correlations file


    Returns
    -------

    out_file : list (nifti files)
        list of z_scores for mask or ROI
    """

    import nibabel as nb
    import numpy as np
    import os

    if isinstance(timeseries_one_d, basestring):
        if '.1D' in timeseries_one_d or '.csv' in timeseries_one_d:
            timeseries_file = timeseries_one_d

    else:
        for timeseries in timeseries_one_d:
            if '.1D' in timeseries or '.csv' in timeseries:
                timeseries_file = timeseries

    # get the specific roi number
    filename = correlation_file.split("/")[-1]
    filename = filename.replace(".nii", "")
    if ".gz" in filename:
        filename = filename.replace(".gz", "")

    corr_img = nb.load(correlation_file)
    corr_data = corr_img.get_data()

    hdr = corr_img.get_header()

    # calculate the Fisher r-to-z transformation
    corr_data = np.log((1 + corr_data) / (1 - corr_data)) / 2.0

    z_score_img = nb.Nifti1Image(corr_data, header=hdr,
                                 affine=corr_img.get_affine())

    out_file = os.path.join(os.getcwd(), filename + '_fisher_zstd.nii.gz')

    z_score_img.to_filename(out_file)

    return out_file


def get_operand_string(mean, std_dev):
    """
    Method to get operand string for Fsl Maths

    Parameters
    ----------
    mean : string
        path to img containing mean
    std_dev : string
        path to img containing standard deviation

    Returns
    ------
    op_string : string
        operand string
    """

    str1 = "-sub %f -div %f" % (float(mean), float(std_dev))
    op_string = str1 + " -mas %s"
    return op_string


def get_roi_num_list(timeseries_file, prefix=None):
    # extracts the ROI labels from the 3dROIstats output CSV file
    with open(timeseries_file, "r") as f:
        roi_file_lines = f.read().splitlines()

    roi_err = "\n\n[!] The output of 3dROIstats, used in extracting the " \
              "timeseries, is either empty, or not in the expected " \
              "format.\n\nROI output file: {0}\n\nIf there are no rows " \
              "in the output file, double-check your ROI/mask selection." \
              "\n\n".format(str(timeseries_file))

    for line in roi_file_lines:
        if "Mean_" in line:
            try:
                roi_list = line.split("\t")
                # clear out any blank strings/non ROI labels in the list
                roi_list = [x for x in roi_list if "Mean" in x]
                # rename labels
                roi_list = [x.replace("Mean", "ROI").replace(" ", "") \
                            for x in roi_list]
            except:
                raise Exception(roi_err)
            break
    else:
        raise Exception(roi_err)

    if prefix:
        temp_rois = []
        for roi in roi_list:
            roi = prefix + "_" + str(roi)
            temp_rois.append(roi)
        roi_list = temp_rois

    return roi_list


def safe_shape(*vol_data):
    """
    Checks if the volume (first three dimensions) of multiple ndarrays
    are the same shape.

    Parameters
    ----------
    vol_data0, vol_data1, ..., vol_datan : ndarray
        Volumes to check

    Returns
    -------
    same_volume : bool
        True only if all volumes have the same shape.
    """
    same_volume = True

    first_vol_shape = vol_data[0].shape[:3]
    for vol in vol_data[1:]:
        same_volume &= (first_vol_shape == vol.shape[:3])

    return same_volume


def extract_one_d(list_timeseries):
    if isinstance(list_timeseries, basestring):
        if '.1D' in list_timeseries or '.csv' in list_timeseries:
            return list_timeseries

    for timeseries in list_timeseries:
        if '.1D' in timeseries or '.csv' in timeseries:
            return timeseries

    raise Exception("Unable to retrieve roi timeseries 1D or csv" \
                    " file. Files found:" + list_timeseries)


def extract_txt(list_timeseries):
    """
    Method to extract txt file containing
    roi timeseries required for dual regression
    """
    if isinstance(list_timeseries, basestring):
        if list_timeseries.endswith('.txt'):
            return list_timeseries

    out_file = None
    for timeseries in list_timeseries:
        if timeseries.endswith('.txt'):
            out_file = timeseries

    if not out_file:
        raise Exception("Unable to retrieve roi timeseries txt" \
                        " file required for dual regression." \
                        " Existing files are:%s" % (list_timeseries))

    return out_file


def set_gauss(fwhm):
    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string


def check(params_dct, subject, scan, val, throw_exception):

    if val not in params_dct:

        if throw_exception:
            raise Exception("Missing Value for {0} for subject "
                            "{1}".format(val, subject))

        return None

    if isinstance(params_dct[val], dict):
        ret_val = params_dct[val][scan]
    else:
        ret_val = params_dct[val]

    if ret_val == 'None':
        if throw_exception:
            raise Exception("None Parameter Value for {0} for subject "
                            "{1}".format(val, subject))
        else:
            ret_val = None

    if ret_val == '' and throw_exception:
        raise Exception("Missing Value for {0} for subject "
                        "{1}".format(val, subject))

    return ret_val


def try_fetch_parameter(scan_parameters, subject, scan, keys):
    
    scan_parameters = dict(
        (k.lower(), v)
        for k, v in scan_parameters.iteritems()
    )

    for key in keys:

        key = key.lower()
        
        if key not in scan_parameters:
            continue

        if isinstance(scan_parameters[key], dict):
            value = scan_parameters[key][scan]
        else:
            value = scan_parameters[key]

        # Explicit none value
        if value == 'None':
            return None

        if value is not None:
            return value

    return None
    #raise Exception("Missing Value for {0} for subject "
    #                "{1}".format(' or '.join(keys), subject))


def get_scan_params(subject_id, scan, pipeconfig_tr, pipeconfig_tpattern,
                    pipeconfig_start_indx, pipeconfig_stop_indx,
                    data_config_scan_params=None):
    """
    Method to extract slice timing correction parameters
    and scan parameters.

    Parameters
    ----------
    subject_id: a string
        subject id
    scan : a string
        scan id
    subject_map : a dictionary
        subject map containing all subject information
    start_indx : an integer
        starting volume index
    stop_indx : an integer
        ending volume index

    Returns
    -------
    TR : a string
        TR value
    pattern : a string
        slice aquisition pattern string or file path
    ref_slice : an integer
        reference slice which is used to allign all other slices
    first_tr : an integer
        starting TR or starting volume index
    last_tr : an integer
        ending TR or ending volume index
    """

    import os
    import json
    import warnings

    check2 = lambda val: val if val == None or val == '' or \
                                isinstance(val, str) else int(val)

    # initialize vars to empty
    TR = ''
    TE = None
    pattern = ''
    ref_slice = ''
    first_tr = ''
    last_tr = ''
    unit = 's'

    if isinstance(pipeconfig_tpattern, list) or isinstance(pipeconfig_tpattern, str):
        if "None" in pipeconfig_tpattern:
            pipeconfig_tpattern = None

    if isinstance(pipeconfig_tr, str):
        if "None" in pipeconfig_tr or "none" in pipeconfig_tr:
            pipeconfig_tr = None

    if isinstance(pipeconfig_stop_indx, str):
        if "End" in pipeconfig_stop_indx or "end" in pipeconfig_stop_indx:
            pipeconfig_stop_indx = None

    if data_config_scan_params:

        if ".json" in data_config_scan_params:

            if not os.path.exists(data_config_scan_params):
                err = "\n[!] WARNING: Scan parameters JSON file listed in " \
                      "your data configuration file does not exist:\n{0}" \
                      "\n\n".format(data_config_scan_params)
                raise Exception(err)

            with open(data_config_scan_params, "r") as f:
                params_dct = json.load(f)

            # get details from the configuration
            # if this is a JSON file, the key values are the BIDS format
            # standard
            # TODO: better handling of errant key values!!!
            if "RepetitionTime" in params_dct.keys():
                TR = float(check(params_dct, subject_id, scan,
                                 'RepetitionTime', False))
            if "SliceTiming" in params_dct.keys():
                pattern = str(check(params_dct, subject_id, scan,
                                    'SliceTiming', False))
            elif "SliceAcquisitionOrder" in params_dct.keys():
                pattern = str(check(params_dct, subject_id, scan,
                                    'SliceAcquisitionOrder', False))

        elif len(data_config_scan_params) > 0 and \
                isinstance(data_config_scan_params, dict):

            try:
                params_dct = data_config_scan_params
            except:
                err = "\n[!] Could not parse the scan parameter information "\
                      "included in your data configuration file for " \
                      "participant: {0}\n\n".format(subject_id)
                raise Exception(err)

            # TODO: better handling of errant key values!!!
            # TODO: use schema validator to deal with it
            # get details from the configuration
            TR = float(
                try_fetch_parameter(
                    params_dct,
                    subject_id, 
                    scan, 
                    ['TR', 'RepetitionTime']
                )
            )

            pattern = str(
                try_fetch_parameter(
                    params_dct,
                    subject_id,
                    scan,
                    ['acquisition', 'SliceTiming', 'SliceAcquisitionOrder']
                )
            )
            
            ref_slice = check(params_dct, subject_id, scan, 'reference', False)
            if ref_slice:
                ref_slice = int(ref_slice)

            first_tr = check(params_dct, subject_id, scan, 'first_TR', False)
            if first_tr:
                first_tr = check2(first_tr)

            last_tr = check(params_dct, subject_id, scan, 'last_TR', False)
            if last_tr:
                last_tr = check2(last_tr)

        else:
            err = "\n\n[!] Could not read the format of the scan parameters "\
                  "information included in the data configuration file for " \
                  "the participant {0}.\n\n".format(subject_id)
            raise Exception(err)

    # if values are still empty, override with GUI config
    if TR == '':
        if pipeconfig_tr:
            TR = float(pipeconfig_tr)
        else:
            TR = None

    if first_tr == '':
        first_tr = pipeconfig_start_indx

    if last_tr == '':
        last_tr = pipeconfig_stop_indx

    unit = 's'

    if 'None' in pattern or 'none' in pattern:
        pattern = None

    if not pattern:
        if pipeconfig_tpattern:
            if "Use NIFTI Header" in pipeconfig_tpattern:
                pattern = ''
            else:
                pattern = pipeconfig_tpattern

    # pattern can be one of a few keywords, a filename, or blank which
    # indicates that the images header information should be used
    tpattern_file = None

    valid_patterns = ['alt+z', 'altplus', 'alt+z2', 'alt-z', 'altminus',
                      'alt-z2', 'seq+z', 'seqplus', 'seq-z', 'seqminus']

    if pattern and pattern != '' and pattern not in valid_patterns:

        if isinstance(pattern, list) or \
                ("[" in pattern and "]" in pattern and "," in pattern):
            # if we got the slice timing as a list, from a BIDS-format scan
            # parameters JSON file

            if not isinstance(pattern, list):
                pattern = pattern.replace("[", "").replace("]", "").split(",")

            slice_timings = [float(x) for x in pattern]

            # write out a tpattern file for AFNI 3dTShift
            tpattern_file = os.path.join(os.getcwd(), "tpattern.txt")
            try:
                with open(tpattern_file, "wt") as f:
                    for time in slice_timings:
                        f.write("{0}\n".format(time).replace(" ", ""))
            except:
                err = "\n[!] Could not write the slice timing file meant as "\
                      "an input for AFNI 3dTshift (slice timing correction):"\
                      "\n{0}\n\n".format(tpattern_file)
                raise Exception(err)

        elif ".txt" in pattern and not os.path.exists(pattern):
            # if the user provided an acquisition pattern text file for
            # 3dTshift
            raise Exception("Invalid Pattern file path {0}, Please provide "
                            "the correct path".format(pattern))
        elif ".txt" in pattern:
            with open(pattern, "r") as f:
                lines = f.readlines()
            if len(lines) < 2:
                raise Exception('Invalid slice timing file format. The file '
                                'should contain only one value per row. Use '
                                'new line char as delimiter')
            tpattern_file = pattern
            slice_timings = [float(l.rstrip('\r\n')) for l in lines]
        else:
            # this only happens if there is a non-path string set in the data
            # config dictionary for acquisition pattern (like "alt+z"), except
            # the pattern is not listed in that list
            err = "\n[!] The slice timing acquisition pattern provided is " \
                  "not supported by AFNI 3dTshift:\n" \
                  "{0}\n".format(str(pattern))
            raise Exception(err)

        pattern = tpattern_file

        slice_timings.sort()
        max_slice_offset = slice_timings[-1]

        # checking if the unit of TR and slice timing match or not
        # if slice timing in ms convert TR to ms as well
        if TR and max_slice_offset > TR:
            warnings.warn("TR is in seconds and slice timings are in "
                          "milliseconds. Converting TR into milliseconds")
            TR = TR * 1000
            print("New TR value {0} ms".format(TR))
            unit = 'ms'

    else:
        # check to see, if TR is in milliseconds, convert it into seconds
        if TR and TR > 10:
            warnings.warn('TR is in milliseconds, Converting it into seconds')
            TR = TR / 1000.0
            print("New TR value {0} s".format(TR))
            unit = 's'

    print("scan_parameters -> {0} {1} {2} {3} {4} "
          "{5} {6}".format(subject_id, scan, str(TR) + unit, pattern,
                           ref_slice, first_tr, last_tr))

    # swap back in
    if TR:
        tr = "{0}{1}".format(str(TR), unit)
    else:
        tr = ""
    tpattern = pattern
    start_indx = first_tr
    stop_indx = last_tr

    return tr, tpattern, ref_slice, start_indx, stop_indx


def get_tr(tr):
    """
    Method to return TR in seconds
    """
    import re
    if tr:
        tr = re.search("\d+.\d+", str(tr)).group(0)
        tr = float(tr)
        if tr > 10:
            tr = tr / 1000.0
    else:
        tr = ""
    return tr


def check_tr(tr, in_file):
    # imageData would have to be the image data from the funcFlow workflow,
    # funcFlow outputspec.subject
    import nibabel as nib
    img = nib.load(in_file)

    # get header from image data, then extract TR information, TR is fourth
    # item in list returned by get_zooms()
    imageHeader = img.get_header()
    imageZooms = imageHeader.get_zooms()
    header_tr = imageZooms[3]

    # If the TR information from header_tr (funcFlow) and convert_tr node
    # (TR from config file) do not match, prepare to update the TR information
    # from either convert_tr or header_tr using afni 3drefit, then append to
    # func_to_mni
    if header_tr != tr:
        if tr != None and tr != "":
            TR = tr
        else:
            TR = header_tr

        import warnings
        warnings.warn(
            'Warning: The TR information does not match between the config and subject list files.')

    return TR


def add_afni_prefix(tpattern):
    if ".txt" in tpattern:
        tpattern = "@{0}".format(tpattern)
    return tpattern


def write_to_log(workflow, log_dir, index, inputs, scan_id):
    """
    Method to write into log file the status of the workflow run.
    """

    import os
    import time
    import datetime
    
    from CPAC import __version__
    from nipype import logging

    iflogger = logging.getLogger('nipype.interface')

    version = __version__
    subject_id = os.path.basename(log_dir)

    if scan_id is None:
        scan_id = "scan_anat"

    strategy = ""
    ts = time.time()
    stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    try:
        if workflow != 'DONE':
            wf_path = \
                os.path.dirname((os.getcwd()).split(workflow)[1]).strip("/")

            if wf_path and wf_path != "":
                if '/' in wf_path:
                    scan_id, strategy = wf_path.split('/', 1)
                    scan_id = scan_id.strip('_')
                    strategy = strategy.replace("/", "")
                else:
                    scan_id = wf_path.strip('_')

            file_path = os.path.join(log_dir, scan_id, workflow)

            try:
                os.makedirs(file_path)
            except Exception:
                iflogger.info(
                    "filepath already exist, filepath- {0}, "
                    "curr_dir - {1}".format(file_path, os.getcwd()))

        else:
            file_path = os.path.join(log_dir, scan_id)
    except Exception:
        print("ERROR in write log")
        raise

    try:
        os.makedirs(file_path)
    except Exception:
        iflogger.info(
            "filepath already exist, "
            "filepath: {0}, "
            "curr_dir: {1}".format(file_path, os.getcwd())
        )

    out_file = os.path.join(file_path, 'log_{0}.yaml'.format(strategy))

    iflogger.info("CPAC custom log:")

    if isinstance(inputs, list):
        inputs = inputs[0]

    if os.path.exists(inputs):
        status_msg = "wf_status: DONE"
        iflogger.info(
            "version: {0}, "
            "timestamp: {1}, "
            "subject_id: {2}, "
            "scan_id: {3}, "
            "strategy: {4}, "
            "workflow: {5}, "
            "status: COMPLETED".format(
                str(version), str(stamp), subject_id,
                scan_id, strategy, workflow
            )
        )
    else:
        status_msg = "wf_status: ERROR"
        iflogger.info(
            "version: {0}, "
            "timestamp: {1}, "
            "subject_id: {2}, "
            "scan_id: {3}, "
            "strategy: {4}, "
            "workflow: {5}, "
            "status: ERROR".format(
                str(version), str(stamp), subject_id,
                scan_id, strategy, workflow
            )
        )

    with open(out_file, 'w') as f:
        f.write("version: {0}\n".format(str(version)))
        f.write("timestamp: {0}\n".format(str(stamp)))
        f.write("pipeline_index: {0}\n".format(index))
        f.write("subject_id: {0}\n".format(subject_id))
        f.write("scan_id: {0}\n".format(scan_id))
        f.write("strategy: {0}\n".format(strategy))
        f.write("workflow_name: {0}\n".format(workflow))
        f.write(status_msg)

    return out_file


def create_log(wf_name="log", scan_id=None):
    """
    Workflow to create log
    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import CPAC.utils.function as function

    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['workflow',
                                                        'log_dir',
                                                        'index',
                                                        'inputs']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['out_file']),
                          name='outputspec')

    write_log = pe.Node(function.Function(input_names=['workflow',
                                                       'log_dir',
                                                       'index',
                                                       'inputs',
                                                       'scan_id'],
                                          output_names=['out_file'],
                                          function=write_to_log,
                                          as_module=True),
                        name='write_log')

    write_log.inputs.scan_id = scan_id

    wf.connect([
        (
            input_node, write_log, [
                ('workflow', 'workflow'),
                ('log_dir', 'log_dir'),
                ('index', 'index'),
                ('inputs', 'inputs')
            ]
        ),
        (
            write_log, output_node, [
                ('out_file', 'out_file')
            ]
        )
    ])

    return wf


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def extract_output_mean(in_file, output_name):
    '''
    function takes 'in_file', which should be an intermediary 1D file
    from individual-level analysis, containing the mean of the output across
    all voxels

    it then parses this value and writes it to a .csv file named
    output_means.csv located in the subject's output directory
    '''

    if os.path.exists(in_file):

        with open(in_file, 'r') as f:
            line = f.readline()

        line = line.split('[')[0].strip(' ')

        # get filename of input maskave 1D file
        filename = in_file.split("/")[-1]
        filename = filename[0:-3]

        split_fullpath = in_file.split("/")

        if "_mask_" in in_file and \
           ("sca_roi" in in_file or "sca_tempreg" in in_file):

            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname

            filename = split_fullpath[-1]

            if ".1D" in filename:
                filename = filename.replace(".1D", "")

            resource_name = "{0}_{1}_{2}".format(output_name, maskname,
                                                 filename)

        elif "_spatial_map_" in in_file and "dr_tempreg" in in_file:

            for dirname in split_fullpath:
                if "_spatial_map_" in dirname:
                    mapname = dirname

            filename = split_fullpath[-1]

            if ".1D" in filename:
                filename = filename.replace(".1D", "")

            resource_name = "{0}_{1}_{2}".format(output_name, mapname,
                                                 filename)

        elif "_mask_" in in_file and "centrality" in in_file:

            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname

            filename = split_fullpath[-1]

            if ".1D" in filename:
                filename = filename.replace(".1D", "")

            resource_name = "{0}_{1}_{2}".format(output_name, maskname,
                                                 filename)

        else:
            resource_name = output_name

        output_means_file = os.path.join(os.getcwd(),
                                         'mean_{0}.txt'.format(resource_name))

        with open(output_means_file, 'w') as f:
            f.write(line)

    return output_means_file


def create_output_mean_csv(subject_dir):
    '''
    this function finds all of the mean_{output}.txt files in the subject's
    output directory, collects the data and organizes them into one .csv
    file in the subject directory
    '''

    import os
    import csv

    output_vals = {}

    subID = subject_dir.split('/')[len(subject_dir.split('/')) - 1]
    means_dir = os.path.join(subject_dir, 'output_means')

    # extract the mean values
    for root, _, files in os.walk(means_dir):

        for filename in files:

            if 'mean_' in filename:

                output = filename.replace('mean_', '')
                output = output.replace('.txt', '')

                filepath = os.path.join(root, filename)

                if os.path.exists(filepath):
                    try:
                        mean_file = open(filepath, 'rU')
                        val = mean_file.readline()
                        val = val.strip('\n')
                    except:
                        print '\n\n[!] CPAC says: Could not open the output ' \
                              'mean text file.\n'
                        print 'Path: ', filepath, '\n\n'
                        raise Exception

                else:
                    print '\n\n[!] CPAC says: Could not find the output mean ' \
                          'text file.\n'
                    print 'Path not found: ', filepath, '\n\n'
                    raise Exception

                output_vals[output] = val

    # now take the extracted mean values and write them into the .csv file!
    csv_file_path = os.path.join(subject_dir, 'output_means_%s.csv' % subID)
    with open(csv_file_path, 'wt') as csv_file:

        output_items = list(output_vals.items())

        deriv_string = ','.join(v for v, _ in output_items)
        val_string = ','.join(v for _, v in output_items)

        csv_file.write(deriv_string + '\n')
        csv_file.write(val_string + '\n')


# Setup log file
def setup_logger(logger_name, file_path, level, to_screen=False):
    '''
    Function to initialize and configure a logger that can write to file
    and (optionally) the screen.

    Parameters
    ----------
    logger_name : string
        name of the logger
    file_path : string
        file path to the log file on disk
    level : integer
        indicates the level at which the logger should log; this is
        controlled by integers that come with the python logging
        package. (e.g. logging.INFO=20, logging.DEBUG=10)
    to_screen : boolean (optional)
        flag to indicate whether to enable logging to the screen

    Returns
    -------
    logger : logging.Logger object
        Python logging.Logger object which is capable of logging run-
        time information about the program to file and/or screen
    '''

    # Import packages
    import logging

    # Init logger, formatter, filehandler, streamhandler
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    formatter = logging.Formatter('%(asctime)s : %(message)s')

    # Write logs to file
    fileHandler = logging.FileHandler(file_path)
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    # Write to screen, if desired
    if to_screen:
        streamHandler = logging.StreamHandler()
        streamHandler.setFormatter(formatter)
        logger.addHandler(streamHandler)

    # Return the logger
    return logger


def check_system_deps(check_ants=False, check_ica_aroma=False):
    '''
    Function to check system for neuroimaging tools AFNI, C3D, FSL,
    and (optionally) ANTs
    '''

    # Import packages
    import os

    # Init variables
    missing_install = []

    # Check AFNI
    if os.system("3dcalc >/dev/null 2>&1") == 32512:
        missing_install.append("AFNI")
    # Check FSL
    if os.system("fslmaths >/dev/null 2>&1") == 32512:
        missing_install.append("FSL")
    # Check ANTs/C3D
    if check_ants:
        if os.system("c3d_affine_tool >/dev/null 2>&1") == 32512:
            missing_install.append("C3D")
        if os.system("antsRegistration >/dev/null 2>&1") == 32512:
            missing_install.append("ANTS")
    # Check ICA-AROMA
    if check_ica_aroma:
        if os.system("ICA_AROMA.py >/dev/null 2>&1") == 32512:
            missing_install.append("ICA-AROMA")

    # If we're missing deps, raise Exception
    if len(missing_install) > 0:
        missing_string = ""
        for string in missing_install:
            missing_string = missing_string + string + "\n"
        err = "\n\n[!] CPAC says: It appears the following software " \
              "packages are not installed or configured properly:\n\n%s\n" \
              "Consult the CPAC Installation Guide for instructions.\n\n" \
              % missing_string
        raise Exception(err)


# Check pipeline config againts computer resources
def check_config_resources(c):
    '''
    docstring
    '''

    # Import packages
    import psutil
    from multiprocessing import cpu_count

    # Init variables
    sys_virt_mem = psutil.virtual_memory()
    num_cores = cpu_count()

    # Check for pipeline memory for subject
    if c.maximumMemoryPerParticipant is None:
        # Get system memory and numSubsAtOnce
        sys_mem_gb = sys_virt_mem.total / (1024.0 ** 3)
        sub_mem_gb = sys_mem_gb / c.numParticipantsAtOnce
    else:
        sub_mem_gb = c.maximumMemoryPerParticipant

    # If centrality is enabled, check to mem_sub >= mem_centrality
    if c.runNetworkCentrality[0]:
        if sub_mem_gb < c.memoryAllocatedForDegreeCentrality:
            err_msg = 'Memory allocated for subject: %d needs to be greater ' \
                      'than the memory allocated for centrality: %d. Fix ' \
                      'and try again.' % (c.maximumMemoryPerParticipant,
                                          c.memoryAllocatedForDegreeCentrality)
            raise Exception(err_msg)

    # Check for pipeline threads
    # Check if user specified cores
    if c.maxCoresPerParticipant:
        total_user_cores = c.numParticipantsAtOnce * c.maxCoresPerParticipant
        if total_user_cores > num_cores:
            err_msg = 'Config file specifies more subjects running in ' \
                      'parallel than number of threads available. Change ' \
                      'this and try again'
            raise Exception(err_msg)
        else:
            num_cores_per_sub = c.maxCoresPerParticipant
    else:
        num_cores_per_sub = num_cores / c.maxCoresPerParticipant

    # Now check ANTS
    if 'ANTS' in c.regOption:
        if c.num_ants_threads is None:
            num_ants_cores = num_cores_per_sub
        elif c.num_ants_threads > c.maxCoresPerParticipant:
            err_msg = 'Number of threads for ANTS: %d is greater than the ' \
                      'number of threads per subject: %d. Change this and ' \
                      'try again.' % (c.num_ants_threads,
                                      c.maxCoresPerParticipant)
            raise Exception(err_msg)
        else:
            num_ants_cores = c.num_ants_threads
    else:
        num_ants_cores = 1

    # Return memory and cores
    return sub_mem_gb, num_cores_per_sub, num_ants_cores
