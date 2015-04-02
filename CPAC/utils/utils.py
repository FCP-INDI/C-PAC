from inspect import currentframe, getframeinfo

import threading
global_lock = threading.Lock()


files_folders_wf = {
    'anatomical_brain': 'anat',
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
    'slice_time_corrected': 'func',
    'slice_timing_corrected': 'func',
    'movement_parameters': 'parameters',
    'max_displacement':'parameters',
    'xform_matrix' : 'parameters',
    'output_means' : 'parameters',
    'preprocessed':'func',
    'functional_brain_mask':'func',
    'motion_correct':'func',
    'motion_correct_smooth':'func',
    'itk_func_anat_affine_motion_correct_to_standard':'func',
    'itk_collected_warps_motion_correct_to_standard':'func',
    'motion_correct_to_standard':'func',
    'motion_correct_to_standard_smooth':'func',
    'mean_functional_in_anat' : 'func',
    'coordinate_transformation' : 'func',
    'raw_functional' : 'func',
    'selected_func_volume' : 'func',
    'anatomical_wm_edge' : 'registration',
    'anatomical_to_functional_xfm':'registration',
    'inverse_anatomical_to_functional_xfm':'registration',
    'functional_gm_mask':'segmentation',
    'functional_wm_mask':'segmentation',
    'functional_csf_mask':'segmentation',
    'frame_wise_displacement':'parameters',
    'functional_nuisance_residuals':'func',
    'functional_median_angle_corrected':'func',
    'power_spectrum_distribution':'alff',
    'functional_freq_filtered':'func',
    'scrubbing_movement_parameters':'parameters',
    'scrubbing_frames_included':'parameters',
    'scrubbing_frames_excluded':'parameters',
    'motion_params':'parameters',
    'power_params':'parameters',
    'scrubbed_preprocessed':'func',
    'itk_func_anat_affine_functional_mni':'func',
    'itk_func_anat_affine_functional_brain_mask_to_standard':'func',
    'itk_func_anat_affine_mean_functional_in_mni' : 'func',
    'itk_collected_warps_functional_mni':'func',
    'itk_collected_warps_functional_brain_mask_to_standard':'func',
    'itk_collected_warps_mean_functional_in_mni' : 'func', 
    'functional_mni':'func',
    'functional_brain_mask_to_standard':'func',
    'mean_functional_in_mni' : 'func',
    'functional_to_anat_linear_xfm':'registration',
    'functional_to_mni_linear_xfm':'registration',
    'mni_to_functional_linear_xfm':'registration',
    'ants_symmetric_initial_xfm':'registration',
    'ants_symmetric_rigid_xfm':'registration',
    'ants_symmetric_affine_xfm':'registration',
    'anatomical_to_symmetric_mni_nonlinear_xfm':'registration',
    'symmetric_mni_to_anatomical_nonlinear_xfm':'registration',
    'symmetric_mni_to_anatomical_linear_xfm':'registration',
    'anat_to_symmetric_mni_ants_composite_xfm':'registration',
    'symmetric_mni_normalized_anatomical':'registration',
    'anatomical_to_symmetric_mni_linear_xfm':'registration',
    'mni_normalized_anatomical':'anat',
    'vmhc_raw_score':'vmhc',
    'vmhc_fisher_zstd':'vmhc',
    'vmhc_fisher_zstd_zstat_map':'vmhc',
    'alff_img':'alff',
    'falff_img':'alff',
    'alff_smooth':'alff',
    'falff_smooth':'alff',
    'alff_to_standard':'alff',
    'falff_to_standard':'alff',
    'alff_to_standard_smooth':'alff',
    'falff_to_standard_smooth':'alff',
    'alff_to_standard_zstd':'alff',
    'falff_to_standard_zstd':'alff',
    'alff_to_standard_smooth_zstd':'alff',
    'falff_to_standard_smooth_zstd':'alff',
    'raw_reho_map':'reho',
    'reho_smooth':'reho',
    'reho_to_standard':'reho',
    'reho_to_standard_smooth':'reho',
    'reho_to_standard_zstd':'reho',
    'reho_to_standard_smooth_zstd':'reho',
    'voxel_timeseries':'timeseries',
    'voxel_timeseries_for_SCA':'timeseries',
    'roi_timeseries':'timeseries',
    'roi_timeseries_for_SCA':'timeseries',
    'sca_seed_correlation_files':'sca_mask',
    'sca_seed_smooth':'sca_mask',
    'sca_seed_to_standard':'sca_mask',
    'sca_seed_to_standard_smooth':'sca_mask',
    'sca_seed_to_standard_fisher_zstd':'sca_mask',
    'sca_seed_to_standard_smooth_fisher_zstd':'sca_mask',
    'sca_roi_correlation_stack':'sca_roi',
    'sca_roi_correlation_files':'sca_roi',
    'sca_roi_stack_smooth':'sca_roi',
    'sca_roi_files_smooth':'sca_roi',
    'sca_roi_stack_to_standard':'sca_roi',
    'sca_roi_files_to_standard':'sca_roi',
    'sca_roi_stack_to_standard_smooth':'sca_roi',
    'sca_roi_files_to_standard_smooth':'sca_roi',
    'sca_roi_stack_to_standard_fisher_zstd':'sca_roi',
    'sca_roi_stack_to_standard_smooth_fisher_zstd':'sca_roi',
    'sca_roi_files_to_standard_fisher_zstd':'sca_roi',
    'sca_roi_files_to_standard_smooth_fisher_zstd':'sca_roi',
    'bbregister_registration': 'surface_registration',
    'left_hemisphere_surface': 'surface_registration',
    'right_hemisphere_surface': 'surface_registration',
    'vertices_timeseries': 'timeseries',
    'centrality_outputs':'centrality',
    'centrality_outputs_smoothed':'centrality',
    'centrality_outputs_zstd':'centrality',
    'centrality_outputs_smoothed_zstd': 'centrality',
    'centrality_graphs':'centrality',
    'seg_probability_maps': 'anat',
    'seg_mixeltype': 'anat',
    'seg_partial_volume_map': 'anat',
    'seg_partial_volume_files': 'anat',
    'spatial_map_timeseries': 'timeseries',
    'dr_tempreg_maps_stack': 'spatial_regression',
    'dr_tempreg_maps_files': 'spatial_regression',
    'dr_tempreg_maps_zstat_stack': 'spatial_regression',
    'dr_tempreg_maps_zstat_files': 'spatial_regression',
    'dr_tempreg_maps_stack_to_standard': 'spatial_regression',
    'dr_tempreg_maps_files_to_standard': 'spatial_regression',
    'dr_tempreg_maps_zstat_stack_to_standard': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_to_standard': 'spatial_regression',
    'dr_tempreg_maps_stack_to_standard_smooth': 'spatial_regression',
    'dr_tempreg_maps_files_to_standard_smooth': 'spatial_regression',
    'dr_tempreg_maps_zstat_stack_to_standard_smooth': 'spatial_regression',
    'dr_tempreg_maps_zstat_files_to_standard_smooth':'spatial_regression',
    'sca_tempreg_maps_stack': 'sca_roi',
    'sca_tempreg_maps_files': 'sca_roi',
    'sca_tempreg_maps_files_smooth': 'sca_roi',
    'sca_tempreg_maps_zstat_stack': 'sca_roi',
    'sca_tempreg_maps_zstat_files': 'sca_roi',
    'sca_tempreg_maps_stack_smooth': 'sca_roi',
    'sca_tempreg_maps_zstat_stack_smooth': 'sca_roi',
    'sca_tempreg_maps_zstat_files_smooth': 'sca_roi',
}


def get_zscore(input_name, wf_name = 'z_score'):
    
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
    
    wflow = pe.Workflow(name = wf_name)
    
    inputNode = pe.Node(util.IdentityInterface(fields=['input_file',
                                                       'mask_file']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['z_score_img']),
                          name='outputspec')

    mean = pe.Node(interface=fsl.ImageStats(),
                   name='mean')
    mean.inputs.op_string = '-k %s -m'    
    wflow.connect(inputNode, 'input_file',
                  mean, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  mean, 'mask_file')


    standard_deviation = pe.Node(interface=fsl.ImageStats(),
                                 name='standard_deviation')
    standard_deviation.inputs.op_string = '-k %s -s'
    wflow.connect(inputNode, 'input_file',
                  standard_deviation, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  standard_deviation, 'mask_file')
    
    
    op_string = pe.Node(util.Function(input_names=['mean','std_dev'],
                                      output_names=['op_string'],
                                      function=get_operand_string),
                        name='op_string')
    wflow.connect(mean, 'out_stat',
                  op_string, 'mean')
    wflow.connect(standard_deviation, 'out_stat',
                  op_string, 'std_dev')
    
    
    z_score = pe.Node(interface=fsl.MultiImageMaths(),
                        name='z_score')

    z_score.inputs.out_file = input_name + '_zstd.nii.gz'

    wflow.connect(op_string, 'op_string',
                  z_score, 'op_string')
    wflow.connect(inputNode, 'input_file',
                  z_score, 'in_file')
    wflow.connect(inputNode, 'mask_file',
                  z_score, 'operand_files')
  
    wflow.connect(z_score, 'out_file', outputNode, 'z_score_img')
    
    return wflow



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

    tsfile = open(timeseries_file, "rb")

    roi_list = tsfile.readlines()[0].strip("\r\n").replace("#","").split("\t")

    if prefix != None:

        temp_rois = []

        for roi in roi_list:
            roi = prefix + "_" + roi
            temp_rois.append(roi)

        roi_list = temp_rois


    return roi_list



def get_fisher_zscore(input_name, map_node, wf_name = 'fisher_z_score'):

    """
    Runs the compute_fisher_z_score function as part of a one-node workflow.
    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl
    
    wflow = pe.Workflow(name = wf_name)
    
    inputNode = pe.Node(util.IdentityInterface(fields=['correlation_file',
                                                       'timeseries_one_d']),
                        name='inputspec')

    outputNode = pe.Node(util.IdentityInterface(fields=['fisher_z_score_img']),
                          name='outputspec')


    if map_node == 0:

        fisher_z_score = pe.Node(util.Function(input_names=['correlation_file', 'timeseries_one_d', 'input_name'],
                                   output_names=['out_file'],
                     function=compute_fisher_z_score), name='fisher_z_score')

    else:

        # node to separate out 

        fisher_z_score = pe.MapNode(util.Function(input_names=['correlation_file', 'timeseries_one_d', 'input_name'],
                                   output_names=['out_file'],
                     function=compute_fisher_z_score), name='fisher_z_score',
                     iterfield=['correlation_file'])


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

    for timeseries_file_string in timeseries_one_d:
        if ".1D" in timeseries_file_string:
            timeseries_file = timeseries_file_string

    roi_numbers = []
    if '#' in open(timeseries_file, 'r').readline().rstrip('\r\n'):
        roi_numbers = open(timeseries_file, 'r').readline().rstrip('\r\n').replace('#', '').split('\t')


    # get the specific roi number
    filename = correlation_file.split("/")[-1]
    filename = filename.replace(".nii.gz","")


    corr_img = nb.load(correlation_file)
    corr_data = corr_img.get_data()

    hdr = corr_img.get_header()

    # calculate the Fisher r-to-z transformation
    corr_data = np.log((1 + corr_data) / (1 - corr_data)) / 2.0

    dims = corr_data.shape

    out_file = []

    # dims = tuple of dimensions of correlation NIFTI file
    # roi_numbers = list of label numbers for each ROI; length of this will be
    #               how many ROIs you have

    # I think the point of this check is to see if there are multiple volumes
    # in the correlation file (i.e. is a stack), or is a file with ROIs, and if
    # so, to deal with it appropriately
    if len(dims) == 5 or len(roi_numbers) > 0:

        if len(dims) == 5:

            x, y, z, one, roi_number = dims

            corr_data = np.reshape(corr_data, (x * y * z, roi_number), order='F')


        #for i in range(0, len(roi_numbers)):

        sub_data = corr_data
        #if len(dims) == 5:
        #    sub_data = np.reshape(corr_data[:, i], (x, y, z), order='F')

        sub_img = nb.Nifti1Image(sub_data, header=corr_img.get_header(), affine=corr_img.get_affine())

        sub_z_score_file = os.path.join(os.getcwd(), (filename + '_fisher_zstd.nii.gz'))

        sub_img.to_filename(sub_z_score_file)

        out_file.append(sub_z_score_file)


    # if the correlation file is a single volume image
    else:

        z_score_img = nb.Nifti1Image(corr_data, header=hdr, affine=corr_img.get_affine())

        z_score_file = os.path.join(os.getcwd(), input_name + '_fisher_zstd.nii.gz')

        z_score_img.to_filename(z_score_file)

        out_file.append(z_score_file)


    return out_file





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


    for timeseries in list_timeseries:

        if '1D' in timeseries:


            return timeseries

        else:

            print "Error : ROI/Voxel TimeSeries 1D file not found"

            return None


def extract_txt(list_timeseries):
    """
    Method to extract txt file containing 
    roi timeseries required for dual regression
    """

    out_file = None
    for timeseries in list_timeseries:
        if timeseries.endswith('.txt'):
            out_file = timeseries

    if not out_file:
        raise Exception("Unable to retrieve roi timeseries txt"\
                          " file required for dual regression")

    return out_file


def set_gauss(fwhm):

    fwhm = float(fwhm)

    sigma = float(fwhm / 2.3548)

    op = "-kernel gauss %f -fmean -mas " % (sigma) + "%s"
    op_string = op

    return op_string



def get_path_score(path, entry):

    import os

    parent_dir = os.path.dirname(path)

    dirs = parent_dir.split('/')
    dirs.remove('')


    score = 0

    for element in entry:

        if element in dirs:

            score += 1


    return score


def get_strategies_for_path(path, strategies):

    from CPAC.utils.utils import get_path_score
    max_score = 0
    score_dict = {}
    for strategy in strategies:

        score = get_path_score(path, strategy)

        if score >= max_score:
            max_score = score

            if str(max_score) in score_dict:

                lst = score_dict[str(max_score)]
                lst.append(strategy)

            else:

                score_dict[str(max_score)] = [strategy]

    return score_dict[str(max_score)]



def get_workflow(remainder_path):

    # this iterates over the hard-coded list at the top of this file
    # (utils.py) and matches workflow output paths provided to the function to
    # more user-friendly labels. for example, a path to the output of
    # 'functional_preprocessed_mask' is matched with the label 'func'

    global files_folders_wf
    lst = remainder_path.split('/')

    lst = [x for x in lst if not ('' == x) ]

    return lst[0], files_folders_wf[lst[0]], remainder_path.split(lst[0])[1]



def get_session(remainder_path):

    session = 'scan_'


    lst = remainder_path.split('/')

    lst = [x for x in lst if not ('' == x) ]


    for element in lst:

        if 'scan_' in element:

            session += element.split('scan_')[1] + '_'


    if session.endswith('_'):

        session = session.rstrip('_')

    return session



def get_hplpfwhmseed_(parameter, remainder_path):

    # this function extracts the filtering and smoothing parameters info from
    # a path leading to an output of individual level analysis, and then
    # returns these parameter values so that descriptive sym-link directories
    # can be generated 

    partial_parameter_value = remainder_path.split(parameter)[1]

    value = partial_parameter_value.split('/')[0]

    return parameter.lstrip('/_') + value



def create_seeds_(seedOutputLocation, seed_specification_file, FSLDIR):

    import commands
    import os
    import re

    seed_specifications = [line.rstrip('\r\n') for line in open(seed_specification_file, 'r').readlines() if (not line.startswith('#') and not (line == '\n')) ]

    seed_resolutions = {}

    for specification in seed_specifications:

        seed_label, x, y, z, radius, resolution = re.split(r'[\t| |,]+', specification)

        if resolution not in seed_resolutions.keys():
            seed_resolutions[resolution] = []
        seed_resolutions[resolution].append((seed_label, x, y, z, radius, resolution))

    return_roi_files = []
    for resolution_key in seed_resolutions:


        index = 0
        seed_files = []
        for roi_set in seed_resolutions[resolution_key]:


            seed_label, x, y, z, radius, resolution = roi_set
            if not os.path.exists(seedOutputLocation):
                os.makedirs(seedOutputLocation)

            print 'checking if file exists ', '%s/data/standard/MNI152_T1_%s_brain.nii.gz' % (FSLDIR, resolution)
            assert(os.path.exists('%s/data/standard/MNI152_T1_%s_brain.nii.gz' % (FSLDIR, resolution)))
            cmd = "echo %s %s %s | 3dUndump -prefix %s.nii.gz -master %s/data/standard/MNI152_T1_%s_brain.nii.gz \
-srad %s -orient LPI -xyz -" % (x, y, z, os.path.join(seedOutputLocation, str(index) + '_' + seed_label + '_' + resolution), FSLDIR, resolution, radius)

            print cmd
            try:
                commands.getoutput(cmd)
                seed_files.append((os.path.join(seedOutputLocation, '%s.nii.gz' % (str(index) + '_' + seed_label + '_' + resolution)), seed_label))
                print seed_files
            except:
                raise

            index += 1

        print index, ' ', seed_files
        seed_str = ' '
        intensities = ''
        for idx in range(0, index):

            seed, intensity = seed_files[idx]

            intensities += intensity + '_'

            cmd = "3dcalc -a %s -expr 'a*%s' -prefix %s" % (seed, intensity, os.path.join(seedOutputLocation, 'ic_' + os.path.basename(seed)))
            print cmd
            try:
                commands.getoutput(cmd)

                seed_str += "%s " % os.path.join(seedOutputLocation, 'ic_' + os.path.basename(seed))
            except:
                raise


        cmd = '3dMean  -prefix %s.nii.gz -sum %s' % (os.path.join(seedOutputLocation, 'rois_' + resolution), seed_str)
        print cmd
        try:
            commands.getoutput(cmd)
        except:
            raise

        try:
            cmd = 'rm -f %s' % seed_str
            print cmd
            commands.getoutput(cmd)
            for seed, intensity in seed_files:
                try:
                    cmd = 'rm -f  ' + seed
                    print cmd
                    commands.getoutput(cmd)
                except:
                    raise
        except:
            raise
        return_roi_files.append(os.path.join(seedOutputLocation, 'rois_' + resolution + '.nii.gz'))

    print return_roi_files
    return return_roi_files


def create_paths_and_links(pipeline_id, relevant_strategies, path, subject_id, create_sym_links):

    import os
    import commands
    from CPAC.utils.utils import get_workflow, get_session, \
                     get_hplpfwhmseed_

    # path (one of the inputs of this function) is a path to a file output by
    # individual-level analysis, and this function runs once per output file

    # relevant_strategies is a list of lists, where each list contains all of
    # the nuisance correction selections per strategy
    for strategy in relevant_strategies:

        base_path, remainder_path = path.split(subject_id, 1)

        sym_path = path.split(pipeline_id)[0]

        # create 'file_path', the path to a subject's output folder
        # example: {path}/{output folder}/{pipeline id}/{subject id}
        #
        # this is used later to generate the QC and path_files directory paths
        file_path = os.path.join(sym_path, pipeline_id)
        file_path = os.path.join(file_path, subject_id)

        # create the sym-link directory paths
        sym_path = os.path.join(sym_path, 'sym_links')
        sym_path = os.path.join(sym_path, pipeline_id)


        if create_sym_links == True:

            try:
                os.makedirs(sym_path)
            except:
                # don't raise an exception here because multiple runs of the
                # same os.makedirs are expected
                pass
        

        strategy_identifier = None


        try:

            short_names = {'_threshold':'SCRUB_', '_csf_threshold':'CSF_',
                    '_gm_threshold':'GM_',
                    '_compcor_':'compcor',
                    '_target_angle_deg':'MEDIANangle_', '_wm_threshold':'WM_'}

            strategy_identifier = ''

            # iterate over each attribute of the nuisance correction strategy,
            # then process these into the strategy_identifier string which is
            # is used to name the sym-link folders under the pipeline ID
            # folder and above the subject ID folders in the sym-links
            # directory, describing the user's correction selections
            for el in strategy:

                key, value = el.rsplit('_', 1)

                if '_compcor_'in key:

                    if 'compcor0' in value:
                        strategy_identifier += (value + '_')
                        continue
                    else:
                        val1 = key.split('_selector')[0]
                        strategy_identifier += val1 + '_' + value + '_'
                        continue

                # the pipeline name is included as one of the instances of
                # 'el' in this iteration, but we don't want this as part of
                # the strategy identifier, so exclude 'pipeline_{name}'
                if not 'pipeline' in key:
                    strategy_identifier += short_names[key] + value + '_'

            strategy_identifier = strategy_identifier.rsplit('_', 1)[0]

        except:
            print str(strategy), " not in labels_dict"
            raise


        # this removes unused corrections from the strategy_identifier string,
        # keeping in mind that corrections with a 0 appended to the end of the
        # name denotes they were not included
        strategy_identifier = strategy_identifier.replace('pc10.', '')
        strategy_identifier = strategy_identifier.replace('linear0.', '')
        strategy_identifier = strategy_identifier.replace('wm0.', '')
        strategy_identifier = strategy_identifier.replace('global0.', '')
        strategy_identifier = strategy_identifier.replace('motion0.', '')
        strategy_identifier = strategy_identifier.replace('quadratic0.', '')
        strategy_identifier = strategy_identifier.replace('gm0.', '')
        strategy_identifier = strategy_identifier.replace('csf0_', '')
        strategy_identifier = strategy_identifier.replace('compcor0.', '')

        # start making basic sym link directories
        sym_path = os.path.join(sym_path, strategy_identifier)
        new_sub_path = os.path.join(sym_path, subject_id)

        # 'get_workflow' iterates over the hard-coded list at the top of this
        # file (utils.py) and matches workflow output paths provided to the
        # function to more user-friendly labels. for example, a path to the
        # output of 'functional_preprocessed_mask' is matched with the label
        # 'func'
        file_name, wf, remainder_path = get_workflow(remainder_path)
        session = get_session(remainder_path)

        new_session_path = os.path.join(new_sub_path, session)
        new_wf_path = os.path.join(new_session_path, wf)
        new_path = new_wf_path

        # new_path is now the full path to a workflow output's symlink folder
        # for one subject and one nuisance correction strategy
        # example: {path to output folder}/sym_links/{pipeline id}/
        #              {correction selections}/{subject id}/{scan id}/
        #                  {workflow label, such as alff}



        # now bring into use the tier 2 iterables for recursive directory
        # structure

        # the 'get_hplpfwhmseed_' function extracts the filtering and
        # smoothing parameters info from a path leading to an output of
        # individual level analysis, and then returns these parameter values
        # so that descriptive sym-link directories can be generated

        # example: /{path to output folder}/sym_links/{pipeline id}/
        #              {correction selections}/{subject id}/{scan id}/
        #                  alff/hp_0.01/lp_0.1
        #
        # with 'hp_0.01' and 'lp_0.1' being the outputs of the
        # 'get_hplpfwhmseed_' function for the remainder_path that included
        # the alff output 

        scan_info = ''
        if '/_scan_' in path:

            scan_info = get_hplpfwhmseed_('/_scan_', path)


        if '/_mask_' in remainder_path:

            new_path = os.path.join(new_path, \
                get_hplpfwhmseed_('/_mask_', remainder_path))


        if '/_roi_' in remainder_path:

            new_path = os.path.join(new_path, \
                get_hplpfwhmseed_('/_roi_', remainder_path))


        if '/_sca_roi_' in remainder_path:

            new_path = os.path.join(new_path, \
                get_hplpfwhmseed_('/_sca_roi_', remainder_path))


        hp_str = ''
        if ('_hp_'  in remainder_path):

            hp_str = get_hplpfwhmseed_('/_hp_', remainder_path)
            new_path = os.path.join(new_path, hp_str)


        lp_str = ''
        if ('_lp_' in remainder_path):

            lp_str = get_hplpfwhmseed_('/_lp_', remainder_path)
            new_path = os.path.join(new_path, lp_str)


        bp_freq = ''
        if ('_bandpass_freqs_' in remainder_path):

            bp_freq = get_hplpfwhmseed_('/_bandpass_freqs_', remainder_path)
            new_path = os.path.join(new_path, bp_freq)


        spatial_map = ''
        if('_spatial_map_' in remainder_path):

            spatial_map = get_hplpfwhmseed_('/_spatial_map_', remainder_path)
            new_path = os.path.join(new_path, spatial_map)


        fwhm_str = ''
        if ('_fwhm_' in remainder_path):

            fwhm_str = get_hplpfwhmseed_('/_fwhm_', remainder_path)
            new_path = os.path.join(new_path, fwhm_str)


        if create_sym_links == True:

            # create the final symlink path for the output, if not created
            # already
            try:
                os.makedirs(new_path)
            except:
                # don't raise an exception here because multiple runs of
                # os.makedirs are expected
                pass


        # prepare paths and filenames for QC text files and output paths_file
        # text files that are written to the output directory

        try:

            if wf == 'qc':
                # if the output file is QC related, send it over to
                # 'qc_files_here'. these files are often the .png images for
                # the QC dashboard
                new_f_path = os.path.join(file_path, 'qc_files_here')
                os.makedirs(new_f_path)
            else:
                # if the output file is not QC related, send it over to
                # 'path_files_here'
                new_f_path = os.path.join(file_path, 'path_files_here')
                os.makedirs(new_f_path)

        except:
            # don't raise an exception here because multiple runs of
            # os.makedirs are expected
            pass


        try:

            global global_lock
            global_lock.acquire()

            f_n = None

            # use the strategy identifier strings and the filtering and
            # smoothing parameter values extracted above to name the qc and
            # output paths_file text files descriptively

            if wf == 'qc':
                f_n = os.path.join(new_f_path, 'qc_%s.txt') % (scan_info + '_' + strategy_identifier + '_' + bp_freq + '_' + hp_str + '_' + lp_str + '_' + fwhm_str)
            else:
                f_n = os.path.join(new_f_path, 'paths_file_%s.txt') % (scan_info + '_' + strategy_identifier + '_' + bp_freq + '_' + hp_str + '_' + lp_str + '_' + fwhm_str)

            f = open(f_n, 'a')
            print >>f, path

            global_lock.release()

        except:

            print 'trouble acquiring locks or opening file skipping :', os.path.join(new_f_path, 'paths_file_%s.txt') % new_path.replace('/', '_')
            raise



        if create_sym_links == True:

            # create the actual sym-links now

            # fname is the filename of the current individual level output
            # file
            fname = os.path.basename(path)

            # ext is the extension of the current individual level output file
            ext = fname.split('.', 1)[1]
            ext = '.' + (ext)


            # special case for ROI , need the ROI number
            if '_ROI_' in fname and 'sca_' in path:

                # extracts the ROI number from the output file's path in the
                # output folder and appends it to the symlink file's name
                import re
                roi_number = re.findall(r'\d+', fname)[0]
                file_name += '_' + roi_number


            dont_change_fname = ['vertices_timeseries',
            'centrality_outputs_smoothed',
            'centrality_outputs_zscore',
            'centrality_outputs',
            'centrality_graphs'
            'voxel_timeseries',
            'roi_timeseries',
            'seg_probability_maps',
            'seg_mixeltype',
            'seg_partial_volume_map',
            'seg_partial_volume_files',
            'dr_tempreg_maps_zstat_files',
            'dr_tempreg_maps_zstat_files_smooth',
            'sca_tempreg_maps_zstat_files',
            'sca_tempreg_maps_zstat_files_smooth']


            # split up the files between QC or not, and also between the ones
            # where the filename has changed for the symlink or not

            # if wf == 'qc', those are files used in the QC pages, such as
            # .png files and other images

            if (file_name in dont_change_fname) or (wf == 'qc'):

                cmd = 'ln -s %s %s' % (path, os.path.join(new_path, fname))
                print cmd
                commands.getoutput(cmd)

            else:

                cmd = 'ln -s %s %s' % (path, os.path.join(new_path, file_name + ext))

                try:
                    f1 = open(os.path.join(new_path, file_name + ext))
                except:
                    print cmd
                    commands.getoutput(cmd)



def prepare_gp_links(in_file, resource):

    import os
    import re

    in_file = os.path.abspath(in_file)
    def get_param_val_(splittext, in_file):

        param = in_file.split(splittext, 1)[1]

        return param.split('/')[0]


    in_file = os.path.abspath(in_file)

    sink_dir, pipeline_id = in_file.split('/pipeline_')

    pipeline_id = ''.join(['pipeline_', pipeline_id.split('/')[0]])

    # sym link directory
    sink_dir = os.path.join(sink_dir, 'sym_links')

    # pipeline directory
    sink_dir = os.path.join(sink_dir, pipeline_id)

    # group directory
    sink_dir = os.path.join(sink_dir, 'group_analysis_results')

    strategy_identifier = ''
    if '_selector_' in in_file:

        prepath, strategy_identifier = in_file.split('_selector_')

        strategy_identifier = strategy_identifier.split('/', 1)[0]

        strategy_identifier = strategy_identifier.replace('pc10.', '')
        strategy_identifier = strategy_identifier.replace('linear0.', '')
        strategy_identifier = strategy_identifier.replace('wm0.', '')
        strategy_identifier = strategy_identifier.replace('global0.', '')
        strategy_identifier = strategy_identifier.replace('motion0.', '')
        strategy_identifier = strategy_identifier.replace('quadratic0.', '')
        strategy_identifier = strategy_identifier.replace('gm0.', '')
        strategy_identifier = strategy_identifier.replace('csf0', '')
        strategy_identifier = strategy_identifier.replace('compcor0.', '')

        ncomponents = ''
        if 'compcor1' in strategy_identifier:

            ncomponents = prepath.split('_ncomponents_')[1]


        strategy_identifier = strategy_identifier.replace('pc11.', 'pc1.')
        strategy_identifier = strategy_identifier.replace('linear1.', 'lin.')
        strategy_identifier = strategy_identifier.replace('wm1.', 'wm.')
        strategy_identifier = strategy_identifier.replace('global1.', 'gl.')
        strategy_identifier = strategy_identifier.replace('motion1.', 'motion.')
        strategy_identifier = strategy_identifier.replace('quadratic1.', 'quad.')
        strategy_identifier = strategy_identifier.replace('gm0.', 'gm.')
        strategy_identifier = strategy_identifier.replace('csf1', 'csf')
        strategy_identifier = strategy_identifier.replace('compcor1.', 'compcor_nc_%s.' % ncomponents)

    scan_info = ''

    if '/_scan_' in in_file:

        scan_info = get_param_val_('/_scan_', in_file)

    scrub = ''
    if '/_threshold_' in in_file:
        scrub = get_param_val_('/_threshold_', in_file)
        strategy_identifier = 'SCRUB_%s' % scrub + '_' + strategy_identifier

    if '/_wm_threshold' in in_file:
        strategy_identifier += '_wmT_' + get_param_val_('/_wm_threshold_', in_file)

    if '/_csf_threshold' in in_file:
        strategy_identifier += '_csfT_' + get_param_val_('/_csf_threshold_', in_file)

    if '/_gm_threshold' in in_file:
        strategy_identifier += '_gmT_' + get_param_val_('/_gm_threshold_', in_file)


    if not scan_info == '':
        strategy_identifier = scan_info + '_' + strategy_identifier


    sink_dir = os.path.join(sink_dir, strategy_identifier)


    second_tier = ''
    if '/_bandpass_freqs_' in in_file:
        second_tier = 'bp_freqs_' + get_param_val_('/_bandpass_freqs_', in_file)

    if '/_hp_' in in_file:
        second_tier += '_hp_' + get_param_val_('/_hp_', in_file)

    if '/_lp_' in in_file:
        second_tier += '_lp_' + get_param_val_('/_lp_', in_file)

    if '/_fwhm_' in in_file:
        second_tier += '_fwhm_' + get_param_val_('/_fwhm_', in_file)


    sink_dir = os.path.join(sink_dir, second_tier)
    gp_dir = str(sink_dir)

    if '/_grp_model_' in in_file:

        model_info = get_param_val_('/_grp_model_', in_file)
        sink_dir = os.path.join(sink_dir, model_info)



    third_tier = ''

    fourth_tier = ''

    if 'sca_roi' in resource and '/_roi_' in in_file:

        third_tier = resource + '_' + get_param_val_('/_roi_', in_file)

        roi_number = ''.join(['ROI_', get_param_val_('/ROI_number_', in_file)])
        third_tier = third_tier + '/' + roi_number

    elif ('dr_tempreg_maps_zstat_files' in resource and '/temp_reg_map_z_' in in_file):

        third_tier = resource + '_' + get_param_val_('/_spatial_map_', in_file)
        third_tier = third_tier + '/' + get_param_val_('/temp_reg_map_z_', in_file)

    elif ('sca_tempreg_maps_zstat_files' in resource and '/sca_tempreg_z_maps_roi_' in in_file):
        third_tier = resource + '_' + get_param_val_('/_roi_', in_file)
        roi_number = ''.join(['ROI_', get_param_val_('/sca_tempreg_z_maps_roi_', in_file)])
        third_tier = third_tier + '/' + roi_number


    elif ('sca_seed' in resource or 'centrality_outputs' in resource)  and '/_mask_' in in_file:

        third_tier = resource + '_' + get_param_val_('/_mask_', in_file)

        if 'centrality' in resource:

            if 'degree_centrality_binarize' in in_file:
                centrality_type = 'degree_centrality_binarize'
            elif 'degree_centrality_weighted' in in_file:
                centrality_type = 'degree_centrality_weighted'
            elif 'eigenvector_centrality_binarize' in in_file:
                centrality_type = 'eigenvector_centrality_binarize'
            elif 'eigenvector_centrality_weighted' in in_file:
                centrality_type = 'eigenvector_centrality_weighted'

            else:
                raise ValueError('centrality type not in degree_centrality_binarize, \
                        degree_centrality_weighted eigenvector_centrality_binarize, \
                        eigenvector_centrality_weighted')

            third_tier = third_tier + '/' + centrality_type


    else:

        third_tier = resource


    sink_dir = os.path.join(sink_dir, third_tier)

    grp = None
    if 'model_files' in in_file:
        grp = re.search(r'model_files(.)*', in_file)

    if 'merged' in in_file:

        grp = re.search(r'merged(.)*', in_file)

    if 'rendered' in in_file:

        grp = re.search(r'rendered(.)*', in_file)

    if 'stats' in in_file:

        grp = re.search(r'stats(.)*', in_file)

    tier_4 = os.path.dirname(re.sub(r'_grp_model_(.)+/', '', grp.group()))


    sink_dir = os.path.join(sink_dir, tier_4)

    if 'merged' in in_file:
        sink_dir = os.path.join(sink_dir, 'merged.nii.gz')
    else:
        sink_dir = os.path.join(sink_dir, os.path.basename(in_file))

    try:
        os.makedirs(os.path.dirname(sink_dir))

    except Exception, e:

        print '.'


    import commands

    cmd = 'ln -s %s %s' % (in_file, sink_dir)
    print cmd
    commands.getoutput(cmd)



def clean_strategy(strategies, helper):

# ##
# ## If segmentation or scrubbing or nuisance or median is turned off
# ## in the pipeline then remove them from the strategy tag list

    new_strat = []


    for strat in strategies:

        tmpstrat = []
        for el in strat:

            key = el.rsplit('_', 1)[0]

            print key, ': ', el
            if not ('compcor' in key):

                if 'pipeline' in key:

                    tmpstrat.append(el)
                    continue

                try:
                    todos = helper[key]

                    if not (todos == 0):

                        tmpstrat.append(el)

                except:

                    print 'key ', key, 'from ', el, ' not in ', helper
                    raise

            else:

                if not helper['nuisance'] == 0:

                    tmpstrat.append(el)

        new_strat.append(tmpstrat)


    return new_strat



def process_outputs(in_file, strategies, subject_id, pipeline_id, helper, create_sym_links):

    from CPAC.utils.utils import get_strategies_for_path, create_paths_and_links, clean_strategy
        
    for path in in_file:

        for strategy in strategies:

            strategy.append(pipeline_id)

        relevant_strategies = get_strategies_for_path(path, strategies)
        
        cleaned_strategies = clean_strategy(relevant_strategies, helper)
        
        create_paths_and_links(pipeline_id, cleaned_strategies, path, subject_id, create_sym_links)



def modify_model(input_sublist, output_sublist, mat_file, grp_file):
    """
    Method to modify .grp and .mat fsl group analysis model files
     
    Parameters
    ----------
    input_sublist : string (list)
         Path to group analysis input subject list containing all the subjects 
         for which CPAC is run
    output_sublist : string (list)
        Path to subject list for that were successfully run for a particular 
        derivative
    mat_file : string (fsl mat file)
        path to mat file containing  matrix for design
    grp_file : string (fsl grp file)
         path to file containing matrix specifying 
         the groups the covariance is split into 
         
    Returns
    -------
    new_grp_file : string (grp file)
        modified covariance group matrix model file
    new_mat_file : string (mat file)
        modified design matrix file
    new_sub_file : string (txt file)
        new model subject list 
    """

    import os
    def read_model_file(file):
        """
        Method to read model file and
        create a map out of it
        """
        dict1 = {}
        for line in open(file, 'r'):
            if line.strip() != '':
                if 'NumWaves' in line:
                    dict1['/NumWaves'] = line.split()[1]
                elif 'NumPoints' in line:
                    dict1['/NumPoints'] = line.split()[1]
                elif 'PPheights' in line:
                    dict1['/PPHeights'] = line.split()[1:]
                elif 'Matrix' in line:
                    dict1['/Matrix'] = []
                else:
                    dict1.get('/Matrix').append(line)
        return dict1

    def write_model_file(model_map, model_file, remove_index):

        out_file, ext = os.path.splitext(os.path.basename(model_file))
        out_file = out_file + '_new' + ext
        out_file = os.path.join(os.getcwd(), out_file)

        # create an index of all subjects for a derivative for which
        # CPAC did not run successfully

        f = open(out_file, 'wb')
        print >> f, '/NumWaves\t' + model_map['/NumWaves']

        num_points = int(model_map['/NumPoints']) - len(remove_index)
        print >> f, '/NumPoints\t' + str(num_points)
        if ext == ".mat":
            f.write('/PPHeights\t\t')
            for val in model_map['/PPHeights']:
                f.write(val + '\t')
            f.write('\n')
        print >> f, ''
        print >> f, '/Matrix'
        count = 0
        for values in model_map['/Matrix']:
            # remove the row form matrix for all unsuccessful subjects
            if count not in remove_index:
                f.write(values)
            count += 1

        f.close()

        return out_file

    # get new subject list
    new_sub_file = os.path.join(os.getcwd(), 'model_subject_list.txt')
    f = open(new_sub_file, 'wb')
    remove_index = []
    for subject in input_sublist:
         if subject not in output_sublist:
              print "Derivative output not found for subject %s " % (subject)
              remove_index.append(input_sublist.index(subject))
         else:
              print >> f, subject

    f.close()

    print "removing subject at the indices", remove_index
    print "modifying the mat and grp files"

    model_map = read_model_file(mat_file)
    new_mat_file = write_model_file(model_map, mat_file, remove_index)

    model_map = read_model_file(grp_file)
    new_grp_file = write_model_file(model_map, grp_file, remove_index)

    return new_grp_file, new_mat_file, new_sub_file



def select_model_files(model, ftest, model_name):

    """
    Method to select model files
    """

    import os
    import glob

    files = glob.glob(os.path.join(model, '*'))

    if len(files) == 0:
        raise Exception("No files found inside directory %s" % model)

    fts_file = ''

    for filename in files:
        if (model_name + '.mat') in filename:
            mat_file = filename
        elif (model_name + '.grp') in filename:
            grp_file = filename
        elif ((model_name + '.fts') in filename) and ftest:
            fts_file = filename
        elif (model_name + '.con') in filename:
            con_file = filename

    if ftest == True and fts_file == '':

        errmsg = "\n[!] CPAC says: You have f-tests included in your group " \
                 "analysis model '%s', but no .fts files were found in the " \
                 "output folder specified for group analysis: %s.\n\nThe " \
                 ".fts file is automatically generated by CPAC, and if you " \
                 "are seeing this error, it is because something went wrong " \
                 "with the generation of this file, or it has been moved." \
                 "\n\n" % (model_name, model)

        raise Exception(errmsg)


    return fts_file, con_file, grp_file, mat_file



def get_scan_params(subject, scan, subject_map, start_indx, stop_indx, tr, tpattern):

    """
    Method to extract slice timing correction parameters
    and scan parameters.
    
    Parameters
    ----------
    subject: a string
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
    import warnings

    def check(val, throw_exception):

        if isinstance(subject_map['scan_parameters'][val], dict):
            ret_val = subject_map['scan_parameters'][val][scan]
        else:
            ret_val = subject_map['scan_parameters'][val]

        if ret_val == 'None':
            if throw_exception:
                raise Exception("None Parameter Value for %s for subject %s" % (val, subject))
            else:
                ret_val = None

        if ret_val == '' and throw_exception:
            raise Exception("Missing Value for %s for subject %s" % (val, subject))

        return ret_val

    check2 = lambda val : val if val == None or val == '' else int(val)


    # initialize vars to empty
    TR=''
    pattern=''
    ref_slice=''
    first_tr=''
    last_tr=''

    if 'scan_parameters' in subject_map.keys():
        # get details from the configuration
        TR = float(check('tr', False))
        pattern = str(check('acquisition', False))
        ref_slice = int(check('reference', False))
        first_tr = check2(check('first_tr', False))
        last_tr = check2(check('last_tr', False))


    # if values are still empty, override with GUI config
    if TR == '':
        if tr:
            TR = float(tr)
        else:
            TR = None

    if first_tr == '':
        first_tr = start_indx

    if last_tr == '':
        last_tr = stop_indx

    unit = 's'

    # if the user has mandated that we the timining informaiton in the header,
    # that takes precedence
    if "Use NIFTI Header" in tpattern:
        pattern = ''
    else:
    # otherwise he slice acquisition pattern in the subject file takes precedence, but if it 
    # isn't set we use the value in the configuration file
        if pattern == '':
            pattern = tpattern

    # pattern can be one of a few keywords, a filename, or blank which indicates that the 
    # images header information should be used
    if pattern and pattern not in ['alt+z', 'altplus', 'alt+z2', 'alt-z', 'altminus',
                   'alt-z2', 'seq+z', 'seqplus', 'seq-z', 'seqminus']:
        if not os.path.exists(pattern):
            raise Exception ("Invalid Pattern file path %s , Please provide the correct path" % pattern)
        else:
            lines = open(pattern, 'r').readlines()
            if len(lines) < 2:
                raise Exception('Invalid slice timing file format. The file should contain '\
                                'only one value per row. Use new line char as delimiter')
            pattern = '@' + pattern

            slice_timings = [float(l.rstrip('\r\n')) for l in lines]
            slice_timings.sort()
            max_slice_offset = slice_timings[-1]
            # checking if the unit of TR and slice timing match or not
            # if slice timing in ms convert TR to ms as well
            if  TR and max_slice_offset > TR:
                warnings.warn("TR is in seconds and slice timings are in milliseconds."\
                              "Converting TR into milliseconds")
                TR = TR * 1000
                print "New TR value %.2f ms" % TR
                unit = 'ms'

    else:
        # check to see, if TR is in milliseconds, convert it into seconds
        if TR and TR > 10:
            warnings.warn('TR is in milliseconds, Converting it into seconds')
            TR = TR / 1000.0
            print "New TR value %.2f s" % TR
            unit = 's'

    print "scan_parameters -> ", subject, scan, str(TR) + unit, pattern, ref_slice, first_tr, last_tr

    return str(TR) + unit, pattern, ref_slice, first_tr, last_tr


def get_tr (tr):
    """
    Method to return TR in seconds
    """
    import re
    if 'None' in tr:
        tr = None
    if tr:
        tr = re.search("\d+.\d+", str(tr)).group(0)
        tr = float(tr)
        if tr > 10:
            tr = tr / 1000.0
    return tr


def check_tr(tr, in_file):

    # imageData would have to be the image data from the funcFlow workflow, funcFlow outputspec.subject
    import nibabel as nib
    img = nib.load(in_file)
    
    # get header from image data, then extract TR information, TR is fourth item in list returned by get_zooms()
    imageHeader = img.get_header()
    imageZooms = imageHeader.get_zooms()
    header_tr = imageZooms[3]
    
                
    # If the TR information from header_tr (funcFlow) and convert_tr node (TR from config file)
    # do not match, prepare to update the TR information from either convert_tr or header_tr using
    # afni 3drefit, then append to func_to_mni
    if header_tr != tr:
        
        if tr != None and tr != "":
            TR = tr
        else:
            TR = header_tr
            
        import warnings
        warnings.warn('Warning: The TR information does not match between the config and subject list files.')
    
    return TR



def write_to_log(workflow, log_dir, index, inputs, scan_id ):
    """
    Method to write into log file the status of the workflow run.
    """
    
    import os
    import CPAC
    from nipype import logging
    iflogger = logging.getLogger('interface')
     
    version = CPAC.__version__
         
    subject_id = os.path.basename(log_dir)
    
    if scan_id == None:
        scan_id = "scan_anat"
    
    strategy = ""
    
    import time
    import datetime
    ts = time.time()
    
    stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    try:
        if workflow!= 'DONE':
            wf_path = os.path.dirname((os.getcwd()).split(workflow)[1]).strip("/")
            
            if wf_path and wf_path != "":
                if '/' in wf_path:
                    scan_id, strategy = wf_path.split('/',1)
                    scan_id = scan_id.strip('_')
                    strategy = strategy.replace("/","")
                else:
                    scan_id = wf_path.strip('_')
        
            file_path = os.path.join(log_dir, scan_id, workflow)
            
            try:
                os.makedirs(file_path)
            except Exception:
                iflogger.info("filepath already exist, filepath- %s, curr_dir - %s"%(file_path, os.getcwd()))

        else:
            file_path = os.path.join(log_dir, scan_id)
    except Exception:
        print "ERROR in write log"
        raise

    try:
        os.makedirs(file_path)
    except Exception:
        iflogger.info("filepath already exist, filepath- %s, curr_dir - %s"%(file_path, os.getcwd()))
        
    out_file = os.path.join(file_path, 'log_%s.yaml'%strategy)
    
    f = open(out_file, 'w')
    
    
    print >>f, "version : %s"%(str(version))
    print >>f, "timestamp: %s"%(str(stamp))
    print >>f, "pipeline_index: %d"%(index) 
    print >>f, "subject_id: %s"%(subject_id)
    print >>f, "scan_id: %s"%(scan_id)
    print >>f, "strategy: %s"%(strategy)
    print >>f, "workflow_name: %s"%(workflow)
        
        

    iflogger.info("CPAC custom log :")
    
    if isinstance(inputs, list):
        inputs = inputs[0]
        
    if os.path.exists(inputs):

        print >>f,  "wf_status: DONE"
  
        iflogger.info(" version - %s, timestamp -%s, subject_id -%s, scan_id - %s, strategy -%s, workflow - %s, status -%s"\
                      %(str(version), str(stamp), subject_id, scan_id,strategy,workflow,'COMPLETED') )
  
    else:
        
        iflogger.info(" version - %s, timestamp -%s, subject_id -%s, scan_id - %s, strategy -%s, workflow - %s, status -%s"\
                      %(str(version), str(stamp), subject_id, scan_id,strategy,workflow,'ERROR') )
    
        print>>f, "wf_status: ERROR"
    
    f.close()
    
    #os.system("/home2/haipan/tmp/C-PAC/scripts/log_py2js.py %s %s"%(out_file, log_dir))   ###
    
    return out_file


def create_log( wf_name = "log", 
                scan_id = None):
    
    """
    Workflow to create log 
    
    """
    
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    
    wf = pe.Workflow(name=wf_name)
    
    inputNode = pe.Node(util.IdentityInterface(fields=['workflow',
                                                       'log_dir',
                                                       'index',
                                                       'inputs'
                                                    ]),
                        name='inputspec')


    outputNode = pe.Node(util.IdentityInterface(fields=['out_file']),
                        name='outputspec')

    write_log = pe.Node(util.Function(input_names=[ 'workflow',
                                                    'log_dir',
                                                    'index',
                                                    'inputs',
                                                    'scan_id'],
                                               output_names=['out_file'],
                                               function=write_to_log),
                                 name='write_log')


    wf.connect(inputNode, 'workflow',
               write_log, 'workflow')
    wf.connect(inputNode, 'log_dir',
               write_log, 'log_dir')
    wf.connect(inputNode, 'index',
               write_log, 'index')
    wf.connect(inputNode, 'inputs',
               write_log, 'inputs')
    
    write_log.inputs.scan_id = scan_id

    wf.connect(write_log, 'out_file',
               outputNode, 'out_file')
    
    return wf
    

def create_log_template(pip_ids, wf_list, scan_ids, subject_id, log_dir):
   
    import datetime, os
    from os import path as op
    from jinja2 import Template
    import pkg_resources as p
    import CPAC
    import itertools    
    
    now = datetime.datetime.now()
    
    chain = itertools.chain(*wf_list)
    wf_keys = list(chain)
    wf_keys = list(set(wf_keys))
    
    
    tvars = {}
    tvars['subject_id']  = subject_id
    tvars['scans']       = scan_ids
    tvars['pipelines']   = pip_ids
    tvars['wf_list']     = "%s" % wf_list
    tvars['wf_keys']     = "%s" % wf_keys
    tvars['pipeline_indices'] = range(len(tvars['pipelines']))
    tvars['resources'] = os.path.join(CPAC.__path__[0], 'resources')
    tvars['gui_resources'] = os.path.join(CPAC.__path__[0], 'GUI', 'resources')
    
    
    reportdir = op.join(log_dir, "reports")
    if not op.exists(reportdir):
        os.mkdir(reportdir)
    
    for scan in scan_ids:
        jsfile   = op.join(reportdir, "%s.js" % scan)
        open(jsfile, 'w').close()
        
        tvars['cur_scan']    = scan
        tvars['logfile']     = jsfile
        tvars['timestamp']   = now.strftime("%Y-%m-%d %H:%M:%S")
        
        fname    = p.resource_filename('CPAC','resources/templates/cpac_runner.html')
        tfile    = open(fname, 'r')
        raw_text = tfile.read()
        tfile.close()
        
        template = Template(raw_text)
        text     = template.render(**tvars)
        
        htmlfile = op.join(reportdir, "%s.html" % scan)
        html     = open(htmlfile, 'w')
        html.write(text)
        html.close()
    
    # Index File
    fname       = p.resource_filename('CPAC','resources/templates/logger_subject_index.html')
    tfile       = open(fname, 'r')
    raw_text    = tfile.read()
    tfile.close()
    
    template    = Template(raw_text)
    text        = template.render(**tvars)
    
    htmlfile    = op.join(reportdir, "index.html")
    html        = open(htmlfile, 'w')
    html.write(text)
    html.close()
    
    
    return


def create_group_log_template(subject_scan_map, log_dir):
    
    import os
    from os import path as op
    from jinja2 import Template
    import pkg_resources as p
    import CPAC

    tvars = {}
    tvars['subject_ids'] = subject_scan_map.keys()
    tvars['scan_ids']    = subject_scan_map
    tvars['resources']   = op.join(CPAC.__path__[0], 'resources')
    tvars['log_dir']     = log_dir
    
    reportdir = op.join(log_dir, "reports")
    if not op.exists(reportdir):
        os.makedirs(reportdir)
    
    fname    = p.resource_filename('CPAC','resources/templates/logger_group_index.html')
    tfile    = open(fname, 'r')
    raw_text = tfile.read()
    tfile.close()
    
    template = Template(raw_text)
    text     = template.render(**tvars)
    
    htmlfile = op.join(reportdir, "index.html")
    html     = open(htmlfile, 'w')
    html.write(text)
    html.close()
    
    return



def extract_output_mean(in_file, output_name):

    '''
    function takes 'in_file', which should be an intermediary 1D file
    from individual-level analysis, containing the mean of the output across
    all voxels

    it then parses this value and writes it to a .csv file named
    output_means.csv located in the subject's output directory
    '''

    import os

    if os.path.exists(in_file):

        mean_oned_file = open(in_file, 'rU')
        line = mean_oned_file.readline()
        mean_oned_file.close()

        line = line.split('[')[0].strip(' ')
        
        # get filename of input maskave 1D file
        filename = in_file.split("/")[-1]
        filename = filename[0:-3]
        
        
        split_fullpath = in_file.split("/")
        
        if ("_mask_" in in_file) and (("sca_roi" in in_file) or \
            ("sca_tempreg" in in_file)):
            
            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname
                    
            filename = split_fullpath[-1]
            
            if ".1D" in filename:
                filename = filename.replace(".1D","")
            
            resource_name = output_name + "_%s_%s" % (maskname, filename)

            
        elif ("_spatial_map_" in in_file) and \
            ("dr_tempreg" in in_file):
            
            for dirname in split_fullpath:
                if "_spatial_map_" in dirname:
                    mapname = dirname
                    
            filename = split_fullpath[-1]
            
            if ".1D" in filename:
                filename = filename.replace(".1D","")
            
            resource_name = output_name + "_%s_%s" % (mapname, filename)
            
            
        elif ("_mask_" in in_file) and ("centrality" in in_file):
            
            for dirname in split_fullpath:
                if "_mask_" in dirname:
                    maskname = dirname
                    
            filename = split_fullpath[-1]
            
            if ".1D" in filename:
                filename = filename.replace(".1D","")
            
            resource_name = output_name + "_%s_%s" % (maskname, filename)
            
            
        else:
        
            resource_name = output_name
        

        output_means_file = os.path.join(os.getcwd(), 'mean_%s.txt' % resource_name)
        output_means = open(output_means_file, 'wb')

        print >>output_means, line

        output_means.close()


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

    subID = subject_dir.split('/')[len(subject_dir.split('/'))-1]
    means_dir = os.path.join(subject_dir, 'output_means')

    # extract the mean values
    for root, dirs, files in os.walk(means_dir):

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
    csv_file = open(csv_file_path, 'wt')

    deriv_string = ''
    val_string = ''

    for deriv in output_vals.keys():
        if deriv_string == '':
            deriv_string = deriv
            val_string = output_vals[deriv]
        else:
            deriv_string = deriv_string + ',' + deriv
            val_string = val_string + ',' + output_vals[deriv]

    print >>csv_file, deriv_string
    print >>csv_file, val_string


    csv_file.close()

    return


def dbg_file_lineno():
    cf=currentframe()
    return cf.f_back.f_code.co_filename, cf.f_back.f_lineno


