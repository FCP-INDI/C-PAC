import threading
global_lock = threading.Lock()



files_folders_wf = {
    'anatomical_brain': 'anat',
    'anatomical_reorient': 'anat',
    'anatomical_to_mni_linear_xfm': 'anat',
    'mni_to_anatomical_linear_xfm': 'anat',
    'anatomical_to_mni_nonlinear_xfm': 'anat',
    'anatomical_gm_mask': 'anat',
    'anatomical_csf_mask': 'anat',
    'anatomical_wm_mask': 'anat',
    'mean_functional': 'func',
    'functional_preprocessed_mask': 'func',
    'functional_to_spatial_map': 'func',
    'functional_mask_to_spatial_map': 'func',
    'slice_timing_corrected': 'func',
    'movement_parameters': 'parameters',
    'max_displacement':'parameters',
    'preprocessed':'func',
    'functional_brain_mask':'func',
    'motion_correct':'func',
    'mean_functional_in_mni' : 'func',
    'mean_functional_in_anat' : 'func',
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
    'alff_img':'alff',
    'falff_img':'alff',
    'alff_Z_img':'alff',
    'falff_Z_img':'alff',
    'functional_freq_filtered':'func',
    'scrubbing_movement_parameters':'parameters',
    'scrubbing_frames_included':'parameters',
    'scrubbing_frames_excluded':'parameters',
    'motion_params':'parameters',
    'power_params':'parameters',
    'scrubbed_preprocessed':'func',
    'functional_mni':'func',
    'functional_brain_mask_to_standard':'func',
    'functional_to_anat_linear_xfm':'registration',
    'functional_to_mni_linear_xfm':'registration',
    'mni_to_functional_linear_xfm':'registration',
    'mni_normalized_anatomical':'anat',
    'vmhc_raw_score':'vmhc',
    'vmhc_z_score':'vmhc',
    'vmhc_z_score_stat_map':'vmhc',
    'raw_reho_map':'reho',
    'reho_Z_img':'reho',
    'alff_Z_to_standard':'alff',
    'falff_Z_to_standard':'alff',
    'alff_Z_smooth':'alff',
    'falff_Z_smooth':'alff',
    'alff_Z_to_standard_smooth':'alff',
    'falff_Z_to_standard_smooth':'alff',
    'reho_Z_to_standard':'reho',
    'reho_Z_smooth':'reho',
    'reho_Z_to_standard_smooth':'reho',
    'voxel_timeseries':'timeseries',
    'roi_timeseries':'timeseries',
    'sca_roi_correlations':'sca_roi',
    'sca_roi_Z':'sca_roi',
    'sca_seed_correlations':'sca_mask',
    'sca_seed_Z':'sca_mask',
    'sca_seed_Z_to_standard':'sca_mask',
    'sca_roi_Z_to_standard':'sca_roi',
    'sca_seed_Z_smooth':'sca_mask',
    'sca_seed_Z_to_standard_smooth':'sca_mask',
    'sca_roi_Z_smooth':'sca_roi',
    'sca_roi_Z_to_standard_smooth':'sca_roi',
    'bbregister_registration': 'surface_registration',
    'left_hemisphere_surface': 'surface_registration',
    'right_hemisphere_surface': 'surface_registration',
    'vertices_timeseries': 'timeseries',
    'centrality_outputs_smoothed':'centrality',
    'centrality_outputs_zscore':'centrality',
    'centrality_outputs':'centrality',
    'centrality_graphs':'centrality',
    'seg_probability_maps': 'anat',
    'seg_mixeltype': 'anat',
    'seg_partial_volume_map': 'anat',
    'seg_partial_volume_files': 'anat',
    'spatial_map_timeseries': 'timeseries',
    'dr_tempreg_maps_stack': 'spatial_regression',
    'dr_tempreg_maps_z_stack': 'spatial_regression',
    'dr_tempreg_maps_z_files': 'spatial_regression',
    'dr_tempreg_maps_stack_smooth': 'spatial_regression',
    'dr_tempreg_maps_z_stack_smooth': 'spatial_regression',
    'dr_tempreg_maps_z_files_smooth':'spatial_regression',
    'sca_tempreg_maps_stack': 'sca_roi',
    'sca_tempreg_maps_z_stack': 'sca_roi',
    'sca_tempreg_maps_z_files' : 'sca_roi',
    'sca_tempreg_maps_stack_smooth': 'sca_roi',
    'sca_tempreg_maps_z_stack_smooth': 'sca_roi',
    'sca_tempreg_maps_z_files_smooth': 'sca_roi',
}

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

    op_string = ""

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


    partial_parameter_value = remainder_path.split(parameter)[1]

    # print partial_parameter_value, ' ~~~~', parameter

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




def create_symbolic_links(pipeline_id, relevant_strategies, path, subject_id):

    import os
    import commands
    from CPAC.utils.utils import get_workflow, get_session, \
                     get_hplpfwhmseed_

    # from sink import global_lock


    for strategy in relevant_strategies:

        base_path, remainder_path = path.split(subject_id, 1)


        sym_path = path.split(pipeline_id)[0]

        file_path = os.path.join(sym_path, pipeline_id)
        file_path = os.path.join(file_path, subject_id)
        sym_path = os.path.join(sym_path, 'sym_links')

        sym_path = os.path.join(sym_path, pipeline_id)

        try:
            os.makedirs(sym_path)
        except:
            print '.'

        strategy_identifier = None

        try:

            short_names = {'_threshold':'SCRUB_', '_csf_threshold':'CSF_',
                    '_gm_threshold':'GM_',
                    '_compcor_':'compcor',
                    '_target_angle_deg':'MEDIANangle_', '_wm_threshold':'WM_'}

            strategy_identifier = ''

            for el in strategy:
                key, value = el.rsplit('_', 1)

                print key, ' ----------> ', value

                if '_compcor_'in key:

                    if 'compcor0' in value:
                        strategy_identifier += (value + '_')

                        continue
                    else:
                        val1 = key.split('_selector')[0]
                        strategy_identifier += val1 + '_' + value + '_'
                        continue

                if not 'pipeline' in key:
                    strategy_identifier += short_names[key] + value + '_'

            strategy_identifier = strategy_identifier.rsplit('_', 1)[0]



        except:
            print str(strategy), " not in labels_dict"
            raise


        # remove unused corrections
        strategy_identifier = strategy_identifier.replace('pc10.', '')
        strategy_identifier = strategy_identifier.replace('linear0.', '')
        strategy_identifier = strategy_identifier.replace('wm0.', '')
        strategy_identifier = strategy_identifier.replace('global0.', '')
        strategy_identifier = strategy_identifier.replace('motion0.', '')
        strategy_identifier = strategy_identifier.replace('quadratic0.', '')
        strategy_identifier = strategy_identifier.replace('gm0.', '')
        strategy_identifier = strategy_identifier.replace('csf0_', '')
        strategy_identifier = strategy_identifier.replace('compcor0.', '')

#        strategy_identifier = 'regressors.' + strategy_identifier

        # start making basic sym link directories
        sym_path = os.path.join(sym_path, strategy_identifier)
        new_sub_path = os.path.join(sym_path, subject_id)
        file_name, wf, remainder_path = get_workflow(remainder_path)
        session = get_session(remainder_path)
        new_session_path = os.path.join(new_sub_path, session)
        new_wf_path = os.path.join(new_session_path, wf)
        new_path = new_wf_path


        # bring into use the tier 2 iterables for recursive directory structure

        scan_info = '~~~'
        if '/_scan_' in path:

            scan_info = get_hplpfwhmseed_('/_scan_', path)

        print scan_info, ' ~~~', remainder_path

        if '/_mask_' in remainder_path:

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_mask_', remainder_path))

        if '/_roi_' in remainder_path:

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_roi_', remainder_path))

#        if '/_sca_roi_' in remainder_path:

#            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_sca_roi_', remainder_path))

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


        fwhm_str = ''


        spatial_map = ''
        if('_spatial_map_' in remainder_path):
            spatial_map = get_hplpfwhmseed_('/_spatial_map_', remainder_path)
            new_path = os.path.join(new_path, spatial_map)

        if ('_fwhm_' in remainder_path):


            fwhm_str = get_hplpfwhmseed_('/_fwhm_', remainder_path)
            new_path = os.path.join(new_path, fwhm_str)


        try:
            os.makedirs(new_path)
        except:
            print '.'


        try:
            new_f_path = os.path.join(file_path, 'path_files_here')
            os.makedirs(new_f_path)
        except:
            print '.'


        try:


            global global_lock
            global_lock.acquire()

            if not '~~~' in scan_info:
                f = open(os.path.join(new_f_path, 'paths_file_%s.txt') % (scan_info + '_' + strategy_identifier + '_' + bp_freq + '_' + hp_str + '_' + lp_str + '_' + fwhm_str), 'a')

                print >> f, path
            global_lock.release()

        except:
            print 'trouble acquiring locks or opening file skipping :', os.path.join(new_f_path, 'paths_file_%s.txt') % new_path.replace('/', '_')
            raise

        fname = os.path.basename(path)

        ext = fname.split('.', 1)[1]
        ext = '.' + (ext)

        # special case for ROI , need the ROI number

        if '_ROI_' in fname and 'sca_' in path:
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
        'dr_tempreg_maps_z_files',
        'dr_tempreg_maps_z_files_smooth',
        'sca_tempreg_maps_z_files',
        'sca_tempreg_maps_z_files_smooth']

        if file_name in dont_change_fname:

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

    if 'sca_roi_Z' in resource and '/_roi_' in in_file:

        third_tier = resource + '_' + get_param_val_('/_roi_', in_file)

        roi_number = ''.join(['ROI_', get_param_val_('/ROI_number_', in_file)])
        third_tier = third_tier + '/' + roi_number

    elif ('dr_tempreg_maps_z_files' in resource and '/temp_reg_map_z_' in in_file):

        third_tier = resource + '_' + get_param_val_('/_spatial_map_', in_file)
        third_tier = third_tier + '/' + get_param_val_('/temp_reg_map_z_', in_file)

    elif ('sca_tempreg_maps_z_files' in resource and '/sca_tempreg_z_maps_roi_' in in_file):
        third_tier = resource + '_' + get_param_val_('/_roi_', in_file)
        roi_number = ''.join(['ROI_', get_param_val_('/sca_tempreg_z_maps_roi_', in_file)])
        third_tier = third_tier + '/' + roi_number


    elif ('sca_seed_Z' in resource or 'centrality_outputs' in resource)  and '/_mask_' in in_file:

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

            print '~~~~~~ ', key, ' ~~~ ', el
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


def prepare_symbolic_links(in_file, strategies, subject_id, pipeline_id, helper):

    from  CPAC.utils.utils import get_strategies_for_path, create_symbolic_links, clean_strategy


    for path in in_file:

        for strategy in strategies:

            strategy.append(pipeline_id)

        relevant_strategies = get_strategies_for_path(path, strategies)

        cleaned_strategies = clean_strategy(relevant_strategies, helper)

        create_symbolic_links(pipeline_id, cleaned_strategies, path, subject_id)


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


def select_model_files(model, ftest):
    """
    Method to select model files
    """
    import os
    import glob

    try:
        files = glob.glob(os.path.join(model, '*'))

        if len(files) == 0:
            raise Exception("No files foudn inside model %s" % model)

        fts_file = ''
        for file in files:
            if file.endswith('.mat'):
                mat_file = file
            elif file.endswith('.grp'):
                grp_file = file
            elif file.endswith('.fts') and ftest:
                 fts_file = file
            elif file.endswith('.con'):
                 con_file = file

    except Exception:
        print "All the model files are not present. Please check the model folder %s" % model
        raise

    return fts_file, con_file, grp_file, mat_file



def get_scan_params(subject, scan, subject_map, start_indx, stop_indx):

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

    TR = float(check('tr', True))
    pattern = str(check('acquisition', True))
    ref_slice = int(check('reference', True))
    first_tr = check2(check('first_tr', False))
    last_tr = check2(check('last_tr', False))
    unit = 's'
    # if empty override with config information
    if first_tr == '':
        first_tr = start_indx

    if last_tr == '':
        last_tr = stop_indx

    if pattern not in ['alt+z', 'altplus', 'alt+z2', 'alt-z', 'altminus',
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
            if  max_slice_offset > TR:
                warnings.warn("TR is in seconds and slice timings are in milliseconds."\
                              "Converting TR into milliseconds")
                TR = TR * 1000
                print "New TR value %.2f ms" % TR
                unit = 'ms'

    else:
        # check to see, if TR is in milliseconds, convert it into seconds
        if TR > 10:
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
    if tr != None:
       tr = re.search("\d+.\d+", str(tr)).group(0)
       tr = float(tr)
       if tr > 10:
           tr = tr / 1000.0
    return tr

