
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
    'movement_parameters': 'parameters',
    'max_displacement':'parameters',
    'preprocessed':'func',
    'functional_brain_mask':'func',
    'motion_correct':'func',
    'anatomical_to_functional_xfm':'registration',
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
    'sca_correlations_roi_based':'sca_roi',
    'sca_z_trans_correlations_roi_based':'sca_roi',
    'sca_correlations_seed_based':'sca_mask',
    'sca_z_trans_correlations_seed_based':'sca_mask',
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
    'seg_partial_volume_files': 'anat'
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

    session = 'scan_func_'

    if 'anat_1' in remainder_path:
        print '----------R----------------->', remainder_path

    lst = remainder_path.split('/')

    lst = [x for x in lst if not ('' == x) ]


    for element in lst:

        if 'scan_func' in element:

            session += element.split('scan_func_')[1] + '_'


    if session.endswith('_'):

        session = session.rstrip('_')

    return session


def get_hplpfwhmseed_(parameter, remainder_path):


    partial_parameter_value = remainder_path.split(parameter)[1]

    print partial_parameter_value, ' ~~~~~~~~ ', parameter

    value = partial_parameter_value.split('/')[0]

    return parameter.lstrip('/_')+value



def create_symbolic_links(pipeline_id, relevant_strategies, path, subject_id):


    import os
    import commands
    from CPAC.utils.utils import get_workflow, get_session, \
                     get_hplpfwhmseed_

    #from sink import global_lock


    for strategy in relevant_strategies:

        base_path, remainder_path = path.split(subject_id)

        sym_path = path.split(pipeline_id)[0]

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
                    '_compcor_':'',
                    '_target_angle_deg':'MEDIANangle_', '_wm_threshold':'WM_'}

            strategy_identifier = ''

            for el in strategy:
                key, value = el.rsplit('_', 1)

                print key, ' ----------> ', value

                if '_compcor_'in key:
                    key = '_compcor_'

                if not 'pipeline' in key:
                    strategy_identifier += short_names[key] + value + '_'

            strategy_identifier = strategy_identifier.rsplit('_', 1)[0]

#            global_lock.acquire()

            f = open(os.path.join(os.path.dirname(base_path), 'paths_file.txt'), 'a')

            print >>f, path
#            global_lock.release()


        except:
            print str(strategy), " not in labels_dict"
            raise


        #remove unused corrections
        strategy_identifier = strategy_identifier.replace('pc10.', '')
        strategy_identifier = strategy_identifier.replace('linear0.', '')
        strategy_identifier = strategy_identifier.replace('wm0.', '')
        strategy_identifier = strategy_identifier.replace('global0.', '')
        strategy_identifier = strategy_identifier.replace('motion0.', '')
        strategy_identifier = strategy_identifier.replace('quadratic0.', '')
        strategy_identifier = strategy_identifier.replace('gm0.', '')
        strategy_identifier = strategy_identifier.replace('compcor0.', '')


#        strategy_identifier = 'regressors.' + strategy_identifier

        #start making basic sym link directories
        sym_path = os.path.join(sym_path, strategy_identifier)
        new_sub_path = os.path.join(sym_path, subject_id)
        file_name, wf, remainder_path = get_workflow(remainder_path)
        new_wf_path = os.path.join(new_sub_path, wf)
        session = get_session(remainder_path)
        new_session_path = os.path.join(new_wf_path, session)

        new_path = new_session_path


        #bring into use the tier 2 iterables for recursive directory structure

        if '/_mask_' in remainder_path:

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_mask_', remainder_path))

        if '/_roi_' in remainder_path:

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_roi_', remainder_path))


        if ('_hp_'  in remainder_path):

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_hp_', remainder_path))


        if ('_lp_' in remainder_path):

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_lp_', remainder_path))

        if ('_bandpass_freqs_' in remainder_path):

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_bandpass_freqs_', remainder_path))


        if ('_fwhm_' in remainder_path):

            new_path = os.path.join(new_path, get_hplpfwhmseed_('/_fwhm_', remainder_path))


        try:
            os.makedirs(new_path)
        except:
            print '.'

        fname = os.path.basename(path)

        ext = fname.split('.', 1)[1]
        ext = '.' + (ext)

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
        'seg_partial_volume_files']

        if file_name in dont_change_fname:

            cmd = 'ln -s %s %s' % (path, os.path.join(new_path, fname))
            print cmd
            commands.getoutput(cmd)
        else:

            cmd = 'ln -s %s %s' % (path, os.path.join(new_path, file_name+ ext))
            try:

                f1 = open(os.path.join(new_path, file_name+ ext))

                print >>f, 'Error-----------> ',file_name+ ext, ' already exists'
                print >>f, '~~~~~~~~~~', fname
            except:
                print cmd
                commands.getoutput(cmd)




def prepare_symbolic_links(in_file, strategies, subject_id, pipeline_id):

    from  CPAC.utils.utils import get_strategies_for_path, create_symbolic_links


    for path in in_file:

        for strategy in strategies:

            strategy.append(pipeline_id)

        relevant_strategies = get_strategies_for_path(path, strategies)

        create_symbolic_links(pipeline_id, relevant_strategies, path, subject_id)


def modify_model_files(model_files, sublist):

    def read_model_file(file):

        dict1 ={}
        f = open(file, 'r')
        for line in open(file, 'r'):
            if line.strip() !='':
                if 'NumWaves' in line: 
                    dict1['/NumWaves']= line.split()[1]
                elif 'NumPoints' in line:
                    dict1['/NumPoints'] = line.split()[1]
                elif 'PPheights' in line:
                    dict1['/PPHeights'] = line.split()[1:]
                elif 'Matrix' in line:
                    dict1['/Matrix'] = []
                else:
                    dict1.get('/Matrix').append(line)
        return dict1

