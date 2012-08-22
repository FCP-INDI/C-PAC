def test_nuisance():
    from CPAC.nuisance import create_nuisance
    cn = create_nuisance()
    subjects_list = open('/home/data/Projects/ADHD200/adhd200_sublist_for_basc_withGSR').readlines()
    subjects_list = [ subject.strip() for subject in subjects_list ]
    
    subject = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/funcpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/func_scale/mapflow/_func_scale0/lfo_3dc_RPI_3dv_3dc_maths.nii.gz'
    csf_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/segpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_csf_threshold_0.4/seg_mask/mapflow/_seg_mask0/segment_prob_0_flirt_maths_maths_maths.nii.gz'
    wm_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/segpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_wm_threshold_0.66/seg_mask1/mapflow/_seg_mask10/segment_prob_2_flirt_maths_maths_maths.nii.gz'
    gm_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/segpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_gm_threshold_0.2/seg_mask2/mapflow/_seg_mask20/segment_prob_1_flirt_maths_maths_maths.nii.gz'
    motion_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/funcpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/func_volreg_1/mapflow/_func_volreg_10/lfo_3dc_RPI_3dv1D.1D'
    
    stor = {'compcor' : True,
            'wm' : True,
            'csf' : True,
            'gm' : True,
            'global' : True,
            'pc1' : True,
            'motion' : True,
            'linear' : True,
            'quadratic' : True}
    
    stor2 = {'compcor' : True,
            'wm' : True,
            'csf' : True,
            'gm' : True,
            'global' : False,
            'pc1' : True,
            'motion' : False,
            'linear' : True,
            'quadratic' : False}
      
    cn.inputs.inputspec.subject = subject
    cn.inputs.inputspec.gm_mask = gm_file
    cn.inputs.inputspec.wm_mask = wm_file
    cn.inputs.inputspec.csf_mask = csf_file
    cn.inputs.inputspec.motion_components = motion_file
    cn.get_node('residuals').iterables = ('selector',[stor, stor2])
    cn.inputs.inputspec.compcor_ncomponents = 5
    
    cn.run(plugin='MultiProc', plugin_args={'n_procs' : 2})