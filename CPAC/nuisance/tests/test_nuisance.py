def test_calc_residuals():
    import numpy as np
    from CPAC.nuisance import calc_residuals
    from scipy.io import loadmat
    from scipy.stats import pearsonr
    import nibabel as nb
    
    def normalize(X):
        Xc = X - X.mean(0)
        return Xc/np.sqrt( (Xc**2).sum(0) )
    
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

    calc_residuals(subject, stor, wm_file, csf_file, gm_file, motion_file, compcor_ncomponents=3)
    X = loadmat('nuisance_regressors.mat')
    
    nii = nb.load('residual.nii.gz')
    data = nii.get_data().astype(np.float64)
    global_mask = (data != 0).sum(-1) != 0
    Y = data[global_mask].T

    r_glb = normalize(Y).T.dot(X['global'])
    r_csf = normalize(Y).T.dot(X['csf'])
    r_wm = normalize(Y).T.dot(X['wm'])
    

def test_nuisance():
#    from CPAC.nuisance import create_nuisance
#    cn = create_nuisance()
#    subjects_list = open('/home/data/Projects/ADHD200/adhd200_sublist_for_basc_withGSR').readlines()
#    subjects_list = [ subject.strip() for subject in subjects_list ]
#    
#    subject = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/funcpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/func_scale/mapflow/_func_scale0/lfo_3dc_RPI_3dv_3dc_maths.nii.gz'
#    csf_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/segpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_csf_threshold_0.4/seg_mask/mapflow/_seg_mask0/segment_prob_0_flirt_maths_maths_maths.nii.gz'
#    wm_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/segpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_wm_threshold_0.66/seg_mask1/mapflow/_seg_mask10/segment_prob_2_flirt_maths_maths_maths.nii.gz'
#    gm_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/segpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_gm_threshold_0.2/seg_mask2/mapflow/_seg_mask20/segment_prob_1_flirt_maths_maths_maths.nii.gz'
#    motion_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/funcpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/func_volreg_1/mapflow/_func_volreg_10/lfo_3dc_RPI_3dv1D.1D'
#    
#    stor = {'compcor' : True,
#            'wm' : True,
#            'csf' : True,
#            'gm' : True,
#            'global' : True,
#            'pc1' : True,
#            'motion' : True,
#            'linear' : True,
#            'quadratic' : True}
#    
#    stor2 = {'compcor' : True,
#            'wm' : True,
#            'csf' : True,
#            'gm' : True,
#            'global' : False,
#            'pc1' : True,
#            'motion' : False,
#            'linear' : True,
#            'quadratic' : False}
#      
#    cn.inputs.inputspec.subject = subject
#    cn.inputs.inputspec.gm_mask = gm_file
#    cn.inputs.inputspec.wm_mask = wm_file
#    cn.inputs.inputspec.csf_mask = csf_file
#    cn.inputs.inputspec.motion_components = motion_file
#    cn.get_node('residuals').iterables = ('selector',[stor, stor2])
#    cn.inputs.inputspec.compcor_ncomponents = 5
#    
#    cn.run(plugin='MultiProc', plugin_args={'n_procs' : 2})
    
    from nipype import config
    config.update_config({'execution': {'remove_unnecessary_outputs':False}})
    from CPAC.nuisance import create_nuisance
    cn = create_nuisance()
    stor = {'compcor' : True,
            'wm' : True,
            'csf' : True,
            'gm' : True,
            'global' : True,
            'pc1' : True,
            'motion' : True,
            'linear' : True,
            'quadratic' : True}
    cn.inputs.inputspec.selector = stor
    cn.inputs.inputspec.compcor_ncomponents = 5
    cn.inputs.inputspec.motion_components = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/movement_parameters/_scan_rest_1_rest/rest_3dc_RPI_3dv1D.1D'
    cn.inputs.inputspec.func_to_anat_linear_xfm = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/functional_to_anat_linear_xfm/_scan_rest_1_rest/rest_3dc_RPI_3dv_3dc_3dT_flirt.mat'
    cn.inputs.inputspec.mni_to_anat_linear_xfm = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/mni_to_anatomical_linear_xfm/mprage_RPI_3dc_flirt_inv.mat'
    cn.inputs.inputspec.gm_mask = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/anatomical_gm_mask/_gm_threshold_0.2/segment_prob_1_maths_maths_maths.nii.gz'
    cn.inputs.inputspec.csf_mask = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/anatomical_csf_mask/_csf_threshold_0.4/segment_prob_0_maths_maths_maths.nii.gz'
    cn.inputs.inputspec.wm_mask = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/anatomical_wm_mask/_wm_threshold_0.66/segment_prob_2_maths_maths_maths.nii.gz'
    cn.inputs.inputspec.harvard_oxford_mask = '/usr/share/fsl/4.1/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz'
    cn.inputs.inputspec.subject = '/home/data/PreProc/ABIDE_CPAC_test_1/pipeline_0/0050102_session_1/preprocessed/_scan_rest_1_rest/rest_3dc_RPI_3dv_3dc_maths.nii.gz'
    cn.base_dir = '/home/bcheung/cn_run'
