def test_median_angle_correct():
    from CPAC.median_angle import median_angle_correct
    import numpy as np
    import nibabel as nb
    
    def getY(filepath):
        nii = nb.load(filepath)
        data = nii.get_data().astype(np.float64)
        mask = (data != 0).sum(-1) != 0
        
        return data[mask].T
    
    def normalize(X):
        Xc = X - X.mean(0)
        return Xc/np.sqrt( (Xc**2).sum(0) )
    
    subject = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/funcpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/func_scale/mapflow/_func_scale0/lfo_3dc_RPI_3dv_3dc_maths.nii.gz'
    target_angle = 88.0
    
    Y_orig = normalize(getY(subject))
    U_orig, S, Vh = np.linalg.svd(Y_orig, full_matrices=False)
    
    corrected_file, angles_file = median_angle_correct(target_angle, subject)
    
    Y_corr = normalize(getY(corrected_file))
    
    median_angle_orig = np.median(np.arccos(U_orig[:,0].T.dot(Y_orig)))
    median_angle_corr = np.median(np.arccos(U_orig[:,0].T.dot(Y_corr)))
    
    print(median_angle_orig*180.0/np.pi, median_angle_corr*180.0/np.pi)
    
    