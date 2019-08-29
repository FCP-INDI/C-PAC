
def test_nonlinear_register():
    from ..registration import create_nonlinear_register
    
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    
    ## necessary inputs
    ## -input_brain
    ## -input_skull
    ## -reference_brain
    ## -reference_skull
    ## -fnirt_config
    ## -fnirt_warp_res
    
    ## input_brain
    anat_bet_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/anatpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/anat_skullstrip/mprage_anonymized_RPI_3dT.nii.gz'
    
    ## input_skull
    
    ## reference_brain
    mni_file = '/usr/share/fsl/4.1/data/standard/MNI152_T1_3mm_brain.nii.gz'
    
    ## reference_skull
    
    ## fnirt_config
    fnirt_config = 'T1_2_MNI152_3mm'
    
    ## fnirt_warp_res
    fnirt_warp_res = None
    
    #?? what is this for?:
    func_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/nuisance_preproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_csf_threshold_0.4/_gm_threshold_0.2/_wm_threshold_0.66/_run_scrubbing_False/_nc_5/_selector_6.7/regress_nuisance/mapflow/_regress_nuisance0/residual.nii.gz'
    
    
    mni_workflow = pe.Workflow(name='mni_workflow')

    linear_reg = pe.Node(interface=fsl.FLIRT(),
                         name='linear_reg_0')
    linear_reg.inputs.cost = 'corratio'
    linear_reg.inputs.dof = 6
    linear_reg.inputs.interp = 'nearestneighbour'
    
    linear_reg.inputs.in_file = func_file
    linear_reg.inputs.reference = anat_bet_file
    
    #T1 to MNI Node
    c = create_nonlinear_register()
    c.inputs.inputspec.input = anat_bet_file
    c.inputs.inputspec.reference = '/usr/share/fsl/4.1/data/standard/MNI152_T1_3mm_brain.nii.gz'
    c.inputs.inputspec.fnirt_config = 'T1_2_MNI152_3mm'
    
    #EPI to MNI warp Node
    mni_warp = pe.Node(interface=fsl.ApplyWarp(),
                       name='mni_warp')
    mni_warp.inputs.ref_file = '/usr/share/fsl/4.1/data/standard/MNI152_T1_3mm_brain.nii.gz'
    mni_warp.inputs.in_file = func_file

    mni_workflow.connect(c, 'outputspec.nonlinear_xfm',
                         mni_warp, 'field_file')
    mni_workflow.connect(linear_reg, 'out_matrix_file',
                         mni_warp, 'premat')
    
    mni_workflow.base_dir = './'
    mni_workflow.run()    
    
def test_registration():
    from ..registration import create_nonlinear_register
    from ..registration import create_register_func_to_mni
    
    import nipype.pipeline.engine as pe
    import nipype.interfaces.fsl as fsl
    
    func_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/nuisance_preproc/_session_id_NYU_TRT_session1_subject_id_sub05676/_csf_threshold_0.4/_gm_threshold_0.2/_wm_threshold_0.66/_run_scrubbing_False/_nc_5/_selector_6.7/regress_nuisance/mapflow/_regress_nuisance0/residual.nii.gz'
    anat_skull_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/anatpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/anat_reorient/mprage_anonymized_RPI.nii.gz'
    anat_bet_file = '/home/data/Projects/nuisance_reliability_paper/working_dir_CPAC_order/resting_preproc/anatpreproc/_session_id_NYU_TRT_session1_subject_id_sub05676/anat_skullstrip/mprage_anonymized_RPI_3dT.nii.gz'
    mni_brain_file = '/usr/share/fsl/4.1/data/standard/MNI152_T1_3mm_brain.nii.gz'
    mni_skull_file = '/usr/share/fsl/4.1/data/standard/MNI152_T1_3mm.nii.gz'

    
    mni_workflow = pe.Workflow(name='mni_workflow')
    
    nr = create_nonlinear_register()
    nr.inputs.inputspec.input_brain = anat_bet_file
    nr.inputs.inputspec.input_skull = anat_skull_file
    nr.inputs.inputspec.reference_brain = mni_brain_file
    nr.inputs.inputspec.reference_skull = mni_skull_file
    nr.inputs.inputspec.fnirt_config = '/usr/share/fsl/4.1/etc/flirtsch/T1_2_MNI152_3mm.cnf'
    func2mni = create_register_func_to_mni()
    func2mni.inputs.inputspec.func = func_file
    func2mni.inputs.inputspec.mni = mni_brain_file
    func2mni.inputs.inputspec.anat = anat_bet_file
    
    mni_workflow.connect(nr, 'outputspec.nonlinear_xfm',
                         func2mni, 'inputspec.anat_to_mni_xfm')
    mni_workflow.base_dir = './mni_05676_3'
    mni_workflow.run()


def test_registration_lesion():
    import os
    import nipype.pipeline.engine as pe
    from ..registration import create_wf_calculate_ants_warp
    from CPAC.anat_preproc.anat_preproc import create_anat_preproc
    from CPAC.anat_preproc.lesion_preproc import create_lesion_preproc

    # Skull stripped anat image
    anat_file = '/bids_dataset/sub-0027228/ses-1/anat/sub-0027228_ses-1_run-1_T1w.nii.gz'
    lesion_file = '/bids_dataset/sub-0027228/ses-1/anat/sub-0027228_ses-1_run-1_T1w_lesion-mask.nii.gz'
    mni_brain_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain.nii.gz'

    if not os.path.exists(anat_file):
        raise IOError(anat_file + ' not found')
    if not os.path.exists(lesion_file):
        raise IOError(lesion_file + ' not found')
    if not os.path.exists(mni_brain_file):
        raise IOError(mni_brain_file + ' not found')

    wf = pe.Workflow(name='test_reg_lesion')

    anat_preproc = create_anat_preproc(method='mask',
                                       already_skullstripped=True,
                                       wf_name='anat_preproc')

    anat_preproc.inputs.inputspec.anat = anat_file

    lesion_preproc = create_lesion_preproc(
        wf_name='lesion_preproc'
    )

    lesion_preproc.inputs.inputspec.lesion = lesion_file

    ants_reg_anat_mni = \
        create_wf_calculate_ants_warp(
            'anat_mni_ants_register',
            0,
            num_threads=4
        )

    # pass the reference file
    ants_reg_anat_mni.inputs.inputspec.reference_brain = mni_brain_file

    wf.connect(
        anat_preproc, 'outputspec.reorient',
        ants_reg_anat_mni, 'inputspec.anatomical_brain'
    )

    wf.connect(
        lesion_preproc, 'outputspec.reorient',
        ants_reg_anat_mni, 'inputspec.fixed_image_mask'
    )

    ants_reg_anat_mni.inputs.inputspec.set(
        dimension=3,
        use_histogram_matching=True,
        winsorize_lower_quantile=0.01,
        winsorize_upper_quantile=0.99,
        metric=['MI', 'MI', 'CC'],
        metric_weight=[1, 1, 1],
        radius_or_number_of_bins=[32, 32, 4],
        sampling_strategy=['Regular', 'Regular', None],
        sampling_percentage=[0.25, 0.25, None],
        number_of_iterations=[
            [1000, 500, 250, 100],
            [1000, 500, 250, 100],
            [100, 100, 70, 20]
        ],
        convergence_threshold=[1e-8, 1e-8, 1e-9],
        convergence_window_size=[10, 10, 15],
        transforms=['Rigid', 'Affine', 'SyN'],
        transform_parameters=[[0.1], [0.1], [0.1, 3, 0]],
        shrink_factors=[
            [8, 4, 2, 1],
            [8, 4, 2, 1],
            [6, 4, 2, 1]
        ],
        smoothing_sigmas=[
            [3, 2, 1, 0],
            [3, 2, 1, 0],
            [3, 2, 1, 0]
        ]
    )

    wf.run()
