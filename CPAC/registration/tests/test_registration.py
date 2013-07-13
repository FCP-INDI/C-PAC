
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
