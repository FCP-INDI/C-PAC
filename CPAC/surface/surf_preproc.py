import os
import nipype.interfaces.utility as util
from CPAC.utils.interfaces.function import Function
from CPAC.pipeline import nipype_pipeline_engine as pe


def run_surface(post_freesurfer_folder,
                freesurfer_folder,
                subject,
                t1w_restore_image,
                atlas_space_t1w_image,
                atlas_transform,
                inverse_atlas_transform,
                atlas_space_bold,
                scout_bold,
                surf_atlas_dir,
                gray_ordinates_dir,
                gray_ordinates_res,
                high_res_mesh,
                low_res_mesh,
                subcortical_gray_labels,
                freesurfer_labels,
                fmri_res,
                smooth_fwhm):
    """
    Returns
    -------
    dtseries : str
        Path to the dtseries file.

    desikan_killiany : str
        Path to the Desikan-Killiany parcellation file.

    destrieux : str
        Path to the Destrieux parcellation file.
    """

    import os
    import subprocess

    freesurfer_folder = os.path.join(freesurfer_folder, 'recon_all')

    # DCAN-HCP PostFreeSurfer
    # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/run.sh', '--post_freesurfer_folder', post_freesurfer_folder, \
        '--freesurfer_folder', freesurfer_folder, '--subject', subject, \
        '--t1w_restore', t1w_restore_image, '--atlas_t1w', atlas_space_t1w_image, \
        '--atlas_transform', atlas_transform, '--inverse_atlas_transform', inverse_atlas_transform, \
        '--surfatlasdir', surf_atlas_dir, '--grayordinatesdir', gray_ordinates_dir, '--grayordinatesres', gray_ordinates_res, \
        '--hiresmesh', high_res_mesh, '--lowresmesh', low_res_mesh, \
        '--subcortgraylabels', subcortical_gray_labels, '--freesurferlabels', freesurfer_labels]
    subprocess.check_output(cmd)

    # DCAN-HCP fMRISurface
    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh
    cmd = ['bash', '/code/CPAC/surface/fMRISurface/run.sh',
           '--post_freesurfer_folder', post_freesurfer_folder,
           '--subject', subject, '--fmri', atlas_space_bold, '--scout',
           scout_bold, '--lowresmesh', low_res_mesh, '--grayordinatesres',
           gray_ordinates_res, '--fmrires', fmri_res, '--smoothingFWHM',
           smooth_fwhm]
    subprocess.check_output(cmd)

    dtseries = os.path.join(post_freesurfer_folder,
                            'MNINonLinear/Results/task-rest01/'
                            'task-rest01_Atlas.dtseries.nii')
    aparc = {'desikan_killiany': {
        164: os.path.join(post_freesurfer_folder, 'MNINonLinear'
                          f'{subject}.aparc.164k_fs_LR.dlabel.nii'),
        32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                         'fsaverage_LR32k',
                         f'{subject}.aparc.32k_fs_LR.dlabel.nii')},
             'destrieux': {
        164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.aparc.a2009s.164k_fs_LR.dlabel.nii'),
        32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                         'fsaverage_LR32k',
                         f'{subject}.aparc.a2009s.32k_fs_LR.dlabel.nii')}}

    return (dtseries, aparc['desikan_killiany'][164], aparc['destrieux'][164],
            aparc['desikan_killiany'][32], aparc['destrieux'][32])

def run_get_mri_info(mri_info):
    
        import os
        import subprocess
         
        cmd = ['mri_info', mri_info]
        mri_info_file = subprocess.check_output(cmd)
        
        with open('mri_info.txt', 'w') as f:
            for v in mri_info_file:

                f.write(str(mri_info_file))




        out_file = os.path.join(os.getcwd(), 'mri_info.txt')

        for line in open(out_file, 'r'):


            cr_val = line.strip().split('c_r =')
            cr = cr_val[1].strip().split('\\n')

            ca_val = line.strip().split('c_a =')
            ca = ca_val[1].strip().split('\\n')

            cs_val = line.strip().split('c_s =')
            cs = cs_val[1].strip().split('\\n')

        import re

        cr = re.sub(r"[\n\t\s]*", "", cr[0])
        cr = float(cr)

        ca = re.sub(r"[\n\t\s]*", "", ca[0])
        ca = float(ca)

        cs = re.sub(r"[\n\t\s]*", "", cs[0])
        cs = float(cs)


        cr_matrix = np.array([1, 0, 0, cr])
        ca_matrix = np.array([0, 1, 0, ca])
        cs_matrix = np.array([1, 0, 0, cs])
        id_mastrix = np.array([0, 0, 0, 1])

        final_mat = np.concatenate((cr_matrix, ca_matrix, cs_matrix, id_mastrix), axis=0)
        final_mat = final_mat.reshape(4,4)


        file_path = '/home/tgeorge/pfreesurfer-bash-runs/c_ras.mat'
        np.savetxt('file_path ',final_mat)
        
        return mri_info_file



def create_resmat(wf, cfg, strat_pool, pipe_num, opt=None):

    '''
    {"name": "create_resmat",
     "config": ["surface_analysis", "post_freesurfer"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["wmparc.mgz"],
     "outputs": ["final_mat"]}
    '''
    
    

    get_mri_info_imports = ['import os', 'import subprocess']
    get_mri_info = pe.Node(util.Function(input_names=['mri_info'],
                                               output_names=['mri_info_file'],
                                               function=run_get_mri_info,
                                               imports=get_mri_info_imports),
                                 name=f'get_mri_info_{pipe_num}')


    
    node, out = strat_pool.get_data('wmparc.mgz')
    wf.connect(node, out, get_mri_info, 'mri_info')
    
    

    



 
# def post_freesurfer_run(wf, cfg, strat_pool, pipe_num, opt=None):

    # '''
    # {"name": "surface_preproc",
     # "config": ["surface_analysis", "post_freesurfer"],
     # "switch": ["run"],
     # "option_key": "None",
     # "option_val": "None",
     # "inputs": ["freesurfer-subject-dir",
                # "desc-restore_T1w",
                # "space-template_desc-head_T1w",
                # "from-T1w_to-template_mode-image_xfm",
                # "from-template_to-T1w_mode-image_xfm",
                # "space-template_desc-brain_bold",
                # "space-template_desc-scout_bold"],
     # "outputs": ["space-fsLR_den-32k_bold-dtseries"]}
    # '''
    # image_list = ["wmparc", "aparc.a2009s+aseg",  "aparc+aseg"]
    # #Convert FreeSurfer Volumes
    # for Image in image_list:

	    # if os.path.exists(os.path.join(os.getcwd(), '$FreeSurferFolder"/mri', '$Image".mgz')):
		   # #mri_convert -rt nearest -rl "$T1wFolder"/"$T1wImage".nii.gz "$FreeSurferFolder"/mri/"$Image".mgz "$T1wFolder"/"$Image"_1mm.nii.gz

           # mri_convert = pe.Node(interface=freesurfer.preprocess.MRIConvert(),name='mri_convert')
           # mri_convert.inputs.in_file = "$FreeSurferFolder"/mri/"$Image".mgz
           # mri_convert.out_file = "$T1wFolder"/"$Image"_1mm.nii.gz
           # mri_convert.resample_type = 'nearest'
           # mri_convert.reslice_like = "$T1wFolder"/"$T1wImage".nii.gz


		   # #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$Image".nii.gz

           # apply_warp_1 = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_1')
           # apply_warp_1.inputs.interp = 'nn'
           # apply_warp_1.inputs.ref_file = "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage"
           # apply_warp_1.out_file = "$T1wFolder"/"$Image".nii.gz
           # apply_warp_1.relwarp = 'TRUE'
           # apply_warp_1.premat = '$FSLDIR/etc/flirtsch/ident.mat'
           # wf.connect(mri_convert ,'out_file',apply_warp_1,'in_file')



		   # #applywarp --rel --interp=nn -i "$T1wFolder"/"$Image"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$Image".nii.gz

           # apply_warp_2 = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_2')
           # apply_warp_2.inputs.interp = 'nn'
           # apply_warp_2.inputs.ref_file = "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage"
           # apply_warp_2.out_file = "$AtlasSpaceFolder"/"$Image".nii.gz
           # apply_warp_2.relwarp = 'TRUE'
           # apply_warp_2.field_file =  "$AtlasTransform"
           # wf.connect(apply_warp_1 ,'out_file',apply_warp_2,'in_file')



		   # #${CARET7DIR}/wb_command -volume-label-import "$T1wFolder"/"$Image".nii.gz "$FreeSurferLabels" "$T1wFolder"/"$Image".nii.gz -drop-unused-labels

           # wb_volume_label_1 = pe.Node(interface=base.CommandLine(),name='volume_label')
           # wb_volume_label_1.command='wb_command -volume-label-import', environ={'DISPLAY': ':1'})
           # volume_label_1.inputs.args = "in_file" "$FreeSurferLabels" "$T1wFolder"/"$Image".nii.gz "-drop-unused-labels"
           # wf.connect(apply_warp_1 ,'out_file',wb_volume_label_1,'in_file')



		   # #${CARET7DIR}/wb_command -volume-label-import "$AtlasSpaceFolder"/"$Image".nii.gz "$FreeSurferLabels" "$AtlasSpaceFolder"/"$Image".nii.gz -drop-unused-labels
           # wb_volume_label_2 = pe.Node(interface=base.CommandLine(),name='volume_label')
           # wb_volume_label_2.command="wb_command -volume-label-import" , environ={'DISPLAY': ':1'})
           # volume_label_2.inputs.args = "in_file" "$FreeSurferLabels" "$AtlasSpaceFolder"/"$Image".nii.gz '-drop-unused-labels'
           # wf.connect(apply_warp_2 ,'out_file',wb_volume_label_2,'in_file')



    # #Create FreeSurfer Brain Mask (Now done in PostFreeSurfer.sh so brainmask_fs.nii.gz exists for ANTs Registration)


    # #fslmaths "$T1wFolder"/wmparc_1mm.nii.gz -bin -dilD -dilD -dilD -ero -ero "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz


    # fsl_maths_1 = pe.Node(interface=fsl.maths.UnaryMaths(),name='fsl_maths_1')
    # fsl_maths_1.out_file = "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
    # fsl_maths_1.inputs.args = '-dilD -dilD -dilD -ero -ero'
    # fsl_maths_1.inputs.operation = 'bin'
    # fsl_maths_1.inputs.in_file = "$T1wFolder"/wmparc_1mm.nii.gz


    # #${CARET7DIR}/wb_command -volume-fill-holes "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz

    # volume-fill-holes = pe.Node(interface=base.CommandLine(),name='volume_label')
    # volume-fill-holes.command = "wb_command -volume-fill-holes" , environ={'DISPLAY': ':1'})
    # volume-fill-holes.inputs.args = "in_file" "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
    # wf.connect(fsl_maths_1 ,'out_file',volume-fill-holes,'in_file')

    # #fslmaths "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -bin "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
    # fsl_maths_2 = pe.Node(interface=fsl.maths.UnaryMaths(),name='fsl_maths_2')
    # fsl_maths_2.out_file = "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz
    # fsl_maths_2.inputs.operation = 'bin'
    # wf.connect(volume-fill-holes,'out_file',fsl_maths_2,'in_file')

    # #applywarp --rel --interp=nn -i "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" --premat=$FSLDIR/etc/flirtsch/ident.mat -o "$T1wFolder"/"$T1wImageBrainMask".nii.gz
    # apply_warp_1 = pe.Node(interface=fsl.ApplyWarp(),name='apply_warp_1')
    # apply_warp_1.inputs.interp = 'nn'
    # apply_warp_1.inputs.ref_file = "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage"
    # apply_warp_1.out_file = "$T1wFolder"/"$T1wImageBrainMask".nii.gz
    # apply_warp_1.relwarp = 'TRUE'
    # apply_warp_1.premat = '$FSLDIR/etc/flirtsch/ident.mat'
    # wf.connect(fsl_maths_2 ,'out_file',apply_warp_1,'in_file')

    # #applywarp --rel --interp=nn -i "$T1wFolder"/"$T1wImageBrainMask"_1mm.nii.gz -r "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage" -w "$AtlasTransform" -o "$AtlasSpaceFolder"/"$T1wImageBrainMask".nii.gz
    # apply_warp_2 = pe.Node(interface=fsl.ApplyWarp(),name="apply_warp_2")
    # apply_warp_2.inputs.interp = 'nn'
    # apply_warp_1.inputs.ref_file = "$AtlasSpaceFolder"/"$AtlasSpaceT1wImage"
    # apply_warp_2.out_file = "$AtlasSpaceFolder"/"$T1wImageBrainMask".nii.gz
    # apply_warp_2.field_file =  "$AtlasTransform"
    # apply_warp_2.relwarp = 'TRUE'
    # wf.connect(fsl_maths_2 ,'out_file',apply_warp_2,'in_file')




def surface_connector(wf, cfg, strat_pool, pipe_num, opt):
    
    surf = pe.Node(util.Function(input_names=['post_freesurfer_folder',
                                              'freesurfer_folder',
                                              'subject',
                                              't1w_restore_image',
                                              'atlas_space_t1w_image',
                                              'atlas_transform',
                                              'inverse_atlas_transform',
                                              'atlas_space_bold',
                                              'scout_bold',
                                              'surf_atlas_dir',
                                              'gray_ordinates_dir',
                                              'gray_ordinates_res',
                                              'high_res_mesh',
                                              'low_res_mesh',
                                              'subcortical_gray_labels',
                                              'freesurfer_labels',
                                              'fmri_res',
                                              'smooth_fwhm'],
                                 output_names=['dtseries',
                                               'desikan_killiany_164',
                                               'destrieux_164',
                                               'desikan_killiany_32',
                                               'destrieux_32'],
                                 function=run_surface),
                   name=f'post_freesurfer_{pipe_num}')

    surf.inputs.subject = cfg['subject_id']

    surf.inputs.post_freesurfer_folder = os.path.join(cfg.pipeline_setup['working_directory']['path'],
        'cpac_'+cfg['subject_id'],
        f'post_freesurfer_{pipe_num}')

    surf.inputs.surf_atlas_dir = cfg.surface_analysis['post_freesurfer']['surf_atlas_dir']
    surf.inputs.gray_ordinates_dir = cfg.surface_analysis['post_freesurfer']['gray_ordinates_dir']
    surf.inputs.subcortical_gray_labels = cfg.surface_analysis['post_freesurfer']['subcortical_gray_labels']
    surf.inputs.freesurfer_labels = cfg.surface_analysis['post_freesurfer']['freesurfer_labels']

    # convert integers to strings as subprocess requires string inputs
    surf.inputs.gray_ordinates_res = str(cfg.surface_analysis['post_freesurfer']['gray_ordinates_res'])
    surf.inputs.high_res_mesh = str(cfg.surface_analysis['post_freesurfer']['high_res_mesh'])
    surf.inputs.low_res_mesh = str(cfg.surface_analysis['post_freesurfer']['low_res_mesh'])
    surf.inputs.fmri_res = str(cfg.surface_analysis['post_freesurfer']['fmri_res'])
    surf.inputs.smooth_fwhm = str(cfg.surface_analysis['post_freesurfer']['smooth_fwhm'])

    node, out = strat_pool.get_data('freesurfer-subject-dir')
    wf.connect(node, out, surf, 'freesurfer_folder')

    node, out = strat_pool.get_data('desc-restore_T1w')
    wf.connect(node, out, surf, 't1w_restore_image')

    node, out = strat_pool.get_data('space-template_desc-head_T1w')
    wf.connect(node, out, surf, 'atlas_space_t1w_image')

    node, out = strat_pool.get_data('from-T1w_to-template_mode-image_xfm')
    wf.connect(node, out, surf, 'atlas_transform')

    node, out = strat_pool.get_data('from-template_to-T1w_mode-image_xfm')
    wf.connect(node, out, surf, 'inverse_atlas_transform')

    node, out = strat_pool.get_data('space-template_desc-brain_bold')
    wf.connect(node, out, surf, 'atlas_space_bold')

    node, out = strat_pool.get_data('space-template_desc-scout_bold')
    wf.connect(node, out, surf, 'scout_bold')

    outputs = {
        'atlas-DesikanKilliany_space-fsLR_den-32k_dlabel': (surf,
                                                            'desikan_'
                                                            'killiany_32'),
        'atlas-Destrieux_space-fsLR_den-32k_dlabel': (surf, 'destrieux_32'),
        'atlas-DesikanKilliany_space-fsLR_den-164k_dlabel': (surf,
                                                             'desikan_'
                                                             'killiany_164'),
        'atlas-Destrieux_space-fsLR_den-164k_dlabel': (surf, 'destrieux_164'),
        'space-fsLR_den-32k_bold-dtseries': (surf, 'dtseries')
    }

    return wf, outputs


def surface_preproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_preproc",
     "config": ["surface_analysis", "post_freesurfer"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["freesurfer-subject-dir",
                "desc-restore_T1w",
                "space-template_desc-head_T1w",
                "from-T1w_to-template_mode-image_xfm",
                "from-template_to-T1w_mode-image_xfm",
                "space-template_desc-brain_bold",
                "space-template_desc-scout_bold",
                "wmparc.mgz"],
     "outputs": ["atlas-DesikanKilliany_space-fsLR_den-32k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-32k_dlabel",
                 "atlas-DesikanKilliany_space-fsLR_den-164k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-164k_dlabel",
                 "space-fsLR_den-32k_bold-dtseries"]}
    '''
    
    
    wf, outputs = surface_connector(wf, cfg, strat_pool, pipe_num, opt)
    #raise Exception("Entered the node")
    wf, outputs = create_resmat(wf, cfg, strat_pool, pipe_num, opt=None)
    return (wf, outputs)
