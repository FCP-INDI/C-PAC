from logging import raiseExceptions
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

    # # DCAN-HCP PostFreeSurfer
    # # Ref: https://github.com/DCAN-Labs/DCAN-HCP/blob/master/PostFreeSurfer/PostFreeSurferPipeline.sh
    # cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/run.sh', '--post_freesurfer_folder', post_freesurfer_folder, \
    #     '--freesurfer_folder', freesurfer_folder, '--subject', subject, \
    #     '--t1w_restore', t1w_restore_image, '--atlas_t1w', atlas_space_t1w_image, \
    #     '--atlas_transform', atlas_transform, '--inverse_atlas_transform', inverse_atlas_transform, \
    #     '--surfatlasdir', surf_atlas_dir, '--grayordinatesdir', gray_ordinates_dir, '--grayordinatesres', gray_ordinates_res, \
    #     '--hiresmesh', high_res_mesh, '--lowresmesh', low_res_mesh, \
    #     '--subcortgraylabels', subcortical_gray_labels, '--freesurferlabels', freesurfer_labels]

    # subprocess.check_output(cmd)

        # DCAN-HCP PostFreeSurfer Block1
   
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block1.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]


    subprocess.check_output(cmd)

       # DCAN-HCP PostFreeSurfer Block2

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block2.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    subprocess.check_output(cmd)

      # DCAN-HCP PostFreeSurfer Block3

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block3.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    subprocess.check_output(cmd)

    #DCAN-HCP PostFreeSurfer Block3.5
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block3-5.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    subprocess.check_output(cmd)

    # DCAN-HCP PostFreeSurfer Block4

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block4.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    subprocess.check_output(cmd)

  # DCAN-HCP PostFreeSurfer Block5

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block5.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    subprocess.check_output(cmd)

   # DCAN-HCP PostFreeSurfer Block6
      
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block6.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

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
        164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
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
    
    subcortical_atlas = os.path.join(post_freesurfer_folder,
                            'MNINonLinear/Results/task-rest01/'
                            'task-rest01_AtlasSubcortical_s2.nii.gz')

    good_voxels = os.path.join(post_freesurfer_folder,
                            'MNINonLinear/Results/task-rest01/RibbonVolumeToSurfaceMapping/'
                            'goodvoxels.nii.gz')

    ribbon_only = os.path.join(post_freesurfer_folder,
                            'MNINonLinear/Results/task-rest01/RibbonVolumeToSurfaceMapping/'
                            'ribbon_only.nii.gz')

    atlas_roi = {
        'func': {
            'L': os.path.join(post_freesurfer_folder, 
                            'MNINonLinear/Results/task-rest01/',
                            'task-rest01.L.atlasroi.32k_fs_LR.func.gii'),

            'R': os.path.join(post_freesurfer_folder, 
                            'MNINonLinear/Results/task-rest01/',
                            'task-rest01.R.atlasroi.32k_fs_LR.func.gii')},
        'shape': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                     'fsaverage_LR32k',
                            f'{subject}.L.atlasroi.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                     'fsaverage_LR32k',
                            f'{subject}.R.atlasroi.32k_fs_LR.shape.gii')}}
    task_rest = {
        'L': os.path.join(post_freesurfer_folder, 
                            'MNINonLinear/Results/task-rest01/',
                            'task-rest01.L.native.func.gii'),
        'R': os.path.join(post_freesurfer_folder, 
                            'MNINonLinear/Results/task-rest01/',
                            'task-rest01.R.native.func.gii')}
    spec = {
        'LR_32k': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                            f'{subject}.32k_fs_LR.wb.spec'),
        'native': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.native.wb.spec')}
    
    areal_distortion = {
        'FS': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.L.ArealDistortion_FS.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.R.ArealDistortion_FS.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.ArealDistortion_FS.32k_fs_LR.dscalar.nii')},                   
                        
        'MSMSulc': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.L.ArealDistortion_MSMSulc.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.R.ArealDistortion_MSMSulc.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.ArealDistortion_MSMSulc.32k_fs_LR.dscalar.nii')}}
    edge_distortion = {
        'FS': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.L.EdgeDistortion_FS.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.R.EdgeDistortion_FS.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.EdgeDistortion_FS.32k_fs_LR.dscalar.nii')},                   
                        
        'MSMSulc': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.L.EdgeDistortion_MSMSulc.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.R.EdgeDistortion_MSMSulc.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                        'fsaverage_LR32k',
                          f'{subject}.EdgeDistortion_MSMSulc.32k_fs_LR.dscalar.nii')}}
    curvature = {
        'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.curvature.32k_fs_LR.shape.gii'),
        'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.curvature.32k_fs_LR.shape.gii'),
        'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.curvature.32k_fs_LR.dscalar.nii')}

    flat = {
        'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.flat.32k_fs_LR.surf.gii'),
        'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.flat.32k_fs_LR.surf.gii')}  

    inflate_32k = {
        'inflated': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.inflated.32k_fs_LR.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.inflated.32k_fs_LR.surf.gii')},
        'very_inflated': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.very_inflated.32k_fs_LR.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.very_inflated.32k_fs_LR.surf.gii')}}
    inflate = {
        'inflated': {
            'L': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.L.inflated.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.R.inflated.native.surf.gii')},
        'very_inflated': {
            'L': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.L.very_inflated.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.R.very_inflated.native.surf.gii')}}
    midthickness = {
        'L': {
            164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.L.midthickness.164k_fs_LR.surf.gii'),
            32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.midthickness.32k_fs_LR.surf.gii'),
            'native': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.L.midthickness.native.surf.gii')},
        'R': {
            164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.R.midthickness.164k_fs_LR.surf.gii'),
            32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.midthickness.32k_fs_LR.surf.gii'),
            'native': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.R.midthickness.native.surf.gii')}}
    pial = {
        'L': {
            32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.pial.32k_fs_LR.surf.gii'),
            'native': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.L.pial.native.surf.gii')},
        'R': {
            32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.pial.32k_fs_LR.surf.gii'), 
            'native': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.R.pial.native.surf.gii')}}
    sphere = {
        '32k_fs_LR': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.sphere.32k_fs_LR.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.sphere.32k_fs_LR.surf.gii')},
        'MSMSulc': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.L.sphere.MSMSulc.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.R.sphere.MSMSulc.native.surf.gii')},
        'native': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.L.sphere.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.R.sphere.native.surf.gii')},
        'reg': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.L.sphere.reg.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.R.sphere.reg.native.surf.gii')},
        'reg_reg_LR': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.L.sphere.reg.reg_LR.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.R.sphere.reg.reg_LR.native.surf.gii')},
        'rot': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.L.sphere.rot.native.surf.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear/Native/',
                          f'{subject}.R.sphere.rot.native.surf.gii')}}
    StrainJ = {
        'FS': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.StrainJ_FS.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.StrainJ_FS.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.StrainJ_FS.32k_fs_LR.dscalar.nii')},
        'MSMSulc': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.StrainJ_MSMSulc.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.StrainJ_MSMSulc.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.StrainJ_MSMSulc.32k_fs_LR.dscalar.nii')}}
    StrainR = {
        'FS': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.StrainR_FS.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.StrainR_FS.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.StrainR_FS.32k_fs_LR.dscalar.nii')},
        'MSMSulc': {
            'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.StrainR_MSMSulc.32k_fs_LR.shape.gii'),
            'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.StrainR_MSMSulc.32k_fs_LR.shape.gii'),
            'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.StrainR_MSMSulc.32k_fs_LR.dscalar.nii')}}
    
    sulc = {
        'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.sulc.32k_fs_LR.shape.gii'),
        'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.sulc.32k_fs_LR.shape.gii'),
        'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.sulc.32k_fs_LR.dscalar.nii')} 
    thickness = {
        'L': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.thickness.32k_fs_LR.shape.gii'),
        'R': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.thickness.32k_fs_LR.shape.gii'),
        'dscalar': os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.thickness.32k_fs_LR.dscalar.nii')}
    white = {
        'L': {
            164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.L.white.164k_fs_LR.surf.gii'),
            32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.L.white.32k_fs_LR.surf.gii'),
            'native': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.L.white.native.surf.gii')},
        'R': {
            164: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                          f'{subject}.R.white.164k_fs_LR.surf.gii'),
            32: os.path.join(post_freesurfer_folder, 'MNINonLinear',
                                                    'fsaverage_LR32k',
                          f'{subject}.R.white.32k_fs_LR.surf.gii'),
            'native': os.path.join(post_freesurfer_folder, 'T1w/Native/',
                          f'{subject}.R.white.native.surf.gii')}}

    return (dtseries, aparc['desikan_killiany'][164], aparc['destrieux'][164], aparc['desikan_killiany'][32], 
            aparc['destrieux'][32], subcortical_atlas, good_voxels, ribbon_only, atlas_roi['func']['L'], atlas_roi['func']['R'], 
            atlas_roi['shape']['L'], atlas_roi['shape']['R'], task_rest['L'], task_rest['R'], spec['LR_32k'], spec['native'],
            areal_distortion['FS']['L'], areal_distortion['FS']['R'], areal_distortion['FS']['dscalar'], areal_distortion['MSMSulc']['L'], 
            areal_distortion['MSMSulc']['R'], areal_distortion['MSMSulc']['dscalar'], edge_distortion['FS']['L'], 
            edge_distortion['FS']['R'], edge_distortion['FS']['dscalar'], edge_distortion['MSMSulc']['L'], edge_distortion['MSMSulc']['R'], 
            edge_distortion['MSMSulc']['dscalar'], curvature['L'], curvature['R'], curvature['dscalar'], flat['L'], flat['R'], 
            inflate_32k['inflated']['L'], inflate_32k['inflated']['R'], inflate_32k['very_inflated']['L'], inflate_32k['very_inflated']['R'], 
            inflate['inflated']['L'], inflate['inflated']['R'], inflate['very_inflated']['L'], inflate['very_inflated']['R'], 
            midthickness['L'][164], midthickness['L'][32], midthickness['L']['native'], midthickness['R'][164], midthickness['R'][32], 
            midthickness['R']['native'], pial['L'][32], pial['L']['native'], pial['R'][32], pial['R']['native'], sphere['32k_fs_LR']['L'], 
            sphere['32k_fs_LR']['R'], sphere['MSMSulc']['L'], sphere['MSMSulc']['R'], sphere['native']['L'], sphere['native']['R'], 
            sphere['reg']['L'], sphere['reg']['R'], sphere['reg_reg_LR']['L'], sphere['reg_reg_LR']['R'], sphere['rot']['L'], sphere['rot']['R'],
            StrainJ['FS']['L'], StrainJ['FS']['R'], StrainJ['FS']['dscalar'], StrainJ['MSMSulc']['L'], StrainJ['MSMSulc']['R'], 
            StrainJ['MSMSulc']['dscalar'], StrainR['FS']['L'], StrainR['FS']['R'], StrainR['FS']['dscalar'], StrainR['MSMSulc']['L'], 
            StrainR['MSMSulc']['R'], StrainR['MSMSulc']['dscalar'], sulc['L'], sulc['R'], sulc['dscalar'], thickness['L'], thickness['R'], 
            thickness['dscalar'], white['L'][164], white['L'][32], white['L']['native'], white['R'][164], white['R'][32], white['R']['native'])

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
                                               'destrieux_32',
                                               'subcortical_atlas',
                                               'good_voxels',
                                               'ribbon_only',
                                               'atlas_roi_func_L', 
                                               'atlas_roi_func_R', 
                                               'atlas_roi_shape_L', 
                                               'atlas_roi_shape_R', 
                                               'native_L', 
                                               'native_R', 
                                               'spec_LR_32k',
                                               'spec_native', 
                                               'areal_distortion_FS_L', 
                                               'areal_distortion_FS_R', 
                                               'areal_distortion_FS_dscalar', 
                                               'areal_distortion_MSMSulc_L', 
                                               'areal_distortion_MSMSulc_R', 
                                               'areal_distortion_MSMSulc_dscalar',
                                               'edge_distortion_FS_L',
                                               'edge_distortion_FS_R', 
                                               'edge_distortion_FS_dscalar', 
                                               'edge_distortion_MSMSulc_L', 
                                               'edge_distortion_MSMSulc_R', 
                                               'edge_distortion_MSMSulc_dscalar',
                                               'curvature_L', 
                                               'curvature_R', 
                                               'curvature_dscalar', 
                                               'flat_L', 
                                               'flat_R',
                                               '32k_inflated_L', 
                                               '32k_inflated_R', 
                                               '32k_very_inflated_L',
                                               '32k_very_inflated_R', 
                                               'inflated_L', 
                                               'inflated_R', 
                                               'very_inflated_L', 
                                               'very_inflated_R',
                                               'midthickness_L_164',
                                               'midthickness_L_32', 
                                               'midthickness_L_native', 
                                               'midthickness_R_164', 
                                               'midthickness_R_32',
                                               'midthickness_R_native',
                                               'pial_L_32', 
                                               'pial_L_native', 
                                               'pial_R_32', 
                                               'pial_R_native',
                                               'sphere_32k_fs_LR_L', 
                                               'sphere_32k_fs_LR_R', 
                                               'sphere_MSMSulc_L', 
                                               'sphere_MSMSulc_R', 
                                               'sphere_native_L', 
                                               'sphere_native_R', 
                                               'sphere_reg_L', 
                                               'sphere_reg_R', 
                                               'sphere_reg_reg_LR_L', 
                                               'sphere_reg_reg_LR_R', 
                                               'sphere_rot_L', 
                                               'sphere_rot_R',
                                               'StrainJ_FS_L',
                                               'StrainJ_FS_R', 
                                               'StrainJ_FS_dscalar', 
                                               'StrainJ_MSMSulc_L', 
                                               'StrainJ_MSMSulc_R', 
                                               'StrainJ_MSMSulc_dscalar', 
                                               'StrainR_FS_L', 
                                               'StrainR_FS_R', 
                                               'StrainR_FS_dscalar', 
                                               'StrainR_MSMSulc_L', 
                                               'StrainR_MSMSulc_R', 
                                               'StrainR_MSMSulc_dscalar', 
                                               'sulc_L', 
                                               'sulc_R', 
                                               'sulc_dscalar',
                                               'thickness_L', 
                                               'thickness_R', 
                                               'thickness_dscalar', 
                                               'white_L_164', 
                                               'white_L_32', 
                                               'white_L_native', 
                                               'white_R_164', 
                                               'white_R_32', 
                                               'white_R_native'
                                            ],
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

    restore =           ["desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w",
                         "space-longitudinal_desc-reorient_T1w"]
    space_temp =        ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask",]
    atlas_xfm =         ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"]
    atlas_xfm_inv =     ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"]
    atlas_space_bold =  ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"]
    scout_bold =        ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold",
                        "space-template_desc-preproc_bold", "space-template_desc-motion_bold", "space-template_bold"]


    node, out = strat_pool.get_data('freesurfer-subject-dir')
    wf.connect(node, out, surf, 'freesurfer_folder')

    node, out = strat_pool.get_data(restore) 
    wf.connect(node, out, surf, 't1w_restore_image')

    node, out = strat_pool.get_data(space_temp) 
    wf.connect(node, out, surf, 'atlas_space_t1w_image')

    node, out = strat_pool.get_data(atlas_xfm) 
    wf.connect(node, out, surf, 'atlas_transform')

    node, out = strat_pool.get_data(atlas_xfm_inv) 
    wf.connect(node, out, surf, 'inverse_atlas_transform')

    node, out = strat_pool.get_data(atlas_space_bold) 
    wf.connect(node, out, surf, 'atlas_space_bold')

    node, out = strat_pool.get_data(scout_bold)
    wf.connect(node, out, surf, 'scout_bold')

    outputs = {
        'space-fsLR_den-32k_bold-dtseries': (surf, 'dtseries'),
        'atlas-DesikanKilliany_space-fsLR_den-32k_dlabel': (surf,
                                                            'desikan_'
                                                            'killiany_32'),
        'atlas-Destrieux_space-fsLR_den-32k_dlabel': (surf, 'destrieux_32'),
        'atlas-DesikanKilliany_space-fsLR_den-164k_dlabel': (surf,
                                                             'desikan_'
                                                             'killiany_164'),
        'atlas-Destrieux_space-fsLR_den-164k_dlabel': (surf, 'destrieux_164'),
        'AtlasSubcortical_s2': (surf, 'subcortical_atlas'), 
        'goodvoxels': (surf, 'good_voxels'),
        'ribbon_only': (surf, 'ribbon_only'),
        'atlasroi_hemi-L_space-fsLR_den-32k_func': (surf, 
                                    'atlas_roi_func_'
                                    'L'),
        'atlasroi_hemi-R_space-fsLR_den-32k_func': (surf, 
                                    'atlas_roi_func_'
                                    'R'),
        'atlasroi_hemi-L_space-fsLR_den-32k_shape': (surf, 
                                    'atlas_roi_shape_'
                                    'L'),
        'atlasroi_hemi-R_space-fsLR_den-32k_shape': (surf, 
                                    'atlas_roi_shape_'
                                    'R'),                            
        'space-native_hemi-L_func': (surf, 'native_' 'L'),
        'space-native_hemi-R_func': (surf, 'native_' 'R'),
        'space-fsLR_den-32k_wb-spec': (surf, 'spec_' 'LR_32k'), 
        'space-native_wb-spec': (surf, 'spec_' 'native'),
        'arealdistortion-FS_hemi-L_space-fsLR_den-32k_shape': (surf, 'areal_distortion_'
                                                                        'FS_'
                                                                        'L'),
        'arealdistortion-FS_hemi-R_space-fsLR_den-32k_shape': (surf, 'areal_distortion_'
                                                                        'FS_'
                                                                        'R'),
        'arealdistortion-FS_space-fsLR_den-32k_dscalar': (surf, 'areal_distortion_'
                                                                        'FS_'
                                                                        'dscalar'),
        'arealdistortion-MSMSulc_hemi-L_space-fsLR_den-32k_shape': (surf, 'areal_distortion_'
                                                                        'MSMSulc_'
                                                                        'L'),
        'arealdistortion-MSMSulc_hemi-R_space-fsLR_den-32k_shape': (surf, 'areal_distortion_'
                                                                        'MSMSulc_'
                                                                        'R'),
        'arealdistortion-MSMSulc_space-fsLR_den-32k_dscalar': (surf, 'areal_distortion_'
                                                                        'MSMSulc_'
                                                                        'dscalar'),
        'edgedistortion-FS_hemi-L_space-fsLR_den-32k_shape': (surf, 'edge_distortion_'
                                                                        'FS_'
                                                                        'L'),
        'edgedistortion-FS_hemi-R_space-fsLR_den-32k_shape': (surf, 'edge_distortion_'
                                                                        'FS_'
                                                                        'R'),
        'edgedistortion-FS_space-fsLR_den-32k_dscalar': (surf, 'edge_distortion_'
                                                                        'FS_'
                                                                        'dscalar'),
        'edgedistortion-MSMSulc_hemi-L_space-fsLR_den-32k_shape': (surf, 'edge_distortion_'
                                                                        'MSMSulc_'
                                                                        'L'),
        'edgedistortion-MSMSulc_hemi-R_space-fsLR_den-32k_shape': (surf, 'edge_distortion_'
                                                                        'MSMSulc_'
                                                                        'R'),
        'edgedistortion-MSMSulc_space-fsLR_den-32k_dscalar': (surf, 'edge_distortion_'
                                                                        'MSMSulc_'
                                                                        'dscalar'),
        'hemi-L_space-fsLR_den-32k_curv_shape': (surf, 'curvature_' 'L'),
        'hemi-R_space-fsLR_den-32k_curv_shape': (surf, 'curvature_' 'R'),
        'space-fsLR_den-32k_curv_dscalar': (surf, 'curvature_' 'dscalar'),
        'hemi-L_space-fsLR_den-32k_flat_surf': (surf, 'flat_L'),
        'hemi-R_space-fsLR_den-32k_flat_surf': (surf, 'flat_R'),
        'hemi-L_space-fsLR_den-32k_inflated_surf': (surf, '32k_inflated_L'),
        'hemi-R_space-fsLR_den-32k_inflated_surf': (surf, '32k_inflated_R'),
        'hemi-L_space-fsLR_den-32k_very-inflated_surf': (surf, '32k_very_inflated_L'),
        'hemi-R_space-fsLR_den-32k_very-inflated_surf': (surf, '32k_very_inflated_R'),
        'hemi-L_space-native_inflated_surf': (surf, 'inflated_L'),
        'hemi-R_space-native_inflated_surf': (surf, 'inflated_R'),
        'hemi-L_space-native_very-inflated_surf': (surf, 'very_inflated_L'),
        'hemi-R_space-native_very-inflated_surf': (surf, 'very_inflated_R'),
        'hemi-L_space-fsLR_den-164k_midthickness_surf': (surf, 'midthickness_L_164'),
        'hemi-L_space-fsLR_den-32k_midthickness_surf': (surf, 'midthickness_L_32'),
        'hemi-L_space-native_midthickness_surf': (surf, 'midthickness_L_native'),
        'hemi-R_space-fsLR_den-164k_midthickness_surf': (surf, 'midthickness_R_164'),
        'hemi-R_space-fsLR_den-32k_midthickness_surf': (surf, 'midthickness_R_32'),
        'hemi-R_space-native_midthickness_surf': (surf, 'midthickness_R_native'),
        'hemi-L_space-fsLR_den-32k_pial_surf': (surf, 'pial_L_32'),
        'hemi-L_space-native_den-32k_pial_surf': (surf, 'pial_L_native'),
        'hemi-R_space-fsLR_den-32k_pial_surf': (surf, 'pial_R_32'),
        'hemi-R_space-native_den-32k_pial_surf': (surf, 'pial_R_native'),
        'hemi-L_space-fsLR_den-32k_sphere_surf': (surf, 'sphere_32k_fs_LR_L'),
        'hemi-R_space-fsLR_den-32k_sphere_surf': (surf, 'sphere_32k_fs_LR_R'),                            
        'hemi-L_MSMSulc_space-native_sphere_surf': (surf, 'sphere_MSMSulc_L'),
        'hemi-R_MSMSulc_space-native_sphere_surf': (surf, 'sphere_MSMSulc_R'),
        'hemi-L_space-native_sphere_surf': (surf, 'sphere_native_L'),
        'hemi-R_space-native_sphere_surf': (surf, 'sphere_native_R'),
        'hemi-L_space-native_sphere-reg_surf': (surf, 'sphere_reg_L'),
        'hemi-R_space-native_sphere-reg_surf': (surf, 'sphere_reg_R'),
        'hemi-L_space-native_sphere-reg-reg_surf': (surf, 'sphere_reg_reg_LR_L'),
        'hemi-R_space-native_sphere-reg-reg_surf': (surf, 'sphere_reg_reg_LR_R'),
        'hemi-L_space-native_sphere-rot_surf': (surf, 'sphere_rot_L'),
        'hemi-R_space-native_sphere-rot_surf': (surf, 'sphere_rot_R'),
        'hemi-L_strainJ-FS_space-fsLR_den-32k_shape': (surf, 'StrainJ_FS_L'),
        'hemi-R_strainJ-FS_space-fsLR_den-32k_shape': (surf, 'StrainJ_FS_R'),
        'strainJ-FS_space-fsLR_den-32k_dscalar': (surf, 'StrainJ_FS_dscalar'),
        'hemi-L_strainJ-MSMSulc_space-fsLR_den-32k_shape': (surf, 'StrainJ_MSMSulc_L'),
        'hemi-R_strainJ-MSMSulc_space-fsLR_den-32k_shape': (surf, 'StrainJ_MSMSulc_R'),
        'strainJ-MSMSulc_space-fsLR_den-32k_dscalar': (surf, 'StrainJ_MSMSulc_dscalar'),
        'hemi-L_strainR-FS_space-fsLR_den-32k_shape': (surf, 'StrainR_FS_L'),
        'hemi-R_strainR-FS_space-fsLR_den-32k_shape': (surf, 'StrainR_FS_R'),
        'strainR-FS_space-fsLR_den-32k_dscalar': (surf, 'StrainR_FS_dscalar'),
        'hemi-L_strainR-MSMSulc_space-fsLR_den-32k_shape': (surf, 'StrainR_MSMSulc_L'),
        'hemi-R_strainR-MSMSulc_space-fsLR_den-32k_shape': (surf, 'StrainR_MSMSulc_R'),
        'strainR-MSMSulc_space-fsLR_den-32k_dscalar': (surf, 'StrainR_MSMSulc_dscalar'),
        'hemi-L_space-fsLR_den-32k_sulc_shape': (surf, 'sulc_L'),
        'hemi-R_space-fsLR_den-32k_sulc_shape': (surf, 'sulc_R'),
        'space-fsLR_den-32k_sulc_dscalar': (surf, 'sulc_dscalar'),                                   
        'hemi-L_space-fsLR_den-32k_thickness_shape': (surf, 'thickness_L'),
        'hemi-R_space-fsLR_den-32k_thickness_shape': (surf, 'thickness_R'),
        'space-fsLR_den-32k_thickness_dscalar': (surf, 'thickness_dscalar'),                                   
        'hemi-L_space-fsLR_den-164k_white_surf': (surf, 'white_L_164'),
        'hemi-L_space-fsLR_den-32k_white_surf': (surf, 'white_L_32'),
        'hemi-L_space-native_white_surf': (surf, 'white_L_native'),
        'hemi-R_space-fsLR_den-164k_white_surf': (surf, 'white_R_164'),
        'hemi-R_space-fsLR_den-32k_white_surf': (surf, 'white_R_32'),
        'hemi-R_space-native_white_surf': (surf, 'white_R_native')
    }

    return wf, outputs

def surface_postproc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_postproc",
     "config": ["surface_analysis", "post_freesurfer"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["freesurfer-subject-dir",
                ["desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w", 
                "space-longitudinal_desc-reorient_T1w"],
                ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask"],
                ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"],
                ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"],
                ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"],
                ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold", 
                "space-template_desc-preproc_bold", "space-template_desc-motion_bold", "space-template_bold"]],
     "outputs": ["space-fsLR_den-32k_bold-dtseries",
                 "atlas-DesikanKilliany_space-fsLR_den-32k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-32k_dlabel",
                 "atlas-DesikanKilliany_space-fsLR_den-164k_dlabel",
                 "atlas-Destrieux_space-fsLR_den-164k_dlabel",
                 "AtlasSubcortical_s2",
                 "goodvoxels",
                 "ribbon_only",
                 "atlasroi_hemi-L_space-fsLR_den-32k_func",
                 "atlasroi_hemi-R_space-fsLR_den-32k_func",
                 "atlasroi_hemi-L_space-fsLR_den-32k_shape",
                 "atlasroi_hemi-R_space-fsLR_den-32k_shape",
                 "space-native_hemi-L_func",
                 "space-native_hemi-R_func",
                 "space-fsLR_den-32k_wb-spec",
                 "space-native_wb-spec",
                 "arealdistortion-FS_hemi-L_space-fsLR_den-32k_shape",
                 "arealdistortion-FS_hemi-R_space-fsLR_den-32k_shape",
                 "arealdistortion-FS_space-fsLR_den-32k_dscalar",
                 "arealdistortion-MSMSulc_hemi-L_space-fsLR_den-32k_shape",
                 "arealdistortion-MSMSulc_hemi-R_space-fsLR_den-32k_shape",
                 "arealdistortion-MSMSulc_space-fsLR_den-32k_dscalar",
                 "edgedistortion-FS_hemi-L_space-fsLR_den-32k_shape",
                 "edgedistortion-FS_hemi-R_space-fsLR_den-32k_shape",
                 "edgedistortion-FS_space-fsLR_den-32k_dscalar",
                 "edgedistortion-MSMSulc_hemi-L_space-fsLR_den-32k_shape",
                 "edgedistortion-MSMSulc_hemi-R_space-fsLR_den-32k_shape",
                 "edgedistortion-MSMSulc_space-fsLR_den-32k_dscalar",
                 "hemi-L_space-fsLR_den-32k_curv_shape",
                 "hemi-R_space-fsLR_den-32k_curv_shape",
                 "space-fsLR_den-32k_curv_dscalar",
                 "hemi-L_space-fsLR_den-32k_flat_surf",
                 "hemi-R_space-fsLR_den-32k_flat_surf",
                 "hemi-L_space-fsLR_den-32k_inflated_surf",
                 "hemi-R_space-fsLR_den-32k_inflated_surf",
                 "hemi-L_space-fsLR_den-32k_very-inflated_surf",
                 "hemi-R_space-fsLR_den-32k_very-inflated_surf",
                 "hemi-L_space-native_inflated_surf",
                 "hemi-R_space-native_inflated_surf",
                 "hemi-L_space-native_very-inflated_surf",
                 "hemi-R_space-native_very-inflated_surf",
                 "hemi-L_space-fsLR_den-164k_midthickness_surf",
                 "hemi-R_space-fsLR_den-164k_midthickness_surf",
                 "hemi-L_space-fsLR_den-32k_midthickness_surf",
                 "hemi-R_space-fsLR_den-32k_midthickness_surf",
                 "hemi-L_space-native_midthickness_surf",
                 "hemi-R_space-native_midthickness_surf",
                 "hemi-L_space-fsLR_den-32k_pial_surf",
                 "hemi-R_space-fsLR_den-32k_pial_surf",
                 "hemi-L_space-native_den-32k_pial_surf",
                 "hemi-R_space-native_den-32k_pial_surf",
                 "hemi-L_space-fsLR_den-32k_sphere_surf",
                 "hemi-R_space-fsLR_den-32k_sphere_surf",
                 "hemi-L_MSMSulc_space-native_sphere_surf",
                 "hemi-R_MSMSulc_space-native_sphere_surf",
                 "hemi-L_space-native_sphere_surf",
                 "hemi-R_space-native_sphere_surf",
                 "hemi-L_space-native_sphere-reg_surf",
                 "hemi-R_space-native_sphere-reg_surf",
                 "hemi-L_space-native_sphere-reg-reg_surf",
                 "hemi-R_space-native_sphere-reg-reg_surf",
                 "hemi-L_space-native_sphere-rot_surf",
                 "hemi-R_space-native_sphere-rot_surf",
                 "hemi-L_strainJ-FS_space-fsLR_den-32k_shape",
                 "hemi-R_strainJ-FS_space-fsLR_den-32k_shape",
                 "strainJ-FS_space-fsLR_den-32k_dscalar",
                 "hemi-L_strainJ-MSMSulc_space-fsLR_den-32k_shape",
                 "hemi-R_strainJ-MSMSulc_space-fsLR_den-32k_shape",
                 "strainJ-MSMSulc_space-fsLR_den-32k_dscalar",
                 "hemi-L_strainR-FS_space-fsLR_den-32k_shape",
                 "hemi-R_strainR-FS_space-fsLR_den-32k_shape",
                 "strainR-FS_space-fsLR_den-32k_dscalar",
                 "hemi-L_strainR-MSMSulc_space-fsLR_den-32k_shape",
                 "hemi-R_strainR-MSMSulc_space-fsLR_den-32k_shape",
                 "strainR-MSMSulc_space-fsLR_den-32k_dscalar",
                 "hemi-L_space-fsLR_den-32k_sulc_shape",
                 "hemi-R_space-fsLR_den-32k_sulc_shape",
                 "space-fsLR_den-32k_sulc_dscalar",
                 "hemi-L_space-fsLR_den-32k_thickness_shape",
                 "hemi-R_space-fsLR_den-32k_thickness_shape",
                 "space-fsLR_den-32k_thickness_dscalar",
                 "hemi-L_space-fsLR_den-164k_white_surf",
                 "hemi-R_space-fsLR_den-164k_white_surf",
                 "hemi-L_space-fsLR_den-32k_white_surf",
                 "hemi-R_space-fsLR_den-32k_white_surf",
                 "hemi-L_space-native_white_surf",
                 "hemi-R_space-native_white_surf"]}
    '''
    wf, outputs = surface_connector(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)

def cal_surface_falff(wf, cfg, strat_pool, pipe_num, opt):

    falff = pe.Node(util.Function(input_names=['subject','dtseries'], 
                                output_names=['surf_falff'],
                                function=run_surf_falff),
                   name=f'surf_falff_{pipe_num}')
    
    falff.inputs.subject = cfg['subject_id']
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold-dtseries') 
    wf.connect(node, out, falff, 'dtseries')
    
    outputs = {
        'surf_falff': (falff,'surf_falff')}
    return wf, outputs

def cal_surface_alff(wf, cfg, strat_pool, pipe_num, opt):

    alff = pe.Node(util.Function(input_names=['subject','dtseries'], 
                             output_names=['surf_alff'],
                               function=run_surf_alff),
                 name=f'surf_alff_{pipe_num}')
    
    alff.inputs.subject = cfg['subject_id']
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold-dtseries') 
    wf.connect(node, out,alff, 'dtseries')
    outputs = {
        'surf_alff': (alff,'surf_alff')}
    return wf, outputs

# def cal_reho(wf, cfg, strat_pool, pipe_num, opt):

#     L_cortex_file = pe.Node(util.Function(input_names=['dtseries', 'structure', 'cortex_filename'], 
#                                 output_names=['L_cortex_file'],
#                                 function=run_get_cortex),
#                   name=f'L_surf_cortex_{pipe_num}')

#     L_cortex_file.inputs.structure = "LEFT"
#     L_cortex_file.inputs.cortex_filename = "L_cortex.func.gii"
#     node, out = strat_pool.get_data(space-fsLR_den-32k_bold-dtseries) 
#     wf.connect(node, out, L_cortex_file, 'dtseries')

#     R_cortex_file = pe.Node(util.Function(input_names=['dtseries', 'structure', 'cortex_filename'], 
#                                 output_names=['R_cortex_file'],
#                                 function=run_get_cortex),
#                   name=f'R_surf_cortex_{pipe_num}')


#     R_cortex_file.inputs.structure = "RIGHT"
#     R_cortex_file.inputs.cortex_filename = "R_cortex.func.gii"
#     wf.connect(node, out,R_cortex_file, 'dtseries')


#     mean_timeseries = pe.Node(util.Function(input_names=['dtseries'], 
#                                 output_names=['mean_timeseries'],
#                                 function=run_mean_timeseries),
#                   name=f'mean_timeseries_{pipe_num}')


#     wf.connect(node, out, mean_timeseries, 'dtseries')

#     L_reho = pe.Node(util.Function(input_names=['dtseries', 'mask' , 'cortex_file', 'surface_file', 'mean_timeseries','reho_filename'], 
#                                 output_names=['L_reho'],
#                                 function=run_surf_reho),
#                     name=f'surf_reho{pipe_num}')

#     wf.connect(L_cortex_file, 'L_cortex_file', L_reho, 'cortex_file')
#     wf.connect(surf, 'midthickness_L_32', L_reho, 'surface_file')
#     wf.connect(surf,'atlas_roi_shape_L', L_reho, 'mask')
#     wf.connect(mean_timeseries, 'mean_timeseries', L_reho, 'mean_timeseries') 
#     L_reho.inputs.reho_filename = L_surf_reho.dscalar.nii
#     wf.connect(node, out, L_reho, 'dtseries')

#     R_reho = pe.Node(util.Function(input_names=['dtseries', 'mask' , 'cortex_file', 'surface_file', 'mean_timeseries', 'reho_filename'], 
#                                 output_names=['R_reho'],
#                                 function=run_surf_reho),
#                     name=f'surf_reho{pipe_num}')

#     wf.connect(R_cortex_file, 'R_cortex_file', R_reho, 'cortex_file')
#     wf.connect(surf, 'midthickness_R_32', R_reho, 'surface_file')
#     wf.connect(surf,'atlas_roi_shape_R', L_reho, 'mask')
#     wf.connect(mean_timeseries, 'mean_timeseries', R_reho, 'mean_timeseries')  
#     R_reho.inputs.reho_filename = R_surf_reho.dscalar.nii
#     wf.connect(node, out, R_reho, 'dtseries')

#     outputs = {
#         'surf-L_reho': (L_reho,'L_reho'),
#         'surf-R_reho': (R_reho,'R_reho')}

#     return wf, outputs

# def cal_connectivity_matrix(wf, cfg, strat_pool, pipe_num, opt):
        
#     connectivity_parcellation = pe.Node(util.Function(input_names=['dtseries', 'surf_atlaslabel'], 
#                                 output_names=['parcellation_file'],
#                                 function=run_ciftiparcellate),
#                             name=f'connectivity_parcellation{pipe_num}')
                
#     wf.connect(node, out, connectivity, 'dtseries') 
#     connectivity_parcellation.inputs.surf_atlaslabel = ## path to the label file

#     correlation_matrix = pe.Node(util.Function(input_names=['ptseries'], 
#                                 output_names=['correlation_matrix'],
#                                 function=run_cifticorrelation),
#                             name=f'correlation_matrix{pipe_num}')
                
    
#     wf.connect(connectivity_parcellation, 'parcellation_file', correlation_matrix 'ptseries') 
    
#      outputs = {
#         'surf-correlation_matrix': (correlation_matrix,'correlation_matrix')}

#     return wf, outputs


def surface_falff(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_falff",
     "config": ["amplitude_low_frequency_fluctuation"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["space-fsLR_den-32k_bold-dtseries"],
     "outputs": ["surf_falff"]}
    '''
    wf, outputs = cal_surface_falff(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)

def surface_alff(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "surface_alff",
     "config": ["amplitude_low_frequency_fluctuation"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": ["space-fsLR_den-32k_bold-dtseries"],
     "outputs": ["surf_alff"]}
    '''
    wf, outputs = cal_surface_alff(wf, cfg, strat_pool, pipe_num, opt)

    return (wf, outputs)

# def surface_reho(wf, cfg, strat_pool, pipe_num, opt=None):
#     '''
#     {"name": "ReHo",
#      "config": ["regional_homogeneity"],
#      "switch": ["run"],
#      "option_key": "surface_reho",
#      "option_val": "None",
#      "inputs": ["space-fsLR_den-32k_bold-dtseries"],
#      "outputs": ["surf-L_reho", "surf-R_reho"]}
#     '''
#     wf, outputs = cal_reho(wf, cfg, strat_pool, pipe_num, opt)

#     return (wf, outputs)

# def surface_connectivity_matrix(wf, cfg, strat_pool, pipe_num, opt=None):
#     '''
#     {"name": "surface_connectivity_matrix",
#      "config": ["surface_analysis", "post_freesurfer"],
#      "switch": ["run"],
#      "option_key": "None",
#      "option_val": "None",
#      "inputs": ["space-fsLR_den-32k_bold-dtseries",
#      "outputs": ["surf-correlation_matrix"}
#     '''
#     wf, outputs = cal_connectivity_matrix(wf, cfg, strat_pool, pipe_num, opt)

#     return (wf, outputs)

def run_surf_falff(subject,dtseries):
    import os
    import subprocess
    falff = os.path.join(os.getcwd(), f'{subject}_falff_surf.dscalar.nii.gz')
    cmd = ['ciftify_falff', dtseries, falff, '--min-low-freq', '0.01', '--max-low-freq' , '0.1']
    subprocess.check_output(cmd)
    return falff

def run_surf_alff(subject,dtseries):
   import os
   import subprocess
   alff = os.path.join(os.getcwd(), f'{subject}_alff_surf.dscalar.nii.gz')
   cmd = ['ciftify_falff', dtseries, alff, '--min-low-freq', '0.01', '--max-low-freq' , '0.1' , '--calc-alff']
   subprocess.check_output(cmd)
   return alff
    
# def run_get_cortex(dtseries, structure, cortex_filename):
#     import os
#     import subprocess
#     cortex_file = os.path.join(os.getcwd(), cortex_filename)
#     cmd = ['wb_command', '-cifti-separate', dtseries , 'COLUMN', '-label', structure, cortex_file]
#     subprocess.check_output(cmd)
#     return cortex_file

# def run_mean_timeseries(dtseries,post_freesurfer_folder):

#     import os
#     import subprocess
#     mean_timeseries = os.path.join(os.getcwd(), mean.dscalar.nii)
#     cmd = ['wb_command', '-cifti-reduce', dtseries, 'MEAN', mean_timeseries]
#     subprocess.check_output(cmd)
#     return mean_timeseries
    
# def run_surf_reho(dtseries, mask, cortex_file, surface_file,mean_timeseries,post_freesurfer_folder, reho_filename):

#     import os
#     import subprocess
#     surf_reho = os.path.join(os.getcwd(), reho_filename)
#     cmd = ['python', '/code/CPAC/surface/PostFreeSurfer/surf_reho.py', dtseries, mask, cortex_file, surface_file, mean_timeseries, surf_reho]
#     subprocess.check_output(cmd)
#     return surf_reho
    
# def run_ciftiparcellate(dtseries, surf_atlaslabel):
#     import os  
#     import subprocess
#     parcellation_file = os.path.join(os.getcwd(), 'parcellation.ptseries.nii')
#     cmd = ['wb_command', '-cifti-parcellate', dtseries , surf_atlaslabel, 'COLUMN', parcellation_file ]
#     subprocess.check_output(cmd)
   
#     return parcellation_file

# def run_cifticorrelation(ptseries):
#     import os  
#     import subprocess
#     correlation_matrix = os.path.join(os.getcwd(), 'cifti_corr.pconn.nii')
#     cmd = ['wb_command', '-cifti-correlation ', ptseries , correlation_matrix]
#     subprocess.check_output(cmd)
#     return correlation_matrix
