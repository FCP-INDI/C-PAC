import os
from CPAC.pipeline.nodeblock import nodeblock
import nipype.interfaces.utility as util
from CPAC.utils.interfaces.function import Function
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.surface.PostFreeSurfer.surf_reho import run_surf_reho


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
    
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    import os
    
    recon_all_path = os.path.join(freesurfer_folder, 'recon_all')

    if os.path.isdir(recon_all_path):
        freesurfer_folder = recon_all_path

        # DCAN-HCP PostFreeSurfer Block1
   
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block1.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]


    log_subprocess(cmd)

       # DCAN-HCP PostFreeSurfer Block2

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block2.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    log_subprocess(cmd)

      # DCAN-HCP PostFreeSurfer Block3

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block3.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    log_subprocess(cmd)

    #DCAN-HCP PostFreeSurfer Block3.5
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block4.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    log_subprocess(cmd)

    # DCAN-HCP PostFreeSurfer Block4

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block5.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    log_subprocess(cmd)

  # DCAN-HCP PostFreeSurfer Block5

    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block6.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    log_subprocess(cmd)

   # DCAN-HCP PostFreeSurfer Block6
      
    cmd = ['bash', '/code/CPAC/surface/PostFreeSurfer/block7.sh', post_freesurfer_folder, \
          freesurfer_folder, subject, \
          t1w_restore_image, atlas_space_t1w_image, \
          atlas_transform, inverse_atlas_transform, \
          surf_atlas_dir, gray_ordinates_dir, gray_ordinates_res, \
          high_res_mesh, low_res_mesh, \
          subcortical_gray_labels, freesurfer_labels]

    log_subprocess(cmd)
    # DCAN-HCP fMRISurface
    # https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRISurface/GenericfMRISurfaceProcessingPipeline.sh
    cmd = ['bash', '/code/CPAC/surface/fMRISurface/run.sh',
           '--post_freesurfer_folder', post_freesurfer_folder,
           '--subject', subject, '--fmri', atlas_space_bold, '--scout',
           scout_bold, '--lowresmesh', low_res_mesh, '--grayordinatesres',
           gray_ordinates_res, '--fmrires', fmri_res, '--smoothingFWHM',
           smooth_fwhm]
    log_subprocess(cmd)

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

@nodeblock(
    name="surface_postproc",
    config=["surface_analysis", "post_freesurfer"],
    switch=["run"],
    inputs=[
        (
            "freesurfer-subject-dir",
            [
                "pipeline-fs_desc-restore_T1w",
                "desc-preproc_T1w",
                "desc-reorient_T1w",
                "T1w",
                "space-longitudinal_desc-reorient_T1w",
            ],
            [
                "space-template_desc-head_T1w",
                "space-template_desc-brain_T1w",
                "space-template_desc-T1w_mask",
            ],
            [
                "from-T1w_to-template_mode-image_xfm",
                "from-T1w_to-template_mode-image_desc-linear_xfm",
            ],
            [
                "from-template_to-T1w_mode-image_xfm",
                "from-template_to-T1w_mode-image_desc-linear_xfm",
            ],
            ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"],
            [
                "space-template_desc-scout_bold",
                "space-template_desc-cleaned_bold",
                "space-template_desc-brain_bold",
                "space-template_desc-motion_bold",
                "space-template_bold",
            ],
        )
    ],
    outputs=[
                "space-fsLR_den-32k_bold",
                 "atlas-DesikanKilliany_space-fsLR_den-32k",
                 "atlas-Destrieux_space-fsLR_den-32k",
                 "atlas-DesikanKilliany_space-fsLR_den-164k",
                 "atlas-Destrieux_space-fsLR_den-164k",
                 "space-fsLR_den-32k_bold-dtseries",
                 "AtlasSubcortical-s2",
                 "goodvoxels",
                 "ribbon-only",
                 "hemi-L_space-fsLR_den-32k_desc-atlasroi_bold",
                 "hemi-R_space-fsLR_den-32k_desc-atlasroi_bold",
                 "hemi-L_space-fsLR_den-32k_desc-atlasroi_mask",
                 "hemi-R_space-fsLR_den-32k_desc-atlasroi_mask",
                 "hemi-L_space-native_bold",
                 "hemi-R_space-native_bold",
                 "space-fsLR_den-32k_wb-spec",
                 "space-native_wb-spec",
                 "hemi-L_space-fsLR_den-32k_desc-FS_arealdistortion",
                 "hemi-R_space-fsLR_den-32k_desc-FS_arealdistortion",
                 "space-fsLR_den-32k_desc-FS_arealdistortion",
                 "hemi-L_space-fsLR_den-32k_desc-MSMSulc_arealdistortion",
                 "hemi-R_space-fsLR_den-32k_desc-MSMSulc_arealdistortion",
                 "space-fsLR_den-32k_desc-MSMSulc_arealdistortion",
                 "hemi-L_space-fsLR_den-32k_desc-FS_edgedistortion",
                 "hemi-R_space-fsLR_den-32k_desc-FS_edgedistortion",
                 "space-fsLR_den-32k_desc-FS_edgedistortion",
                 "hemi-L_space-fsLR_den-32k_desc-MSMSulc_edgedistortion",
                 "hemi-R_space-fsLR_den-32k_desc-MSMSulc_edgedistortion",
                 "space-fsLR_den-32k_desc-MSMSulc_edgedistortion",
                 "hemi-L_space-fsLR_den-32k_curv",
                 "hemi-R_space-fsLR_den-32k_curv",
                 "space-fsLR_den-32k_curv",
                 "hemi-L_space-fsLR_den-32k_flat",
                 "hemi-R_space-fsLR_den-32k_flat",
                 "hemi-L_space-fsLR_den-32k_inflated",
                 "hemi-R_space-fsLR_den-32k_inflated",
                 "hemi-L_space-fsLR_den-32k_veryinflated",
                 "hemi-R_space-fsLR_den-32k_veryinflated",
                 "hemi-L_space-native_inflated",
                 "hemi-R_space-native_inflated",
                 "hemi-L_space-native_veryinflated",
                 "hemi-R_space-native_veryinflated",
                 "hemi-L_space-fsLR_den-164k_midthickness",
                 "hemi-R_space-fsLR_den-164k_midthickness",
                 "hemi-L_space-fsLR_den-32k_midthickness",
                 "hemi-R_space-fsLR_den-32k_midthickness",
                 "hemi-L_space-native_midthickness",
                 "hemi-R_space-native_midthickness",
                 "hemi-L_space-fsLR_den-32k_pial",
                 "hemi-R_space-fsLR_den-32k_pial",
                 "hemi-L_space-native_den-32k_pial",
                 "hemi-R_space-native_den-32k_pial",
                 "hemi-L_space-fsLR_den-32k_sphere",
                 "hemi-R_space-fsLR_den-32k_sphere",
                 "hemi-L_space-native_desc-MSMSulc_sphere",
                 "hemi-R_space-native_desc-MSMSulc_sphere",
                 "hemi-L_space-native_sphere",
                 "hemi-R_space-native_sphere",
                 "hemi-L_space-native_desc-reg_sphere",
                 "hemi-R_space-native_desc-reg_sphere",
                 "hemi-L_space-native_desc-reg-reg_sphere",
                 "hemi-R_space-native_desc-reg-reg_sphere",
                 "hemi-L_space-native_desc-rot_sphere",
                 "hemi-R_space-native_desc-rot_sphere",
                 "hemi-L_space-fsLR_den-32k_desc-FS_strainJ",
                 "hemi-R_space-fsLR_den-32k_desc-FS_strainJ",
                 "space-fsLR_den-32k_desc-FS_strainJ",
                 "hemi-L_space-fsLR_den-32k_desc-MSMSulc_strainJ",
                 "hemi-R_space-fsLR_den-32k_desc-MSMSulc_strainJ",
                 "space-fsLR_den-32k_desc-MSMSulc_strainJ",
                 "hemi-L_space-fsLR_den-32k_desc-FS_strainR",
                 "hemi-R_space-fsLR_den-32k_desc-FS_strainR",
                 "space-fsLR_den-32k_desc-FS_strainR",
                 "hemi-L_space-fsLR_den-32k_desc-MSMSulc_strainR",
                 "hemi-R_space-fsLR_den-32k_desc-MSMSulc_strainR",
                 "space-fsLR_den-32k_desc-MSMSulc_strainR",
                 "hemi-L_space-fsLR_den-32k_sulc",
                 "hemi-R_space-fsLR_den-32k_sulc",
                 "space-fsLR_den-32k_sulc",
                 "hemi-L_space-fsLR_den-32k_thickness",
                 "hemi-R_space-fsLR_den-32k_thickness",
                 "space-fsLR_den-32k_thickness",
                 "hemi-L_space-fsLR_den-164k_white",
                 "hemi-R_space-fsLR_den-164k_white",
                 "hemi-L_space-fsLR_den-32k_white",
                 "hemi-R_space-fsLR_den-32k_white",
                 "hemi-L_space-native_white",
                 "hemi-R_space-native_white",
    ],
)

def surface_postproc(wf, cfg, strat_pool, pipe_num, opt=None):

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
        'cpac_' + cfg['subject_id'],
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

    restore = ["pipeline-fs_desc-restore_T1w", "desc-preproc_T1w", "desc-reorient_T1w", "T1w",
                  "space-longitudinal_desc-reorient_T1w"]
    space_temp = ["space-template_desc-head_T1w", "space-template_desc-brain_T1w", "space-template_desc-T1w_mask"]
    atlas_xfm = ["from-T1w_to-template_mode-image_xfm", "from-T1w_to-template_mode-image_desc-linear_xfm"]
    atlas_xfm_inv = ["from-template_to-T1w_mode-image_xfm", "from-template_to-T1w_mode-image_desc-linear_xfm"]
    atlas_space_bold = ["space-template_desc-brain_bold", "space-template_desc-preproc_bold"]
    scout_bold = ["space-template_desc-scout_bold", "space-template_desc-cleaned_bold", "space-template_desc-brain_bold",
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
        'space-fsLR_den-32k_bold': (surf, 'dtseries'),
        'atlas-DesikanKilliany_space-fsLR_den-32k': (surf,
                                                            'desikan_'
                                                            'killiany_32'),
        'atlas-Destrieux_space-fsLR_den-32k': (surf, 'destrieux_32'),
        'atlas-DesikanKilliany_space-fsLR_den-164k': (surf,
                                                             'desikan_'
                                                             'killiany_164'),
        'atlas-Destrieux_space-fsLR_den-164k': (surf, 'destrieux_164'),
        'AtlasSubcortical-s2': (surf, 'subcortical_atlas'), 
        'goodvoxels': (surf, 'good_voxels'),
        'ribbon-only': (surf, 'ribbon_only'),
        'hemi-L_space-fsLR_den-32k_desc-atlasroi_bold': (surf, 
                                    'atlas_roi_func_'
                                    'L'),
        'hemi-R_space-fsLR_den-32k_desc-atlasroi_bold': (surf, 
                                    'atlas_roi_func_'
                                    'R'),
        'hemi-L_space-fsLR_den-32k_desc-atlasroi_mask': (surf, 
                                    'atlas_roi_shape_'
                                    'L'),
        'hemi-R_space-fsLR_den-32k_desc-atlasroi_mask': (surf, 
                                    'atlas_roi_shape_'
                                    'R'),                            
        'hemi-L_space-native_bold': (surf, 'native_' 'L'),
        'hemi-R_space-native_bold': (surf, 'native_' 'R'),
        'space-fsLR_den-32k_wb-spec': (surf, 'spec_' 'LR_32k'), 
        'space-native_wb-spec': (surf, 'spec_' 'native'),
        'hemi-L_space-fsLR_den-32k_desc-FS_arealdistortion': (surf, 'areal_distortion_'
                                                                        'FS_'
                                                                        'L'),
        'hemi-R_space-fsLR_den-32k_desc-FS_arealdistortion': (surf, 'areal_distortion_'
                                                                        'FS_'
                                                                        'R'),
        'space-fsLR_den-32k_desc-FS_arealdistortion': (surf, 'areal_distortion_'
                                                                        'FS_'
                                                                        'dscalar'),
        'hemi-L_space-fsLR_den-32k_desc-MSMSulc_arealdistortion': (surf, 'areal_distortion_'
                                                                        'MSMSulc_'
                                                                        'L'),
        'hemi-R_space-fsLR_den-32k_desc-MSMSulc_arealdistortion': (surf, 'areal_distortion_'
                                                                        'MSMSulc_'
                                                                        'R'),
        'space-fsLR_den-32k_desc-MSMSulc_arealdistortion': (surf, 'areal_distortion_'
                                                                        'MSMSulc_'
                                                                        'dscalar'),
        'hemi-L_space-fsLR_den-32k_desc-FS_edgedistortion': (surf, 'edge_distortion_'
                                                                        'FS_'
                                                                        'L'),
        'hemi-R_space-fsLR_den-32k_desc-FS_edgedistortion': (surf, 'edge_distortion_'
                                                                        'FS_'
                                                                        'R'),
        'space-fsLR_den-32k_desc-FS_edgedistortion': (surf, 'edge_distortion_'
                                                                        'FS_'
                                                                        'dscalar'),
        'hemi-L_space-fsLR_den-32k_desc-MSMSulc_edgedistortion': (surf, 'edge_distortion_'
                                                                        'MSMSulc_'
                                                                        'L'),
        'hemi-R_space-fsLR_den-32k_desc-MSMSulc_edgedistortion': (surf, 'edge_distortion_'
                                                                        'MSMSulc_'
                                                                        'R'),
        'space-fsLR_den-32k_desc-MSMSulc_edgedistortion': (surf, 'edge_distortion_'
                                                                        'MSMSulc_'
                                                                        'dscalar'),
        'hemi-L_space-fsLR_den-32k_curv': (surf, 'curvature_' 'L'),
        'hemi-R_space-fsLR_den-32k_curv': (surf, 'curvature_' 'R'),
        'space-fsLR_den-32k_curv': (surf, 'curvature_' 'dscalar'),
        'hemi-L_space-fsLR_den-32k_flat': (surf, 'flat_L'),
        'hemi-R_space-fsLR_den-32k_flat': (surf, 'flat_R'),
        'hemi-L_space-fsLR_den-32k_inflated': (surf, '32k_inflated_L'),
        'hemi-R_space-fsLR_den-32k_inflated': (surf, '32k_inflated_R'),
        'hemi-L_space-fsLR_den-32k_veryinflated': (surf, '32k_very_inflated_L'),
        'hemi-R_space-fsLR_den-32k_veryinflated': (surf, '32k_very_inflated_R'),
        'hemi-L_space-native_inflated': (surf, 'inflated_L'),
        'hemi-R_space-native_inflated': (surf, 'inflated_R'),
        'hemi-L_space-native_veryinflated': (surf, 'very_inflated_L'),
        'hemi-R_space-native_veryinflated': (surf, 'very_inflated_R'),
        'hemi-L_space-fsLR_den-164k_midthickness': (surf, 'midthickness_L_164'),
        'hemi-L_space-fsLR_den-32k_midthickness': (surf, 'midthickness_L_32'),
        'hemi-L_space-native_midthickness': (surf, 'midthickness_L_native'),
        'hemi-R_space-fsLR_den-164k_midthickness': (surf, 'midthickness_R_164'),
        'hemi-R_space-fsLR_den-32k_midthickness': (surf, 'midthickness_R_32'),
        'hemi-R_space-native_midthickness': (surf, 'midthickness_R_native'),
        'hemi-L_space-fsLR_den-32k_pial': (surf, 'pial_L_32'),
        'hemi-L_space-native_den-32k_pial': (surf, 'pial_L_native'),
        'hemi-R_space-fsLR_den-32k_pial': (surf, 'pial_R_32'),
        'hemi-R_space-native_den-32k_pial': (surf, 'pial_R_native'),
        'hemi-L_space-fsLR_den-32k_sphere': (surf, 'sphere_32k_fs_LR_L'),
        'hemi-R_space-fsLR_den-32k_sphere': (surf, 'sphere_32k_fs_LR_R'),                            
        'hemi-L_space-native_desc-MSMSulc_sphere': (surf, 'sphere_MSMSulc_L'),
        'hemi-R_space-native_desc-MSMSulc_sphere': (surf, 'sphere_MSMSulc_R'),
        'hemi-L_space-native_sphere': (surf, 'sphere_native_L'),
        'hemi-R_space-native_sphere': (surf, 'sphere_native_R'),
        'hemi-L_space-native_desc-reg_sphere': (surf, 'sphere_reg_L'),
        'hemi-R_space-native_desc-reg_sphere': (surf, 'sphere_reg_R'),
        'hemi-L_space-native_desc-reg-reg_sphere': (surf, 'sphere_reg_reg_LR_L'),
        'hemi-R_space-native_desc-reg-reg_sphere': (surf, 'sphere_reg_reg_LR_R'),
        'hemi-L_space-native_desc-rot_sphere': (surf, 'sphere_rot_L'),
        'hemi-R_space-native_desc-rot_sphere': (surf, 'sphere_rot_R'),
        'hemi-L_space-fsLR_den-32k_desc-FS_strainJ': (surf, 'StrainJ_FS_L'),
        'hemi-R_space-fsLR_den-32k_desc-FS_strainJ': (surf, 'StrainJ_FS_R'),
        'space-fsLR_den-32k_desc-FS_strainJ': (surf, 'StrainJ_FS_dscalar'),
        'hemi-L_space-fsLR_den-32k_desc-MSMSulc_strainJ': (surf, 'StrainJ_MSMSulc_L'),
        'hemi-R_space-fsLR_den-32k_desc-MSMSulc_strainJ': (surf, 'StrainJ_MSMSulc_R'),
        'space-fsLR_den-32k_desc-FS_strainJ': (surf, 'StrainJ_MSMSulc_dscalar'),
        'hemi-L_space-fsLR_den-32k_desc-FS_strainR': (surf, 'StrainR_FS_L'),
        'hemi-R_space-fsLR_den-32k_desc-FS_strainR': (surf, 'StrainR_FS_R'),
        'space-fsLR_den-32k_desc-FS_strainR': (surf, 'StrainR_FS_dscalar'),
        'hemi-L_space-fsLR_den-32k_desc-MSMSulc_strainR': (surf, 'StrainR_MSMSulc_L'),
        'hemi-R_space-fsLR_den-32k_desc-MSMSulc_strainR': (surf, 'StrainR_MSMSulc_R'),
        'space-fsLR_den-32k_desc-MSMSulc_strainR': (surf, 'StrainR_MSMSulc_dscalar'),
        'hemi-L_space-fsLR_den-32k_sulc': (surf, 'sulc_L'),
        'hemi-R_space-fsLR_den-32k_sulc': (surf, 'sulc_R'),
        'space-fsLR_den-32k_sulc': (surf, 'sulc_dscalar'),                                   
        'hemi-L_space-fsLR_den-32k_thickness': (surf, 'thickness_L'),
        'hemi-R_space-fsLR_den-32k_thickness': (surf, 'thickness_R'),
        'space-fsLR_den-32k_thickness': (surf, 'thickness_dscalar'),                                   
        'hemi-L_space-fsLR_den-164k_white': (surf, 'white_L_164'),
        'hemi-L_space-fsLR_den-32k_white': (surf, 'white_L_32'),
        'hemi-L_space-native_white': (surf, 'white_L_native'),
        'hemi-R_space-fsLR_den-164k_white': (surf, 'white_R_164'),
        'hemi-R_space-fsLR_den-32k_white': (surf, 'white_R_32'),
        'hemi-R_space-native_white': (surf, 'white_R_native')
    }

    return wf, outputs


@nodeblock(
    name="surface_falff",
    config=["surface_analysis", "amplitude_low_frequency_fluctuation"],
    switch=["run"],
    inputs=["space-fsLR_den-32k_bold"],
    outputs=[
        "space-fsLR_den-32k_bold_surf_falff"],
)
def surface_falff(wf, cfg, strat_pool, pipe_num, opt):

    falff = pe.Node(util.Function(input_names=['subject','dtseries'], 
                                output_names=['surf_falff'],
                                function=run_surf_falff),
                   name=f'surf_falff_{pipe_num}')
    
    falff.inputs.subject = cfg['subject_id']
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold') 
    wf.connect(node, out, falff, 'dtseries')
    
    outputs = {
        'space-fsLR_den-32k_bold_surf_falff': (falff,'surf_falff')}
    return wf, outputs


@nodeblock(
    name="surface_alff",
    config=["surface_analysis", "amplitude_low_frequency_fluctuation"],
    switch=["run"],
    inputs=["space-fsLR_den-32k_bold"],
    outputs=[
        "space-fsLR_den-32k_bold_surf_alff"],
)   

def surface_alff(wf, cfg, strat_pool, pipe_num, opt):

    alff = pe.Node(util.Function(input_names=['subject', 'dtseries'], 
                             output_names=['surf_alff'],
                               function=run_surf_alff),
                 name=f'surf_alff_{pipe_num}')
    
    alff.inputs.subject = cfg['subject_id']
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold') 
    wf.connect(node, out, alff, 'dtseries')
    outputs = {
        'space-fsLR_den-32k_bold_surf_alff': (alff, 'surf_alff')}
    return wf, outputs


@nodeblock(
    name="surface_reho",
    config=["surface_analysis", "regional_homogeneity"],
    switch=["run"],
    inputs=["space-fsLR_den-32k_bold",
            "hemi-L_space-fsLR_den-32k_midthickness", 
            "hemi-L_space-fsLR_den-32k_desc-atlasroi_mask",
            "hemi-R_space-fsLR_den-32k_midthickness",   
            "hemi-R_space-fsLR_den-32k_desc-atlasroi_mask"],
    outputs=[
        "space-fsLR_den-32k_bold_surf-L_reho",
        "space-fsLR_den-32k_bold_surf-R_reho"],
)  
def surface_reho(wf, cfg, strat_pool, pipe_num, opt):

    L_cortex_file = pe.Node(util.Function(input_names=['subject', 'dtseries', 'structure', 'cortex_filename'], 
                                output_names=['L_cortex_file'],
                                function=run_get_cortex),
                  name=f'L_surf_cortex_{pipe_num}')
    
    L_cortex_file.inputs.subject = cfg['subject_id']
    L_cortex_file.inputs.structure = "CORTEX_LEFT"
    L_cortex_file.inputs.cortex_filename = "L_cortex.func.gii"
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold') 
    wf.connect(node, out, L_cortex_file, 'dtseries')

    R_cortex_file = pe.Node(util.Function(input_names=['subject', 'dtseries', 'structure', 'cortex_filename'], 
                                output_names=['R_cortex_file'],
                                function=run_get_cortex),
                  name=f'R_surf_cortex_{pipe_num}')

    R_cortex_file.inputs.subject = cfg['subject_id']
    R_cortex_file.inputs.structure = "CORTEX_RIGHT"
    R_cortex_file.inputs.cortex_filename = "R_cortex.func.gii"
    wf.connect(node, out,R_cortex_file, 'dtseries')


    mean_timeseries = pe.Node(util.Function(input_names=['subject', 'dtseries'], 
                                output_names=['mean_timeseries'],
                                function=run_mean_timeseries),
                  name=f'mean_timeseries_{pipe_num}')

    mean_timeseries.inputs.subject = cfg['subject_id']
    wf.connect(node, out, mean_timeseries, 'dtseries')

    L_reho = pe.Node(util.Function(input_names=['subject', 'dtseries', 'mask' , 'cortex_file', 'surface_file', 'mean_timeseries','reho_filename', 'structure_name'], 
                                output_names=['L_reho'],
                                function=run_surf_reho),
                    name=f'L_surf_reho_{pipe_num}')
    
    L_reho.inputs.subject = cfg['subject_id']
    wf.connect(L_cortex_file, 'L_cortex_file', L_reho, 'cortex_file')
    node, out = strat_pool.get_data('hemi-L_space-fsLR_den-32k_midthickness')
    wf.connect(node, out, L_reho, 'surface_file')
    
    node, out = strat_pool.get_data('hemi-L_space-fsLR_den-32k_desc-atlasroi_mask')
    wf.connect(node, out, L_reho, 'mask')
    wf.connect(mean_timeseries, 'mean_timeseries', L_reho, 'mean_timeseries') 
    L_reho.inputs.reho_filename = "L_surf_reho.dscalar.nii"
    L_reho.inputs.structure_name = "CIFTI_STRUCTURE_CORTEX_LEFT"
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold')
    wf.connect(node, out, L_reho, 'dtseries')

    R_reho = pe.Node(util.Function(input_names=['subject','dtseries', 'mask' , 'cortex_file', 'surface_file', 'mean_timeseries', 'reho_filename', 'structure_name'], 
                                output_names=['R_reho'],
                                function=run_surf_reho),
                    name=f'R_surf_reho_{pipe_num}')
    
    R_reho.inputs.subject = cfg['subject_id']
    wf.connect(R_cortex_file, 'R_cortex_file', R_reho, 'cortex_file')
    node, out = strat_pool.get_data('hemi-R_space-fsLR_den-32k_midthickness')
    wf.connect(node, out, R_reho, 'surface_file')
    R_reho.inputs.structure_name = "CIFTI_STRUCTURE_CORTEX_RIGHT"
    node, out = strat_pool.get_data('hemi-R_space-fsLR_den-32k_desc-atlasroi_mask')
    wf.connect(node, out, R_reho, 'mask')
    wf.connect(mean_timeseries, 'mean_timeseries', R_reho, 'mean_timeseries') 
    R_reho.inputs.reho_filename = "R_surf_reho.dscalar.nii"
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold')
    wf.connect(node, out, R_reho, 'dtseries')

    outputs = {
        'space-fsLR_den-32k_bold_surf-L_reho': (L_reho,'L_reho'),
        'space-fsLR_den-32k_bold_surf-R_reho': (R_reho,'R_reho')}

    return wf, outputs

@nodeblock(
    name="surface_connectivity_matrix",
    config=["surface_analysis", "surface_connectivity"],
    switch=["run"],
    inputs=["space-fsLR_den-32k_bold"],
    outputs=["space-fsLR_den-32k_bold_surf-correlation_matrix"],
)  

def surface_connectivity_matrix(wf, cfg, strat_pool, pipe_num, opt):

     
    connectivity_parcellation = pe.Node(util.Function(input_names=['subject', 'dtseries', 'surf_atlaslabel'], 
                                output_names=['parcellation_file'],
                                function=run_ciftiparcellate),
                            name=f'connectivity_parcellation_{pipe_num}')

    connectivity_parcellation.inputs.subject = cfg['subject_id'] 
    node, out = strat_pool.get_data('space-fsLR_den-32k_bold')        
    wf.connect(node, out, connectivity_parcellation, 'dtseries') 
    connectivity_parcellation.inputs.surf_atlaslabel = cfg.surface_analysis['surface_connectivity']['surface_parcellation_template']

    correlation_matrix = pe.Node(util.Function(input_names=['subject','ptseries'], 
                                output_names=['correlation_matrix'],
                                function=run_cifticorrelation),
                            name=f'correlation_matrix_{pipe_num}')
                
    correlation_matrix.inputs.subject = cfg['subject_id'] 
    wf.connect(connectivity_parcellation, 'parcellation_file', correlation_matrix, 'ptseries') 
    
    outputs = {
        'space-fsLR_den-32k_bold_surf-correlation_matrix': (correlation_matrix,'correlation_matrix')}

    return wf, outputs


def run_surf_falff(subject, dtseries):
    import os
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    falff = os.path.join(os.getcwd(), f'{subject}_falff_surf.dscalar.nii')
    cmd = ['ciftify_falff', dtseries, falff, '--min-low-freq', '0.01', '--max-low-freq' , '0.1']
    log_subprocess(cmd)
    return falff

def run_surf_alff(subject, dtseries):
   import os
   from CPAC.utils.monitoring.custom_logging import log_subprocess
   alff = os.path.join(os.getcwd(), f'{subject}_alff_surf.dscalar.nii')
   cmd = ['ciftify_falff', dtseries, alff, '--min-low-freq', '0.01', '--max-low-freq' , '0.1' , '--calc-alff']
   log_subprocess(cmd)
   return alff
    
def run_get_cortex(subject, dtseries, structure, cortex_filename):
    import os
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    cortex_file = os.path.join(os.getcwd(), f'{subject}_{cortex_filename}')
    cmd = ['wb_command', '-cifti-separate', dtseries , 'COLUMN', '-metric', structure, cortex_file]
    log_subprocess(cmd)
    return cortex_file

def run_mean_timeseries(subject, dtseries):

    import os
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    mean_timeseries = os.path.join(os.getcwd(), f'{subject}_mean.dscalar.nii')
    cmd = ['wb_command', '-cifti-reduce', dtseries, 'MEAN', mean_timeseries]
    log_subprocess(cmd)
    return mean_timeseries
    
def run_ciftiparcellate(subject, dtseries, surf_atlaslabel):
    import os  
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    parcellation_file = os.path.join(os.getcwd(), f'{subject}_parcellation.ptseries.nii')
    cmd = ['wb_command', '-cifti-parcellate', dtseries , surf_atlaslabel, 'COLUMN', parcellation_file ]
    log_subprocess(cmd)
    return parcellation_file

def run_cifticorrelation(subject, ptseries):
    import os  
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    correlation_matrix = os.path.join(os.getcwd(), f'{subject}_cifti_corr.pconn.nii')
    cmd = ['wb_command', '-cifti-correlation', ptseries , correlation_matrix]
    log_subprocess(cmd)
    return correlation_matrix
