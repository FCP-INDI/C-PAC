# -*- coding: utf-8 -*-
import nipype.interfaces.utility as util

from CPAC.pipeline import nipype_pipeline_engine as pe
from nibabel import load as nib_load, Nifti1Image
from numpy import zeros

def get_shape(nifti_image):
    return nib_load(nifti_image).shape

def pad(cropped_image_path, target_image_path):
    """
    Pad a cropped image to match the dimensions of a target image along the z-axis, 
    while keeping padded image aligned with target_image.

    Parameters:
    - cropped_image_path (str): The file path to the cropped image (NIfTI format).
    - target_image_path (str): The file path to the target image (NIfTI format).

    Returns:
    - str: The file path to the saved padded image (NIfTI format).

    The function loads cropped and target iamges, calculates the z-dimension shift required for alignment such 
    that the mask generated from padded image will work correctly on the target image. The result padded image is
    saved as an NIfTI file in the working directory/node and file path is returned as output.

    Note: The function assumes that the input images are in NIfTI format and have compatible dimensions. The cropped
    and target image should only differ in z-axis dimension.
    """
    from numpy import asanyarray, zeros_like, ndarray
    from nibabel import load, save, Nifti1Image
    from os import path, getcwd
    from typing import Optional

    cropped_image: Optional[ndarray] = asanyarray(load(cropped_image_path).dataobj)
    target_image: Optional[ndarray] = asanyarray(load(target_image_path).dataobj)

    # Taking 1 slice to calculate the z dimension shift from top
    center_row:int =target_image.shape[0]//2
    center_column:int = target_image.shape[1]//2
    z_slice_cropped_image: Optional[ndarray] = cropped_image[center_row, center_column, :]
    z_slice_target_image: Optional[ndarray] = target_image[center_row, center_column, :]

    for z_shift in range(len(z_slice_target_image) - len(z_slice_cropped_image) + 1):
        if (z_slice_target_image[z_shift:z_shift+len(z_slice_cropped_image)] == z_slice_cropped_image).all():
            break

    padded_image_matrix: Optional[ndarray] = zeros_like(target_image)
    padded_image_matrix[:, :, z_shift:cropped_image.shape[2]+z_shift] = cropped_image
    padded_image_path:str = path.join(getcwd(),"padded_image_T1w.nii.gz")
    cropped_image = load(cropped_image_path)
    save(Nifti1Image(padded_image_matrix, affine=cropped_image.affine), padded_image_path)
    return padded_image_path

def fsl_aff_to_rigid(in_xfm, out_name):

    out_mat = os.path.join(os.getcwd(), out_name)

    #   Script for getting a 6 DOF approx to a 12 DOF standard transformation
    #
    #   Mark Jenkinson
    #   FMRIB Image Analysis Group
    #
    #   Copyright (C) 2012 University of Oxford
    #
    #   Part of FSL - FMRIB's Software Library
    #   http://www.fmrib.ox.ac.uk/fsl
    #   fsl@fmrib.ox.ac.uk
    #
    #   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    #   Imaging of the Brain), Department of Clinical Neurology, Oxford
    #   University, Oxford, UK
    #
    #
    #   LICENCE
    #
    #   FMRIB Software Library, Release 5.0 (c) 2012, The University of
    #   Oxford (the "Software")
    #
    #   The Software remains the property of the University of Oxford ("the
    #   University").
    #
    #   The Software is distributed "AS IS" under this Licence solely for
    #   non-commercial use in the hope that it will be useful, but in order
    #   that the University as a charitable foundation protects its assets for
    #   the benefit of its educational and research purposes, the University
    #   makes clear that no condition is made or to be implied, nor is any
    #   warranty given or to be implied, as to the accuracy of the Software,
    #   or that it will be suitable for any particular purpose or for use
    #   under any specific conditions. Furthermore, the University disclaims
    #   all responsibility for the use which is made of the Software. It
    #   further disclaims any liability for the outcomes arising from using
    #   the Software.
    #
    #   The Licensee agrees to indemnify the University and hold the
    #   University harmless from and against any and all claims, damages and
    #   liabilities asserted by third parties (including claims for
    #   negligence) which arise directly or indirectly from the use of the
    #   Software or the sale of any products based on the Software.
    #
    #   No part of the Software may be reproduced, modified, transmitted or
    #   transferred in any form or by any means, electronic or mechanical,
    #   without the express permission of the University. The permission of
    #   the University is not required if the said reproduction, modification,
    #   transmission or transference is done without financial return, the
    #   conditions of this Licence are imposed upon the receiver of the
    #   product, and all original and amended source code is included in any
    #   transmitted product. You may be held legally responsible for any
    #   copyright infringement that is caused or encouraged by your failure to
    #   abide by these terms and conditions.
    #
    #   You are not permitted under this Licence to use this Software
    #   commercially. Use for which any financial return is received shall be
    #   defined as commercial use, and includes (1) integration of all or part
    #   of the source code or the Software into a product for sale or license
    #   by or on behalf of Licensee to third parties or (2) use of the
    #   Software or any derivative of it for research with the final aim of
    #   developing software products for sale or license to a third party or
    #   (3) use of the Software or any derivative of it for research with the
    #   final aim of developing non-software products for sale or license to a
    #   third party, or (4) use of the Software to provide any service to an
    #   external organisation for which payment is received. If you are
    #   interested in using the Software commercially, please contact Isis
    #   Innovation Limited ("Isis"), the technology transfer company of the
    #   University, to negotiate a licence. Contact details are:
    #   innovation@isis.ox.ac.uk quoting reference DE/9564.

    # Load in the necessary info
    a = loadtxt(in_xfm)
    # set specific AC and PC coordinates in FLIRT convention (x1=AC, x2=PC, x3=point above x1 in the mid-sag plane)
    x1 = matrix([[91], [129], [67], [1]])
    x2 = matrix([[91], [100], [70], [1]])
    x3 = matrix([[91], [129], [117], [1]])

    ainv = linalg.inv(a)

    # vectors v are in MNI space, vectors w are in native space
    v21 = (x2 - x1)
    v31 = (x3 - x1)
    # normalise and force orthogonality
    v21 = v21 / linalg.norm(v21)
    v31 = v31 - multiply(v31.T * v21, v21)
    v31 = v31 / linalg.norm(v31)
    tmp = cross(v21[0:3, 0].T, v31[0:3, 0].T).T
    v41 = mat(zeros((4, 1)))
    v41[0:3, 0] = tmp
    # Map vectors to native space
    w21 = ainv * (v21)
    w31 = ainv * (v31)
    # normalise and force orthogonality
    w21 = w21 / linalg.norm(w21)
    w31 = w31 - multiply(w31.T * w21, w21)
    w31 = w31 / linalg.norm(w31)
    tmp = cross(w21[0:3, 0].T, w31[0:3, 0].T).T
    w41 = mat(zeros((4, 1)))
    w41[0:3, 0] = tmp

    # setup matrix: native to MNI space
    r1 = matrix(eye(4))
    r1[0:4, 0] = w21
    r1[0:4, 1] = w31
    r1[0:4, 2] = w41
    r2 = matrix(eye(4))
    r2[0, 0:4] = v21.T
    r2[1, 0:4] = v31.T
    r2[2, 0:4] = v41.T
    r = r2.T * r1.T

    # Fix the translation (keep AC=x1 in the same place)
    ACmni = x1
    ACnat = ainv * x1
    trans = ACmni - r * ACnat
    r[0:3, 3] = trans[0:3]

    # Save out the result
    savetxt(out_mat, r, fmt='%14.10f')

    return out_mat


def freesurfer_hemispheres(wf, reconall, pipe_num):
    """Function to return various hemisphere-specific FreeSurfer outputs.

    Parameters
    ----------
    wf : nipype.pipeline.engine.workflows.Workflow
        Workflow object.

    reconall : nipype.pipeline.engine.nodes.Node

    pipe_num : int

    Returns
    -------
    wf : nipype.pipeline.engine.workflows.Workflow

    outputs : dict
    """
    def split_hemi(multi_file):
        # pylint: disable=invalid-name
        lh = None
        rh = None
        for filepath in multi_file:
            if 'lh.' in filepath:
                lh = filepath
            if 'rh.' in filepath:
                rh = filepath
        return (lh, rh)

    def split_hemi_interface():
        """Returns a function interface for split_hemi."""
        return util.Function(input_names=['multi_file'],
                             output_names=['lh', 'rh'],
                             function=split_hemi)

    splits = {
        label: pe.Node(split_hemi_interface(),
                       name=f'split_{label}_{pipe_num}') for
        label in ['curv', 'pial', 'smoothwm', 'sphere', 'sulc', 'thickness',
                  'volume', 'white']
    }
    for label in splits:
        wf.connect(reconall, label, splits[label], 'multi_file')
    outputs = {
        'pipeline-fs_hemi-L_desc-surface_curv': (splits['curv'], 'lh'),
        'pipeline-fs_hemi-R_desc-surface_curv': (splits['curv'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMesh_pial': (splits['pial'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMesh_pial': (splits['pial'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMesh_smoothwm': (splits['smoothwm'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMesh_smoothwm': (splits['smoothwm'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMesh_sphere': (splits['sphere'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMesh_sphere': (splits['sphere'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMap_sulc': (splits['sulc'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMap_sulc': (splits['sulc'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMap_thickness': (splits['thickness'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMap_thickness': (splits['thickness'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMap_volume': (splits['volume'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMap_volume': (splits['volume'], 'rh'),
        'pipeline-fs_hemi-L_desc-surfaceMesh_white': (splits['white'], 'lh'),
        'pipeline-fs_hemi-R_desc-surfaceMesh_white': (splits['white'], 'rh')}

    return wf, outputs


def create_3dskullstrip_arg_string(shrink_fac, var_shrink_fac,
                                   shrink_fac_bot_lim, avoid_vent, niter,
                                   pushout, touchup, fill_hole, avoid_eyes,
                                   use_edge, exp_frac, NN_smooth, smooth_final,
                                   push_to_edge, use_skull, perc_int,
                                   max_inter_iter, blur_fwhm, fac, monkey,
                                   mask_vol):
    """
    Method to return option string for 3dSkullStrip

    Parameters
    ----------
    shrink_fac : float
        Parameter controlling the brain vs non-brain intensity threshold (tb)

    var_shrink_fac : boolean
        Vary the shrink factor with the number of iterations

    shrink_fac_bot_lim : float
        Do not allow the varying SF to go below SFBL

    avoid_vent : boolean
        Avoid ventricles

    niter : float
        Number of iterations

    pushout : boolean
        Consider values above each node in addition to values below the node when deciding on expansion

    touchup : boolean
        Perform touchup operations at end to include areas not covered by surface expansion

    fill_hole : float
         Fill small holes that can result from small surface intersections caused by the touchup operation

    avoid_eyes : boolean
        Avoid eyes

    use_edge : boolean
        Use edge detection to reduce leakage into meninges and eyes

    exp_frac : float
        Speed of expansion

    NN_smooth : float
        Perform nearest neighbor coordinate interpolation every few iterations. Default is 72.

    smooth_final : float
        Perform final surface smoothing after all iterations

    push_to_edge : boolean
        Perform aggressive push to edge at the end

    use_skull : boolean
        Use outer skull to limit expansion of surface into the skull due to very strong shading artifacts

    perc_int : float
        Percentage of segments allowed to intersect surface

    max_inter_iter : float
        Number of iteration to remove intersection problems

    blur_fwhm : float
        Blur dset after spatial normalization

    fac : float
         Multiply input dataset by FAC if range of values is too small

    monkey : boolean
        Use monkey option in SkullStripping
    
    mask_vol : boolean
        Output a mask volume instead of a skull-stripped volume.

    Returns
    -------
    opt_str : string
        Command args
    
    """

    expr = ''
    defaults = dict(
        fill_hole=10 if touchup else 0,
        shrink_fac=0.6,
        shrink_fac_bot_lim=0.4 if use_edge else 0.65,
        niter=250,
        exp_frac=0.1,
        NN_smooth=72,
        smooth_final=20,
        perc_int=0,
        max_inter_iter=4,
        blur_fwhm=0,
        fac=1.0,
        monkey=False,
        mask_vol=False
    )

    if float(shrink_fac) != defaults['shrink_fac']:
        expr += ' -shrink_fac {0}'.format(shrink_fac)

    if not var_shrink_fac:
        expr += ' -no_var_shrink_fac'

    if mask_vol:
        expr += ' -mask_vol'

    if monkey:
        expr += ' -monkey'

    if float(shrink_fac_bot_lim) != defaults['shrink_fac_bot_lim']:
        expr += ' -shrink_fac_bot_lim {0}'.format(shrink_fac_bot_lim)

    if not use_edge:
        expr += ' -no_use_edge'

    if not avoid_vent:
        expr += ' -no_avoid_vent'

    if int(niter) != defaults['niter']:
        expr += ' -niter {0}'.format(niter)

    if not pushout:
        expr += ' -no_pushout'

    if not touchup:
        expr += ' -no_touchup'

    if int(fill_hole) != defaults['fill_hole']:
        expr += ' -fill_hole {0}'.format(fill_hole)

    if not avoid_eyes:
        expr += ' -no_avoid_eyes'

    if float(exp_frac) != defaults['exp_frac']:
        expr += ' -exp_frac {0}'.format(exp_frac)

    if int(NN_smooth) != defaults['NN_smooth']:
        expr += ' -NN_smooth {0}'.format(NN_smooth)

    if int(smooth_final) != defaults['smooth_final']:
        expr += ' -smooth_final {0}'.format(smooth_final)

    if push_to_edge:
        expr += ' -push_to_edge'

    if use_skull:
        expr += ' -use_skull'

    if float(perc_int) != defaults['perc_int']:
        expr += ' -perc_int {0}'.format(perc_int)

    if int(max_inter_iter) != defaults['max_inter_iter']:
        expr += ' -max_inter_iter {0}'.format(max_inter_iter)

    if float(blur_fwhm) != defaults['blur_fwhm']:
        expr += ' -blur_fwhm {0}'.format(blur_fwhm)

    if float(fac) != defaults['fac']:
        expr += ' -fac {0}'.format(fac)

    return expr


def mri_convert(in_file, reslice_like=None, out_file=None, args=None):
    """
    Method to convert files from mgz to nifti format

    Parameters
    ----------
    in_file : string
        A path of mgz input file
    args : string
        Arguments of mri_convert
    Returns
    -------
    out_file : string
        A path of nifti output file
    """

    import os

    if out_file is None:
        out_file = in_file.replace('.mgz','.nii.gz')

    cmd = 'mri_convert %s %s' % (in_file, out_file)

    if reslice_like is not None:
        cmd = cmd + ' -rl ' + reslice_like

    if args is not None:
        cmd = cmd + ' ' +args

    os.system(cmd)

    return out_file


def wb_command(in_file):

    import os

    out_file = in_file.replace('.nii.gz','_fill_holes.nii.gz')

    cmd = 'wb_command -volume-fill-holes %s %s' % (in_file, out_file)

    os.system(cmd)

    return out_file


def fslmaths_command(in_file, number, out_file_suffix):

    import os

    out_filename = in_file.replace('.nii.gz', out_file_suffix+'.nii.gz')

    out_file = os.path.join(os.getcwd(), out_filename[out_filename.rindex('/')+1:])

    cmd = 'fslmaths %s -div %f -mul 150 -abs %s' % (in_file, number, out_file)

    os.system(cmd)

    return out_file

def normalize_wmparc(source_file, target_file, xfm, out_file):
    from CPAC.utils.monitoring.custom_logging import log_subprocess
    import os

    cmd = ['mri_vol2vol', '--mov', source_file, \
                '--targ', target_file, '--o', out_file, '--lta', xfm]
    log_subprocess(cmd)
    output = os.path.join(os.getcwd(), out_file)
    return output

"""This module provides interfaces for workbench -volume-remove-islands commands"""
from nipype.interfaces.base import TraitedSpec, File, traits, CommandLineInputSpec
from nipype.interfaces.workbench.base import WBCommand
from nipype import logging

iflogger = logging.getLogger("nipype.interface")


class VolumeRemoveIslandsInputSpec(CommandLineInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        argstr="%s",
        position=0,
        desc="the input ROI volume",
    )
    out_file = File(
        name_source=["in_file"],
        name_template="%s_generated.nii.gz",
        argstr="%s",
        position=1,
        desc="the output ROI volume",
    )



class VolumeRemoveIslandsOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="the output ROI volume")


class VolumeRemoveIslands(WBCommand):
    """
    workbench
    -volume-remove-islands
    REMOVE ISLANDS FROM AN ROI VOLUME
    wb_command -volume-remove-islands
        <volume-in> - the input ROI volume
        <volume-out> - output - the output ROI volume

        Finds all face-connected parts of the ROI, and zeros out all but the
        largest one.

    """

    input_spec = VolumeRemoveIslandsInputSpec
    output_spec = VolumeRemoveIslandsOutputSpec
    _cmd = "wb_command -volume-remove-islands" 
    
