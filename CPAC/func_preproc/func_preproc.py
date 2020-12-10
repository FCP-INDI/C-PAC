from nipype import logging
from nipype.interfaces import ants

logger = logging.getLogger('workflow')

import nipype.pipeline.engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from nipype.interfaces import afni
from nipype.interfaces.afni import preprocess
from nipype.interfaces.afni import utils as afni_utils
from CPAC.func_preproc.utils import add_afni_prefix, nullify, chunk_ts, \
    split_ts_chunks, oned_text_concat, notch_filter_motion
from CPAC.utils.interfaces.function import Function
from CPAC.generate_motion_statistics import motion_power_statistics


def collect_arguments(*args):
    command_args = []
    if args[0]:
        command_args += [args[1]]
    command_args += args[2:]
    return ' '.join(command_args)


def anat_refined_mask(init_bold_mask = True, wf_name='init_bold_mask'):
               
    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['func',
                                                        'anatomical_brain_mask',
                                                        'anat_brain',
                                                        'init_func_brain_mask']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['func_brain_mask']),
                         name='outputspec')

    # 1 Take single volume of func 
    func_single_volume = pe.Node(interface=afni.Calc(),
                            name='func_single_volume')

    func_single_volume.inputs.set(
        expr='a',
        single_idx=1,
        outputtype='NIFTI_GZ'
    )

    wf.connect(input_node, 'func',
                func_single_volume, 'in_file_a')
    
    # 2 get temporary func brain 
    func_tmp_brain = pe.Node(interface=afni_utils.Calc(),
                            name='func_tmp_brain')
    func_tmp_brain.inputs.expr = 'a*b'
    func_tmp_brain.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(func_single_volume, 'out_file', 
                func_tmp_brain, 'in_file_a')

    # 2.1 get a tmp func brain mask
    if init_bold_mask == True :
        # 2.1.1 N4BiasFieldCorrection single volume of raw_func 
        func_single_volume_n4_corrected = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=200), shrink_factor=2, 
                                        name='func_single_volume_n4_corrected')
        func_single_volume_n4_corrected.inputs.args = '-r True'

        wf.connect(func_single_volume, 'out_file', 
                    func_single_volume_n4_corrected, 'input_image')
    
        # 2.1.2 bet n4 corrected image - generate tmp func brain mask
        func_tmp_brain_mask = pe.Node(interface=fsl.BET(),
                                        name='func_tmp_brain_mask_pre')
        func_tmp_brain_mask.inputs.mask = True

        wf.connect(func_single_volume_n4_corrected, 'output_image', 
                    func_tmp_brain_mask, 'in_file')

        # 2.1.3 dilate func tmp brain mask 
        func_tmp_brain_mask_dil = pe.Node(interface=fsl.ImageMaths(), 
                                                name='func_tmp_brain_mask_dil')
        func_tmp_brain_mask_dil.inputs.op_string = '-dilM' 

        wf.connect(func_tmp_brain_mask, 'mask_file',   
                    func_tmp_brain_mask_dil, 'in_file')
        
        wf.connect(func_tmp_brain_mask_dil, 'out_file',
                    func_tmp_brain, 'in_file_b')
    else :
        # 2.1.1 connect dilated init func brain mask 
        wf.connect(input_node, 'init_func_brain_mask', 
                    func_tmp_brain, 'in_file_b')
    
    # 3. get transformation of anat to func
    # 3.1 Register func tmp brain to anat brain to get func2anat matrix
    linear_reg_func_to_anat = pe.Node(interface=fsl.FLIRT(),
                        name='func_to_anat_linear_reg')
    linear_reg_func_to_anat.inputs.cost = 'mutualinfo'
    linear_reg_func_to_anat.inputs.dof = 6
    
    wf.connect(func_tmp_brain, 'out_file', 
                linear_reg_func_to_anat, 'in_file')

    wf.connect(input_node, 'anat_brain', 
                linear_reg_func_to_anat, 'reference')

    # 3.2 Inverse func to anat affine
    inv_func_to_anat_affine = pe.Node(interface=fsl.ConvertXFM(),
                            name='inv_func2anat_affine')
    inv_func_to_anat_affine.inputs.invert_xfm = True

    wf.connect(linear_reg_func_to_anat, 'out_matrix_file',
                inv_func_to_anat_affine, 'in_file')

    # 4. anat mask to func space
    # Transform anatomical mask to functional space to get BOLD mask
    reg_anat_mask_to_func = pe.Node(interface=fsl.FLIRT(),
                        name='reg_anat_mask_to_func')
    reg_anat_mask_to_func.inputs.apply_xfm = True
    reg_anat_mask_to_func.inputs.cost = 'mutualinfo'
    reg_anat_mask_to_func.inputs.dof = 6
    reg_anat_mask_to_func.inputs.interp = 'nearestneighbour'

    wf.connect(input_node, 'anatomical_brain_mask',
                reg_anat_mask_to_func, 'in_file')

    wf.connect(func_tmp_brain, 'out_file', 
                reg_anat_mask_to_func, 'reference')

    wf.connect(inv_func_to_anat_affine, 'out_file',
                reg_anat_mask_to_func, 'in_matrix_file')
                
    # 5. get final func mask: refine func tmp mask with anat_mask_in_func mask
    func_mask = pe.Node(interface=fsl.MultiImageMaths(), name='func_mask')
    func_mask.inputs.op_string = "-mul %s"

    wf.connect(reg_anat_mask_to_func, 'out_file', 
                func_mask, 'operand_files')
    
    if init_bold_mask == True :
        wf.connect(func_tmp_brain_mask_dil, 'out_file', 
                    func_mask, 'in_file')
    else :
        wf.connect(input_node, 'init_func_brain_mask',
                    func_mask, 'in_file')
    
    wf.connect(func_mask, 'out_file', 
                output_node, 'func_brain_mask')

    return wf

def anat_based_mask(wf_name='bold_mask'):
# reference DCAN lab BOLD mask
# https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased.sh
    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['func',
                                                        'anat_brain',
                                                        'anat_head']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['func_brain_mask']),
                         name='outputspec')

    # 0. Take single volume of func 
    func_single_volume = pe.Node(interface=afni.Calc(),
                            name='func_single_volume')

    func_single_volume.inputs.set(
        expr='a',
        single_idx=1,
        outputtype='NIFTI_GZ'
    )

    wf.connect(input_node, 'func',
                func_single_volume, 'in_file_a')

    # 1. Register func head to anat head to get func2anat matrix
    linear_reg_func_to_anat = pe.Node(interface=fsl.FLIRT(),
                        name='func_to_anat_linear_reg')
    linear_reg_func_to_anat.inputs.dof = 6
    linear_reg_func_to_anat.inputs.interp = 'spline'
    linear_reg_func_to_anat.inputs.searchr_x = [30, 30]
    linear_reg_func_to_anat.inputs.searchr_y = [30, 30]
    linear_reg_func_to_anat.inputs.searchr_z = [30, 30]

    wf.connect(func_single_volume, 'out_file', 
                linear_reg_func_to_anat, 'in_file')

    wf.connect(input_node, 'anat_head', 
                linear_reg_func_to_anat, 'reference')

    # 2. Inverse func to anat affine, to get anat-to-func transform 
    inv_func_to_anat_affine = pe.Node(interface=fsl.ConvertXFM(),
                            name='inv_func2anat_affine')
    inv_func_to_anat_affine.inputs.invert_xfm = True

    wf.connect(linear_reg_func_to_anat, 'out_matrix_file',
                inv_func_to_anat_affine, 'in_file')

    # 3. get BOLD mask
    # 3.1 Apply anat-to-func transform to transfer anatomical brain to functional space 
    reg_anat_brain_to_func = pe.Node(interface=fsl.ApplyWarp(),
                            name='reg_anat_brain_to_func')
    reg_anat_brain_to_func.inputs.interp = 'nn'
    reg_anat_brain_to_func.inputs.relwarp = True

    wf.connect(input_node, 'anat_brain', 
                reg_anat_brain_to_func, 'in_file')

    wf.connect(input_node, 'func', 
                reg_anat_brain_to_func, 'ref_file')

    wf.connect(inv_func_to_anat_affine, 'out_file', 
                reg_anat_brain_to_func, 'premat')

    # 3.2 Binarize transfered image
    func_mask_bin = pe.Node(interface=fsl.ImageMaths(),
                                name='func_mask_bin')
    func_mask_bin.inputs.op_string = '-abs -bin'

    wf.connect(reg_anat_brain_to_func, 'out_file', 
                func_mask_bin, 'in_file')

    # 3.3 Fill holes to get BOLD mask
    func_mask_fill_holes = pe.Node(interface=afni.MaskTool(),
                                name='func_mask_fill_holes')
    func_mask_fill_holes.inputs.fill_holes = True
    func_mask_fill_holes.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(func_mask_bin, 'out_file', 
                func_mask_fill_holes, 'in_file')

    wf.connect(func_mask_fill_holes, 'out_file', 
                output_node, 'func_brain_mask')

    return wf


def normalize_motion_parameters(in_file):
    """
    Convert FSL mcflirt motion parameters to AFNI space  
    """
    import os 
    import numpy as np

    motion_params = np.genfromtxt(in_file).T
    motion_params = np.vstack((motion_params[2,:]*180/np.pi,
                                motion_params[0,:]*180/np.pi,
                                -motion_params[1,:]*180/np.pi,
                                motion_params[5,:],
                                motion_params[3,:],
                                -motion_params[4,:]))
    motion_params = np.transpose(motion_params)

    out_file = os.path.join(os.getcwd(), 'motion_params.1D')
    np.savetxt(out_file, motion_params)                                

    return out_file


def skullstrip_functional(skullstrip_tool='afni', config=None, wf_name='skullstrip_functional'):

    skullstrip_tool = skullstrip_tool.lower()
    if skullstrip_tool != 'afni' and skullstrip_tool != 'fsl' and skullstrip_tool != 'fsl_afni' and skullstrip_tool != 'anatomical_refined' and skullstrip_tool != 'anatomical_based':
        raise Exception("\n\n[!] Error: The 'tool' parameter of the "
                        "'skullstrip_functional' workflow must be either "
                        "'afni' or 'fsl' or 'fsl_afni' or 'anatomical_refined' or 'anatomical_based'.\n\nTool input: "
                        "{0}\n\n".format(skullstrip_tool))
               
    wf = pe.Workflow(name=wf_name)

    input_node = pe.Node(util.IdentityInterface(fields=['raw_func',
                                                        'func',
                                                        'anatomical_brain_mask',
                                                        'anat_brain',
                                                        'anat_head']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['func_brain',
                                                         'func_brain_mask']),
                         name='outputspec')

    if skullstrip_tool == 'afni':
        func_get_brain_mask = pe.Node(interface=preprocess.Automask(),
                                      name='func_get_brain_mask_AFNI')
        func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

        wf.connect(input_node, 'func', func_get_brain_mask, 'in_file')

        wf.connect(func_get_brain_mask, 'out_file',
                   output_node, 'func_brain_mask')

    elif skullstrip_tool == 'fsl':
        inputnode_bet = pe.Node(
            util.IdentityInterface(fields=['frac',
                                            'mesh_boolean',
                                            'outline',
                                            'padding',
                                            'radius',
                                            'reduce_bias',
                                            'remove_eyes',
                                            'robust',
                                            'skull',
                                            'surfaces',
                                            'threshold',
                                            'vertical_gradient']),
            name='BET_options')        

        func_get_brain_mask = pe.Node(interface=fsl.BET(),
                                    name='func_get_brain_mask_BET')
        func_get_brain_mask.inputs.output_type = 'NIFTI_GZ'
        func_get_brain_mask.inputs.mask = True

        inputnode_bet.inputs.set(
                frac=config.bold_bet_frac, # 0.3
                mesh_boolean=config.bold_bet_mesh_boolean,
                outline=config.bold_bet_outline,
                padding=config.bold_bet_padding,
                radius=config.bold_bet_radius,
                reduce_bias=config.bold_bet_reduce_bias,
                remove_eyes=config.bold_bet_remove_eyes,
                robust=config.bold_bet_robust,
                skull=config.bold_bet_skull,
                surfaces=config.bold_bet_surfaces,
                threshold=config.bold_bet_threshold,
                vertical_gradient=config.bold_bet_vertical_gradient,
            )

        wf.connect([
            (inputnode_bet, func_get_brain_mask, [
                ('frac', 'frac'),
                ('mesh_boolean', 'mesh'),
                ('outline', 'outline'),
                ('padding', 'padding'),
                ('radius', 'radius'),
                ('reduce_bias', 'reduce_bias'),
                ('remove_eyes', 'remove_eyes'),
                ('robust', 'robust'),
                ('skull', 'skull'),
                ('surfaces', 'surfaces'),
                ('threshold', 'threshold'),
                ('vertical_gradient', 'vertical_gradient'),
            ])
        ])

        if config.bold_bet_functional_mean_boolean : 
            func_skull_mean = pe.Node(interface=afni_utils.TStat(),
                                        name='func_mean_skull_{0}'.format(wf_name))
            func_skull_mean.inputs.options = '-mean'
            func_skull_mean.inputs.outputtype = 'NIFTI_GZ'
            
            wf.connect(input_node, 'func', 
                        func_skull_mean, 'in_file')
            wf.connect(func_skull_mean, 'out_file', 
                        func_get_brain_mask, 'in_file')

        else:
            func_get_brain_mask.inputs.functional = True
            wf.connect(input_node, 'func', 
                        func_get_brain_mask, 'in_file')

        # erode one voxel of functional brian mask
        erode_one_voxel = pe.Node(interface=fsl.ErodeImage(),
                                  name='erode_one_voxel')

        erode_one_voxel.inputs.kernel_shape = 'box'
        erode_one_voxel.inputs.kernel_size = 1.0

        wf.connect(func_get_brain_mask, 'mask_file',
                   erode_one_voxel, 'in_file')

        wf.connect(erode_one_voxel, 'out_file',
                   output_node, 'func_brain_mask')

    elif skullstrip_tool == 'fsl_afni':
        func_skull_mean = pe.Node(interface=afni_utils.TStat(),
                                    name='func_mean_skull')
        func_skull_mean.inputs.options = '-mean'
        func_skull_mean.inputs.outputtype = 'NIFTI_GZ'

        skullstrip_first_pass = pe.Node(fsl.BET(frac=0.2, mask=True, functional=False), name='skullstrip_first_pass')
        bet_dilate = pe.Node(fsl.DilateImage(operation='max', kernel_shape='sphere', kernel_size=6.0, internal_datatype='char'), name='skullstrip_first_dilate')                                                  
        bet_mask = pe.Node(fsl.ApplyMask(), name='skullstrip_first_mask')
        unifize = pe.Node(afni_utils.Unifize(t2=True, outputtype='NIFTI_GZ', args='-clfrac 0.2 -rbt 18.3 65.0 90.0', out_file="uni.nii.gz"), name='unifize')
        skullstrip_second_pass = pe.Node(preprocess.Automask(dilate=1, outputtype='NIFTI_GZ'), name='skullstrip_second_pass')
        combine_masks = pe.Node(fsl.BinaryMaths(operation='mul'), name='combine_masks')

        wf.connect([(input_node, func_skull_mean, [('func', 'in_file')]),
                        (func_skull_mean, skullstrip_first_pass, [('out_file', 'in_file')]),
                        (skullstrip_first_pass, bet_dilate, [('mask_file', 'in_file')]),
                        (bet_dilate, bet_mask, [('out_file', 'mask_file')]),
                        (skullstrip_first_pass, bet_mask, [('out_file' , 'in_file')]),
                        (bet_mask, unifize, [('out_file', 'in_file')]),
                        (unifize, skullstrip_second_pass, [('out_file', 'in_file')]),
                        (skullstrip_first_pass, combine_masks, [('mask_file', 'in_file')]),
                        (skullstrip_second_pass, combine_masks, [('out_file', 'operand_file')]),
                        (combine_masks, output_node, [('out_file', 'func_brain_mask')])])
    
    # Refine functional mask by registering anatomical mask to functional space
    elif skullstrip_tool == 'anatomical_refined':
        
        # binarize anat mask, in case of it is not a binary mask. 
        anat_brain_mask_bin = pe.Node(interface=fsl.ImageMaths(),
                                    name='anat_brain_mask_bin')
        anat_brain_mask_bin.inputs.op_string = '-bin'

        wf.connect(input_node, 'anatomical_brain_mask',
                    anat_brain_mask_bin, 'in_file')

        # fill holes of anat mask 
        anat_mask_filled = pe.Node(interface=afni.MaskTool(),
                        name='anat_brain_mask_filled')
        anat_mask_filled.inputs.fill_holes = True
        anat_mask_filled.inputs.outputtype = 'NIFTI_GZ'

        wf.connect(anat_brain_mask_bin, 'out_file', 
                    anat_mask_filled, 'in_file')

        # init_bold_mask : input raw func
        init_bold_mask = anat_refined_mask(init_bold_mask = True, wf_name='init_bold_mask')

        func_deoblique = pe.Node(interface=afni_utils.Refit(),
                                name='raw_func_deoblique')
        func_deoblique.inputs.deoblique = True

        wf.connect(input_node, 'raw_func',
                    func_deoblique, 'in_file')

        func_reorient = pe.Node(interface=afni_utils.Resample(),
                                name='raw_func_reorient')

        func_reorient.inputs.orientation = 'RPI'
        func_reorient.inputs.outputtype = 'NIFTI_GZ'

        wf.connect(func_deoblique, 'out_file',
                    func_reorient, 'in_file')
                    
        wf.connect(func_reorient, 'out_file',
                    init_bold_mask, 'inputspec.func')

        wf.connect(anat_mask_filled, 'out_file',
                    init_bold_mask, 'inputspec.anatomical_brain_mask')

        wf.connect(input_node, 'anat_brain',
                    init_bold_mask, 'inputspec.anat_brain')
        
        # dilate init func brain mask
        func_tmp_brain_mask = pe.Node(interface=fsl.ImageMaths(),
                        name='func_tmp_brain_mask_dil')
        func_tmp_brain_mask.inputs.op_string = '-dilM' 
        
        wf.connect(init_bold_mask, 'outputspec.func_brain_mask',
                    func_tmp_brain_mask, 'in_file')

        # refined_bold_mask : input motion corrected func 
        refined_bold_mask = anat_refined_mask(init_bold_mask = False, wf_name='refined_bold_mask')

        wf.connect(input_node, 'func',
                    refined_bold_mask, 'inputspec.func')

        wf.connect(input_node, 'anat_brain', 
                    refined_bold_mask, 'inputspec.anat_brain')

        wf.connect(func_tmp_brain_mask, 'out_file', 
                    refined_bold_mask, 'inputspec.init_func_brain_mask')

        # Dialate anatomical mask, if 'anatomical_mask_dilation : True' in config file
        if config.anatomical_mask_dilation :
            anat_mask_dilate = pe.Node(interface=afni.MaskTool(),
                            name='anat_mask_dilate')
            anat_mask_dilate.inputs.dilate_inputs = '1'
            anat_mask_dilate.inputs.outputtype = 'NIFTI_GZ'

            wf.connect(anat_mask_filled, 'out_file', 
                        anat_mask_dilate, 'in_file')
            wf.connect(anat_mask_dilate, 'out_file',
                        refined_bold_mask, 'inputspec.anatomical_brain_mask')

        else : 
            wf.connect(anat_mask_filled, 'out_file',
                        refined_bold_mask, 'inputspec.anatomical_brain_mask')

        # get final func mask
        func_mask_final = pe.Node(interface=fsl.MultiImageMaths(), name='func_mask_final')
        func_mask_final.inputs.op_string = "-mul %s"

        wf.connect(func_tmp_brain_mask, 'out_file', 
                    func_mask_final, 'in_file')

        wf.connect(refined_bold_mask, 'outputspec.func_brain_mask', 
                    func_mask_final, 'operand_files')

        wf.connect(func_mask_final, 'out_file', 
                    output_node, 'func_brain_mask')
    
    # simply transforms the anatomical brain mask to the functional BOLD data.
    elif skullstrip_tool == 'anatomical_based':

        # refined_bold_mask : input motion corrected func 
        bold_mask = anat_based_mask(wf_name='bold_mask')

        wf.connect(input_node, 'func',
                    bold_mask, 'inputspec.func')

        wf.connect(input_node, 'anat_brain', 
                    bold_mask, 'inputspec.anat_brain')

        wf.connect(input_node, 'anat_head', 
                    bold_mask, 'inputspec.anat_head')

        wf.connect(bold_mask, 'outputspec.func_brain_mask', 
                    output_node, 'func_brain_mask')

    func_edge_detect = pe.Node(interface=afni_utils.Calc(),
                               name='func_extract_brain')

    func_edge_detect.inputs.expr = 'a*b'
    func_edge_detect.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(input_node, 'func', func_edge_detect, 'in_file_a')

    if skullstrip_tool == 'afni':
        wf.connect(func_get_brain_mask, 'out_file',
                   func_edge_detect, 'in_file_b')
    elif skullstrip_tool == 'fsl':
        wf.connect(erode_one_voxel, 'out_file',
                   func_edge_detect, 'in_file_b')
    elif skullstrip_tool == 'fsl_afni':
        wf.connect(combine_masks, 'out_file',
                    func_edge_detect, 'in_file_b')
    elif skullstrip_tool == 'anatomical_refined':
        wf.connect(func_mask_final, 'out_file',
                    func_edge_detect, 'in_file_b')
    elif skullstrip_tool == 'anatomical_based':
        wf.connect(bold_mask, 'outputspec.func_brain_mask', 
                    func_edge_detect, 'in_file_b')

    wf.connect(func_edge_detect, 'out_file',  output_node, 'func_brain')

    return wf

# TODO
def create_scale_func_wf(scaling_factor, wf_name='scale_func'):

    """Workflow to scale func data.
    Parameters
    ----------
        scaling_factor : float
            Scale the size of the dataset voxels by the factor. 
        wf_name : string
            name of the workflow
    
    Workflow Inputs::
        inputspec.func : func file or a list of func/rest nifti file
            User input functional(T2*) Image
    Workflow Outputs::
        outputspec.scaled_func : string (nifti file)
            Path to Output image with scaled data
    Order of commands:
    - Scale the size of the dataset voxels by the factor 'fac'. For details see `3dcalc <https://afni.nimh.nih.gov/pub/dist/doc/program_help/3drefit.html>`_::
        3drefit -xyzscale fac rest.nii.gz
    """

    # allocate a workflow object
    preproc = pe.Workflow(name=wf_name)

    # configure the workflow's input spec
    inputNode = pe.Node(util.IdentityInterface(fields=['func']),
                        name='inputspec')

    # configure the workflow's output spec
    outputNode = pe.Node(util.IdentityInterface(fields=['scaled_func']),
                         name='outputspec')

    # allocate a node to edit the functional file
    func_scale = pe.Node(interface=afni_utils.Refit(),
                            name='func_scale')

    func_scale.inputs.xyzscale = scaling_factor

    # wire in the func_get_idx node
    preproc.connect(inputNode, 'func',
                    func_scale, 'in_file')

    # wire the output
    preproc.connect(func_scale, 'out_file',
                    outputNode, 'scaled_func')

    return preproc

    
def create_wf_edit_func(wf_name="edit_func"):
    """Workflow to edit the scan to the proscribed TRs.
    
    Workflow Inputs::

        inputspec.func : func file or a list of func/rest nifti file
            User input functional(T2*) Image

        inputspec.start_idx : string
            Starting volume/slice of the functional image (optional)

        inputspec.stop_idx : string
            Last volume/slice of the functional image (optional)

    Workflow Outputs::

        outputspec.edited_func : string (nifti file)
            Path to Output image with the initial few slices dropped


    Order of commands:

    - Get the start and the end volume index of the functional run. If not defined by the user, return the first and last volume.

        get_idx(in_files, stop_idx, start_idx)

    - Dropping the initial TRs. For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::

        3dcalc -a rest.nii.gz[4..299]
               -expr 'a'
               -prefix rest_3dc.nii.gz

    """

    # allocate a workflow object
    preproc = pe.Workflow(name=wf_name)

    # configure the workflow's input spec
    inputNode = pe.Node(util.IdentityInterface(fields=['func',
                                                       'start_idx',
                                                       'stop_idx']),
                        name='inputspec')

    # configure the workflow's output spec
    outputNode = pe.Node(util.IdentityInterface(fields=['edited_func']),
                         name='outputspec')

    # allocate a node to check that the requested edits are
    # reasonable given the data
    func_get_idx = pe.Node(util.Function(input_names=['in_files',
                                                      'stop_idx',
                                                      'start_idx'],
                                         output_names=['stopidx',
                                                       'startidx'],
                                         function=get_idx),
                           name='func_get_idx')

    # wire in the func_get_idx node
    preproc.connect(inputNode, 'func',
                    func_get_idx, 'in_files')
    preproc.connect(inputNode, 'start_idx',
                    func_get_idx, 'start_idx')
    preproc.connect(inputNode, 'stop_idx',
                    func_get_idx, 'stop_idx')

    # allocate a node to edit the functional file
    func_drop_trs = pe.Node(interface=afni_utils.Calc(),
                            name='func_drop_trs')

    func_drop_trs.inputs.expr = 'a'
    func_drop_trs.inputs.outputtype = 'NIFTI_GZ'

    # wire in the inputs
    preproc.connect(inputNode, 'func',
                    func_drop_trs, 'in_file_a')

    preproc.connect(func_get_idx, 'startidx',
                    func_drop_trs, 'start_idx')

    preproc.connect(func_get_idx, 'stopidx',
                    func_drop_trs, 'stop_idx')

    # wire the output
    preproc.connect(func_drop_trs, 'out_file',
                    outputNode, 'edited_func')

    return preproc


# functional preprocessing
def create_func_preproc(skullstrip_tool, motion_correct_tool,
                        motion_correct_ref, motion_estimate_filter_option=False,
                        config=None, wf_name='func_preproc'):
    """

    The main purpose of this workflow is to process functional data. Raw rest file is deobliqued and reoriented
    into RPI. Then take the mean intensity values over all time points for each voxel and use this image
    to calculate motion parameters. The image is then skullstripped, normalized and a processed mask is
    obtained to use it further in Image analysis.

    Parameters
    ----------

    wf_name : string
        Workflow name

    Returns
    -------
    func_preproc : workflow object
        Functional Preprocessing workflow object

    Notes
    -----

    `Source <https://github.com/FCP-INDI/C-PAC/blob/master/CPAC/func_preproc/func_preproc.py>`_

    Workflow Inputs::

        inputspec.func : func nifti file
            User input functional(T2) Image, in any of the 8 orientations

        inputspec.twopass : boolean
            Perform two-pass on volume registration

    Workflow Outputs::

        outputspec.refit : string (nifti file)
            Path to deobliqued anatomical data

        outputspec.reorient : string (nifti file)
            Path to RPI oriented anatomical data

        outputspec.motion_correct_ref : string (nifti file)
             Path to Mean intensity Motion corrected image
             (base reference image for the second motion correction run)

        outputspec.motion_correct : string (nifti file)
            Path to motion corrected output file

        outputspec.max_displacement : string (Mat file)
            Path to maximum displacement (in mm) for brain voxels in each volume

        outputspec.movement_parameters : string (Mat file)
            Path to 1D file containing six movement/motion parameters(3 Translation, 3 Rotations)
            in different columns (roll pitch yaw dS  dL  dP)

        outputspec.skullstrip : string (nifti file)
            Path to skull stripped Motion Corrected Image

        outputspec.mask : string (nifti file)
            Path to brain-only mask

        outputspec.func_mean : string (nifti file)
            Mean, Skull Stripped, Motion Corrected output T2 Image path
            (Image with mean intensity values across voxels)

        outputpsec.preprocessed : string (nifti file)
            output skull stripped, motion corrected T2 image
            with normalized intensity values

        outputspec.preprocessed_mask : string (nifti file)
           Mask obtained from normalized preprocessed image

    Order of commands:

    - Deobliqing the scans.  For details see `3drefit <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3drefit.html>`_::

        3drefit -deoblique rest_3dc.nii.gz

    - Re-orienting the Image into Right-to-Left Posterior-to-Anterior Inferior-to-Superior (RPI) orientation. For details see `3dresample <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dresample.html>`_::

        3dresample -orient RPI
                   -prefix rest_3dc_RPI.nii.gz
                   -inset rest_3dc.nii.gz

    - Calculate voxel wise statistics. Get the RPI Image with mean intensity values over all timepoints for each voxel. For details see `3dTstat <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_::

        3dTstat -mean
                -prefix rest_3dc_RPI_3dT.nii.gz
                rest_3dc_RPI.nii.gz

    - Motion Correction. For details see `3dvolreg <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dvolreg.html>`_::

        3dvolreg -Fourier
                 -twopass
                 -base rest_3dc_RPI_3dT.nii.gz/
                 -zpad 4
                 -maxdisp1D rest_3dc_RPI_3dvmd1D.1D
                 -1Dfile rest_3dc_RPI_3dv1D.1D
                 -prefix rest_3dc_RPI_3dv.nii.gz
                 rest_3dc_RPI.nii.gz

      The base image or the reference image is the mean intensity RPI image obtained in the above the step.For each volume
      in RPI-oriented T2 image, the command, aligns the image with the base mean image and calculates the motion, displacement
      and movement parameters. It also outputs the aligned 4D volume and movement and displacement parameters for each volume.

    - Calculate voxel wise statistics. Get the motion corrected output Image from the above step, with mean intensity values over all timepoints for each voxel.
      For details see `3dTstat <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_::

        3dTstat -mean
                -prefix rest_3dc_RPI_3dv_3dT.nii.gz
                rest_3dc_RPI_3dv.nii.gz

    - Motion Correction and get motion, movement and displacement parameters. For details see `3dvolreg <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dvolreg.html>`_::

        3dvolreg -Fourier
                 -twopass
                 -base rest_3dc_RPI_3dv_3dT.nii.gz
                 -zpad 4
                 -maxdisp1D rest_3dc_RPI_3dvmd1D.1D
                 -1Dfile rest_3dc_RPI_3dv1D.1D
                 -prefix rest_3dc_RPI_3dv.nii.gz
                 rest_3dc_RPI.nii.gz

      The base image or the reference image is the mean intensity motion corrected image obtained from the above the step (first 3dvolreg run).
      For each volume in RPI-oriented T2 image, the command, aligns the image with the base mean image and calculates the motion, displacement
      and movement parameters. It also outputs the aligned 4D volume and movement and displacement parameters for each volume.

    - Create a brain-only mask. For details see `3dautomask <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dAutomask.html>`_::

        3dAutomask
                   -prefix rest_3dc_RPI_3dv_automask.nii.gz
                   rest_3dc_RPI_3dv.nii.gz

    - Edge Detect(remove skull) and get the brain only. For details see `3dcalc <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dcalc.html>`_::

        3dcalc -a rest_3dc_RPI_3dv.nii.gz
               -b rest_3dc_RPI_3dv_automask.nii.gz
               -expr 'a*b'
               -prefix rest_3dc_RPI_3dv_3dc.nii.gz

    - Normalizing the image intensity values. For details see `fslmaths <http://www.fmrib.ox.ac.uk/fsl/avwutils/index.html>`_::

        fslmaths rest_3dc_RPI_3dv_3dc.nii.gz
                 -ing 10000 rest_3dc_RPI_3dv_3dc_maths.nii.gz
                 -odt float

      Normalized intensity = (TrueValue*10000)/global4Dmean

    - Calculate mean of skull stripped image. For details see `3dTstat <http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTstat.html>`_::

        3dTstat -mean -prefix rest_3dc_RPI_3dv_3dc_3dT.nii.gz rest_3dc_RPI_3dv_3dc.nii.gz

    - Create Mask (Generate mask from Normalized data). For details see `fslmaths <http://www.fmrib.ox.ac.uk/fsl/avwutils/index.html>`_::

        fslmaths rest_3dc_RPI_3dv_3dc_maths.nii.gz
               -Tmin -bin rest_3dc_RPI_3dv_3dc_maths_maths.nii.gz
               -odt char

    .. exec::
        import yaml
        from urllib.request import urlopen
        from CPAC.func_preproc import create_func_preproc
        from CPAC.utils.configuration import Configuration

        wf = create_func_preproc('fsl', '3dvolreg', 'mean',
            config=Configuration(yaml.safe_load(urlopen(
                'https://raw.githubusercontent.com/FCP-INDI/C-PAC/develop/'
                'dev/docker_data/default_pipeline.yml'
            )))
        )
        wf.write_graph(
            graph2use='orig',
            dotfilename='./images/generated/func_preproc.dot'
        )

    High Level Workflow Graph:

    .. image:: ../../images/generated/func_preproc.png
       :width: 1000

    Detailed Workflow Graph:

    .. image:: ../../images/generated/func_preproc_detailed.png
       :width: 1000

    Examples
    --------

    >>> import func_preproc
    >>> preproc = create_func_preproc(bet=True)
    >>> preproc.inputs.inputspec.func='sub1/func/rest.nii.gz'
    >>> preproc.run() #doctest: +SKIP


    >>> import func_preproc
    >>> preproc = create_func_preproc(bet=False)
    >>> preproc.inputs.inputspec.func='sub1/func/rest.nii.gz'
    >>> preproc.run() #doctest: +SKIP

    """

    if motion_correct_ref != 'mean' and motion_correct_ref != 'median' and motion_correct_ref != 'selected_volume':
        raise Exception("\n\n[!] Error: The 'tool' parameter of the "
                        "'motion_correction_reference' workflow must be either "
                        "'mean' or 'median' or 'selected volume'.\n\nTool input: "
                        "{0}\n\n".format(motion_correct_ref))

    if motion_correct_tool != '3dvolreg' and motion_correct_tool != 'mcflirt':
        raise Exception("\n\n[!] Error: The 'tool' parameter of the "
                        "'motion_correction' workflow must be either "
                        "'3dvolreg' or 'mcflirt'.\n\nTool input: "
                        "{0}\n\n".format(motion_correct_tool))

    preproc = pe.Workflow(name=wf_name)
    input_node = pe.Node(util.IdentityInterface(fields=['raw_func',
                                                        'func',
                                                        'twopass',
                                                        'anatomical_brain_mask',
                                                        'anat_brain',
                                                        'anat_head',
                                                        'TR']),
                         name='inputspec')

    output_node = pe.Node(util.IdentityInterface(fields=['refit',
                                                         'reorient',
                                                         'reorient_mean',
                                                         'motion_correct',
                                                         'motion_correct_ref',
                                                         'motion_correct_median',
                                                         'movement_parameters',
                                                         'max_displacement',
                                                         'mask',
                                                         'skullstrip',
                                                         'func_mean',
                                                         'preprocessed',
                                                         'preprocessed_median',
                                                         'preprocessed_mask',
                                                         'slice_time_corrected',
                                                         'transform_matrices',
                                                         'center_of_mass',
                                                         'motion_filter_info',
                                                         'motion_filter_plot']),
                          name='outputspec')

    func_deoblique = pe.Node(interface=afni_utils.Refit(),
                             name='func_deoblique')
    func_deoblique.inputs.deoblique = True

    preproc.connect(input_node, 'func',
                    func_deoblique, 'in_file')

    func_reorient = pe.Node(interface=afni_utils.Resample(),
                            name='func_reorient')

    func_reorient.inputs.orientation = 'RPI'
    func_reorient.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(func_deoblique, 'out_file',
                    func_reorient, 'in_file')

    preproc.connect(func_reorient, 'out_file', output_node, 'reorient')

    if config:
        if int(config.maxCoresPerParticipant) > 1:
            chunk_imports = ['import nibabel as nb']
            chunk = pe.Node(Function(input_names=['func_file',
                                                  'n_cpus'],
                                     output_names=['TR_ranges'],
                                     function=chunk_ts,
                                     imports=chunk_imports),
                            name='chunk')

            chunk.inputs.n_cpus = int(config.maxCoresPerParticipant)
            preproc.connect(func_reorient, 'out_file', chunk, 'func_file')

            split_imports = ['import os', 'import subprocess']
            split = pe.Node(Function(input_names=['func_file',
                                                  'tr_ranges'],
                                     output_names=['split_funcs'],
                                     function=split_ts_chunks,
                                     imports=split_imports),
                            name='split')

            preproc.connect(func_reorient, 'out_file', split, 'func_file')
            preproc.connect(chunk, 'TR_ranges', split, 'tr_ranges')

            out_split_func = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_split_func')

            preproc.connect(split, 'split_funcs', out_split_func, 'out_file')

            func_motion_correct = pe.MapNode(interface=preprocess.Volreg(),
                                             name='func_generate_ref',
                                             iterfield=['in_file'])

            preproc.connect(out_split_func, 'out_file',
                            func_motion_correct, 'in_file')

            func_concat = pe.Node(interface=afni_utils.TCat(),
                                  name='func_concat')
            func_concat.inputs.outputtype = 'NIFTI_GZ'

            preproc.connect(func_motion_correct, 'out_file',
                            func_concat, 'in_files')

            out_motion = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_motion')

            preproc.connect(func_concat, 'out_file', out_motion, 'out_file')

        else:
            out_split_func = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_split_func')
            preproc.connect(func_reorient, 'out_file', out_split_func, 'out_file')

            func_motion_correct = pe.Node(interface=preprocess.Volreg(),
                                          name='func_generate_ref')

            preproc.connect(out_split_func, 'out_file',
                            func_motion_correct, 'in_file')

            out_motion = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_motion')

            preproc.connect(func_motion_correct, 'out_file',
                            out_motion, 'out_file')

    else:
        out_split_func = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_split_func')
        preproc.connect(func_reorient, 'out_file', out_split_func, 'out_file')

        func_motion_correct = pe.Node(interface=preprocess.Volreg(),
                                      name='func_generate_ref')

        preproc.connect(out_split_func, 'out_file',
                        func_motion_correct, 'in_file')

        out_motion = pe.Node(
            interface=util.IdentityInterface(fields=['out_file']),
            name='out_motion')

        preproc.connect(func_motion_correct, 'out_file',
                        out_motion, 'out_file')

    func_motion_correct.inputs.zpad = 4
    func_motion_correct.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect([(input_node, func_motion_correct, [(('twopass',
                                                          collect_arguments,
                                                          '-twopass',
                                                          '-Fourier'),
                                                         'args')]), ])

    # Calculate motion correction reference
    if motion_correct_ref == 'mean': 
        func_get_mean_RPI = pe.Node(interface=afni_utils.TStat(),
                                    name='func_get_mean_RPI')

        func_get_mean_RPI.inputs.options = '-mean'
        func_get_mean_RPI.inputs.outputtype = 'NIFTI_GZ'

        preproc.connect(func_reorient, 'out_file',
                        func_get_mean_RPI, 'in_file')

        preproc.connect(func_get_mean_RPI, 'out_file',
                        func_motion_correct, 'basefile')

        func_motion_correct_ref = func_get_mean_RPI.clone('func_get_mean_motion')

        preproc.connect(out_motion, 'out_file',
                        func_motion_correct_ref, 'in_file')

    elif motion_correct_ref == 'median':
        func_get_median_RPI = pe.Node(interface=afni_utils.TStat(),
                                      name='func_get_median_RPI')

        func_get_median_RPI.inputs.options = '-median'
        func_get_median_RPI.inputs.outputtype = 'NIFTI_GZ'

        preproc.connect(func_reorient, 'out_file',
                        func_get_median_RPI, 'in_file')

        preproc.connect(func_get_median_RPI, 'out_file',
                        func_motion_correct, 'basefile')

        func_motion_correct_ref = func_get_median_RPI.clone('func_get_median_motion')

        preproc.connect(out_motion, 'out_file',
                        func_motion_correct_ref, 'in_file')

    elif motion_correct_ref == 'selected_volume':
        func_get_selected_RPI = pe.Node(interface=afni.Calc(),
                                          name='func_get_selected_RPI')

        func_get_selected_RPI.inputs.set(
            expr='a',
            single_idx=config.motion_correction_reference_volume,
            outputtype='NIFTI_GZ'
        )

        preproc.connect(func_reorient, 'out_file',
                        func_get_selected_RPI, 'in_file_a')

        preproc.connect(func_get_selected_RPI, 'out_file',
                        func_motion_correct, 'basefile')

        func_motion_correct_ref = func_get_selected_RPI.clone('func_get_selected_motion')

        preproc.connect(out_motion, 'out_file',
                        func_motion_correct_ref, 'in_file_a')

    preproc.connect(func_motion_correct_ref, 'out_file',
                    output_node, 'motion_correct_ref')                 

    # Calculate motion parameters
    if motion_correct_tool == '3dvolreg':
        func_motion_correct_A = func_motion_correct.clone('func_motion_correct_3dvolreg')
        func_motion_correct_A.inputs.md1d_file = 'max_displacement.1D'

        preproc.connect([
            (
                input_node, func_motion_correct_A, [
                    (
                        ('twopass', collect_arguments, '-twopass', '-Fourier'),
                        'args'
                    )]
            ),
        ])

        preproc.connect(out_split_func, 'out_file',
                        func_motion_correct_A, 'in_file')
        preproc.connect(func_motion_correct_ref, 'out_file',
                        func_motion_correct_A, 'basefile')

        if config:
            if int(config.maxCoresPerParticipant) > 1:
                motion_concat = pe.Node(interface=afni_utils.TCat(),
                                        name='motion_concat')
                motion_concat.inputs.outputtype = 'NIFTI_GZ'

                preproc.connect(func_motion_correct_A, 'out_file',
                                motion_concat, 'in_files')

                out_motion_A = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_motion_A')

                preproc.connect(motion_concat, 'out_file', 
                                out_motion_A, 'out_file')

                concat_imports = ['import os']
                md1d_concat = pe.Node(Function(input_names=['in_files'],
                                               output_names=['out_file'],
                                               function=oned_text_concat,
                                               imports=concat_imports),
                                      name='md1d_concat')

                preproc.connect(func_motion_correct_A, 'md1d_file',
                                md1d_concat, 'in_files')

                oned_concat = pe.Node(Function(input_names=['in_files'],
                                               output_names=['out_file'],
                                               function=oned_text_concat,
                                               imports=concat_imports),
                                      name='oned_concat')

                preproc.connect(func_motion_correct_A, 'oned_file',
                                oned_concat, 'in_files')

                oned_matrix_concat = pe.Node(Function(input_names=['in_files'],
                                                      output_names=['out_file'],
                                                      function=oned_text_concat,
                                                      imports=concat_imports),
                                             name='oned_matrix_concat')

                preproc.connect(func_motion_correct_A, 'oned_matrix_save',
                                oned_matrix_concat, 'in_files')

                out_md1d = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_md1d')

                preproc.connect(md1d_concat, 'out_file', 
                                out_md1d, 'out_file')

                out_oned = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_oned')

                preproc.connect(oned_concat, 'out_file', 
                                out_oned, 'out_file')

                out_oned_matrix = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_oned_matrix')

                preproc.connect(oned_matrix_concat, 'out_file', 
                                out_oned_matrix, 'out_file')

            else:
                out_motion_A = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_motion_A')

                preproc.connect(func_motion_correct_A, 'out_file', 
                                out_motion_A, 'out_file')

                out_md1d = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_md1d')

                preproc.connect(func_motion_correct_A, 'md1d_file', 
                                out_md1d, 'out_file')

                out_oned = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_oned')

                preproc.connect(func_motion_correct_A, 'oned_file', 
                                out_oned, 'out_file')

                out_oned_matrix = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_oned_matrix')

                preproc.connect(func_motion_correct_A, 'oned_matrix_save', 
                                out_oned_matrix, 'out_file')

        else:
            out_motion_A = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_motion_A')

            preproc.connect(func_motion_correct_A, 'out_file', 
                            out_motion_A, 'out_file')

            out_md1d = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_md1d')

            preproc.connect(func_motion_correct_A, 'md1d_file', 
                            out_md1d, 'out_file')

            out_oned = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_oned')

            preproc.connect(func_motion_correct_A, 'oned_file', 
                            out_oned, 'out_file')

            out_oned_matrix = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name='out_oned_matrix')

            preproc.connect(func_motion_correct_A, 'oned_matrix_save', 
                            out_oned_matrix, 'out_file')

        if config:
            if motion_estimate_filter_option:
                notch_imports = ['import os', 'import numpy as np',
                                 'from scipy.signal import iirnotch, lfilter, firwin, freqz',
                                 'from matplotlib import pyplot as plt',
                                 'from CPAC.func_preproc.utils import degrees_to_mm, mm_to_degrees']
                notch = pe.Node(Function(input_names=['motion_params',
                                                      'filter_type',
                                                      'TR',
                                                      'fc_RR_min',
                                                      'fc_RR_max',
                                                      'center_freq',
                                                      'freq_bw',
                                                      'lowpass_cutoff',
                                                      'filter_order'],
                                         output_names=['filtered_motion_params',
                                                       'filter_info',
                                                       'filter_plot'],
                                         function=notch_filter_motion,
                                         imports=notch_imports),
                                name='filter_motion_params')

                notch.inputs.filter_type = config.motion_estimate_filter['filter_type']
                notch.inputs.fc_RR_min = config.motion_estimate_filter['breathing_rate_min']
                notch.inputs.fc_RR_max = config.motion_estimate_filter['breathing_rate_max']
                notch.inputs.center_freq = config.motion_estimate_filter['center_frequency']
                notch.inputs.freq_bw = config.motion_estimate_filter['filter_bandwidth']
                notch.inputs.lowpass_cutoff = config.motion_estimate_filter['lowpass_cutoff']
                notch.inputs.filter_order = config.motion_estimate_filter['filter_order']

                preproc.connect(notch, 'filter_info',
                                output_node, 'motion_filter_info')

                preproc.connect(notch, 'filter_plot',
                                output_node, 'motion_filter_plot')

                preproc.connect(out_oned, 'out_file', notch, 'motion_params')
                preproc.connect(input_node, 'TR', notch, 'TR')

                out_oned = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name='out_oned_notch')

                preproc.connect(notch, 'filtered_motion_params', out_oned, 'out_file')

        preproc.connect(out_motion_A, 'out_file',
                        output_node, 'motion_correct')
        preproc.connect(out_md1d, 'out_file',
                        output_node, 'max_displacement')
        preproc.connect(out_oned, 'out_file',
                        output_node, 'movement_parameters')
        preproc.connect(out_oned_matrix, 'out_file',
                        output_node, 'transform_matrices')

        skullstrip_func = skullstrip_functional(skullstrip_tool=skullstrip_tool, config=config, 
                                                wf_name="{0}_skullstrip".format(wf_name))
        
        preproc.connect(out_motion_A, 'out_file',
                        skullstrip_func, 'inputspec.func')

    elif motion_correct_tool == 'mcflirt':
        func_motion_correct_A = pe.Node(interface=fsl.MCFLIRT(save_mats=True, save_plots=True),
                                    name='func_motion_correct_mcflirt')
        
        func_motion_correct_A.inputs.save_mats = True
        func_motion_correct_A.inputs.save_plots = True
        func_motion_correct_A.inputs.save_rms = True

        preproc.connect(func_reorient, 'out_file',
                        func_motion_correct_A, 'in_file')

        preproc.connect(func_motion_correct_ref, 'out_file',
                        func_motion_correct_A, 'ref_file')

        preproc.connect(func_motion_correct_A, 'out_file',
                        output_node, 'motion_correct')

        skullstrip_func = skullstrip_functional(skullstrip_tool=skullstrip_tool, config=config,
                                                 wf_name="{0}_skullstrip".format(wf_name))

        preproc.connect(func_motion_correct_A, 'out_file',
                        skullstrip_func, 'inputspec.func')
        
        normalize_motion_params = pe.Node(Function(input_names=['in_file'],
                                     output_names=['out_file'],
                                     function=normalize_motion_parameters),
                            name='norm_motion_params')

        preproc.connect(func_motion_correct_A, 'par_file',
                        normalize_motion_params, 'in_file')

        if config:
            if motion_estimate_filter_option:
                notch_imports = ['import os', 'import numpy as np',
                                 'from scipy.signal import iirnotch, lfilter, firwin, freqz',
                                 'from matplotlib import pyplot as plt',
                                 'from CPAC.func_preproc.utils import degrees_to_mm, mm_to_degrees']
                notch = pe.Node(Function(input_names=['motion_params',
                                                      'filter_type',
                                                      'TR',
                                                      'fc_RR_min',
                                                      'fc_RR_max',
                                                      'center_freq',
                                                      'freq_bw',
                                                      'lowpass_cutoff',
                                                      'filter_order'],
                                         output_names=['filtered_motion_params',
                                                       'filter_info',
                                                       'filter_plot'],
                                         function=notch_filter_motion,
                                         imports=notch_imports),
                                name='filter_motion_params')

                notch.inputs.filter_type = config.motion_estimate_filter['filter_type']
                notch.inputs.fc_RR_min = config.motion_estimate_filter['breathing_rate_min']
                notch.inputs.fc_RR_max = config.motion_estimate_filter['breathing_rate_max']
                notch.inputs.center_freq = config.motion_estimate_filter['center_frequency']
                notch.inputs.freq_bw = config.motion_estimate_filter['filter_bandwidth']
                notch.inputs.lowpass_cutoff = config.motion_estimate_filter['lowpass_cutoff']
                notch.inputs.filter_order = config.motion_estimate_filter['filter_order']

                preproc.connect(notch, 'filter_info',
                                output_node, 'motion_filter_info')

                preproc.connect(notch, 'filter_plot',
                                output_node, 'motion_filter_plot')

                preproc.connect(normalize_motion_params, 'out_file',
                                notch, 'motion_params')
                preproc.connect(input_node, 'TR', notch, 'TR')

                preproc.connect(notch, 'filtered_motion_params', 
                                output_node, 'movement_parameters')

            else:
                preproc.connect(normalize_motion_params, 'out_file',
                                output_node, 'movement_parameters')
        else:
            preproc.connect(normalize_motion_params, 'out_file',
                            output_node, 'movement_parameters')

        preproc.connect(func_motion_correct_A, 'mat_file',
                         output_node, 'transform_matrices')
        preproc.connect(func_motion_correct_A, 'rms_files',
                        output_node, 'max_displacement')

    preproc.connect(input_node, 'raw_func',
                    skullstrip_func, 'inputspec.raw_func')

    preproc.connect(input_node, 'anatomical_brain_mask',
                    skullstrip_func, 'inputspec.anatomical_brain_mask')

    preproc.connect(input_node, 'anat_brain',
                    skullstrip_func, 'inputspec.anat_brain')                

    preproc.connect(input_node, 'anat_head',
                    skullstrip_func, 'inputspec.anat_head')                

    preproc.connect(skullstrip_func, 'outputspec.func_brain',
                    output_node, 'skullstrip')

    preproc.connect(skullstrip_func, 'outputspec.func_brain_mask',
                    output_node, 'mask')

    func_mean = pe.Node(interface=afni_utils.TStat(),
                        name='func_mean')

    func_mean.inputs.options = '-mean'
    func_mean.inputs.outputtype = 'NIFTI_GZ'

    preproc.connect(skullstrip_func, 'outputspec.func_brain', 
                    func_mean, 'in_file')

    if "Selected Functional Volume" in config.func_reg_input:
        get_func_volume = pe.Node(interface=afni.Calc(),
                                    name='get_func_volume')

        get_func_volume.inputs.set(
            expr='a',
            single_idx=config.func_reg_input_volume,
            outputtype='NIFTI_GZ'
        )
        preproc.connect(skullstrip_func, 'outputspec.func_brain', 
                        get_func_volume, 'in_file_a')
    
        if config.n4_correct_func_reg_input :

            get_func_volume_n4_corrected = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=200), shrink_factor=2, 
                                            name='get_func_volume_n4_corrected')
            get_func_volume_n4_corrected.inputs.args = '-r True'
            
            preproc.connect(get_func_volume, 'out_file', 
                            get_func_volume_n4_corrected, 'input_image')
            preproc.connect(get_func_volume_n4_corrected, 'output_image',
                            output_node, 'get_func_volume')

        else:
            preproc.connect(get_func_volume, 'out_file',
                            output_node, 'get_func_volume')

    if config.n4_correct_func_reg_input :

        func_mean_n4_corrected = pe.Node(interface = ants.N4BiasFieldCorrection(dimension=3, copy_header=True, bspline_fitting_distance=200), shrink_factor=2, 
                                        name='func_mean_n4_corrected')
        func_mean_n4_corrected.inputs.args = '-r True'
        
        preproc.connect(func_mean, 'out_file', 
                        func_mean_n4_corrected, 'input_image')
        preproc.connect(func_mean_n4_corrected, 'output_image',
                        output_node, 'func_mean')

    else:
        preproc.connect(func_mean, 'out_file',
                        output_node, 'func_mean')

    func_normalize = pe.Node(interface=fsl.ImageMaths(),
                             name='func_normalize')
    func_normalize.inputs.op_string = '-ing 10000'
    func_normalize.inputs.out_data_type = 'float'

    preproc.connect(skullstrip_func, 'outputspec.func_brain',
                    func_normalize, 'in_file')

    preproc.connect(func_normalize, 'out_file',
                    output_node, 'preprocessed')

    # TODO XL review forking
    if 'func' in config.run_longitudinal:
        # get median brain for longitudinal
        func_get_preprocessed_median = pe.Node(interface=afni_utils.TStat(),
                            name='func_get_preprocessed_median')

        func_get_preprocessed_median.inputs.options = '-median'
        func_get_preprocessed_median.inputs.outputtype = 'NIFTI_GZ'

        preproc.connect(func_normalize, 'out_file',
                        func_get_preprocessed_median, 'in_file')

        preproc.connect(func_get_preprocessed_median, 'out_file',
                    output_node, 'preprocessed_median')   

        # get median skull for longitudinal
        func_get_motion_correct_median = pe.Node(interface=afni_utils.TStat(),
                            name='func_get_motion_correct_median')

        func_get_motion_correct_median.inputs.options = '-median'
        func_get_motion_correct_median.inputs.outputtype = 'NIFTI_GZ'

        preproc.connect(func_motion_correct_A, 'out_file',
                        func_get_motion_correct_median, 'in_file')

        preproc.connect(func_get_motion_correct_median, 'out_file',
                    output_node, 'motion_correct_median')

    func_mask_normalize = pe.Node(interface=fsl.ImageMaths(),
                                  name='func_mask_normalize')
    func_mask_normalize.inputs.op_string = '-Tmin -bin'
    func_mask_normalize.inputs.out_data_type = 'char'

    preproc.connect(func_normalize, 'out_file',
                    func_mask_normalize, 'in_file')

    preproc.connect(func_mask_normalize, 'out_file',
                    output_node, 'preprocessed_mask')

    return preproc


def slice_timing_wf(name='slice_timing'):

    # allocate a workflow object
    wf = pe.Workflow(name=name)

    # configure the workflow's input spec
    inputNode = pe.Node(util.IdentityInterface(fields=['func_ts',
                                                       'tr',
                                                       'tpattern']),
                        name='inputspec')

    # configure the workflow's output spec
    outputNode = pe.Node(util.IdentityInterface(fields=['slice_time_corrected']),
                         name='outputspec')

    # create TShift AFNI node
    func_slice_timing_correction = pe.Node(interface=preprocess.TShift(),
                                           name='slice_timing')
    func_slice_timing_correction.inputs.outputtype = 'NIFTI_GZ'


    wf.connect([
        (
            inputNode,
            func_slice_timing_correction,
            [
                (
                    'func_ts',
                    'in_file'
                ),
                (
                    # add the @ prefix to the tpattern file going into
                    # AFNI 3dTshift - needed this so the tpattern file
                    # output from get_scan_params would be tied downstream
                    # via a connection (to avoid poofing)
                    ('tpattern', nullify, add_afni_prefix),
                    'tpattern'
                ),
                (
                    ('tr', nullify),
                    'tr'
                ),
            ]
        ),
    ])

    wf.connect(func_slice_timing_correction, 'out_file',
               outputNode, 'slice_time_corrected')

    return wf


def get_idx(in_files, stop_idx=None, start_idx=None):
    """
    Method to get the first and the last slice for
    the functional run. It verifies the user specified
    first and last slice. If the values are not valid, it
    calculates and returns the very first and the last slice

    Parameters
    ----------
    in_file : string (nifti file)
       Path to input functional run

    stop_idx : int
        Last volume to be considered, specified by user
        in the configuration file

    stop_idx : int
        First volume to be considered, specified by user
        in the configuration file

    Returns
    -------
    stop_idx :  int
        Value of first slice to consider for the functional run

    start_idx : int
        Value of last slice to consider for the functional run

    """

    # Import packages
    from nibabel import load

    # Init variables
    img = load(in_files)
    hdr = img.get_header()
    shape = hdr.get_data_shape()

    # Check to make sure the input file is 4-dimensional
    if len(shape) != 4:
        raise TypeError('Input nifti file: %s is not a 4D file' % in_files)
    # Grab the number of volumes
    nvols = int(hdr.get_data_shape()[3])

    if (start_idx == None) or (start_idx < 0) or (start_idx > (nvols - 1)):
        startidx = 0
    else:
        startidx = start_idx

    if (stop_idx == None) or (stop_idx > (nvols - 1)):
        stopidx = nvols - 1
    else:
        stopidx = stop_idx

    return stopidx, startidx


def connect_func_init(workflow, strat_list, c, unique_id=None):

    if c.runScaling is True:
        for num_strat, strat in enumerate(strat_list):
            # scale func data based on configuration information
            scale_func_wf = create_scale_func_wf(
                scaling_factor=c.scaling_factor,
                wf_name="scale_func_%d" % (num_strat)
            )

            # connect the functional data from the leaf node into the wf
            node, out_file = strat.get_leaf_properties()
            workflow.connect(node, out_file,
                            scale_func_wf, 'inputspec.func')

            # replace the leaf node with the output from the recently added
            # workflow
            strat.set_leaf_properties(scale_func_wf, 'outputspec.scaled_func')

    for num_strat, strat in enumerate(strat_list):
        # Truncate scan length based on configuration information
        if unique_id is None:
            trunc_wf = create_wf_edit_func(
                wf_name=f"edit_func_{num_strat}"
            )
        else:
            trunc_wf = create_wf_edit_func(
                wf_name=f"edit_func_{unique_id}_{num_strat}"
            )

        # connect the functional data from the leaf node into the wf
        node, out_file = strat.get_leaf_properties()
        workflow.connect(node, out_file,
                         trunc_wf, 'inputspec.func')

        # connect the other input parameters
        node, node_out = strat['start_idx']
        workflow.connect(node, node_out,
                         trunc_wf, 'inputspec.start_idx')

        node, node_out = strat['stop_idx']
        workflow.connect(node, node_out,
                         trunc_wf, 'inputspec.stop_idx')

        # replace the leaf node with the output from the recently added
        # workflow
        strat.set_leaf_properties(trunc_wf, 'outputspec.edited_func')
        strat.update_resource_pool({
            'raw_functional_trunc': (trunc_wf, 'outputspec.edited_func'),
        })

    # Motion Statistics Workflow
    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        if 0 in c.runMotionStatisticsFirst:
            new_strat_list += [strat.fork()]

        if 1 in c.runMotionStatisticsFirst:

            for skullstrip_tool in c.functionalMasking:

                skullstrip_tool = skullstrip_tool.lower()

                for motion_correct_ref in c.motion_correction_reference:

                    motion_correct_ref = motion_correct_ref.lower()

                    if " " in motion_correct_ref:
                        motion_correct_ref = motion_correct_ref.replace(" ", "_")

                    for motion_correct_tool in c.motion_correction:

                        motion_correct_tool = motion_correct_tool.lower()

                        for motion_estimate_filter_option in  c.motion_estimate_filter['run']:

                            if motion_estimate_filter_option:
                                func_preproc_workflow_name=f'func_preproc_before_stc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_motion_filter_{num_strat}'
                                motion_stats_workflow_name=f'gen_motion_stats_before_stc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_motion_filter_{num_strat}'
                            else:
                                func_preproc_workflow_name=f'func_preproc_before_stc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_{num_strat}'
                                motion_stats_workflow_name=f'gen_motion_stats_before_stc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_{num_strat}'

                            new_strat = strat.fork()

                            func_preproc = create_func_preproc(
                                skullstrip_tool=skullstrip_tool,
                                motion_correct_tool=motion_correct_tool,
                                motion_correct_ref=motion_correct_ref,
                                motion_estimate_filter_option=motion_estimate_filter_option,
                                config=c,
                                wf_name=func_preproc_workflow_name
                            )

                            node, out_file = new_strat['raw_functional_trunc']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.raw_func')

                            node, out_file = new_strat.get_leaf_properties()
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.func')

                            node, out_file = new_strat['anatomical_brain']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.anat_brain')

                            node, out_file = new_strat['anatomical_brain_mask']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.anatomical_brain_mask')

                            node, out_file = new_strat['anatomical_skull_leaf']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.anat_head')

                            node, out_file = strat['tr']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.TR')

                            func_preproc.inputs.inputspec.twopass = \
                                getattr(c, 'functional_volreg_twopass', True)

                            new_strat.update_resource_pool({
                                'bold_masking_method': skullstrip_tool,
                                'motion_correction_method': motion_correct_tool,
                                'motion_correction_ref': motion_correct_ref,
                                'movement_parameters': (
                                func_preproc, 'outputspec.movement_parameters'),
                                'max_displacement': (
                                func_preproc, 'outputspec.max_displacement'),
                                'functional_brain_mask_before_stc': (
                                func_preproc, 'outputspec.mask'),
                                'motion_correct_before_stc': (
                                func_preproc, 'outputspec.motion_correct'),
                                'coordinate_transformation': (
                                func_preproc, 'outputspec.transform_matrices'),
                            })

                            gen_motion_stats = motion_power_statistics(
                                name=motion_stats_workflow_name,
                                motion_correct_tool=motion_correct_tool)

                            # Special case where the workflow is not getting outputs from
                            # resource pool but is connected to functional datasource
                            node, out_file = new_strat['subject']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.subject_id')

                            node, out_file = new_strat['scan']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.scan_id')

                            node, out_file = new_strat[
                                'motion_correct_before_stc']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.motion_correct')

                            node, out_file = new_strat['movement_parameters']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.movement_parameters')

                            node, out_file = new_strat['max_displacement']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.max_displacement')

                            node, out_file = new_strat[
                                'functional_brain_mask_before_stc']
                            workflow.connect(node, out_file,
                                            gen_motion_stats, 'inputspec.mask')

                            node, out_file = new_strat[
                                'coordinate_transformation']
                            workflow.connect(node, out_file,
                                            gen_motion_stats,
                                            'inputspec.transformations')

                            new_strat.append_name(gen_motion_stats.name)

                            new_strat.update_resource_pool({
                                'frame_wise_displacement_power': (
                                gen_motion_stats, 'outputspec.FDP_1D'),
                                'frame_wise_displacement_jenkinson': (
                                gen_motion_stats, 'outputspec.FDJ_1D'),
                                'dvars': (
                                gen_motion_stats, 'outputspec.DVARS_1D'),
                                'power_params': (
                                gen_motion_stats, 'outputspec.power_params'),
                                'motion_params': (
                                gen_motion_stats, 'outputspec.motion_params')
                            })

                            new_strat_list.append(new_strat)

    strat_list = new_strat_list

    # Despike Workflow
    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        if 0 in c.runDespike:
            new_strat_list += [strat.fork()]

        if 1 in c.runDespike:
            new_strat = strat.fork()

            despike = pe.Node(interface=preprocess.Despike(),
                              name='func_despiked_{0}'.format(num_strat))
            despike.inputs.outputtype = 'NIFTI_GZ'

            node, out_file = new_strat.get_leaf_properties()
            workflow.connect(node, out_file,
                             despike, 'in_file')

            new_strat.set_leaf_properties(despike, 'out_file')

            new_strat.update_resource_pool({
                'despiked': (despike, 'out_file')
            })

            new_strat.append_name(despike.name)
            new_strat_list.append(new_strat)

    strat_list = new_strat_list

    # Slice Timing Correction Workflow
    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        if 0 in c.slice_timing_correction:
            new_strat_list += [strat.fork()]

        if 1 in c.slice_timing_correction:
            new_strat = strat.fork()

            if unique_id is None:
                slice_time = slice_timing_wf(
                    name=f'func_slice_timing_correction_{num_strat}')
            else:
                slice_time = slice_timing_wf(
                    name=f'func_slice_timing_correction_{unique_id}_{num_strat}')
            
            node, out_file = new_strat.get_leaf_properties()
            workflow.connect(node, out_file, slice_time,
                             'inputspec.func_ts')

            node, node_out = new_strat['tr']
            workflow.connect(node, node_out,
                             slice_time, 'inputspec.tr')

            node, node_out = new_strat['tpattern']
            workflow.connect(node, node_out,
                             slice_time, 'inputspec.tpattern')

            # add the name of the node to the strat name
            new_strat.append_name(slice_time.name)

            # set the leaf node
            new_strat.set_leaf_properties(slice_time,
                                          'outputspec.slice_time_corrected')

            # add the outputs to the resource pool
            new_strat.update_resource_pool({
                'slice_time_corrected': (
                slice_time, 'outputspec.slice_time_corrected')
            })

            new_strat_list.append(new_strat)

    # add new strats (if forked)
    strat_list = new_strat_list

    return (workflow, strat_list)


def connect_func_preproc(workflow, strat_list, c, unique_id=None):

    from CPAC.func_preproc.func_preproc import create_func_preproc
    
    new_strat_list = []

    for num_strat, strat in enumerate(strat_list):

        if 'motion_params' in strat:

            # skullstripping tool
            skullstrip_tool = strat.get('bold_masking_method')

            # motion correction reference
            motion_correct_ref = strat.get('motion_correction_ref')

            # motion correction tool
            motion_correct_tool = strat.get('motion_correction_method')

            # motion estimate filter option
            motion_estimate_filter_option = any('motion_filter' in node_name for node_name in strat.get_nodes_names())

            if motion_estimate_filter_option:
                if unique_id is None:
                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_motion_filter_{num_strat}'
                else:
                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_motion_filter_{unique_id}_{num_strat}'
            else:
                if unique_id is None:
                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_{num_strat}'
                else:
                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_{unique_id}_{num_strat}'

            func_preproc = create_func_preproc(
                skullstrip_tool=skullstrip_tool,
                motion_correct_tool=motion_correct_tool,
                motion_correct_ref=motion_correct_ref,
                motion_estimate_filter_option=motion_estimate_filter_option,
                config=c,
                wf_name=workflow_name
            )

            node, out_file = strat['raw_functional_trunc']
            workflow.connect(node, out_file, func_preproc,
                            'inputspec.raw_func')

            node, out_file = strat.get_leaf_properties()
            workflow.connect(node, out_file, func_preproc,
                            'inputspec.func')

            node, out_file = strat['anatomical_brain']
            workflow.connect(node, out_file, func_preproc,
                            'inputspec.anat_brain')

            node, out_file = strat['anatomical_brain_mask']
            workflow.connect(node, out_file, func_preproc,
                            'inputspec.anatomical_brain_mask')

            node, out_file = strat['anatomical_skull_leaf']
            workflow.connect(node, out_file, func_preproc,
                            'inputspec.anat_head')

            node, out_file = strat['tr']
            workflow.connect(node, out_file, func_preproc,
                             'inputspec.TR')

            func_preproc.inputs.inputspec.twopass = \
                getattr(c, 'functional_volreg_twopass', True)

            strat.append_name(func_preproc.name)
            strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

            strat.update_resource_pool({
                'mean_functional': (func_preproc, 'outputspec.func_mean'),
                'selected_func_volume': (func_preproc, 'outputspec.get_func_volume'),
                'functional_preprocessed_mask': (func_preproc, 'outputspec.preprocessed_mask'),                              
                'functional_preprocessed': (func_preproc, 'outputspec.preprocessed'),
                'functional_brain_mask': (func_preproc, 'outputspec.mask'),
                'motion_correct': (func_preproc, 'outputspec.motion_correct'),
                'motion_estimate_filter_info_design': (func_preproc, 'outputspec.motion_filter_info'),
                'motion_estimate_filter_info_plot': (func_preproc, 'outputspec.motion_filter_plot')
            })

            if 'func' in c.run_longitudinal:
                strat.update_resource_pool({
                    'functional_preprocessed_median': (func_preproc, 'outputspec.preprocessed_median'),
                    'motion_correct_median': (func_preproc, 'outputspec.motion_correct_median'),                                
                })

            new_strat_list.append(strat)

        else:

            for skullstrip_tool in c.functionalMasking:
                
                skullstrip_tool = skullstrip_tool.lower()

                for motion_correct_ref in c.motion_correction_reference:

                    motion_correct_ref = motion_correct_ref.lower()

                    if " " in motion_correct_ref:
                        motion_correct_ref = motion_correct_ref.replace(" ", "_")

                    for motion_correct_tool in c.motion_correction:

                        motion_correct_tool = motion_correct_tool.lower()
                        
                        for motion_estimate_filter_option in  c.motion_estimate_filter['run']:

                            if motion_estimate_filter_option:
                                if unique_id is None:
                                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_motion_filter_{num_strat}'
                                else:
                                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_motion_filter_{unique_id}_{num_strat}'
                            else:
                                if unique_id is None:
                                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_{num_strat}'
                                else:
                                    workflow_name=f'func_preproc_{skullstrip_tool}_{motion_correct_ref}_{motion_correct_tool}_{unique_id}_{num_strat}'

                            new_strat = strat.fork()

                            func_preproc = create_func_preproc(
                                skullstrip_tool=skullstrip_tool,
                                motion_correct_tool=motion_correct_tool,
                                motion_correct_ref=motion_correct_ref,
                                motion_estimate_filter_option=motion_estimate_filter_option,
                                config=c,
                                wf_name=workflow_name
                            )

                            node, out_file = new_strat['raw_functional_trunc']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.raw_func')

                            node, out_file = new_strat.get_leaf_properties()
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.func')
                            
                            if skullstrip_tool == 'anatomical_refined' or skullstrip_tool == 'anatomical_based':
                                
                                node, out_file = new_strat['anatomical_skull_leaf']
                                workflow.connect(node, out_file, func_preproc,
                                                'inputspec.anat_head')

                                node, out_file = new_strat['anatomical_brain']
                                workflow.connect(node, out_file, func_preproc,
                                                'inputspec.anat_brain')

                                node, out_file = new_strat['anatomical_brain_mask']
                                workflow.connect(node, out_file, func_preproc,
                                                'inputspec.anatomical_brain_mask')

                            node, out_file = strat['tr']
                            workflow.connect(node, out_file, func_preproc,
                                            'inputspec.TR')

                            func_preproc.inputs.inputspec.twopass = \
                                getattr(c, 'functional_volreg_twopass', True)

                            new_strat.append_name(func_preproc.name)
                            new_strat.set_leaf_properties(func_preproc, 'outputspec.preprocessed')

                            new_strat.update_resource_pool({
                                'bold_masking_method': skullstrip_tool,
                                'motion_correction_method': motion_correct_tool,
                                'motion_correction_ref': motion_correct_ref,
                                'mean_functional': (func_preproc, 'outputspec.func_mean'),
                                'selected_func_volume': (func_preproc, 'outputspec.get_func_volume'),
                                'functional_preprocessed_mask': (func_preproc, 'outputspec.preprocessed_mask'),
                                'movement_parameters': (func_preproc, 'outputspec.movement_parameters'),
                                'max_displacement': (func_preproc, 'outputspec.max_displacement'),
                                'functional_preprocessed': (func_preproc, 'outputspec.preprocessed'),
                                'functional_brain_mask': (func_preproc, 'outputspec.mask'),
                                'motion_correct': (func_preproc, 'outputspec.motion_correct'),
                                'coordinate_transformation': (func_preproc, 'outputspec.transform_matrices'),
                                'motion_estimate_filter_info_design': (func_preproc, 'outputspec.motion_filter_info'),
                                'motion_estimate_filter_info_plot': (func_preproc, 'outputspec.motion_filter_plot')
                            })

                            if 'func' in c.run_longitudinal:
                                new_strat.update_resource_pool({
                                    'functional_preprocessed_median': (func_preproc, 'outputspec.preprocessed_median'),
                                    'motion_correct_median': (func_preproc, 'outputspec.motion_correct_median'),                                
                                })
                            
                            new_strat_list.append(new_strat)

    return workflow, new_strat_list
