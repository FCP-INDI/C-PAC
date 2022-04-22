from nipype import logging
from nipype.interfaces import ants

logger = logging.getLogger('nipype.workflow')

from CPAC.pipeline import nipype_pipeline_engine as pe
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
from nipype.interfaces import afni
from nipype.interfaces.afni import preprocess
from nipype.interfaces.afni import utils as afni_utils
from CPAC.func_preproc.utils import add_afni_prefix, nullify, chunk_ts, \
    split_ts_chunks, oned_text_concat, notch_filter_motion
from CPAC.utils.interfaces.function import Function
from CPAC.generate_motion_statistics import motion_power_statistics

from CPAC.utils.utils import check_prov_for_motion_tool

# niworkflows
from ..utils.interfaces.ants import AI, \
    CopyImageHeaderInformation  # noqa: E402


def collect_arguments(*args):
    command_args = []
    if args[0]:
        command_args += [args[1]]
    command_args += args[2:]
    return ' '.join(command_args)


def anat_refined_mask(init_bold_mask=True, wf_name='init_bold_mask'):

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

    # TODO add an option to select volume
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
    if init_bold_mask == True:
        # 2.1.1 N4BiasFieldCorrection single volume of raw_func
        func_single_volume_n4_corrected = pe.Node(
            interface=ants.N4BiasFieldCorrection(dimension=3,
                                                 copy_header=True,
                                                 bspline_fitting_distance=200),
            shrink_factor=2,
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
    else:
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

    if init_bold_mask == True:
        wf.connect(func_tmp_brain_mask_dil, 'out_file',
                   func_mask, 'in_file')
    else:
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

    # 3.2 Binarize transfered image and fill holes to get BOLD mask. 
    # Binarize
    func_mask_bin = pe.Node(interface=fsl.ImageMaths(),
                                name='func_mask')
    func_mask_bin.inputs.op_string = '-bin'

    wf.connect(reg_anat_brain_to_func, 'out_file', 
                func_mask_bin, 'in_file')

    wf.connect(func_mask_bin, 'out_file', 
                output_node, 'func_brain_mask')

    return wf


def normalize_motion_parameters(in_file):
    """
    Convert FSL mcflirt motion parameters to AFNI space
    """
    import os
    import numpy as np

    motion_params = np.genfromtxt(in_file).T
    motion_params = np.vstack((motion_params[2, :] * 180 / np.pi,
                               motion_params[0, :] * 180 / np.pi,
                               -motion_params[1, :] * 180 / np.pi,
                               motion_params[5, :],
                               motion_params[3, :],
                               -motion_params[4, :]))
    motion_params = np.transpose(motion_params)

    out_file = os.path.join(os.getcwd(), 'motion_params.1D')
    np.savetxt(out_file, motion_params)

    return out_file


def get_mcflirt_rms_abs(rms_files):
    for path in rms_files:
        if 'abs.rms' in path:
            abs_file = path
        if 'rel.rms' in path:
            rels_file = path
    return (abs_file, rels_file)


def estimate_reference_image(in_file):

    # fMRIPrep-style BOLD reference
    # Ref: https://github.com/nipreps/niworkflows/blob/maint/1.3.x/niworkflows/interfaces/registration.py#L446-L549

    import os
    import numpy as np
    import nibabel as nb

    ref_input = [in_file]
    mc_out_file = 'bold_mc.nii.gz'

    # Build the nibabel spatial image we will work with
    ref_im = []
    for im_i in ref_input:
        max_new_volumes = 50 - len(ref_im)
        if max_new_volumes <= 0:
            break
        nib_i = nb.squeeze_image(nb.load(im_i))
        if nib_i.dataobj.ndim == 3:
            ref_im.append(nib_i)
        elif nib_i.dataobj.ndim == 4:
            ref_im += nb.four_to_three(nib_i.slicer[..., :max_new_volumes])
    ref_im = nb.squeeze_image(nb.concat_images(ref_im))

    out_file = os.path.join(os.getcwd(), "ref_bold.nii.gz")

    # Slicing may induce inconsistencies with shape-dependent values in extensions.
    # For now, remove all. If this turns out to be a mistake, we can select extensions
    # that don't break pipeline stages.
    ref_im.header.extensions.clear()

    if ref_im.shape[-1] > 40:
        ref_im = nb.Nifti1Image(
            ref_im.dataobj[:, :, :, 20:40], ref_im.affine, ref_im.header
        )

    ref_name = os.path.join(os.getcwd(), "slice.nii.gz")
    ref_im.to_filename(ref_name)
    cmd = '3dvolreg -Fourier -twopass -zpad 4 -prefix %s %s'%(mc_out_file, ref_name)
    os.system(cmd)

    mc_slice_nii = nb.load(mc_out_file)

    median_image_data = np.median(mc_slice_nii.get_fdata(), axis=3)

    nb.Nifti1Image(median_image_data, ref_im.affine, ref_im.header).to_filename(
        out_file
    )

    return out_file


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
                            name='func_drop_trs',
                            mem_gb=0.37,
                            mem_x=(739971956005215 / 151115727451828646838272,
                                   'in_file_a'))

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


def slice_timing_wf(name='slice_timing', tpattern=None, tzero=None):
    # allocate a workflow object
    wf = pe.Workflow(name=name)

    # configure the workflow's input spec
    inputNode = pe.Node(util.IdentityInterface(fields=['func_ts',
                                                       'tr',
                                                       'tpattern']),
                        name='inputspec')

    # configure the workflow's output spec
    outputNode = pe.Node(
        util.IdentityInterface(fields=['slice_time_corrected']),
        name='outputspec')

    # create TShift AFNI node
    func_slice_timing_correction = pe.Node(interface=preprocess.TShift(),
                                           name='slice_timing',
                                           mem_gb=0.45,
                                           mem_x=(5247073869855161 /
                                                  604462909807314587353088,
                                                  'in_file'))
    func_slice_timing_correction.inputs.outputtype = 'NIFTI_GZ'

    if tzero is not None:
        func_slice_timing_correction.inputs.tzero = tzero

    wf.connect([
        (
            inputNode,
            func_slice_timing_correction,
            [
                (
                    'func_ts',
                    'in_file'
                ),
                # (
                #     # add the @ prefix to the tpattern file going into
                #     # AFNI 3dTshift - needed this so the tpattern file
                #     # output from get_scan_params would be tied downstream
                #     # via a connection (to avoid poofing)
                #     ('tpattern', nullify, add_afni_prefix),
                #     'tpattern'
                # ),
                (
                    ('tr', nullify),
                    'tr'
                ),
            ]
        ),
    ])

    if tpattern is not None:
        func_slice_timing_correction.inputs.tpattern = tpattern
    else:
        wf.connect(inputNode, ('tpattern', nullify, add_afni_prefix),
               func_slice_timing_correction, 'tpattern')

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

    if (start_idx == None) or (int(start_idx) < 0) or (
        int(start_idx) > (nvols - 1)):
        startidx = 0
    else:
        startidx = int(start_idx)

    if (stop_idx == None) or (int(stop_idx) > (nvols - 1)):
        stopidx = nvols - 1
    else:
        stopidx = int(stop_idx)

    return stopidx, startidx


def motion_correct_connections(wf, cfg, strat_pool, pipe_num, opt):
    if opt != '3dvolreg' and opt != 'mcflirt':
        raise Exception("\n\n[!] Error: The 'tool' parameter of the "
                        "'motion_correction' workflow must be either "
                        "'3dvolreg' or 'mcflirt'.\n\nTool input: "
                        "{0}\n\n".format(opt))

    if cfg:
        if int(cfg.pipeline_setup['system_config'][
                   'max_cores_per_participant']) > 1:
            chunk_imports = ['import nibabel as nb']
            chunk = pe.Node(Function(input_names=['func_file',
                                                  'n_chunks',
                                                  'chunk_size'],
                                     output_names=['TR_ranges'],
                                     function=chunk_ts,
                                     imports=chunk_imports),
                            name=f'chunk_{pipe_num}')

            #chunk.inputs.n_chunks = int(cfg.pipeline_setup['system_config'][
            #                              'max_cores_per_participant'])
            
            # 10-TR sized chunks
            chunk.inputs.chunk_size = 10
            
            node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
            wf.connect(node, out, chunk, 'func_file')

            split_imports = ['import os', 'import subprocess']
            split = pe.Node(Function(input_names=['func_file',
                                                  'tr_ranges'],
                                     output_names=['split_funcs'],
                                     function=split_ts_chunks,
                                     imports=split_imports),
                            name=f'split_{pipe_num}')

            node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
            wf.connect(node, out, split, 'func_file')
            wf.connect(chunk, 'TR_ranges', split, 'tr_ranges')

            out_split_func = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_split_func_{pipe_num}')

            wf.connect(split, 'split_funcs', out_split_func, 'out_file')

            func_motion_correct = pe.MapNode(interface=preprocess.Volreg(),
                                             name='func_generate_'
                                                  f'ref_{pipe_num}',
                                             iterfield=['in_file'])

            wf.connect(out_split_func, 'out_file',
                       func_motion_correct, 'in_file')

            func_concat = pe.Node(interface=afni_utils.TCat(),
                                  name=f'func_concat_{pipe_num}')
            func_concat.inputs.outputtype = 'NIFTI_GZ'

            wf.connect(func_motion_correct, 'out_file',
                       func_concat, 'in_files')

            out_motion = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_motion_{pipe_num}')

            wf.connect(func_concat, 'out_file', out_motion, 'out_file')

        else:
            out_split_func = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_split_func_{pipe_num}')

            node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
            wf.connect(node, out, out_split_func, 'out_file')

            func_motion_correct = pe.Node(interface=preprocess.Volreg(),
                                          name=f'func_generate_ref_{pipe_num}')

            wf.connect(out_split_func, 'out_file',
                       func_motion_correct, 'in_file')

            out_motion = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_motion_{pipe_num}')

            wf.connect(func_motion_correct, 'out_file',
                       out_motion, 'out_file')

    else:
        out_split_func = pe.Node(
            interface=util.IdentityInterface(fields=['out_file']),
            name=f'out_split_func_{pipe_num}')

        node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
        wf.connect(node, out, out_split_func, 'out_file')

        func_motion_correct = pe.Node(interface=preprocess.Volreg(),
                                      name=f'func_generate_ref_{pipe_num}')

        wf.connect(out_split_func, 'out_file',
                   func_motion_correct, 'in_file')

        out_motion = pe.Node(
            interface=util.IdentityInterface(fields=['out_file']),
            name=f'out_motion_{pipe_num}')

        wf.connect(func_motion_correct, 'out_file', out_motion, 'out_file')

    func_motion_correct.inputs.zpad = 4
    func_motion_correct.inputs.outputtype = 'NIFTI_GZ'

    args = f'-Fourier'
    if cfg.functional_preproc['motion_estimates_and_correction'][
        'motion_correction']['AFNI-3dvolreg']['functional_volreg_twopass']:
        args = f'-twopass {args}'

    func_motion_correct.inputs.args = args

    # Calculate motion parameters
    if opt == '3dvolreg':
        func_motion_correct_A = func_motion_correct.clone(
            f'func_motion_correct_3dvolreg_{pipe_num}')
        func_motion_correct_A.inputs.md1d_file = 'max_displacement.1D'
        func_motion_correct_A.inputs.args = args

        wf.connect(out_split_func, 'out_file',
                   func_motion_correct_A, 'in_file')

        node, out = strat_pool.get_data('motion-basefile')
        wf.connect(node, out, func_motion_correct_A, 'basefile')

        if cfg:
            if int(cfg.pipeline_setup['system_config'][
                       'max_cores_per_participant']) > 1:
                motion_concat = pe.Node(interface=afni_utils.TCat(),
                                        name=f'motion_concat_{pipe_num}')
                motion_concat.inputs.outputtype = 'NIFTI_GZ'

                wf.connect(func_motion_correct_A, 'out_file',
                           motion_concat, 'in_files')

                out_motion_A = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_motion_A_{pipe_num}')

                wf.connect(motion_concat, 'out_file',
                           out_motion_A, 'out_file')

                concat_imports = ['import os']
                md1d_concat = pe.Node(Function(input_names=['in_files'],
                                               output_names=['out_file'],
                                               function=oned_text_concat,
                                               imports=concat_imports),
                                      name=f'md1d_concat_{pipe_num}')

                wf.connect(func_motion_correct_A, 'md1d_file',
                           md1d_concat, 'in_files')

                oned_concat = pe.Node(Function(input_names=['in_files'],
                                               output_names=['out_file'],
                                               function=oned_text_concat,
                                               imports=concat_imports),
                                      name=f'oned_concat_{pipe_num}')

                wf.connect(func_motion_correct_A, 'oned_file',
                           oned_concat, 'in_files')

                oned_matrix_concat = pe.Node(
                    Function(input_names=['in_files'],
                             output_names=['out_file'],
                             function=oned_text_concat,
                             imports=concat_imports),
                    name=f'oned_matrix_concat_{pipe_num}')

                wf.connect(func_motion_correct_A, 'oned_matrix_save',
                           oned_matrix_concat, 'in_files')

                out_md1d = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_md1d_{pipe_num}')

                wf.connect(md1d_concat, 'out_file',
                           out_md1d, 'out_file')

                out_oned = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_oned_{pipe_num}')

                wf.connect(oned_concat, 'out_file',
                           out_oned, 'out_file')

                out_oned_matrix = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_oned_matrix_{pipe_num}')

                wf.connect(oned_matrix_concat, 'out_file',
                           out_oned_matrix, 'out_file')

            else:
                out_motion_A = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_motion_A_{pipe_num}')

                wf.connect(func_motion_correct_A, 'out_file',
                           out_motion_A, 'out_file')

                out_md1d = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_md1d_{pipe_num}')

                wf.connect(func_motion_correct_A, 'md1d_file',
                           out_md1d, 'out_file')

                out_oned = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_oned_{pipe_num}')

                wf.connect(func_motion_correct_A, 'oned_file',
                           out_oned, 'out_file')

                out_oned_matrix = pe.Node(
                    interface=util.IdentityInterface(fields=['out_file']),
                    name=f'out_oned_matrix_{pipe_num}')

                wf.connect(func_motion_correct_A, 'oned_matrix_save',
                           out_oned_matrix, 'out_file')

        else:
            out_motion_A = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_motion_A_{pipe_num}')

            wf.connect(func_motion_correct_A, 'out_file',
                       out_motion_A, 'out_file')

            out_md1d = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_md1d_{pipe_num}')

            wf.connect(func_motion_correct_A, 'md1d_file',
                       out_md1d, 'out_file')

            out_oned = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_oned_{pipe_num}')

            wf.connect(func_motion_correct_A, 'oned_file',
                       out_oned, 'out_file')

            out_oned_matrix = pe.Node(
                interface=util.IdentityInterface(fields=['out_file']),
                name=f'out_oned_matrix_{pipe_num}')

            wf.connect(func_motion_correct_A, 'oned_matrix_save',
                       out_oned_matrix, 'out_file')

        outputs = {
            'desc-preproc_bold': (out_motion_A, 'out_file'),
            'desc-motion_bold': (out_motion_A, 'out_file'),
            'max-displacement': (out_md1d, 'out_file'),
            'movement-parameters': (out_oned, 'out_file'),
            'coordinate-transformation': (out_oned_matrix, 'out_file')
        }

    elif opt == 'mcflirt':
        func_motion_correct_A = pe.Node(
            interface=fsl.MCFLIRT(save_mats=True, save_plots=True),
            name=f'func_motion_correct_mcflirt_{pipe_num}', mem_gb=2.5)

        func_motion_correct_A.inputs.save_mats = True
        func_motion_correct_A.inputs.save_plots = True
        func_motion_correct_A.inputs.save_rms = True

        node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
        wf.connect(node, out, func_motion_correct_A, 'in_file')

        node, out = strat_pool.get_data('motion-basefile')
        wf.connect(node, out, func_motion_correct_A, 'ref_file')

        normalize_motion_params = pe.Node(Function(input_names=['in_file'],
                                                   output_names=['out_file'],
                                                   function=normalize_motion_parameters),
                                          name=f'norm_motion_params_{pipe_num}')

        wf.connect(func_motion_correct_A, 'par_file',
                   normalize_motion_params, 'in_file')

        get_rms_abs = pe.Node(Function(input_names=['rms_files'],
                                       output_names=['abs_file',
                                                     'rels_file'],
                                       function=get_mcflirt_rms_abs),
                              name=f'get_mcflirt_rms_abs_{pipe_num}')

        wf.connect(func_motion_correct_A, 'rms_files',
                   get_rms_abs, 'rms_files')

        outputs = {
            'desc-preproc_bold': (func_motion_correct_A, 'out_file'),
            'desc-motion_bold': (func_motion_correct_A, 'out_file'),
            'max-displacement': (get_rms_abs, 'abs_file'),
            'rels-displacement': (get_rms_abs, 'rels_file'),
            'movement-parameters': (normalize_motion_params, 'out_file'),
            'coordinate-transformation': (func_motion_correct_A, 'mat_file')
        }

    return (wf, outputs)


def func_scaling(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_scaling",
     "config": ["functional_preproc", "scaling"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": ["desc-preproc_bold"]}
    '''

    scale_func_wf = create_scale_func_wf(
        scaling_factor=cfg.scaling_factor,
        wf_name=f"scale_func_{pipe_num}"
    )

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, scale_func_wf, 'inputspec.func')

    outputs = {
        'desc-preproc_bold': (scale_func_wf, 'outputspec.scaled_func')
    }

    return (wf, outputs)


def func_truncate(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_truncate",
     "config": ["functional_preproc", "truncation"],
     "switch": "None",
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": {
         "desc-preproc_bold": {
             "Description": "Truncated functional time-series BOLD data."
         }}
    }
    '''

    # if cfg.functional_preproc['truncation']['start_tr'] == 0 and \
    #                cfg.functional_preproc['truncation']['stop_tr'] == None:
    #    data, key = strat_pool.get_data(["desc-preproc_bold", "bold"],
    #                                           True)
    #    outputs = {key: data}
    #    return (wf, outputs)

    trunc_wf = create_wf_edit_func(
        wf_name=f"edit_func_{pipe_num}"
    )
    trunc_wf.inputs.inputspec.start_idx = cfg.functional_preproc[
        'truncation']['start_tr']
    trunc_wf.inputs.inputspec.stop_idx = cfg.functional_preproc['truncation'][
        'stop_tr']

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, trunc_wf, 'inputspec.func')

    outputs = {
        'desc-preproc_bold': (trunc_wf, 'outputspec.edited_func')
    }

    return (wf, outputs)


def func_despike(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_despike",
     "config": ["functional_preproc", "despiking"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": {
         "desc-preproc_bold": {
             "Description": "De-spiked BOLD time-series via AFNI 3dDespike."
         }}
    }
    '''

    despike = pe.Node(interface=preprocess.Despike(),
                      name=f'func_despiked_{pipe_num}',
                      mem_gb=0.66,
                      mem_x=(8251808479088459 / 1208925819614629174706176,
                             'in_file'))
    despike.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, despike, 'in_file')

    outputs = {
        'desc-preproc_bold': (despike, 'out_file')
    }

    return (wf, outputs)


def func_slice_time(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_slice_time",
     "config": ["functional_preproc", "slice_timing_correction"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"],
                "TR",
                "tpattern"],
     "outputs": {
         "desc-preproc_bold": {
             "Description": "Slice-time corrected BOLD time-series via AFNI 3dTShift."
         },
         "desc-stc_bold": {
             "Description": "Slice-time corrected BOLD time-series via AFNI 3dTShift."}}
    }
    '''

    slice_time = slice_timing_wf(name='func_slice_timing_correction_'
                                 f'{pipe_num}',
                                 tpattern=cfg.functional_preproc[
                                 'slice_timing_correction']['tpattern'],
                                 tzero=cfg.functional_preproc[
                                 'slice_timing_correction']['tzero'])

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, slice_time, 'inputspec.func_ts')

    node, out = strat_pool.get_data('TR')
    wf.connect(node, out, slice_time, 'inputspec.tr')

    node, out = strat_pool.get_data('tpattern')
    wf.connect(node, out, slice_time, 'inputspec.tpattern')

    outputs = {
        'desc-preproc_bold': (slice_time, 'outputspec.slice_time_corrected'),
        'desc-stc_bold': (slice_time, 'outputspec.slice_time_corrected')
    }

    return (wf, outputs)


def func_reorient(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_reorient",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": ["desc-preproc_bold", "desc-reorient_bold"]}
    '''

    func_deoblique = pe.Node(interface=afni_utils.Refit(),
                             name=f'func_deoblique_{pipe_num}',
                             mem_gb=0.68,
                             mem_x=(4664065662093477 /
                                    1208925819614629174706176,
                                    'in_file'))
    func_deoblique.inputs.deoblique = True

    node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
    wf.connect(node, out, func_deoblique, 'in_file')

    func_reorient = pe.Node(interface=afni_utils.Resample(),
                            name=f'func_reorient_{pipe_num}',
                            mem_gb=0,
                            mem_x=(0.0115, 'in_file', 't'))

    func_reorient.inputs.orientation = 'RPI'
    func_reorient.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(func_deoblique, 'out_file', func_reorient, 'in_file')

    outputs = {
        'desc-preproc_bold': (func_reorient, 'out_file'),
        'desc-reorient_bold': (func_reorient, 'out_file')
    }

    return (wf, outputs)


def get_motion_ref(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "get_motion_ref",
     "config": "None",
     "switch": ["functional_preproc", "motion_estimates_and_correction", 
                "run"],
     "option_key": ["functional_preproc", "motion_estimates_and_correction",
                    "motion_correction", "motion_correction_reference"],
     "option_val": ["mean", "median", "selected_volume", "fmriprep_reference"],
     "inputs": [["desc-preproc_bold", "bold"],
                 "bold"],
     "outputs": ["motion-basefile"]}
    '''

    if opt != 'mean' and opt != 'median' and opt != 'selected_volume' and opt != 'fmriprep_reference':
        raise Exception("\n\n[!] Error: The 'tool' parameter of the "
                        "'motion_correction_reference' workflow must be either "
                        "'mean' or 'median' or 'selected_volume' or 'fmriprep_reference'.\n\nTool input: "
                        "{0}\n\n".format(opt))

    if opt == 'mean':
        func_get_RPI = pe.Node(interface=afni_utils.TStat(),
                               name=f'func_get_mean_RPI_{pipe_num}',
                               mem_gb=0.48,
                               mem_x=(1435097126797993 /
                                      302231454903657293676544,
                                      'in_file'))

        func_get_RPI.inputs.options = '-mean'
        func_get_RPI.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
        wf.connect(node, out, func_get_RPI, 'in_file')

    elif opt == 'median':
        func_get_RPI = pe.Node(interface=afni_utils.TStat(),
                               name=f'func_get_median_RPI_{pipe_num}')

        func_get_RPI.inputs.options = '-median'
        func_get_RPI.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
        wf.connect(node, out, func_get_RPI, 'in_file')

    elif opt == 'selected_volume':
        func_get_RPI = pe.Node(interface=afni.Calc(),
                               name=f'func_get_selected_RPI_{pipe_num}')

        func_get_RPI.inputs.set(
            expr='a',
            single_idx=cfg.functional_preproc['motion_estimates_and_correction'][
                'motion_correction']['motion_correction_reference_volume'],
            outputtype='NIFTI_GZ'
        )

        node, out = strat_pool.get_data(['desc-preproc_bold', 'bold'])
        wf.connect(node, out, func_get_RPI, 'in_file_a')

    elif opt == 'fmriprep_reference':
        func_get_RPI = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=estimate_reference_image),
                           name=f'func_get_fmriprep_ref_{pipe_num}')

        node, out = strat_pool.get_data('bold')
        wf.connect(node, out, func_get_RPI, 'in_file')

    outputs = {
        'motion-basefile': (func_get_RPI, 'out_file')
    }

    return (wf, outputs)


def func_motion_correct(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "motion_correction",
     "config": "None",
     "switch": ["functional_preproc", "motion_estimates_and_correction", 
                "run"],
     "option_key": ["functional_preproc", "motion_estimates_and_correction",
                    "motion_correction", "using"],
     "option_val": ["3dvolreg", "mcflirt"],
     "inputs": [(["desc-preproc_bold", "bold"],
                 "motion-basefile")],
     "outputs": ["desc-preproc_bold",
                 "desc-motion_bold",
                 "max-displacement",
                 "rels-displacement",
                 "movement-parameters",
                 "coordinate-transformation"]}
    '''

    wf, outputs = motion_correct_connections(wf, cfg, strat_pool, pipe_num,
                                             opt)

    return (wf, outputs)


def func_motion_estimates(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "motion_estimates",
     "config": "None",
     "switch": ["functional_preproc", "motion_estimates_and_correction", 
                "run"],
     "option_key": ["functional_preproc", "motion_estimates_and_correction",
                    "motion_correction", "using"],
     "option_val": ["3dvolreg", "mcflirt"],
     "inputs": [(["desc-preproc_bold", "bold"],
                 "motion-basefile")],
     "outputs": ["max-displacement",
                 "rels-displacement",
                 "movement-parameters",
                 "coordinate-transformation"]}
    '''

    wf, wf_outputs = motion_correct_connections(wf, cfg, strat_pool, pipe_num,
                                                opt)

    outputs = {
        'max-displacement': wf_outputs['max-displacement'],
        'movement-parameters': wf_outputs['movement-parameters']
    }

    if 'coordinate-transformation' in wf_outputs:
        outputs['coordinate-transformation'] = \
            wf_outputs['coordinate-transformation']

    if 'rels-displacement' in wf_outputs:
        outputs['rels-displacement'] = wf_outputs['rels-displacement']

    return (wf, outputs)


def func_motion_correct_only(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "motion_correction_only",
     "config": "None",
     "switch": ["functional_preproc", "motion_estimates_and_correction", 
                "run"],
     "option_key": ["functional_preproc", "motion_estimates_and_correction",
                    "motion_correction", "using"],
     "option_val": ["3dvolreg", "mcflirt"],
     "inputs": [(["desc-preproc_bold", "bold"],
                 "motion-basefile")],
     "outputs": ["desc-preproc_bold",
                 "desc-motion_bold"]}
    '''

    wf, wf_outputs = motion_correct_connections(wf, cfg, strat_pool, pipe_num,
                                                opt)

    outputs = {
        'desc-motion_bold': wf_outputs['desc-motion_bold']
    }

    return (wf, outputs)


def motion_estimate_filter(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "motion_estimate_filter",
     "config": ["functional_preproc", "motion_estimates_and_correction",
                "motion_estimate_filter"],
     "switch": ["run"],
     "option_key": "filter_type",
     "option_val": ["notch", "lowpass"],
     "inputs": ["movement-parameters",
                "TR"],
     "outputs": ["movement-parameters",
                 "motion-filter-info",
                 "motion-filter-plot"]}
    '''

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
                             output_names=[
                                 'filtered_motion_params',
                                 'filter_info',
                                 'filter_plot'],
                             function=notch_filter_motion,
                             imports=notch_imports),
                    name=f'filter_motion_params_{pipe_num}')

    notch.inputs.filter_type = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['filter_type']
    notch.inputs.fc_RR_min = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['breathing_rate_min']
    notch.inputs.fc_RR_max = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['breathing_rate_max']
    notch.inputs.center_freq = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['center_frequency']
    notch.inputs.freq_bw = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['filter_bandwidth']
    notch.inputs.lowpass_cutoff = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['lowpass_cutoff']
    notch.inputs.filter_order = cfg.functional_preproc[
        "motion_estimates_and_correction"][
        "motion_estimate_filter"]['filter_order']

    node, out = strat_pool.get_data('movement-parameters')
    wf.connect(node, out, notch, 'motion_params')

    node, out = strat_pool.get_data('TR')
    wf.connect(node, out, notch, 'TR')

    outputs = {
        'motion-filter-info': (notch, 'filter_info'),
        'motion-filter-plot': (notch, 'filter_plot'),
        'movement-parameters': (notch, 'filtered_motion_params')
    }

    return (wf, outputs)


def calc_motion_stats(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "calc_motion_stats",
     "config": "None",
     "switch": [["functional_preproc", "run"], 
                ["functional_preproc", "motion_estimates_and_correction", 
                "motion_estimates", "calculate_motion_after"]],
     "option_key": "None",
     "option_val": "None",
     "inputs": [("desc-motion_bold",
                 "space-bold_desc-brain_mask",
                 "movement-parameters",
                 "max-displacement",
                 "rels-displacement",
                 "coordinate-transformation"),
                "subject",
                "scan"],
     "outputs": ["framewise-displacement-power",
                 "framewise-displacement-jenkinson",
                 "dvars",
                 "power-params",
                 "motion-params"]}
    '''

    motion_prov = strat_pool.get_cpac_provenance('movement-parameters')
    motion_correct_tool = check_prov_for_motion_tool(motion_prov)

    gen_motion_stats = motion_power_statistics(
        name=f'gen_motion_stats_{pipe_num}',
        motion_correct_tool=motion_correct_tool)

    # Special case where the workflow is not getting outputs from
    # resource pool but is connected to functional datasource
    node, out_file = strat_pool.get_data('subject')
    wf.connect(node, out_file,
               gen_motion_stats, 'inputspec.subject_id')

    node, out_file = strat_pool.get_data('scan')
    wf.connect(node, out_file,
               gen_motion_stats, 'inputspec.scan_id')

    node, out_file = strat_pool.get_data("desc-motion_bold")
    wf.connect(node, out_file,
               gen_motion_stats, 'inputspec.motion_correct')

    node, out_file = strat_pool.get_data('space-bold_desc-brain_mask')
    wf.connect(node, out_file,
               gen_motion_stats, 'inputspec.mask')

    node, out_file = strat_pool.get_data('movement-parameters')
    wf.connect(node, out_file,
               gen_motion_stats,
               'inputspec.movement_parameters')

    node, out_file = strat_pool.get_data('max-displacement')
    wf.connect(node, out_file,
               gen_motion_stats,
               'inputspec.max_displacement')

    if strat_pool.check_rpool('rels-displacement'):
        node, out_file = strat_pool.get_data('rels-displacement')
        wf.connect(node, out_file, gen_motion_stats,
                   'inputspec.rels_displacement')

    if strat_pool.check_rpool('coordinate-transformation'):
        node, out_file = strat_pool.get_data('coordinate-transformation')
        wf.connect(node, out_file, gen_motion_stats,
                   'inputspec.transformations')

    outputs = {
        'framewise-displacement-power':
            (gen_motion_stats, 'outputspec.FDP_1D'),
        'framewise-displacement-jenkinson':
            (gen_motion_stats, 'outputspec.FDJ_1D'),
        'dvars': (gen_motion_stats, 'outputspec.DVARS_1D'),
        'power-params': (gen_motion_stats, 'outputspec.power_params'),
        'motion-params': (gen_motion_stats, 'outputspec.motion_params')
    }

    return (wf, outputs)


def bold_mask_afni(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "bold_mask_afni",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "AFNI",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": {
         "space-bold_desc-brain_mask": {
             "Description": "Binary brain mask of the BOLD functional "
                            "time-series created by AFNI 3dAutomask."}}}
    '''

    func_get_brain_mask = pe.Node(interface=preprocess.Automask(),
                                  name=f'func_get_brain_mask_AFNI_{pipe_num}')
    func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(["desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, func_get_brain_mask, 'in_file')

    outputs = {
        'space-bold_desc-brain_mask': (func_get_brain_mask, 'out_file')
    }

    return (wf, outputs)


def bold_mask_fsl(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "bold_mask_fsl",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "FSL",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": ["space-bold_desc-brain_mask"]}
    '''

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
        name=f'BET_options_{pipe_num}')

    func_get_brain_mask = pe.Node(interface=fsl.BET(),
                                  name=f'func_get_brain_mask_BET_{pipe_num}')
    func_get_brain_mask.inputs.output_type = 'NIFTI_GZ'
    func_get_brain_mask.inputs.mask = True

    inputnode_bet.inputs.set(
        frac=cfg.functional_preproc['func_masking']['FSL-BET']['frac'],
        mesh_boolean=cfg.functional_preproc['func_masking']['FSL-BET'][
            'mesh_boolean'],
        outline=cfg.functional_preproc['func_masking']['FSL-BET'][
            'outline'],
        padding=cfg.functional_preproc['func_masking']['FSL-BET'][
            'padding'],
        radius=cfg.functional_preproc['func_masking']['FSL-BET']['radius'],
        reduce_bias=cfg.functional_preproc['func_masking']['FSL-BET'][
            'reduce_bias'],
        remove_eyes=cfg.functional_preproc['func_masking']['FSL-BET'][
            'remove_eyes'],
        robust=cfg.functional_preproc['func_masking']['FSL-BET']['robust'],
        skull=cfg.functional_preproc['func_masking']['FSL-BET']['skull'],
        surfaces=cfg.functional_preproc['func_masking']['FSL-BET'][
            'surfaces'],
        threshold=cfg.functional_preproc['func_masking']['FSL-BET'][
            'threshold'],
        vertical_gradient=
        cfg.functional_preproc['func_masking']['FSL-BET'][
            'vertical_gradient'],
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

    if cfg.functional_preproc['func_masking']['FSL-BET'][
        'functional_mean_boolean']:
        func_skull_mean = pe.Node(interface=afni_utils.TStat(),
                                  name=f'func_mean_skull_{pipe_num}')
        func_skull_mean.inputs.options = '-mean'
        func_skull_mean.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
        wf.connect(node, out, func_skull_mean, 'in_file')

        out_node, out_file = (func_skull_mean, 'out_file')

        if cfg.functional_preproc['func_masking']['FSL-BET'][
            'functional_mean_thr']['run']:
            # T=$(fslstats ${subject}_tmean.nii.gz -p 98)
            threshold_T = pe.Node(interface=fsl.ImageStats(),
                                    name=f'func_mean_skull_thr_value_{pipe_num}',
                                    iterfield=['in_file'])
            threshold_T.inputs.op_string = "-p %f " % (cfg.functional_preproc['func_masking']['FSL-BET']['functional_mean_thr']['threshold_value'])
            
            wf.connect(func_skull_mean, 'out_file', threshold_T, 'in_file')
            
            # z=$(echo "$T / 10" | bc -l)
            def form_thr_string(thr):
                threshold_z = str(float(thr/10))
                return '-thr %s' % (threshold_z)

            form_thr_string = pe.Node(util.Function(input_names=['thr'],
                                            output_names=['out_str'],
                                            function=form_thr_string),
                                            name=f'form_thr_string_{pipe_num}')

            wf.connect(threshold_T, 'out_stat', form_thr_string, 'thr')

            # fslmaths ${subject}_tmean.nii.gz -thr ${z} ${subject}_tmean_thr.nii.gz
            func_skull_mean_thr = pe.Node(interface=fsl.ImageMaths(),
                                            name=f'func_mean_skull_thr_{pipe_num}')

            wf.connect(func_skull_mean, 'out_file', func_skull_mean_thr, 'in_file')
            wf.connect(form_thr_string, 'out_str', func_skull_mean_thr, 'op_string')
            
            out_node, out_file = (func_skull_mean_thr, 'out_file')

        if cfg.functional_preproc['func_masking']['FSL-BET'][
            'functional_mean_bias_correction']:

            # fast --nopve -B ${subject}_tmean_thr.nii.gz
            func_mean_skull_fast = pe.Node(interface=fsl.FAST(),
                                                    name=f'func_mean_skull_fast_{pipe_num}')
            func_mean_skull_fast.inputs.no_pve = True
            func_mean_skull_fast.inputs.output_biascorrected = True
            
            wf.connect(out_node, out_file, func_mean_skull_fast, 'in_files')
            
            out_node, out_file = (func_mean_skull_fast, 'restored_image')

        wf.connect(out_node, out_file, func_get_brain_mask, 'in_file')

    else:
        func_get_brain_mask.inputs.functional = True

        node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
        wf.connect(node, out, func_get_brain_mask, 'in_file')

    # erode one voxel of functional brian mask
    erode_one_voxel = pe.Node(interface=fsl.ErodeImage(),
                              name=f'erode_one_voxel_{pipe_num}')

    erode_one_voxel.inputs.kernel_shape = 'box'
    erode_one_voxel.inputs.kernel_size = 1.0

    wf.connect(func_get_brain_mask, 'mask_file',
               erode_one_voxel, 'in_file')

    outputs = {
        'space-bold_desc-brain_mask': (erode_one_voxel, 'out_file')
    }

    return (wf, outputs)


def bold_mask_fsl_afni(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "bold_mask_fsl_afni",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "FSL_AFNI",
     "inputs": [["desc-motion_bold", "desc-preproc_bold", "bold"],
                 "motion-basefile"],
     "outputs": ["space-bold_desc-brain_mask",
                 "desc-ref_bold"]}
    '''

    # fMRIPrep-style BOLD mask
    # Ref: https://github.com/nipreps/niworkflows/blob/maint/1.3.x/niworkflows/func/util.py#L246-L514

    # Initialize transforms with antsAI
    init_aff = pe.Node(
        AI(
            metric=("Mattes", 32, "Regular", 0.2),
            transform=("Affine", 0.1),
            search_factor=(20, 0.12),
            principal_axes=False,
            convergence=(10, 1e-6, 10),
            verbose=True,
        ),
        name=f"init_aff_{pipe_num}",
        n_procs=cfg.pipeline_setup['system_config']['num_OMP_threads'],
    )

    init_aff.inputs.fixed_image = cfg.functional_preproc[
        'func_masking']['FSL_AFNI']['bold_ref']

    init_aff.inputs.fixed_image_mask = cfg.functional_preproc[
        'func_masking']['FSL_AFNI']['brain_mask']

    init_aff.inputs.search_grid = (40, (0, 40, 40))

    # Set up spatial normalization
    norm = pe.Node(
        ants.Registration(
            winsorize_upper_quantile=0.98,
            winsorize_lower_quantile=0.05,
            float=True,
            metric=['Mattes'],
            metric_weight=[1],
            radius_or_number_of_bins=[64],
            transforms=['Affine'],
            transform_parameters=[[0.1]],
            number_of_iterations=[[200]],
            convergence_window_size=[10],
            convergence_threshold=[1.e-9],
            sampling_strategy=['Random', 'Random'],
            smoothing_sigmas=[[2]],
            sigma_units=['mm', 'mm', 'mm'],
            shrink_factors=[[2]],
            sampling_percentage=[0.2],
            use_histogram_matching=[True],
            use_estimate_learning_rate_once=[True]
        ),
        name=f"norm_{pipe_num}",
        n_procs=cfg.pipeline_setup['system_config']['num_OMP_threads'],
    )

    norm.inputs.fixed_image = cfg.functional_preproc[
        'func_masking']['FSL_AFNI']['bold_ref']

    map_brainmask = pe.Node(
        ants.ApplyTransforms(
            interpolation="BSpline",
            float=True,
        ),
        name=f"map_brainmask_{pipe_num}",
    )
    
    # Use the higher resolution and probseg for numerical stability in rounding
    map_brainmask.inputs.input_image = cfg.functional_preproc[
        'func_masking']['FSL_AFNI']['brain_probseg']

    binarize_mask = pe.Node(interface=fsl.maths.MathsCommand(),
                            name=f'binarize_mask_{pipe_num}')
    binarize_mask.inputs.args = '-thr 0.85 -bin'

    # Dilate pre_mask
    pre_dilate = pe.Node(
        fsl.DilateImage(
            operation="max",
            kernel_shape="sphere",
            kernel_size=3.0,
            internal_datatype="char",
        ),
        name=f"pre_mask_dilate_{pipe_num}",
    )

    # Fix precision errors
    # https://github.com/ANTsX/ANTs/wiki/Inputs-do-not-occupy-the-same-physical-space#fixing-precision-errors
    restore_header = pe.Node(CopyImageHeaderInformation(),
                             name=f'restore_header_{pipe_num}')

    # Run N4 normally, force num_threads=1 for stability (images are
    # small, no need for >1)
    n4_correct = pe.Node(
        ants.N4BiasFieldCorrection(
            dimension=3, copy_header=True, bspline_fitting_distance=200
        ),
        shrink_factor=2,
        rescale_intensities = True,
        name=f"n4_correct_{pipe_num}",
        n_procs=1,
    )

    skullstrip_first_pass = pe.Node(
        fsl.BET(frac=0.2, mask=True, functional=False),
        name=f'skullstrip_first_pass_{pipe_num}')

    bet_dilate = pe.Node(
        fsl.DilateImage(operation='max', kernel_shape='sphere',
                        kernel_size=6.0, internal_datatype='char'),
        name=f'skullstrip_first_dilate_{pipe_num}')

    bet_mask = pe.Node(fsl.ApplyMask(), name=f'skullstrip_first_mask_'
                                             f'{pipe_num}')

    unifize = pe.Node(afni_utils.Unifize(t2=True, outputtype='NIFTI_GZ',
                                         args='-clfrac 0.2 -rbt 18.3 65.0 90.0',
                                         out_file="uni.nii.gz"),
                      name=f'unifize_{pipe_num}')

    skullstrip_second_pass = pe.Node(
        preprocess.Automask(dilate=1, outputtype='NIFTI_GZ'),
        name=f'skullstrip_second_pass_{pipe_num}')

    combine_masks = pe.Node(fsl.BinaryMaths(operation='mul'),
                            name=f'combine_masks_{pipe_num}')

    apply_mask = pe.Node(fsl.ApplyMask(), 
                         name=f'extract_ref_brain_bold_{pipe_num}')

    node, out = strat_pool.get_data(["motion-basefile"])

    wf.connect([(node, init_aff, [(out, "moving_image")]),
                (node, map_brainmask, [(out, "reference_image")]),
                (node, norm, [(out, "moving_image")]),
                (node, restore_header, [(out, "refimage")]),
                (init_aff, norm, [
                    ("output_transform", "initial_moving_transform")]),
                (norm, map_brainmask, [
                    ("reverse_invert_flags", "invert_transform_flags"),
                    ("reverse_transforms", "transforms"),
                ]),
                (map_brainmask, binarize_mask, [("output_image", "in_file")]),
                (binarize_mask, pre_dilate, [("out_file", "in_file")]),
                (pre_dilate, restore_header, [
                    ("out_file", "imagetocopyrefimageinfoto")]),
                (restore_header, n4_correct, [("imageout", "mask_image")]),
                (node, n4_correct, [(out, "input_image")]),
                (n4_correct, skullstrip_first_pass,
                 [('output_image', 'in_file')]),
                (skullstrip_first_pass, bet_dilate,
                 [('mask_file', 'in_file')]),
                (bet_dilate, bet_mask, [('out_file', 'mask_file')]),
                (skullstrip_first_pass, bet_mask, [('out_file', 'in_file')]),
                (bet_mask, unifize, [('out_file', 'in_file')]),
                (unifize, skullstrip_second_pass, [('out_file', 'in_file')]),
                (skullstrip_first_pass, combine_masks,
                 [('mask_file', 'in_file')]),
                (skullstrip_second_pass, combine_masks,
                 [('out_file', 'operand_file')]),
                (unifize, apply_mask, [('out_file', 'in_file')]),
                (combine_masks, apply_mask, [('out_file', 'mask_file')]),
                ])

    outputs = {
        'space-bold_desc-brain_mask': (combine_masks, 'out_file'),
        'desc-ref_bold': (apply_mask, 'out_file')
    }

    return (wf, outputs)


def bold_mask_anatomical_refined(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "bold_mask_anatomical_refined",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "Anatomical_Refined",
     "inputs": ["bold",
                ["desc-preproc_bold", "bold"],
                ("desc-brain_T1w",
                 ["space-T1w_desc-brain_mask",
                  "space-T1w_desc-acpcbrain_mask"])],
     "outputs": ["space-bold_desc-brain_mask"]}
    '''

    # binarize anat mask, in case it is not a binary mask.
    anat_brain_mask_bin = pe.Node(interface=fsl.ImageMaths(),
                                  name=f'anat_brain_mask_bin_{pipe_num}')
    anat_brain_mask_bin.inputs.op_string = '-bin'

    node, out = strat_pool.get_data(['space-T1w_desc-brain_mask',
                                     'space-T1w_desc-acpcbrain_mask'])
    wf.connect(node, out, anat_brain_mask_bin, 'in_file')

    # fill holes of anat mask
    anat_mask_filled = pe.Node(interface=afni.MaskTool(),
                               name=f'anat_brain_mask_filled_{pipe_num}')
    anat_mask_filled.inputs.fill_holes = True
    anat_mask_filled.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(anat_brain_mask_bin, 'out_file',
               anat_mask_filled, 'in_file')

    # init_bold_mask : input raw func
    init_bold_mask = anat_refined_mask(init_bold_mask=True,
                                       wf_name=f'init_bold_mask_{pipe_num}')

    func_deoblique = pe.Node(interface=afni_utils.Refit(),
                             name=f'raw_func_deoblique_{pipe_num}')
    func_deoblique.inputs.deoblique = True

    node, out = strat_pool.get_data('bold')
    wf.connect(node, out, func_deoblique, 'in_file')

    func_reorient = pe.Node(interface=afni_utils.Resample(),
                            name=f'raw_func_reorient_{pipe_num}',
                            mem_gb=0,
                            mem_x=(0.0115, 'in_file', 't'))

    func_reorient.inputs.orientation = 'RPI'
    func_reorient.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(func_deoblique, 'out_file',
               func_reorient, 'in_file')

    wf.connect(func_reorient, 'out_file',
               init_bold_mask, 'inputspec.func')

    wf.connect(anat_mask_filled, 'out_file',
               init_bold_mask, 'inputspec.anatomical_brain_mask')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, init_bold_mask, 'inputspec.anat_brain')

    # dilate init func brain mask
    func_tmp_brain_mask = pe.Node(interface=fsl.ImageMaths(),
                                  name=f'func_tmp_brain_mask_dil_{pipe_num}')
    func_tmp_brain_mask.inputs.op_string = '-dilM'

    wf.connect(init_bold_mask, 'outputspec.func_brain_mask',
               func_tmp_brain_mask, 'in_file')

    # refined_bold_mask : input motion corrected func
    refined_bold_mask = anat_refined_mask(init_bold_mask=False,
                                          wf_name='refined_bold_mask'
                                                  f'_{pipe_num}')

    node, out = strat_pool.get_data(["desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, refined_bold_mask, 'inputspec.func')

    node, out = strat_pool.get_data('desc-brain_T1w')
    wf.connect(node, out, refined_bold_mask, 'inputspec.anat_brain')

    wf.connect(func_tmp_brain_mask, 'out_file',
               refined_bold_mask, 'inputspec.init_func_brain_mask')

    # dilate anatomical mask
    if cfg.functional_preproc['func_masking']['Anatomical_Refined'][
        'anatomical_mask_dilation']:
        anat_mask_dilate = pe.Node(interface=afni.MaskTool(),
                                   name=f'anat_mask_dilate_{pipe_num}')
        anat_mask_dilate.inputs.dilate_inputs = '1'
        anat_mask_dilate.inputs.outputtype = 'NIFTI_GZ'

        wf.connect(anat_mask_filled, 'out_file',
                   anat_mask_dilate, 'in_file')
        wf.connect(anat_mask_dilate, 'out_file',
                   refined_bold_mask, 'inputspec.anatomical_brain_mask')

    else:
        wf.connect(anat_mask_filled, 'out_file',
                   refined_bold_mask, 'inputspec.anatomical_brain_mask')

    # get final func mask
    func_mask_final = pe.Node(interface=fsl.MultiImageMaths(),
                              name=f'func_mask_final_{pipe_num}')
    func_mask_final.inputs.op_string = "-mul %s"

    wf.connect(func_tmp_brain_mask, 'out_file',
               func_mask_final, 'in_file')

    wf.connect(refined_bold_mask, 'outputspec.func_brain_mask',
               func_mask_final, 'operand_files')

    outputs = {
        'space-bold_desc-brain_mask': (func_mask_final, 'out_file')
    }

    return (wf, outputs)


def bold_mask_anatomical_based(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Generate the BOLD mask by basing it off of the anatomical brain mask.
    Adapted from DCAN Lab's BOLD mask method from the ABCD pipeline.
        https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased.sh

    Node Block:
    {"name": "bold_mask_anatomical_based",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "Anatomical_Based",
     "inputs": [["desc-preproc_bold", "bold"],
                "desc-brain_T1w",
                ["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-bold_desc-brain_mask"]}
    '''

    # 0. Take single volume of func
    func_single_volume = pe.Node(interface=afni.Calc(),
                                 name='func_single_volume')

    func_single_volume.inputs.set(
        expr='a',
        single_idx=1,
        outputtype='NIFTI_GZ'
    )

    node, out = strat_pool.get_data(["desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, func_single_volume, 'in_file_a')

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

    node, out = strat_pool.get_data(["desc-preproc_T1w", "desc-reorient_T1w",
                                     "T1w"])
    wf.connect(node, out, linear_reg_func_to_anat, 'reference')

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

    node, out = strat_pool.get_data("desc-brain_T1w")
    wf.connect(node, out, reg_anat_brain_to_func, 'in_file')

    node, out = strat_pool.get_data(["desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, reg_anat_brain_to_func, 'ref_file')

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

    outputs = {
        'space-bold_desc-brain_mask': (func_mask_fill_holes, 'out_file')
    }

    return (wf, outputs)


def bold_mask_anatomical_resampled(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Resample anatomical brain mask in standard space to get BOLD brain mask in standard space
    Adapted from DCAN Lab's BOLD mask method from the ABCD pipeline.
        https://github.com/DCAN-Labs/DCAN-HCP/blob/master/fMRIVolume/scripts/OneStepResampling.sh#L121-L132

    Node Block:
    {"name": "bold_mask_anatomical_resampled",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "Anatomical_Resampled",
     "inputs": [["desc-preproc_bold", "bold"],
                "T1w-template-funcreg",
                "space-template_desc-brain_T1w",
                "space-template_desc-T1w_mask"],
     "outputs": ["space-template_res-bold_desc-brain_T1w",
                 "space-template_desc-bold_mask",
                 "space-bold_desc-brain_mask"]}
    '''

    # applywarp --rel --interp=spline -i ${T1wImage} -r ${ResampRefIm} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${T1wImageFile}.${FinalfMRIResolution}
    anat_brain_to_func_res = pe.Node(interface=fsl.ApplyWarp(), 
                                     name=f'resample_anat_brain_in_standard_{pipe_num}')
    
    anat_brain_to_func_res.inputs.interp = 'spline'
    anat_brain_to_func_res.inputs.premat = cfg.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    node, out = strat_pool.get_data('space-template_desc-brain_T1w')
    wf.connect(node, out, anat_brain_to_func_res, 'in_file')

    node, out = strat_pool.get_data('T1w-template-funcreg')
    wf.connect(node, out, anat_brain_to_func_res, 'ref_file')

    # Create brain masks in this space from the FreeSurfer output (changing resolution)
    # applywarp --rel --interp=nn -i ${FreeSurferBrainMask}.nii.gz -r ${WD}/${T1wImageFile}.${FinalfMRIResolution} --premat=$FSLDIR/etc/flirtsch/ident.mat -o ${WD}/${FreeSurferBrainMaskFile}.${FinalfMRIResolution}.nii.gz
    anat_brain_mask_to_func_res = pe.Node(interface=fsl.ApplyWarp(),
                                          name=f'resample_anat_brain_mask_in_standard_{pipe_num}')
    
    anat_brain_mask_to_func_res.inputs.interp = 'nn'
    anat_brain_mask_to_func_res.inputs.premat = cfg.registration_workflows[
        'anatomical_registration']['registration']['FSL-FNIRT']['identity_matrix']

    node, out = strat_pool.get_data('space-template_desc-T1w_mask')
    wf.connect(node, out, anat_brain_mask_to_func_res, 'in_file')

    wf.connect(anat_brain_to_func_res, 'out_file',
               anat_brain_mask_to_func_res, 'ref_file')

    # Resample func mask in template space back to native space
    func_mask_template_to_native = pe.Node(
        interface=afni.Resample(),
        name=f'resample_func_mask_to_native_{pipe_num}',
        mem_gb=0,
        mem_x=(0.0115, 'in_file', 't'))
    func_mask_template_to_native.inputs.resample_mode = 'NN'
    func_mask_template_to_native.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(anat_brain_mask_to_func_res, 'out_file',
        func_mask_template_to_native, 'in_file')

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, func_mask_template_to_native, 'master')

    outputs = {
        'space-template_res-bold_desc-brain_T1w': (anat_brain_to_func_res, 'out_file'),
        'space-template_desc-bold_mask': (anat_brain_mask_to_func_res, 'out_file'),
        'space-bold_desc-brain_mask': (func_mask_template_to_native, 'out_file')
    }

    return (wf, outputs)


def bold_mask_ccs(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Generate the BOLD mask by basing it off of the anatomical brain.
    Adapted from the BOLD mask method from the CCS pipeline.
        https://github.com/TingsterX/CCS/blob/master/ccs_01_funcpreproc.sh#L89-L110
    Node Block:
    {"name": "bold_mask_ccs",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": ["func_masking", "using"],
     "option_val": "CCS_Anatomical_Refined",
     "inputs": [["desc-motion_bold", "desc-preproc_bold", "bold"],
                "desc-brain_T1w",
                ["desc-preproc_T1w", "desc-reorient_T1w", "T1w"]],
     "outputs": ["space-bold_desc-brain_mask",
                "desc-ROIbrain_bold"]}
    '''

    # Run 3dAutomask to generate func initial mask
    func_tmp_brain_mask = pe.Node(interface=preprocess.Automask(),
                                  name=f'func_tmp_brain_mask_AFNI_{pipe_num}')
    func_tmp_brain_mask.inputs.dilate = 1
    func_tmp_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(["desc-motion_bold", 
                                     "desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, func_tmp_brain_mask, 'in_file')

    # Extract 8th volume as func ROI
    func_roi = pe.Node(interface=fsl.ExtractROI(),
                       name=f'extract_func_roi_{pipe_num}')
    func_roi.inputs.t_min = 7
    func_roi.inputs.t_size = 1

    node, out = strat_pool.get_data(["desc-motion_bold", 
                                     "desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, func_roi, 'in_file')

    # Apply func initial mask on func ROI volume
    func_tmp_brain = pe.Node(interface=fsl.maths.ApplyMask(),
                             name=f'get_func_tmp_brain_{pipe_num}')

    wf.connect(func_roi, 'roi_file',
               func_tmp_brain, 'in_file')

    wf.connect(func_tmp_brain_mask, 'out_file',
               func_tmp_brain, 'mask_file')

    # Register func tmp brain to anat brain to get func2anat matrix
    reg_func_to_anat = pe.Node(interface=fsl.FLIRT(),
                               name=f'func_to_anat_linear_reg_{pipe_num}')
    reg_func_to_anat.inputs.interp = 'trilinear'
    reg_func_to_anat.inputs.cost = 'corratio'
    reg_func_to_anat.inputs.dof = 6

    wf.connect(func_tmp_brain, 'out_file',
               reg_func_to_anat, 'in_file')

    node, out = strat_pool.get_data("desc-brain_T1w")
    wf.connect(node, out, reg_func_to_anat, 'reference')

    # Inverse func2anat matrix
    inv_func_to_anat_affine = pe.Node(interface=fsl.ConvertXFM(),
                                      name=f'inv_func2anat_affine_{pipe_num}')
    inv_func_to_anat_affine.inputs.invert_xfm = True

    wf.connect(reg_func_to_anat, 'out_matrix_file',
               inv_func_to_anat_affine, 'in_file')

    # Transform anat brain to func space
    reg_anat_brain_to_func = pe.Node(interface=fsl.FLIRT(),
                                     name=f'reg_anat_brain_to_func_{pipe_num}')
    reg_anat_brain_to_func.inputs.apply_xfm = True
    reg_anat_brain_to_func.inputs.interp = 'trilinear'

    node, out = strat_pool.get_data("desc-brain_T1w")
    wf.connect(node, out, reg_anat_brain_to_func, 'in_file')

    wf.connect(func_roi, 'roi_file',
               reg_anat_brain_to_func, 'reference')

    wf.connect(inv_func_to_anat_affine, 'out_file',
               reg_anat_brain_to_func, 'in_matrix_file')

    # Binarize and dilate anat brain in func space
    bin_anat_brain_in_func = pe.Node(interface=fsl.ImageMaths(),
                                     name=f'bin_anat_brain_in_func_{pipe_num}')
    bin_anat_brain_in_func.inputs.op_string = '-bin -dilM'

    wf.connect(reg_anat_brain_to_func, 'out_file', 
               bin_anat_brain_in_func, 'in_file')

    # Binarize detectable func signals
    bin_func = pe.Node(interface=fsl.ImageMaths(),
                                     name=f'bin_func_{pipe_num}')
    bin_func.inputs.op_string = '-Tstd -bin'

    node, out = strat_pool.get_data(["desc-motion_bold", 
                                     "desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, bin_func, 'in_file')

    # Take intersection of masks
    merge_func_mask = pe.Node(util.Merge(2), 
                              name=f'merge_func_mask_{pipe_num}')

    wf.connect(func_tmp_brain_mask, 'out_file',
               merge_func_mask, 'in1')

    wf.connect(bin_anat_brain_in_func, 'out_file',
               merge_func_mask, 'in2')

    intersect_mask = pe.Node(interface=fsl.MultiImageMaths(),
                             name=f'intersect_mask_{pipe_num}')
    intersect_mask.inputs.op_string = '-mul %s -mul %s'
    intersect_mask.inputs.output_datatype = 'char'

    wf.connect(bin_func, 'out_file', 
               intersect_mask, 'in_file')

    wf.connect(merge_func_mask, 'out', 
               intersect_mask, 'operand_files')

    # this is the func input for coreg in ccs
    # TODO evaluate if it's necessary to use this brain
    example_func_brain = pe.Node(interface=fsl.maths.ApplyMask(),
                                 name=f'get_example_func_brain_{pipe_num}')

    wf.connect(func_roi, 'roi_file',
               example_func_brain, 'in_file')

    wf.connect(intersect_mask, 'out_file',
               example_func_brain, 'mask_file')

    outputs = {
        'space-bold_desc-brain_mask': (intersect_mask, 'out_file'),
        'desc-ROIbrain_bold': (example_func_brain, 'out_file')
    }

    return (wf, outputs)


def bold_masking(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "bold_masking",
     "config": None,
     "switch": [["functional_preproc", "run"],
                ["functional_preproc", "func_masking", "apply_func_mask_in_native_space"]],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-preproc_bold", "bold"],
                 "space-bold_desc-brain_mask")],
     "outputs": {
         "desc-preproc_bold": {
             "Description": "The skull-stripped BOLD time-series.",
             "SkullStripped": True},
         "desc-brain_bold": {
             "Description": "The skull-stripped BOLD time-series.",
             "SkullStripped": True}}
    }
    '''

    func_edge_detect = pe.Node(interface=afni_utils.Calc(),
                               name=f'func_extract_brain_{pipe_num}')

    func_edge_detect.inputs.expr = 'a*b'
    func_edge_detect.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(["desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, func_edge_detect, 'in_file_a')

    node, out = strat_pool.get_data("space-bold_desc-brain_mask")
    wf.connect(node, out, func_edge_detect, 'in_file_b')

    outputs = {
        'desc-preproc_bold': (func_edge_detect, 'out_file'),
        'desc-brain_bold': (func_edge_detect, 'out_file')
    }

    return (wf, outputs)


def func_mean(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_mean",
     "config": "None",
     "switch": [["functional_preproc", "run"],
                ["functional_preproc", "generate_func_mean", "run"]],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": ["desc-mean_bold"]
    }
    '''

    func_mean = pe.Node(interface=afni_utils.TStat(),
                        name=f'func_mean_{pipe_num}')

    func_mean.inputs.options = '-mean'
    func_mean.inputs.outputtype = 'NIFTI_GZ'

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, func_mean, 'in_file')

    outputs = {
        'desc-mean_bold': (func_mean, 'out_file')
    }

    return (wf, outputs)


def func_normalize(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_normalize",
     "config": "None",
     "switch": [["functional_preproc", "run"],
                ["functional_preproc", "normalize_func", "run"]],
     "option_key": "None",
     "option_val": "None",
     "inputs": [["desc-preproc_bold", "bold"]],
     "outputs": ["desc-preproc_bold"]}
    '''

    func_normalize = pe.Node(interface=fsl.ImageMaths(),
                             name=f'func_normalize_{pipe_num}',
                             mem_gb=0.7,
                             mem_x=(4538494663498653 /
                                    604462909807314587353088, 'in_file'))
    func_normalize.inputs.op_string = '-ing 10000'
    func_normalize.inputs.out_data_type = 'float'

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, func_normalize, 'in_file')

    outputs = {
        'desc-preproc_bold': (func_normalize, 'out_file')
    }

    return (wf, outputs)
    
    
def func_mask_normalize(wf, cfg, strat_pool, pipe_num, opt=None):
    '''
    {"name": "func_mask_normalize",
     "config": ["functional_preproc"],
     "switch": ["run"],
     "option_key": "None",
     "option_val": "None",
     "inputs": [(["desc-preproc_bold", "bold"],
                 "space-bold_desc-brain_mask")],
     "outputs": ["space-bold_desc-brain_mask"]}
    '''

    func_mask_normalize = pe.Node(interface=fsl.ImageMaths(),
                                  name=f'func_mask_normalize_{pipe_num}',
                                  mem_gb=0.7,
                                  mem_x=(4538494663498653 /
                                         604462909807314587353088, 'in_file'))
    func_mask_normalize.inputs.op_string = '-Tmin -bin'
    func_mask_normalize.inputs.out_data_type = 'char'

    node, out = strat_pool.get_data(["desc-preproc_bold", "bold"])
    wf.connect(node, out, func_mask_normalize, 'in_file')

    outputs = {
        'space-bold_desc-brain_mask': (func_mask_normalize, 'out_file')
    }

    return (wf, outputs)

