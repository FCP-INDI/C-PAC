# Copyright (C) 2012-2023  C-PAC Developers

# This file is part of C-PAC.

# C-PAC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# C-PAC is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
# License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with C-PAC. If not, see <https://www.gnu.org/licenses/>.
"""Functions for calculating motion parameters"""
# pylint: disable=ungrouped-imports,wrong-import-order,wrong-import-position
from CPAC.pipeline import nipype_pipeline_engine as pe
from nipype.interfaces.afni import preprocess
from nipype.interfaces import afni, fsl, utility as util
from nipype.interfaces.afni import utils as afni_utils
from CPAC.func_preproc.utils import chunk_ts, oned_text_concat, \
                                    split_ts_chunks
from CPAC.func_preproc.utils import notch_filter_motion
from CPAC.generate_motion_statistics import affine_file_from_params_file, \
                                            motion_power_statistics
from CPAC.pipeline.nodeblock import nodeblock
from CPAC.pipeline.schema import valid_options
from CPAC.utils.interfaces.function import Function
from CPAC.utils.utils import check_prov_for_motion_tool


@nodeblock(
    name='calc_motion_stats',
    switch=[['functional_preproc', 'run'],
            ['functional_preproc', 'motion_estimates_and_correction', 'run']],
    inputs=[('desc-preproc_bold', 'space-bold_desc-brain_mask',
             'desc-movementParameters_motion', 'max-displacement',
             'rels-displacement', 'filtered-coordinate-transformation',
             'coordinate-transformation'), 'subject', 'scan'],
    outputs=['dvars', 'framewise-displacement-power',
             'framewise-displacement-jenkinson', 'power-params',
             'motion-params', 'motion', 'desc-summary_motion'])
def calc_motion_stats(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Calculate motion statistics for motion parameters.'''
    motion_prov = strat_pool.get_cpac_provenance(
        'desc-movementParameters_motion')
    motion_correct_tool = check_prov_for_motion_tool(motion_prov)
    coordinate_transformation = ['filtered-coordinate-transformation',
                                 'coordinate-transformation']
    gen_motion_stats = motion_power_statistics(
        name=f'gen_motion_stats_{pipe_num}',
        motion_correct_tool=motion_correct_tool,
        filtered=strat_pool.filtered_movement)

    # Special case where the workflow is not getting outputs from
    # resource pool but is connected to functional datasource
    wf.connect(*strat_pool.get_data('desc-preproc_bold'),
               gen_motion_stats, 'inputspec.motion_correct')
    wf.connect(*strat_pool.get_data('space-bold_desc-brain_mask'),
               gen_motion_stats, 'inputspec.mask')
    wf.connect(*strat_pool.get_data('desc-movementParameters_motion'),
               gen_motion_stats, 'inputspec.movement_parameters')
    wf.connect(*strat_pool.get_data('max-displacement'),
               gen_motion_stats, 'inputspec.max_displacement')

    if strat_pool.check_rpool('rels-displacement'):
        wf.connect(*strat_pool.get_data('rels-displacement'),
                   gen_motion_stats, 'inputspec.rels_displacement')
    if strat_pool.check_rpool(coordinate_transformation):
        wf.connect(*strat_pool.get_data(coordinate_transformation),
                   gen_motion_stats, 'inputspec.transformations')

    outputs = {
        'framewise-displacement-power':
            (gen_motion_stats, 'outputspec.FDP_1D'),
        'framewise-displacement-jenkinson':
            (gen_motion_stats, 'outputspec.FDJ_1D'),
        'dvars': (gen_motion_stats, 'outputspec.DVARS_1D'),
        'power-params': (gen_motion_stats, 'outputspec.power_params'),
        'motion-params': (gen_motion_stats, 'outputspec.motion_params'),
        'motion': (gen_motion_stats, 'outputspec.motion'),
        'desc-summary_motion': (gen_motion_stats,
                                'outputspec.desc-summary_motion')}

    return wf, outputs


def estimate_reference_image(in_file):
    """fMRIPrep-style BOLD reference
    Ref: https://github.com/nipreps/niworkflows/blob/maint/1.3.x/niworkflows/interfaces/registration.py#L446-L549
    """
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
    os.system('3dvolreg -Fourier -twopass -zpad 4 '
              f'-prefix {mc_out_file} {ref_name}')

    mc_slice_nii = nb.load(mc_out_file)

    median_image_data = np.median(mc_slice_nii.get_fdata(), axis=3)

    nb.Nifti1Image(median_image_data, ref_im.affine, ref_im.header
                   ).to_filename(out_file)

    return out_file


_MOTION_CORRECTED_OUTPUTS = {
    "desc-preproc_bold": {"Description": "Motion-corrected BOLD time-series."},
    "desc-motion_bold": {"Description": "Motion-corrected BOLD time-series."}}
# the "filtered" outputs here are just for maintaining expecting
# forking and connections and will not be output
_MOTION_PARAM_OUTPUTS = {
    "max-displacement": {},
    "rels-displacement": {},
    "desc-movementParameters_motion": {
        "Description": "Each line contains for one timepoint a 6-DOF "
                       "rigid transform parameters in the format "
                       "defined by AFNI's 3dvolreg: [roll, pitch, yaw, "
                       "superior displacement, left displacement, "
                       "posterior displacement]. Rotation parameters "
                       "are in degrees counterclockwise, and translation "
                       "parameters are in millimeters."
    },
    "filtered-coordinate-transformation": {"Description": "UNFILTERED"},
    "coordinate-transformation": {
        "Description" : "Each row contains for one timepoint the first "
                        "12 values of a 4x4 affine matrix"}}


@nodeblock(
    name="motion_correction",
    switch=["functional_preproc", "motion_estimates_and_correction", "run"],
    option_key=["functional_preproc", "motion_estimates_and_correction",
                "motion_correction", "using"],
    option_val=["3dvolreg", "mcflirt"],
    inputs=[("desc-preproc_bold", "motion-basefile")],
    outputs={**_MOTION_CORRECTED_OUTPUTS, **_MOTION_PARAM_OUTPUTS})
def func_motion_correct(wf, cfg, strat_pool, pipe_num, opt=None):
    wf, outputs = motion_correct_connections(wf, cfg, strat_pool, pipe_num,
                                             opt)

    return wf, outputs


@nodeblock(
    name="motion_correction_only",
    switch=["functional_preproc", "motion_estimates_and_correction", "run"],
    option_key=["functional_preproc", "motion_estimates_and_correction",
                "motion_correction", "using"],
    option_val=["3dvolreg", "mcflirt"],
    inputs=[("desc-preproc_bold", "motion-basefile")],
    outputs=_MOTION_CORRECTED_OUTPUTS)
def func_motion_correct_only(wf, cfg, strat_pool, pipe_num, opt=None):
    wf, wf_outputs = motion_correct_connections(wf, cfg, strat_pool, pipe_num,
                                                opt)

    outputs = {
        'desc-preproc_bold': wf_outputs['desc-motion_bold'],
        'desc-motion_bold': wf_outputs['desc-motion_bold']
    }

    return (wf, outputs)


@nodeblock(
    name="motion_estimates",
    switch=["functional_preproc", "motion_estimates_and_correction", "run"],
    option_key=["functional_preproc", "motion_estimates_and_correction",
                "motion_correction", "using"],
    option_val=["3dvolreg", "mcflirt"],
    inputs=[("desc-preproc_bold", "motion-basefile")],
    outputs=_MOTION_PARAM_OUTPUTS)
def func_motion_estimates(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Calculate motion estimates using 3dVolReg or MCFLIRT.'''
    from CPAC.pipeline.utils import present_outputs
    wf, wf_outputs = motion_correct_connections(wf, cfg, strat_pool, pipe_num,
                                                opt)
    return (wf, present_outputs(wf_outputs,
                                ['coordinate-transformation',
                                 'filtered-coordinate-transformation',
                                 'max-displacement',
                                 'desc-movementParameters_motion',
                                 'rels-displacement']))


def get_mcflirt_rms_abs(rms_files):
    for path in rms_files:
        if 'abs.rms' in path:
            abs_file = path
        if 'rel.rms' in path:
            rels_file = path
    return abs_file, rels_file


@nodeblock(
    name='get_motion_ref',
    switch=['functional_preproc', 'motion_estimates_and_correction', 'run'],
    option_key=['functional_preproc', 'motion_estimates_and_correction',
                'motion_correction', 'motion_correction_reference'],
    option_val=['mean', 'median', 'selected_volume', 'fmriprep_reference'],
    inputs=['desc-preproc_bold', 'desc-reorient_bold'],
    outputs=['motion-basefile']
)
def get_motion_ref(wf, cfg, strat_pool, pipe_num, opt=None):
    if opt not in get_motion_ref.option_val:
        raise ValueError('\n\n[!] Error: The \'motion_correction_reference\' '
                         'parameter of the \'motion_correction\' workflow '
                         'must be one of:\n\t{0}.\n\nTool input: \'{1}\''
                         '\n\n'.format(
                             ' or '.join([f"'{val}'" for val in
                                          get_motion_ref.option_val]),
                             opt))

    if opt == 'mean':
        func_get_RPI = pe.Node(interface=afni_utils.TStat(),
                               name=f'func_get_mean_RPI_{pipe_num}',
                               mem_gb=0.48,
                               mem_x=(1435097126797993 /
                                      302231454903657293676544,
                                      'in_file'))

        func_get_RPI.inputs.options = '-mean'
        func_get_RPI.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data('desc-preproc_bold')
        wf.connect(node, out, func_get_RPI, 'in_file')

    elif opt == 'median':
        func_get_RPI = pe.Node(interface=afni_utils.TStat(),
                               name=f'func_get_median_RPI_{pipe_num}')

        func_get_RPI.inputs.options = '-median'
        func_get_RPI.inputs.outputtype = 'NIFTI_GZ'

        node, out = strat_pool.get_data('desc-preproc_bold')
        wf.connect(node, out, func_get_RPI, 'in_file')

    elif opt == 'selected_volume':
        func_get_RPI = pe.Node(interface=afni.Calc(),
                               name=f'func_get_selected_RPI_{pipe_num}')

        func_get_RPI.inputs.set(
            expr='a',
            single_idx=cfg.functional_preproc[
                'motion_estimates_and_correction'][
                'motion_correction']['motion_correction_reference_volume'],
            outputtype='NIFTI_GZ'
        )

        node, out = strat_pool.get_data('desc-preproc_bold')
        wf.connect(node, out, func_get_RPI, 'in_file_a')

    elif opt == 'fmriprep_reference':
        func_get_RPI = pe.Node(util.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=estimate_reference_image
                                             ),
                           name=f'func_get_fmriprep_ref_{pipe_num}')

        node, out = strat_pool.get_data('desc-reorient_bold')
        wf.connect(node, out, func_get_RPI, 'in_file')

    outputs = {
        'motion-basefile': (func_get_RPI, 'out_file')
    }

    return (wf, outputs)


def motion_correct_3dvolreg(wf, cfg, strat_pool, pipe_num):
    """Calculate motion parameters with 3dvolreg"""
    if int(cfg.pipeline_setup['system_config']['max_cores_per_participant']
           ) > 1:
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

        node, out = strat_pool.get_data("desc-preproc_bold")
        wf.connect(node, out, chunk, 'func_file')

        split_imports = ['import os', 'import subprocess']
        split = pe.Node(Function(input_names=['func_file',
                                              'tr_ranges'],
                                  output_names=['split_funcs'],
                                  function=split_ts_chunks,
                                  imports=split_imports),
                        name=f'split_{pipe_num}')

        node, out = strat_pool.get_data("desc-preproc_bold")
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

        node, out = strat_pool.get_data('desc-preproc_bold')
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

    func_motion_correct.inputs.zpad = 4
    func_motion_correct.inputs.outputtype = 'NIFTI_GZ'

    args = '-Fourier'
    if cfg.functional_preproc['motion_estimates_and_correction'][
        'motion_correction']['AFNI-3dvolreg']['functional_volreg_twopass']:
        args = f'-twopass {args}'

    func_motion_correct.inputs.args = args

    # Calculate motion parameters
    func_motion_correct_A = func_motion_correct.clone(
        f'func_motion_correct_3dvolreg_{pipe_num}')
    func_motion_correct_A.inputs.md1d_file = 'max_displacement.1D'
    func_motion_correct_A.inputs.args = args

    wf.connect(out_split_func, 'out_file',
                func_motion_correct_A, 'in_file')

    node, out = strat_pool.get_data('motion-basefile')
    wf.connect(node, out, func_motion_correct_A, 'basefile')

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

    outputs = {
        'desc-preproc_bold': (out_motion_A, 'out_file'),
        'desc-motion_bold': (out_motion_A, 'out_file'),
        'max-displacement': (out_md1d, 'out_file'),
        'desc-movementParameters_motion': (out_oned, 'out_file'),
        'coordinate-transformation': (out_oned_matrix, 'out_file'),
        'filtered-coordinate-transformation': (out_oned_matrix, 'out_file')
    }

    return wf, outputs


def motion_correct_mcflirt(wf, cfg, strat_pool, pipe_num):
    """Calculate motion parameters with MCFLIRT"""
    func_motion_correct_A = pe.Node(
        interface=fsl.MCFLIRT(save_mats=True, save_plots=True),
        name=f'func_motion_correct_mcflirt_{pipe_num}', mem_gb=2.5)

    func_motion_correct_A.inputs.save_mats = True
    func_motion_correct_A.inputs.save_plots = True
    func_motion_correct_A.inputs.save_rms = True

    node, out = strat_pool.get_data('desc-preproc_bold')
    wf.connect(node, out, func_motion_correct_A, 'in_file')

    node, out = strat_pool.get_data('motion-basefile')
    wf.connect(node, out, func_motion_correct_A, 'ref_file')

    normalize_motion_params = pe.Node(Function(
        input_names=['in_file'], output_names=['out_file'],
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
        'desc-movementParameters_motion': (normalize_motion_params,
                                           'out_file'),
        'coordinate-transformation': (func_motion_correct_A, 'mat_file'),
        'filtered-coordinate-transformation': (func_motion_correct_A,
                                               'mat_file')
    }

    return wf, outputs


motion_correct = {'3dvolreg': motion_correct_3dvolreg,
                  'mcflirt': motion_correct_mcflirt}


def motion_correct_connections(wf, cfg, strat_pool, pipe_num, opt):
    """Check opt for valid option, then connect that option."""
    motion_correct_options = valid_options['motion_correction']
    if opt not in motion_correct_options:
        raise KeyError("\n\n[!] Error: The 'tool' parameter of the "
                       "'motion_correction' workflow must be one of "
                       f"{str(motion_correct_options).strip('[{()}]')}"
                       f".\n\nTool input: {opt}\n\n")
    return motion_correct[opt](wf, cfg, strat_pool, pipe_num)


@nodeblock(
    name="motion_estimate_filter",
    config=["functional_preproc", "motion_estimates_and_correction",
            "motion_estimate_filter"],
    switch=["run"],
    option_key="filters",
    option_val="USER-DEFINED",
    inputs=[("desc-preproc_bold",
             "space-bold_desc-brain_mask",
             "max-displacement",
             "rels-displacement",
             "coordinate-transformation",
             "desc-movementParameters_motion"),
            "TR"],
     outputs={"filtered-coordinate-transformation": {
                  "Description": "Affine matrix regenerated from"
                                  " filtered motion parameters. Note:"
                                  " the translation vector does not"
                                  " account for recentering inherent"
                                  " in rotation; this omission does"
                                  " not seem to affect framewise"
                                  " displacement calculation, for which"
                                  " this matrix is used."},
              "desc-movementParameters_motion": {
                  "Description": "Filtered movement parameters"
                                 " (3 rotation, 3 translation)."},
              "desc-movementParametersUnfiltered_motion": {
                  "Description": "Unfiltered movement parameters"
                                 " (3 rotation, 3 translation)."},
              "motion-filter-info": {},
              "motion-filter-plot": {}})
def motion_estimate_filter(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Filter motion parameters.

    .. versionchanged:: 1.8.6
       Beginning with version 1.8.6, C-PAC outputs both the unfiltered
       and the filtered motion parameters and uses the unfiltered
       parameters in QC. Previous versions only reported the filtered
       parameters and used the filtered parameters for QC.
    '''
    notch_imports = ['import os', 'import numpy as np',
                     'from scipy.signal import iirnotch, filtfilt, firwin, '
                     'freqz',
                     'from matplotlib import pyplot as plt',
                     'from CPAC.func_preproc.utils import degrees_to_mm, '
                     'mm_to_degrees']
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
                    name=f'filter_motion_params_{opt["Name"]}_{pipe_num}')

    notch.inputs.filter_type = opt.get('filter_type')
    notch.inputs.fc_RR_min = opt.get('breathing_rate_min')
    notch.inputs.fc_RR_max = opt.get('breathing_rate_max')
    notch.inputs.center_freq = opt.get('center_frequency')
    notch.inputs.freq_bw = opt.get('filter_bandwidth')
    notch.inputs.lowpass_cutoff = opt.get('lowpass_cutoff')
    notch.inputs.filter_order = opt.get('filter_order')

    movement_parameters = strat_pool.node_data(
        'desc-movementParameters_motion')
    wf.connect(movement_parameters.node, movement_parameters.out,
               notch, 'motion_params')

    node, out = strat_pool.get_data('TR')
    wf.connect(node, out, notch, 'TR')

    affine = pe.Node(Function(input_names=['params_file'],
                              output_names=['affine_file'],
                              function=affine_file_from_params_file),
                     name='affine_from_filtered_params_'
                          f'{opt["Name"]}_{pipe_num}')
    wf.connect(notch, 'filtered_motion_params', affine, 'params_file')

    outputs = {
        'filtered-coordinate-transformation': (affine, 'affine_file'),
        'motion-filter-info': (notch, 'filter_info'),
        'motion-filter-plot': (notch, 'filter_plot'),
        'desc-movementParameters_motion': (notch, 'filtered_motion_params')
    }

    if not cfg.switch_is_off(["functional_preproc",
                              "motion_estimates_and_correction",
                              "motion_estimate_filter", "run"]):
        outputs['desc-movementParametersUnfiltered_motion'] = (
            movement_parameters.node, movement_parameters.out)

    return (wf, outputs)


def normalize_motion_parameters(in_file):
    """Convert FSL mcflirt motion parameters to AFNI space"""
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

    out_file = os.path.join(os.getcwd(), 'motion_params.tsv')
    np.savetxt(out_file, motion_params)

    return out_file
