import os
from CPAC.pipeline import nipype_pipeline_engine as pe
from CPAC.pipeline.nodeblock import nodeblock
import nipype.algorithms.rapidart as ra
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from .utils import *
from CPAC.vmhc import *
from nipype.interfaces.afni import preprocess
from CPAC.registration.registration import apply_transform
from CPAC.image_utils import spatial_smoothing

from CPAC.utils.utils import check_prov_for_regtool


@nodeblock(
    name="smooth_func_vmhc",
    switch=["voxel_mirrored_homotopic_connectivity", "run"],
    option_key=["post_processing", "spatial_smoothing", "smoothing_method"],
    option_val=["AFNI", "FSL"],
    inputs=[
        ["desc-cleaned_bold", "desc-brain_bold", "desc-preproc_bold", "bold"],
        "space-bold_desc-brain_mask",
    ],
    outputs=["desc-sm_bold", "fwhm"],
)
def smooth_func_vmhc(wf, cfg, strat_pool, pipe_num, opt=None):
    fwhm = cfg.post_processing['spatial_smoothing']['fwhm']

    smooth = spatial_smoothing(f'smooth_symmetric_{pipe_num}',
                               fwhm, opt=opt)

    node, out = strat_pool.get_data(["desc-cleaned_bold",
                                     "desc-brain_bold",
                                     "desc-preproc_bold",
                                     "bold"])
    wf.connect(node, out, smooth, 'inputspec.in_file')

    node, out = strat_pool.get_data("space-bold_desc-brain_mask")
    wf.connect(node, out, smooth, 'inputspec.mask')

    # 'fwhm' output for iterable
    outputs = {
        "desc-sm_bold": (smooth, 'outputspec.out_file'),
        "fwhm": (smooth, 'fwhm_input.fwhm')
    }

    return (wf, outputs)


@nodeblock(
    name="transform_timeseries_to_sym_template",
    config=["voxel_mirrored_homotopic_connectivity"],
    switch=["run"],
    inputs=[
        [
            "desc-cleaned-sm_bold",
            "desc-brain-sm_bold",
            "desc-preproc-sm_bold",
            "desc-sm_bold",
        ],
        "from-bold_to-symtemplate_mode-image_xfm",
        "T1w-brain-template-symmetric",
    ],
    outputs=["space-symtemplate_desc-sm_bold"],
)
def warp_timeseries_to_sym_template(wf, cfg, strat_pool, pipe_num, opt=None):

    xfm_prov = strat_pool.get_cpac_provenance(
        'from-bold_to-symtemplate_mode-image_xfm')
    reg_tool = check_prov_for_regtool(xfm_prov)

    num_cpus = cfg.pipeline_setup['system_config'][
        'max_cores_per_participant']

    num_ants_cores = cfg.pipeline_setup['system_config']['num_ants_threads']

    apply_xfm = apply_transform(f'warp_ts_to_sym_template_{pipe_num}',
                                reg_tool, time_series=True, num_cpus=num_cpus,
                                num_ants_cores=num_ants_cores)

    if reg_tool == 'ants':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'ANTs_pipelines']['interpolation']
    elif reg_tool == 'fsl':
        apply_xfm.inputs.inputspec.interpolation = cfg.registration_workflows[
            'functional_registration']['func_registration_to_template'][
            'FNIRT_pipelines']['interpolation']

    # smoothed BOLD
    connect, resource = strat_pool.get_data(["desc-cleaned-sm_bold",
                                             "desc-brain-sm_bold",
                                             "desc-preproc-sm_bold",
                                             "desc-sm_bold"],
                                            report_fetched=True)
    node, out = connect
    wf.connect(node, out, apply_xfm, 'inputspec.input_image')

    node, out = strat_pool.get_data("T1w-brain-template-symmetric")
    wf.connect(node, out, apply_xfm, 'inputspec.reference')

    node, out = strat_pool.get_data("from-bold_to-symtemplate_mode-image_xfm")
    wf.connect(node, out, apply_xfm, 'inputspec.transform')

    outputs = {
        f'space-symtemplate_{resource}':
            (apply_xfm, 'outputspec.output_image')
    }

    return (wf, outputs)


@nodeblock(
    name="vmhc",
    config=["voxel_mirrored_homotopic_connectivity"],
    switch=["run"],
    inputs=[
        [
            "space-symtemplate_desc-cleaned-sm_bold",
            "space-symtemplate_desc-brain-sm_bold",
            "space-symtemplate_desc-preproc-sm_bold",
            "space-symtemplate_desc-sm_bold",
        ]
    ],
    outputs=["vmhc"],
)
def vmhc(wf, cfg, strat_pool, pipe_num, opt=None):
    '''Compute Voxel-Mirrored Homotopic Connectivity.

    VMHC is the map of brain functional homotopy, the high degree of
    synchrony in spontaneous activity between geometrically corresponding
    interhemispheric (i.e., homotopic) regions.
    '''

    # write out a swapped version of the file
    # copy and L/R swap file
    copy_and_L_R_swap = pe.Node(interface=fsl.SwapDimensions(),
                                name=f'copy_and_L_R_swap_{pipe_num}',
                                mem_gb=3.0)

    copy_and_L_R_swap.inputs.new_dims = ('-x', 'y', 'z')

    node, out = strat_pool.get_data(["space-symtemplate_desc-cleaned-sm_bold",
                                     "space-symtemplate_desc-brain-sm_bold",
                                     "space-symtemplate_desc-preproc-sm_bold",
                                     "space-symtemplate_desc-sm_bold"])
    wf.connect(node, out, copy_and_L_R_swap, 'in_file')

    # calculate correlation between original and swapped images
    pearson_correlation = pe.Node(interface=preprocess.TCorrelate(),
                                  name=f'pearson_correlation_{pipe_num}',
                                  mem_gb=3.0)

    pearson_correlation.inputs.pearson = True
    pearson_correlation.inputs.polort = -1
    pearson_correlation.inputs.outputtype = 'NIFTI_GZ'

    wf.connect(node, out, pearson_correlation, 'xset')

    wf.connect(copy_and_L_R_swap, 'out_file',
               pearson_correlation, 'yset')

    outputs = {
        'vmhc': (pearson_correlation, 'out_file')
    }

    return (wf, outputs)
