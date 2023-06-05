import os
import sys
import subprocess
import numpy as np
import nipype
import nibabel
from nipype.pipeline import engine as pe
from nipype.interfaces import fsl
from nipype.interfaces import utility
import nipype.interfaces.utility as util
from nipype.interfaces.fsl import utils as fsl_utils
from nipype.interfaces.fsl import maths as fsl_maths
from nipype.interfaces.fsl import ExtractROI
from nipype.interfaces.fsl import maths


def run_HCP_gradient_unwarp(phase_vol, input_coeffs):

    import os
    import subprocess

    trilinear = os.path.join(os.getcwd(), "trilinear.nii.gz")

    cmd = [
        "../gradunwarp/core/gradient_unwarp.py",
        phase_vol,
        trilinear,
        "siemens",
        "-g",
        input_coeffs,
        "-n",
    ]

    retcode = subprocess.check_output(cmd)
    abs_fullWarp = os.path.join(os.getcwd(), "fullWarp_abs.nii.gz")

    return trilinear, abs_fullWarp


def run_convertwarp(cw_trilinear, cw_fullWarp_abs):

    import os
    import subprocess

    out_file = os.path.join(os.getcwd(), "out_file.nii.gz")
    out_jac = os.path.join(os.getcwd(), "out_jac.nii.gz")

    cmd = [
        "convertwarp",
        "--abs",
        f"--ref={cw_trilinear}",
        f"--warp1={cw_fullWarp_abs}",
        "--relout",
        f"--out={out_file}",
        f"--j={jac_out}",
    ]
    retcode = subprocess.check_output(cmd)

    return out_file, out_jac


def gradient_distortion_correction(wf, inp_image, name):

    extract_ROI = pe.Node(interface=fsl.ExtractROI(), name=f"extract_ROI_{name}")
    extract_ROI.inputs.in_file = inp_image
    extract_ROI.inputs.t_min = 0
    extract_ROI.inputs.t_size = 1

    grad_unwarp_imports = ["import os", "import subprocess"]
    grad_unwarp = pe.Node(
        util.Function(
            input_names=["phase_vol", "input_coeffs"],
            output_names=["trilinear", "abs_fullWarp"],
            function=run_HCP_gradient_unwarp,
            imports=grad_unwarp_imports,
        ),
        name=f"grad_unwarp_{name}",
    )

    wf.connect(extract_ROI, "roi_file", grad_unwarp, "phase_vol")
    # wf.connect(',grad_unwarp,'input_coeffs')

    convertwarp_imports = ["import os", "import subprocess"]
    convert_warp = pe.Node(
        util.Function(
            input_names=["cw_trilinear", "cw_fullWarp_abs"],
            output_names=["out_file_cw", "out_jac_cw"],
            function=run_convertwarp,
            imports=convertwarp_imports,
        ),
        name=f"convert_warp_{name}",
    )

    wf.connect(grad_unwarp, "trilinear", convert_warp, "cw_trilinear")
    wf.connect(grad_unwarp, "abs_fullWarp", convert_warp, "cw_fullWarp_abs")

    fsl_maths = pe.Node(interface=fsl.maths.UnaryMaths(), name=f"fsl_maths_{name}")
    fsl_maths.inputs.args = "-Tmean"

    wf.connect(convert_warp, "out_jac_cw", fsl_maths, "in_file")

    apply_warp = pe.Node(interface=fsl.ApplyWarp(), name=f"apply_warp_{name}")
    apply_warp.inputs.in_file = inp_image
    apply_warp.inputs.interp = "spline"
    apply_warp.inputs.relwarp = True

    wf.connect(extract_ROI, "roi_file", apply_warp, "ref_file")
    wf.connect(convert_warp, "out_file_cw", apply_warp, "field_file")
    # applyjacobian ='True'
    # if applyjacobian == "True":

    # Apply Jacobian
    apply_jacobian = pe.Node(
        interface=fsl.maths.BinaryMaths(), name=f"apply_jacobian_{name}"
    )
    apply_jacobian.inputs.operation = "mul"

    wf.connect(apply_warp, "out_file", apply_jacobian, "in_file")
    wf.connect(fsl_maths, "out_file", apply_jacobian, "operand_file")

    # Make a dilated mask in the distortion corrected space
    # Apply_mask
    apply_mask = pe.Node(interface=fsl.maths.UnaryMaths(), name=f"apply_mask_{name}")
    apply_mask.inputs.in_file = inp_image
    apply_mask.inputs.operation = "abs"
    apply_mask.inputs.args = "-bin -dilD"

    # ApplyWarp
    apply_warp_mask = pe.Node(interface=fsl.ApplyWarp(), name=f"apply_warp_mask_{name}")
    apply_warp_mask.inputs.interp = "nn"
    apply_warp_mask.inputs.relwarp = True

    wf.connect(apply_mask, "out_file", apply_warp_mask, "in_file")
    wf.connect(apply_mask, "out_file", apply_warp_mask, "ref_file")
    wf.connect(convert_warp, "out_file_cw", apply_warp_mask, "field_file")

    out_warpmask = (apply_warp_mask, "out_file")
    out_applywarp = (apply_warp, "out_file")

    return (wf, out_warpmask, out_applywarp)


def phase_encode(unwarp_dir, phase_one, phase_two, dwell_time_one=None, 
                 dwell_time_two=None, ro_time_one=None, ro_time_two=None):
    """Calculate readout time and populate parameter file

    Parameters
    __________

    unwarp_dir
        phase encoding direction
    phase_one
        $WD/PhaseOne
    phase_two
        $WD/PhaseTwo
    dwell_time_one
        echo spacing of phase one
    dwell_time_two
        echo spacing of phase two
    ro_time_one
        total readout time of phase one
    ro_time_two
        total readout time of phase two

    Returns
    _______

    acq_params
        readout parameters

    """

    meta_data = [dwell_time_one, dwell_time_two,
                 ro_time_one, ro_time_two]
    if not any(meta_data):
        raise Exception("\n[!] Blip-FSL-TOPUP workflow: neither "
                        "TotalReadoutTime nor DwellTime is present in the "
                        "epi field map meta-data.")

    # create text file
    acq_params = os.path.join(os.getcwd(), "acqparams.txt")

    if isinstance(unwarp_dir, bytes):
        unwarp_dir = unwarp_dir.decode()

    if unwarp_dir in ["x", "x-", "-x", "i", "-i", "i-"]:
        if dwell_time_one and dwell_time_two:
            dim = nibabel.load(phase_one).shape[0]
            n_PE_steps = dim - 1
            ro_time_one = np.round(dwell_time_one * n_PE_steps, 6)
            ro_time_two = np.round(dwell_time_two * n_PE_steps, 6)
        if ro_time_one and ro_time_two:
            ro_times = [f"-1 0 0 {ro_time_one}", f"1 0 0 {ro_time_two}"]
        else:
            raise Exception("[!] No dwell time or total readout time "
                            "present for the acq-fMRI EPI field maps.")
    elif unwarp_dir in ["y", "y-", "-y", "j", "-j", "j-"]:
        if dwell_time_one and dwell_time_two:
            dim = nibabel.load(phase_one).shape[1]
            n_PE_steps = dim - 1
            ro_time_one = np.round(dwell_time_one * n_PE_steps, 6)
            ro_time_two = np.round(dwell_time_two * n_PE_steps, 6)
        if ro_time_one and ro_time_two:
            ro_times = [f"0 -1 0 {ro_time_one}", f"0 1 0 {ro_time_two}"]
        else:
            raise Exception("[!] No dwell time or total readout time "
                            "present for the acq-fMRI EPI field maps.")
    else:
        raise Exception(f"unwarp_dir={unwarp_dir} is unsupported.")

    # get number of volumes
    dims = [
        int(subprocess.check_output([f"fslval", phase_one, "dim4"]).decode(sys.stdout.encoding)),
        int(subprocess.check_output([f"fslval", phase_two, "dim4"]).decode(sys.stdout.encoding))
    ]

    # add read out times to text file
    for d, ro in zip(dims, ro_times):
        i = 1
        while i <= d:
            # write to text file
            with open(acq_params, "a") as fp:
                fp.write(f"{ro}\n")
            i += 1

    return acq_params


def z_pad(name="z_pad"):
    """Pad in Z by one slice if odd so that topup does not complain
    (slice consists of zeros that will be dilated by following step)
    """

    wf = pe.Workflow(name=name)

    inputspec = pe.Node(util.IdentityInterface(fields=['input_image']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['output_image']),
                         name='outputspec')



    fslroi = pe.Node(interface=fsl_utils.ExtractROI(),
                     name=f"fslroi_{name}")
    fslroi.inputs.args = "0 -1 0 -1 0 1 0 -1"

    wf.connect(inputspec, 'input_image', fslroi, 'in_file')

    fslmaths = pe.Node(interface=fsl_maths.BinaryMaths(),
                       name=f"fslmaths_{name}")
    fslmaths.inputs.operand_value = 0
    fslmaths.inputs.operation = "mul"

    wf.connect(fslroi, 'roi_file', fslmaths, 'in_file')
    wf.connect(fslroi, 'roi_file', fslmaths, 'operand_file')

    files_list = pe.Node(interface=utility.Merge(2),
                         name=f"create_list_{name}")

    wf.connect(inputspec, 'input_image', files_list, 'in1')
    wf.connect(fslmaths, 'out_file', files_list, 'in2')

    fslmerge = pe.Node(interface=fsl_utils.Merge(),
                       name=f"fslmerge_{name}")
    fslmerge.inputs.dimension = 'z'

    wf.connect(files_list, 'out', fslmerge, 'in_files')
    wf.connect(fslmerge, 'merged_file', outputspec, 'output_image')

    return wf


def choose_phase_image(phase_imgs, unwarp_dir):

    if isinstance(unwarp_dir, bytes):
        unwarp_dir = unwarp_dir.decode()

    if unwarp_dir in ["x","i","y","j"]:
        # select the first volume from PhaseTwo
        VolumeNumber = 1 + 1
        vnum = str(VolumeNumber).zfill(2)
        out_phase_image = phase_imgs[1]

    elif unwarp_dir in ["-x","-i","-y","-j","x-","i-","y-","j-"]:
        # select the first volume from PhaseOne
        VolumeNumber = 0 + 1
        vnum = str(VolumeNumber).zfill(2)
        out_phase_image = phase_imgs[0]

    return (out_phase_image, vnum)


def run_fsl_topup(merged_file, acqparams):

    out_basename = "topup_out"
    jac_basename = "topup_jac"
    xfm_basename = "topup_xfm"
    warp_basename = "topup_warpfield"

    out_fieldcoef = os.path.join(os.getcwd(), 
                                 f"{out_basename}_fieldcoef.nii.gz")
                                 
    out_movpar = os.path.join(os.getcwd(), 
                              f"{out_basename}_movpar.txt")

    corrected_outfile = os.path.join(os.getcwd(), 
                                     f"{out_basename}_corrected.nii.gz")
                                     
    field_out = os.path.join(os.getcwd(), 
                             f"{out_basename}_field.nii.gz")
                             
    out_jacs = [os.path.join(os.getcwd(), f"{jac_basename}_01.nii.gz"),
                os.path.join(os.getcwd(), f"{jac_basename}_02.nii.gz")]
                
    log_out = os.path.join(os.getcwd(), 
                           f"{out_basename}_log.log")
                           
    out_xfms = [os.path.join(os.getcwd(), f"{xfm_basename}_01.mat"),
                os.path.join(os.getcwd(), f"{xfm_basename}_02.mat")]

    out_warps = [os.path.join(os.getcwd(), f"{warp_basename}_01.nii.gz"),
                 os.path.join(os.getcwd(), f"{warp_basename}_02.nii.gz")]

    cmd = ["topup", f"--datain={acqparams}", f"--imain={merged_file}", 
           f"--out={out_basename}", f"--iout={corrected_outfile}", 
           f"--fout={field_out}", f"--jacout={jac_basename}",
           f"--logout={log_out}", f"--rbmout={xfm_basename}",
           f"--dfout={warp_basename}"]
    
    retcode = subprocess.check_output(cmd)
    
    return (out_fieldcoef, out_movpar, corrected_outfile, field_out, out_jacs,
            log_out, out_xfms, out_warps)


def find_vnum_base(vnum, motion_mat_list, jac_matrix_list, warp_field_list):

    for i in motion_mat_list:
        if f'topup_xfm_{vnum}.mat' in i:
            out_motion_mat = i
            break

    for i in jac_matrix_list:
        if f'topup_jac_{vnum}' in i:
            out_jacobian = i
            break

    for i in warp_field_list:
        if f'topup_warpfield_{vnum}' in i:
            out_warp_field = i
            break
    
    return(out_motion_mat, out_jacobian, out_warp_field)
    
