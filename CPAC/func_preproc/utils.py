
import numpy as np
from scipy.signal import iirnotch, filtfilt
import nibabel as nb
import subprocess

def add_afni_prefix(tpattern):
    if tpattern:
        if ".txt" in tpattern:
            tpattern = "@{0}".format(tpattern)
    return tpattern


def nullify(value, function=None):
    from traits.trait_base import Undefined
    if value is None:
        return Undefined
    if function:
        return function(value)
    return value


def chunk_ts(func_file, n_cpus):
    func_img = nb.load(func_file)
    trs = func_img.shape[3]
    chunk = trs/n_cpus
    TR_ranges = []

    for chunk_idx in range(0, n_cpus):
        if chunk_idx == n_cpus - 1:
            TR_ranges.append((int(chunk_idx*chunk), int(trs - 1)))
        else:
            TR_ranges.append((int(chunk_idx*chunk), int((chunk_idx+1)*chunk -1)))
    return TR_ranges


def split_ts_chunks(func_file, tr_ranges):
    if '.nii' in func_file:
        ext = '.nii'
    if '.nii.gz' in func_file:
        ext = '.nii.gz'

    split_funcs = []
    for chunk_idx, tr_range in enumerate(tr_ranges):
        out_file = os.path.join(os.getcwd(), os.path.basename(func_file).replace(ext, "_{0}{1}".format(chunk_idx, ext)))
        in_file = "{0}[{1}..{2}]".format(func_file, tr_range[0], tr_range[1])

        cmd = ["3dcalc", "-a", in_file, "-expr", "a", "-prefix", out_file]

        retcode = subprocess.check_output(cmd)

        split_funcs.append(out_file)

    return split_funcs


def oned_text_concat(in_files):
    out_file = os.path.join(os.getcwd(), os.path.basename(in_files[0].replace("_0", "")))

    out_txt = []
    for txt in in_files:
        with open(txt, 'r') as f:
            txt_lines = f.readlines()
        if not out_txt:
            out_txt = [x for x in txt_lines]
        else:
            for line in txt_lines:
                if "#" in line:
                    continue
                out_txt.append(line)

    with open(out_file, 'wt') as f:
        for line in out_txt:
            f.write(line)

    return out_file


def notch_filter_motion(motion_params, fc_RR_min, fc_RR_max, TR, 
                        filter_order=4):

    # Adapted from DCAN Labs:
    #   https://github.com/DCAN-Labs/dcan_bold_processing/blob/master/
    #       ...matlab_code/filtered_movement_regressors.m

    TR = float(TR.replace("s", ""))
    params_data = np.loadtxt(motion_params)

    fc_RR_bw = [fc_RR_min, fc_RR_max]

    # Respiratory Rate
    rr = [float(fc_RR_min) / float(60),
          float(fc_RR_max) / float(60)]

    # Sampling frequency
    fs = 1 / TR

    # Nyquist frequency
    fNy = fs / 2;

    rr_fNy = [rr[0] + fNy, rr[1] + fNy]

    fa = abs(rr - np.floor(np.divide(rr_fNy, fs)) * fs)

    W_notch = np.divide(fa, fNy)
    Wn = np.mean(W_notch)
    bw = np.diff(W_notch)
    [b_filt, a_filt] = iirnotch(Wn, bw)
    num_f_apply = np.floor(filter_order / 2)

    filtered_params = filtfilt(b_filt, a_filt, params_data.T)

    filtered_motion_params = os.path.join(os.getcwd(),
                                          "{0}_notch-filtered.1D".format(os.path.basename(motion_params)))
    np.savetxt(filtered_motion_params, filtered_params.T, fmt='%f')

    return filtered_motion_params


