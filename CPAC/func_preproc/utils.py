
import numpy as np
from scipy.signal import iirnotch, firwin, filtfilt, freqz
from matplotlib import pyplot as plt
import nibabel as nb
import subprocess
import math


def nullify(value, function=None):
    from traits.trait_base import Undefined
    if value is None:
        return Undefined
    if function:
        return function(value)
    return value


def chunk_ts(func_file, n_chunks=None, chunk_size=None):
    func_img = nb.load(func_file)
    trs = func_img.shape[3]
    TR_ranges = []

    if n_chunks:
        chunk_size = trs/n_chunks
    elif chunk_size:
        n_chunks = int(trs/chunk_size)
    else:
        raise Exception("\n[!] Dev error: Either 'n_chunks' or 'chunk_size' "
                        "arguments must be passed to 'chunk_ts' function.\n")

    for chunk_idx in range(0, n_chunks):
        if chunk_idx == n_chunks - 1:
            TR_ranges.append((int(chunk_idx*chunk_size), int(trs - 1)))
        else:
            TR_ranges.append((int(chunk_idx*chunk_size), int((chunk_idx+1)*chunk_size - 1)))
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


def degrees_to_mm(degrees, head_radius):
    # function to convert degrees of motion to mm
    mm = 2*math.pi*head_radius*(degrees/360)
    return mm


def mm_to_degrees(mm, head_radius):
    # function to convert mm of motion to degrees
    degrees = 360*mm/(2*math.pi*head_radius)
    return degrees

def degrees_to_mm(degrees, head_radius):
    # function to convert degrees of motion to mm
    mm = 2*math.pi*head_radius*(degrees/360)
    return mm


def mm_to_degrees(mm, head_radius):
    # function to convert mm of motion to degrees
    degrees = 360*mm/(2*math.pi*head_radius)
    return degrees

def degrees_to_mm(degrees, head_radius):
    # function to convert degrees of motion to mm
    mm = 2*math.pi*head_radius*(degrees/360)
    return mm


def mm_to_degrees(mm, head_radius):
    # function to convert mm of motion to degrees
    degrees = 360*mm/(2*math.pi*head_radius)
    return degrees


def notch_filter_motion(motion_params, filter_type, TR, fc_RR_min=None,
                        fc_RR_max=None, center_freq=None, freq_bw=None,
                        lowpass_cutoff=None, filter_order=4):
    # Adapted from DCAN Labs:
    #   https://github.com/DCAN-Labs/dcan_bold_processing/blob/master/
    #       ...matlab_code/filtered_movement_regressors.m

    if "ms" in TR:
        TR = float(TR.replace("ms", ""))/1000
    elif "ms" not in TR and "s" in TR:
        TR = float(TR.replace("s", ""))

    params_data = np.loadtxt(motion_params)

    # Sampling frequency
    fs = 1 / TR

    # Nyquist frequency
    fNy = fs / 2

    if filter_type == "notch":

        # Respiratory Rate
        if fc_RR_min and fc_RR_max:
            rr = [float(fc_RR_min) / float(60),
                  float(fc_RR_max) / float(60)]

            rr_fNy = [rr[0] + fNy, rr[1] + fNy]
            fa = abs(rr - np.floor(np.divide(rr_fNy, fs)) * fs)

        elif center_freq and freq_bw:
            tail = float(freq_bw)/float(2)
            fa = [center_freq-tail, center_freq+tail]

        W_notch = np.divide(fa, fNy)

        Wn = np.mean(W_notch)
        bw = np.diff(W_notch)

        # for filter info
        center_freq = Wn * fNy
        bandwidth = fa[1] - fa[0]

        Q = Wn/bw
        [b_filt, a_filt] = iirnotch(Wn, Q)
        num_f_apply = np.floor(filter_order / 2)

        filter_info = f"Motion estimate filter information\n\nType: Notch\n" \
                      f"\nCenter freq: {center_freq}\nBandwidth: {bandwidth}\n\n" \
                      f"Wn: {Wn}\nQ: {Q}\n\n" \
                      f"Based on:\nSampling freq: {fs}\nNyquist freq: {fNy}"

    elif filter_type == "lowpass":

        if fc_RR_min:
            rr = float(fc_RR_min) / float(60)
            rr_fNy = rr + fNy
            fa = abs(rr - np.floor(np.divide(rr_fNy, fs)) * fs)

        elif lowpass_cutoff:
            fa = lowpass_cutoff

        Wn = fa/fNy

        if filter_order:
            b_filt = firwin(filter_order+1, Wn)
            a_filt = 1

        num_f_apply = 0

        filter_info = f"Motion estimate filter information\n\nType: Lowpass" \
                      f"\n\nCutoff freq: {fa}\nWn: {Wn}\n\n" \
                      f"Based on:\nSampling freq: {fs}\nNyquist freq: {fNy}"

    filter_design = os.path.join(os.getcwd(),
                                 "motion_estimate_filter_design.txt")
    filter_plot = os.path.join(os.getcwd(),
                               "motion_estimate_filter_freq-response.png")

    # plot frequency response for user info
    w, h = freqz(b_filt, a_filt, fs=fs)

    fig, ax1 = plt.subplots()
    ax1.set_title('Motion estimate filter frequency response')

    ax1.plot(w, 20 * np.log10(abs(h)), 'b')
    ax1.set_ylabel('Amplitude [dB]', color='b')
    ax1.set_xlabel('Frequency [Hz]')

    plt.savefig(filter_plot)

    with open(filter_design, 'wt') as f:
        f.write(filter_info)

    # convert rotation params from degrees to mm
    params_data[:, 0:3] = degrees_to_mm(params_data[:, 0:3], head_radius=50)

    filtered_params = filtfilt(b_filt, a_filt, params_data.T)

    # back rotation params to degrees
    filtered_params[0:3,:] = mm_to_degrees(filtered_params[0:3,:], head_radius = 50)

    # back rotation params to degrees
    filtered_params[0:3,:] = mm_to_degrees(filtered_params[0:3,:], head_radius = 50)

    filtered_motion_params = os.path.join(os.getcwd(),
                                          "{0}_filtered.1D".format(os.path.basename(motion_params)))
    np.savetxt(filtered_motion_params, filtered_params.T, fmt='%f')

    return (filtered_motion_params, filter_design, filter_plot)
