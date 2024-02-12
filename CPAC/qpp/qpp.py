import numpy as np
import nibabel as nib
import os
from numpy import ndarray
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from nilearn.masking import apply_mask
from nilearn.masking import compute_epi_mask
import numpy.matlib

from CPAC.utils import check_random_state, correlation


def smooth(x):
    """
    Temporary moving average
    """
    return np.array(
        [x[0]] +
        [np.mean(x[0:3])] +
        (np.convolve(x, np.ones(5), 'valid') / 5).tolist() +
        [np.mean(x[-3:])] +
        [x[-1]]
    )


def flattened_segment(data, window_length, pos):
    return data[:, pos:pos + window_length].flatten(order='F')


def normalize_segment(segment, df):
    segment -= np.sum(segment) / df
    segment = segment / np.sqrt(np.dot(segment, segment))
    return segment


def detect_qpp(data, num_scans, window_length,
               permutations, correlation_threshold, 
               iterations, convergence_iterations=1,
               random_state=None):
    """
    This code is adapted from the paper "Quasi-periodic patterns (QP): Large-
    scale dynamics in resting state fMRI that correlate with local infraslow
    electrical activity", Shella Keilholz et al. NeuroImage, 2014.
    """

    random_state = check_random_state(random_state)

    voxels, trs = data.shape

    iterations = int(max(1, iterations))
    convergence_iterations = int(max(1, convergence_iterations))

    if callable(correlation_threshold):
        correlation_thresholds = [correlation_threshold(i) for i in range(iterations)]
    else:
        correlation_thresholds = [correlation_threshold for _ in range(iterations)]

    trs_per_scan = int(trs / num_scans)
    inpectable_trs = np.arange(trs) % trs_per_scan
    inpectable_trs = np.where(inpectable_trs < trs_per_scan - window_length + 1)[0]

    df = voxels * window_length

    initial_trs = random_state.choice(inpectable_trs, permutations)

    permutation_result = [{} for _ in range(permutations)]
    for perm in range(permutations):

        template_holder = np.zeros(trs)
        random_initial_window = normalize_segment(flattened_segment(data, window_length, initial_trs[perm]), df)
        for tr in inpectable_trs:
            scan_window = normalize_segment(flattened_segment(data, window_length, tr), df)
            template_holder[tr] = np.dot(random_initial_window, scan_window)

        template_holder_convergence = np.zeros((convergence_iterations, trs))

        for iteration in range(iterations):

            peak_threshold = correlation_thresholds[iteration]

            peaks, _ = find_peaks(template_holder, height=peak_threshold, distance=window_length)
            peaks = np.delete(peaks, np.where(~np.isin(peaks, inpectable_trs))[0])

            template_holder = smooth(template_holder)

            found_peaks = np.size(peaks)
            if found_peaks < 1:
                break

            peaks_segments = flattened_segment(data, window_length, peaks[0])
            for peak in peaks[1:]:
                peaks_segments = peaks_segments + flattened_segment(data, window_length, peak)

            peaks_segments = peaks_segments / found_peaks
            peaks_segments = normalize_segment(peaks_segments, df)

            for tr in inpectable_trs:
                scan_window = normalize_segment(flattened_segment(data, window_length, tr), df)
                template_holder[tr] = np.dot(peaks_segments, scan_window)

            if np.all(correlation(template_holder, template_holder_convergence) > 0.9999):
                break

            if convergence_iterations > 1:
                template_holder_convergence[1:] = template_holder_convergence[0:-1]
            template_holder_convergence[0] = template_holder

        if found_peaks > 1:
            permutation_result[perm] = {
                'template': template_holder,
                'peaks': peaks,
                'final_iteration': iteration,
                'correlation_score': np.sum(template_holder[peaks]),
            }

    # Retrieve max correlation of template from permutations
    correlation_scores = np.array([
        r['correlation_score'] if r else 0.0 for r in permutation_result
    ])
    if not np.any(correlation_scores):
        raise Exception("C-PAC could not find QPP in your data. "
                        "Please lower your correlation threshold and try again.")

    max_correlation = np.argsort(correlation_scores)[-1]
    best_template = permutation_result[max_correlation]['template']
    best_selected_peaks = permutation_result[max_correlation]['peaks']

    best_template_metrics = [
        np.median(best_template[best_selected_peaks]),
        np.median(np.diff(best_selected_peaks)),
        len(best_selected_peaks),
    ]

    window_length_start = round(window_length / 2)
    window_length_end = window_length_start - window_length % 2

    best_template_segment = np.zeros((voxels, window_length))

    for best_peak in best_selected_peaks:
        start_tr = int(best_peak - np.ceil(window_length / 2.))
        end_tr = int(best_peak + np.floor(window_length / 2.))

        start_segment = np.zeros((voxels, 0))
        if start_tr <= 0:
            start_segment = np.zeros((voxels, abs(start_tr)))
            start_tr = 0

        end_segment = np.zeros((voxels, 0))
        if end_tr > trs:
            end_segment = np.zeros((voxels, end_tr - trs))
            end_tr = trs

        data_segment = data[:, start_tr:end_tr]

        best_template_segment += np.concatenate([
            start_segment,
            data_segment,
            end_segment,
        ], axis=1)

    best_template_segment /= len(best_selected_peaks)

    return best_template_segment, best_selected_peaks, best_template_metrics
