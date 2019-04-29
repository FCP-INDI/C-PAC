import scipy.io
import numpy as np
import nibabel as nib
import os
from scipy.sparse import lil_matrix
from scipy.ndimage.filters import gaussian_filter
from numpy import ndarray
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
from nilearn.masking import apply_mask
from nilearn.masking import compute_epi_mask
import numpy.matlib

from CPAC.isc.utils import correlation


def qpp_wf(img, num_scans, window_length,
           permutations, correlation_threshold, 
           max_iterations, convergence_iterations,
           out_dir, x, y):
    """
    This code is adapted from the paper "Quasi-periodic patterns (QP): Large-
    scale dynamics in resting state fMRI that correlate with local infraslow
    electrical activity", Shella Keilholz et al. NeuroImage, 2014.

    Input:
    ------
    img: 2D nifti image
    mask: mask of the 2D nifti image
    num_scans: number of subjects*number of runs per subject
    window_length: window length
    permutations: number of repetitions
    correlation_threshold: threshold
    n_itr_th: number of iterations
    mx_itr: maximum number of repetitions
    pfs: path to save the template, FTP, ITP and iter files


    Returns:
    -------
    time_course_file: 2D array of time points where QPP is detected in .npy format

    """
    voxels = np.prod(img.shape[0:3])
    trs = img.shape[3]

    trs_per_scan = int(trs / num_scans) 
    inpectable_windows = trs_per_scan - window_length + 1
    df = voxels * window_length

    convergence_iterations = int(max(1, convergence_iterations))

    if callable(correlation_threshold):
        correlation_thresholds = [correlation_threshold(i) for i in range(max_iterations)]
    else:
        correlation_thresholds = [correlation_threshold for _ in range(max_iterations)]

    def flattened_segment(img, window_length, pos):
        return img[:,:,:, pos:pos + window_length].flatten()

    def normalized_flattened_segment(img, window_length, pos, df):
        segment = flattened_segment(img, window_length, pos)
        segment -= np.sum(segment) / df
        segment = segment / np.sqrt(np.dot(segment, segment))
        return segment

    random_selections = np.zeros((num_scans, window_length - 1))

    for i in range(1, num_scans + 1):
        random_selections[i - 1, :] = range(i * trs_per_scan - window_length + 2, i * trs_per_scan + 1)

    initial_trs = np.arange(1, trs + 1)
    random_selections = random_selections.flatten()
    initial_trs = np.delete(initial_trs, random_selections - 1, 0)

    initial_trs = np.random.permutation(initial_trs)
    initial_trs = initial_trs[0:permutations]

    permutation_result = [{} for _ in range(permutations)]
    for perm in range(permutations):

        template_holder = np.zeros(trs)
        for scan in range(num_scans):
            for sliding_window in range(inpectable_windows):
                random_initial_window = normalized_flattened_segment(img, window_length, initial_trs[perm], df)
                scan_window = normalized_flattened_segment(img, window_length, (scan * trs_per_scan) + sliding_window, df)
                template_holder[scan * trs_per_scan + sliding_window] = np.dot(random_initial_window, scan_window)

        peaks = detect_peaks(template_holder, mph=correlation_thresholds[0], mpd=window_length)

        for i in range(num_scans):
            initial_tr = i * trs_per_scan
            if initial_tr in peaks:
                peaks = np.delete(peaks, np.where(peaks == initial_tr))

            final_tr = i * trs_per_scan + inpectable_windows
            if final_tr in peaks:
                peaks = np.delete(peaks, np.where(peaks == final_tr))

        template_holder_convergence = np.array([
            template_holder,
            template_holder,
            template_holder,
        ])

        for iteration in range(1, max_iterations):

            # Smoothing
            template_holder = gaussian_filter(template_holder, 0.5)

            peak_threshold = correlation_thresholds[iteration]
            found_peaks = np.size(peaks)
            
            if found_peaks < 1:
                break

            peaks_segments = flattened_segment(img, window_length, peaks[0])
            for peak in peaks[1:]:
                peaks_segments = peaks_segments + flattened_segment(img, window_length, peak)
            peaks_segments = peaks_segments / found_peaks
            peaks_segments = peaks_segments - np.sum(peaks_segments) / df
            peaks_segments = peaks_segments / np.sqrt(np.dot(peaks_segments, peaks_segments))

            for scan in range(num_scans):
                for window in range(inpectable_windows):
                    final_tr = scan * trs_per_scan + window
                    template_holder[final_tr] = np.dot(
                        peaks_segments,
                        flattened_segment(img, window_length, final_tr)
                    )

            peaks = detect_peaks(template_holder, mph=peak_threshold, mpd=window_length)

            for scan in range(num_scans):
                initial_tr = scan * trs_per_scan
                if initial_tr in peaks:
                    peaks = np.delete(peaks, np.where(peaks == initial_tr))

                final_tr = scan * trs_per_scan + inpectable_windows
                if final_tr in peaks:
                    peaks = np.delete(peaks, np.where(peaks == final_tr))

            if np.all(correlation(peaks_segments, template_holder_convergence) > 0.9999):
                break

            template_holder_convergence[0:2, :] = template_holder_convergence[1:, :]
            template_holder_convergence[2, :] = template_holder

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

    metrics_best_template = [
        np.median(best_template[best_selected_peaks]),
        np.median(np.diff(best_selected_peaks)),
        len(best_selected_peaks),
    ]

    #% building T1, by averaging segments of B with length 2*WL starting at FTP-WL/2; 
    # extra WL/2 at each end is primarily to have die-off effect; it
    #is also used in fine-phase-matching two QPPs when comparing them
    window_length_start = round(window_length / 2)
    window_length_end = window_length_start - window_length % 2

    best_template = []

    for best_peak in best_selected_peaks:
        start_tr = best_peak - window_length_start
        end_tr = best_peak + window_length - 1 + window_length_end

        start_segment = np.zeros((voxels, 0))
        if start_tr <= 0:
            start_segment = np.zeros((voxels, abs(start_tr)))
            start_tr = 0

        end_segment, = np.zeros((voxels, 0))
        if end_tr > trs:
            end_segment, = np.zeros((voxels, end_tr - trs))
            end_tr = trs

        img_segment = np.array(img[:, :, :, start_tr:end_tr])

        best_template += [
            np.concatenate([
                start_segment,
                img_segment,
                end_segment,
            ], axis=1)
        ]

    best_template = np.concatenate(best_template)
    best_template = best_template / len(best_selected_peaks)

    return best_template, best_selected_peaks, metrics_best_template
