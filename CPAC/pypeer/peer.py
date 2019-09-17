import os
import csv
import glob
import nibabel as nb
import numpy as np

# check if they have PyPEER installed
try:
    import PyPEER
    from PyPEER.peer_func import \
        global_signal_regression, \
        prepare_data_for_svr, \
        train_model, \
        save_model, \
        load_model, \
        predict_fixations, \
        save_fixations, \
        estimate_em, \
        load_data
except ImportError:
    raise ImportError("\n\n[!] PyPEER is not installed. Please double-"
                      "check your Python environment and ensure that the "
                      "PyPEER package is available.")


def pypeer_eye_masking(data_path, eye_mask_path):
    eye_mask = nb.load(eye_mask_path).get_data()

    data = load_data(data_path)

    for vol in range(data.shape[3]):
        output = np.multiply(eye_mask, data[:, :, :, vol])
        data[:, :, :, vol] = output

    return data


def pypeer_zscore(data):
    volumes = data.shape[3]

    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            for z in range(data.shape[2]):
                vmean = np.mean(np.array(data[x, y, z, :]))
                vstdv = np.std(np.array(data[x, y, z, :]))

                for time in range(volumes):
                    if vstdv != 0:
                        data[x, y, z, time] = (float(
                            data[x, y, z, time]) - float(vmean)) / vstdv
                    else:
                        data[x, y, z, time] = float(
                            data[x, y, z, time]) - float(vmean)

    return data


def pypeer_ravel_data(data):
    raveled_data = [data[:, :, :, vol].ravel() for vol in
                    np.arange(data.shape[3])]
    return raveled_data


def motion_scrub(_ms_filename, _motion_threshold):
    """
    Determines volumes with high motion artifact
    Parameters
    ----------
    _ms_filename : string
        Pathname of the CSV file containing the framewise displacement per time point for a given fMRI scan
    _motion_threshold  : float
        Threshold for high motion (framewise displacement, defined by Power et al. 2012)
    Returns
    -------
    _removed_indices : int
        List of volumes to remove for motion scrubbing

    Taken from:
        https://github.com/ChildMindInstitute/PyPEER
        Jake Son

    Adapted to accept a direct file path to the mean FD 1D file.
    """

    with open(_ms_filename, 'r') as f:
        reader = csv.reader(f)
        censor_pre = [x[0] for x in reader]

    nuissance_vector = [float(x) for x in censor_pre]

    _removed_indices = [i for i, x in enumerate(nuissance_vector) if x >= float(_motion_threshold)]

    return _removed_indices


def prep_for_pypeer(peer_scan_names, data_scan_names, eye_mask_path,
                    output_dir, sub_id, stim_path, gsr=False, scrub=False,
                    scrub_thresh=None):

    print("\n\n=== C-PAC now executing PyPEER. ===\n\n")
    print("PEER scans for training model:\n{0}"
          "\n\n".format(str(peer_scan_names)))
    print("Data scans to estimate eye movements for:\n{0}"
          "\n\n".format(str(data_scan_names)))

    # note these are non-nuisance-regression strategy paths
    func_standard_paths = glob.glob(os.path.join(output_dir,
                                                 "pipeline_*",
                                                 sub_id,
                                                 "functional_to_standard",
                                                 "_scan_*", "*.nii*"))

    if not func_standard_paths:
        raise Exception("\n\n[!] Could not find any 'functional_to_standard' "
                        "file paths in your output directory - did your "
                        "C-PAC run complete successfully?\n\n")

    pypeer_outdir = func_standard_paths[0].split("functional_to_standard")[0]

    peer_scans = {}
    data_scans = {}

    for func_path in func_standard_paths:
        scan_label = func_path.split("/")[-2].replace("_scan_", "")

        if scan_label in peer_scan_names or scan_label in data_scan_names:
            print("Eye-masking and z-score standardizing "
                  "{0}..".format(scan_label))
            masked_data = pypeer_eye_masking(func_path, eye_mask_path)
            data = pypeer_zscore(masked_data)

        if gsr:
            print("Global signal regression for {0}..".format(scan_label))
            data = global_signal_regression(data, eye_mask_path)

        removed_indices = None
        if scrub and scan_label in peer_scan_names:
            print("Motion scrubbing (Power 2012) for "
                  "{0}..".format(scan_label))
            fd_path = func_path.replace("functional_to_standard",
                                        "frame_wise_displacement_power")
            fd_path = fd_path.replace(fd_path.split("/")[-1], "FD.1D")

            if not os.path.isfile(fd_path):
                raise Exception("\n\n[!] Could not find the mean framewise "
                                "displacement 1D file in your C-PAC output "
                                "directory.")

            removed_indices = motion_scrub(fd_path, scrub_thresh)

        if scan_label in peer_scan_names:
            peer_scans[func_path] = [data, removed_indices]
        elif scan_label in data_scan_names:
            raveled_data = pypeer_ravel_data(data)
            data_scans[func_path] = raveled_data

    for peer_scan_path in peer_scans.keys():
        print("Training the eye estimation model using:\n{0}"
              "\n\n".format(peer_scan_path))
        data = peer_scans[peer_scan_path][0]
        removed_indices = peer_scans[peer_scan_path][1]
        data_for_training, calibration_points_removed = prepare_data_for_svr(data,
                                                                             removed_indices,
                                                                             eye_mask_path)
        xmodel, ymodel = train_model(data_for_training,
                                     calibration_points_removed,
                                     stim_path)

        save_model(xmodel, ymodel, peer_scan_path, scrub, gsr, pypeer_outdir)

    for data_scan_path in data_scans.keys():
        print("Estimating eye movements for:\n{0}"
              "\n\n".format(data_scan_path))
        data = data_scans[data_scan_path]

        xmodel, ymodel, xname, yname = load_model(pypeer_outdir)
        xfix, yfix = predict_fixations(xmodel, ymodel, data)
        fix_xname, fix_yname = save_fixations(xfix, yfix, xname, yname,
                                              pypeer_outdir)
        estimate_em(xfix, yfix, fix_xname, fix_yname, pypeer_outdir)
