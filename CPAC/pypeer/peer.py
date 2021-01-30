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


def make_pypeer_dir(dirpath):
    try:
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
    except:
        raise Exception("\n\n[!] Could not create the output directory for "
                        "PyPEER. Double-check your permissions?\n\nAttempted "
                        "directory path:\n{0}\n".format(dirpath))


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
    _motion_threshold : float
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
                    output_dir, sub_id, pipeline_ids, stim_path, gsr=False, 
                    scrub=False, scrub_thresh=None):

    print("\n\n=== C-PAC now executing PyPEER for {0}. ===\n\n".format(sub_id))
    print("PEER scans for training model:\n{0}"
          "\n\n".format(str(peer_scan_names)))
    print("Data scans to estimate eye movements for:\n{0}"
          "\n\n".format(str(data_scan_names)))

    # note these are non-nuisance-regression strategy paths
    cpac_func_standard_paths = os.path.join(output_dir,
                                            "pipeline_*",
                                            sub_id,
                                            "functional_to_standard",
                                            "_scan_*", "*.nii*")
    func_standard_paths = glob.glob(cpac_func_standard_paths)

    if not func_standard_paths:
        raise Exception("\n\n[!] Could not find any 'functional_to_standard' "
                        "file paths in your output directory - did your "
                        "C-PAC run complete successfully?\n\n")

    eye_mask_glob = os.path.join(output_dir,
                                 "pipeline_*",
                                 sub_id,
                                 "template_eye_mask",
                                 "*")
    eye_mask_path = glob.glob(eye_mask_glob)[0]

    if not os.path.isfile(eye_mask_path):
        raise Exception("\n\n[!] Could not find the template eye mask "
                        "file path in your output directory - did your "
                        "C-PAC run complete successfully?\n\n")

    print("Found input files:\n{0}\n".format(func_standard_paths))

    pypeer_outdir = func_standard_paths[0].split("functional_to_standard")[0]
    pypeer_outdir = os.path.join(pypeer_outdir, "PyPEER")

    make_pypeer_dir(pypeer_outdir)

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
            peer_scans[func_path] = [data, scan_label, removed_indices]
        elif scan_label in data_scan_names:
            raveled_data = pypeer_ravel_data(data)
            data_scans[func_path] = [raveled_data, scan_label]

    for peer_scan_path in peer_scans.keys():
        print("Training the eye estimation model using:\n{0}"
              "\n\n".format(peer_scan_path))
        data = peer_scans[peer_scan_path][0]
        peername = peer_scans[peer_scan_path][1]
        removed_indices = peer_scans[peer_scan_path][2]
        data_for_training, calibration_points_removed = prepare_data_for_svr(data,
                                                                             removed_indices,
                                                                             eye_mask_path)
        xmodel, ymodel = train_model(data_for_training,
                                     calibration_points_removed,
                                     stim_path)

        model_dir = os.path.join(pypeer_outdir, "peer_model-{0}".format(peername))
        make_pypeer_dir(model_dir)

        save_model(xmodel, ymodel, os.path.basename(peer_scan_path), str(scrub), 
                   str(gsr), model_dir)

        for data_scan_path in data_scans.keys():
            print("Estimating eye movements for:\n{0}"
                  "\n\n".format(data_scan_path))
            data = data_scans[data_scan_path][0]
            name = data_scans[data_scan_path][1]

            xmodel, ymodel, xname, yname = load_model(model_dir)
            xfix, yfix = predict_fixations(xmodel, ymodel, data)

            estimate_dir = os.path.join(pypeer_outdir, 
                                        "estimations-{0}_model-{1}".format(name, peername))
            make_pypeer_dir(estimate_dir)

            fix_xname, fix_yname = save_fixations(xfix, yfix, xname, yname,
                                                  estimate_dir)
            estimate_em(xfix, yfix, fix_xname, fix_yname, estimate_dir)
