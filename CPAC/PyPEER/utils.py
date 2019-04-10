#!/usr/bin/env python
"""
Functions to predict eye movements

Authors:
    - Jake Son, 2017-2018  (jake.son@childmind.org)

"""

import os
import sys
import csv
import json
import numpy as np
import pandas as pd
import nibabel as nib
from sklearn.svm import SVR
from sklearn.externals import joblib


def scaffolding():
    """
    Creates the project folder and file hierarchy and returns pathnames

    Returns
    -------
    _project_dir : string
        Pathname of the highest-level project directory
    _data_dir : string
        Pathname of the directory containing data
    _output_dir : string
        Pathname of the output directory
    _stimulus_path : string
        Pathname of the PEER calibration scan stimuli

    """

    _project_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

    _data_dir = os.path.abspath(os.path.join(_project_dir, 'data'))

    _stimulus_path = os.path.abspath(os.path.join(_project_dir, 'peer/stim_vals.csv'))

    if not os.path.exists(_data_dir):

        os.makedirs(_data_dir)

        print('\nThe "data" directory needs to contain at least one subdirectory that contains data before proceeding. Each '
              'subdirectory should contain at least one PEER file '
              'and one file used to estimate eye movements, such as an fMRI scan acquired during movie viewing or rest.\n')

        sys.exit()

    if not os.listdir(_data_dir):

        print('\nThe "data" directory is empty. The "data" directory needs to contain at least one subdirectory that '
              'contains data. Each subdirectory should contain at least two PEER files.\n')

        sys.exit()

    else:

        for dataset in [x for x in os.listdir(_data_dir) if not x.startswith('.')]:

            dataset_path = os.path.abspath(os.path.join(_data_dir, dataset))
            print(dataset_path)

            if 'outputs' not in os.listdir(dataset_path):

                _dataset_dir = os.path.abspath(os.path.join(_data_dir, dataset))

                _output_dir = os.path.abspath(os.path.join(_dataset_dir, 'outputs'))

                os.makedirs(_output_dir)

            else:

                _dataset_dir = os.path.abspath(os.path.join(_data_dir, dataset))

                _output_dir = os.path.abspath(os.path.join(_dataset_dir, 'outputs'))

    return _project_dir, _data_dir, _stimulus_path


def set_parameters(_configs, new=False):
    """
    Sets configuration parameters

    Parameters
    ----------
    _configs :
        Dictionary containing configuration options from the config file (config.json)
    new : bool
        Do you want to start from a new file?

    Returns
    -------
    _configs :
        Updated dictionary containing configuration options from the config file (config.json)

    """

    if new:
        _configs = {x: "NA" for x in _configs}

    print('*Do not include single or double quotes*\n')

    if _configs['eye_mask_path'] == 'NA':
        _eye_mask_path = input('Add the full eye mask filepath: ')
        _configs['eye_mask_path'] = _eye_mask_path

    if _configs['train_file'] == 'NA':
        _train_file = input('Add the name of the file used for training [peer1.nii.gz]: ')
        if not _train_file:
            _configs['train_file'] = 'peer1.nii.gz'
        else:
            _configs['train_file'] = _train_file

    if _configs['test_file'] == 'NA':
        _test_file = input('Which file would you like to predict eye movements from? [movie.nii.gz]: ')
        if not _test_file:
            _configs['test_file'] = 'movie.nii.gz'
        else:
            _configs['test_file'] = _test_file

    if _configs['use_gsr'] == 'NA':

        _use_gsr = input('Use global signal regression? (y/n) [n]: ')

        if (not _use_gsr) or (_use_gsr == 'n'):
            _configs['use_gsr'] = "0"
        else:
            _configs['use_gsr'] = "1"

    if _configs['motion_scrub'] == 'NA':
        _use_ms = input('Use motion scrubbing? (y/n) [n]: ')

        if (not _use_ms) or (_use_ms == 'n'):
            _configs['use_ms'] = "0"
            _configs['motion_threshold'] = "0"
            _configs['motion_scrub'] = "Not implemented"
        elif _use_ms == 'y':
            _configs['use_ms'] = "1"
            _motion_scrub_filename = input('Add the filename of the CSV that contains the framewise displacement \
                                            time series [motion_ts.csv]: ')
            if not _motion_scrub_filename:
                _configs['motion_scrub'] = 'motion_ts.csv'
            else:
                _configs['motion_scrub'] = _motion_scrub_filename
            _motion_threshold = input('Add a motion threshold for motion scrubbing [.2]: ')
            if not _motion_threshold:
                _configs['motion_threshold'] = ".2"
            else:
                _configs['motion_threshold'] = _motion_threshold

    with open('config.json', 'w') as f:
        json.dump(_configs, f)

    return _configs


def load_config():
    """
    Loads configuration parameters as a dictionary from config.json

    Returns
    -------
    _configs :
        Dictionary containing configuration options from the config file (config.json)

    """

    with open('config.json', 'r') as f:
        _configs = json.load(f)

    if any('NA' in x for x in list(_configs.values())):

        print('\nOne or more of the necessary parameters are missing. You can either:\n'
              '     1) add the content using the command line, or\n'
              '     2) exit and manually edit the config.json file\n')

        _configs = set_parameters(_configs, new=False)

    return _configs


def load_data(_filepath):
    """
    Loads fMRI data

    Parameters
    ----------
    _filepath : string
        Pathname of the NIfTI file used to train a model or predict eye movements

    Returns
    -------
    _data : float
        4D numpy array containing fMRI data

    """

    nib_format = nib.load(_filepath)
    _data = nib_format.get_data()

    print('Training data Loaded')

    return _data


def global_signal_regression(_data, _eye_mask_path):
    """
    Performs global signal regression

    Parameters
    ----------
    _data : float
        Data from an fMRI scan as a 4D numpy array
    _eye_mask_path :
        Pathname for the eye mask NIfTI file (the standard MNI152 2mm FSL template is used for the linked preprint)
    Returns
    -------
    _data :
        4D numpy array containing fMRI data after global signal regression

    """

    eye_mask = nib.load(_eye_mask_path).get_data()

    global_mask = np.array(eye_mask, dtype=bool)

    regressor_map = {'constant': np.ones((_data.shape[3], 1))}
    regressor_map['global'] = _data[global_mask].mean(0)

    X = np.zeros((_data.shape[3], 1))
    csv_filename = ''

    for rname, rval in regressor_map.items():
        X = np.hstack((X, rval.reshape(rval.shape[0], -1)))
        csv_filename += '_' + rname

    X = X[:, 1:]

    Y = _data[global_mask].T
    B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)

    Y_res = Y - X.dot(B)

    _data[global_mask] = Y_res.T

    print('GSR completed.')

    return _data


def motion_scrub(_ms_filename, _data_dir, _motion_threshold):
    """
    Determines volumes with high motion artifact

    Parameters
    ----------
    _ms_filename : string
        Pathname of the CSV file containing the framewise displacement per time point for a given fMRI scan
    _data_dir : string
        Pathname of the directory containing data
    _motion_threshold  : float
        Threshold for high motion (framewise displacement, defined by Power et al. 2012)

    Returns
    -------
    _removed_indices : int
        List of volumes to remove for motion scrubbing

    """
    ##nipype this!


    file_path = os.path.abspath(os.path.join(_data_dir, _ms_filename))

    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        censor_pre = [x for x in reader]

    nuissance_vector = [float(x) for x in censor_pre[0]]

    _removed_indices = [i for i, x in enumerate(nuissance_vector) if x >= float(_motion_threshold)]

    return _removed_indices


def prepare_data_for_svr(_data, _removed_time_points, _eye_mask_path):
    """
    Preprocess fMRI data prior to SVR model generation

    Parameters
    ----------
    _data : float
        4D numpy array containing fMRI data after global signal regression
    _removed_time_points : int
        List of volumes to remove for motion scrubbing
    Returns
    -------
    _processed_data : float
        List of numpy arrays, where each array contains the averaged intensity values for each calibration point
    _calibration_points_removed : int
        List of calibration points removed if all volumes for a given calibration point were high motion

    """

    if _removed_time_points is not None:
        print(str('The {}th volume(s) were removed.').format(_removed_time_points))
    else:
        _removed_time_points = []

    _processed_data = []
    _calibration_points_removed = []

    for num in range(int(_data.shape[3]/5)):

        vol_set = [x for x in np.arange(num * 5, (num + 1) * 5) if x not in _removed_time_points]

        if len(vol_set) != 0:

            _processed_data.append(np.average(_data[:, :, :, vol_set], axis=3).ravel())

        else:

            _calibration_points_removed.append(num)

    if  (_calibration_points_removed) and (_removed_time_points):
        print(str('The {}th calibration point(s) were removed.').format(_calibration_points_removed))
    elif (not _calibration_points_removed) and (_removed_time_points):
        print(str('No calibration points were removed.'))

    return _processed_data, _calibration_points_removed


def train_model(_data, _calibration_points_removed, _stimulus_path):
    """
    Trains the SVR model used in the PEER method

    Parameters
    ----------
    _data : float
        List of numpy arrays, where each array contains the averaged intensity values for each calibration point
    _calibration_points_removed : int
        List of calibration points removed if all volumes for a given calibration point were high motion
    _stimulus_path : string
        Pathname of the PEER calibration scan stimuli

    Returns
    -------
    _xmodel :
        SVR model to estimate eye movements in the x-direction
    _ymodel :
        SVR model to estimate eye movements in the y-direction

    """

    monitor_width = 1680
    monitor_height = 1050

    fixations = pd.read_csv(_stimulus_path)
    x_targets = np.repeat(np.array(fixations['pos_x']), 1) * monitor_width / 2
    y_targets = np.repeat(np.array(fixations['pos_y']), 1) * monitor_height / 2

    x_targets = list(np.delete(np.array(x_targets), _calibration_points_removed))
    y_targets = list(np.delete(np.array(y_targets), _calibration_points_removed))

    _xmodel = SVR(kernel='linear', C=100, epsilon=.01, verbose=2)
    _xmodel.fit(_data, x_targets)
    print(x_targets)

    _ymodel = SVR(kernel='linear', C=100, epsilon=.01, verbose=2)
    _ymodel.fit(_data, y_targets)

    return _xmodel, _ymodel


def save_model(_xmodel, _ymodel, _train_file, _ms, _gsr, _output_dir):
    """
    Saves the SVR models used in the PEER method

    Parameters
    ----------
    _xmodel :
        SVR model to estimate eye movements in the x-direction
    _ymodel :
        SVR model to estimate eye movements in the y-direction
    _train_file : string
        Pathname of the NIfTI file used to train the SVR model
    _ms : bool
        Whether or not to use motion scrubbing
    _gsr : bool
        Whether or not to use global signal regression
    _output_dir :
        Pathname of the output directory

    """

    x_name = os.path.abspath(os.path.join(_output_dir,
                                          str('xmodel_' + _train_file.strip('.nii.gz') + '_ms' + _ms + '_gsr' + _gsr + '.pkl')))
    y_name = os.path.abspath(os.path.join(_output_dir,
                                          str('ymodel_' + _train_file.strip('.nii.gz') + '_ms' + _ms + '_gsr' + _gsr + '.pkl')))

    joblib.dump(_xmodel, x_name)
    joblib.dump(_ymodel, y_name)

    print('SVR Models saved. PEER can now be applied to new data.')


def load_model(_output_dir):
    """
    Loads the SVR models used to estimate eye movements

    Parameters
    ----------
    _output_dir : string
        Pathname of the output directory

    Returns
    -------
    _xmodel :
        SVR model to estimate eye movements in the x-direction
    _ymodel :
        SVR model to estimate eye movements in the y-direction
    _xname : string
        Filename of the model used to estimate eye movements in the x-direction
    _yname : stringf
        Filename of the model used to estimate eye movements in the y-direction

    """

    model_selection = [x for x in os.listdir(_output_dir) if ('pkl' in x) and x.startswith('xmodel')]

    if len(model_selection) > 1:

        options = []
        options_index = list(np.arange(len(model_selection)))

        print("List of available models:\n")

        for i, model in enumerate(model_selection):

            model_option = str('    {}: {}').format(str(i), model.replace('xmodel_', ''))
            options.append(model.replace('xmodel', ''))
            print(model_option)

        print('\n')

        selected_model_index = int(input(str('Which model type? ({}): ').format(options_index)))
        selected_model = str('xmodel' + options[selected_model_index])
        selected_model_path = os.path.abspath(os.path.join(_output_dir, selected_model))

        _xname = selected_model.replace('pkl', 'csv')
        _yname = selected_model.replace('x', 'y').replace('pkl', 'csv')

        _xmodel = joblib.load(selected_model_path)
        _ymodel = joblib.load(selected_model_path.replace('x', 'y'))

    else:

        _xname = model_selection[0]
        _yname = model_selection[0].replace('x', 'y')

        x_selected_model_path = os.path.abspath(os.path.join(_output_dir, _xname))
        y_selected_model_path = os.path.abspath(os.path.join(_output_dir, _yname))

        _xmodel = joblib.load(x_selected_model_path)
        _ymodel = joblib.load(y_selected_model_path)

        _xname = model_selection[0].replace('pkl', 'csv')
        _yname = model_selection[0].replace('pkl', 'csv').replace('x', 'y')

    return _xmodel, _ymodel, _xname, _yname


def predict_fixations(_xmodel, _ymodel, _data):
    """
    Predict fixations

    Parameters
    ----------
    _xmodel :
        SVR model to estimate eye movements in the x-direction
    _ymodel :
        SVR model to estimate eye movements in the y-direction
    _data :
        4D numpy array containing fMRI data used to predict eye movements (e.g., movie data)

    Returns
    -------
    _x_fix : float
        List of predicted fixations in the x-direction
    _y_fix : float
        List of predicted fixations in the y-direction

    """

    _x_fix = _xmodel.predict(_data)
    _y_fix = _ymodel.predict(_data)

    return _x_fix, _y_fix


def save_fixations(_x_fix, _y_fix, _xname, _yname, _output_dir):
    """
    Save predicted fixations

    Parameters
    ----------
    _x_fix : float
        List of predicted fixations in the x-direction
    _y_fix : float
        List of predicted fixations in the y-direction
    _xname : string
        Filename of the model used to estimate eye movements in the x-direction
    _yname : string
        Filename of the model used to estimate eye movements in the y-direction
    _output_dir : string
        Pathname of the output directory

    Returns
    -------
    _fix_xname : string
        Filename of the CSV containing fixation predictions in the x-direction
    _fix_yname : string
        Filename of the CSV containing fixation predictions in the y-direction

    """

    _fix_xname = str('xfixations_' + _xname)
    _fix_yname = str('yfixations_' + _yname)

    x_path = os.path.abspath(os.path.join(_output_dir, _fix_xname))
    y_path = os.path.abspath(os.path.join(_output_dir, _fix_yname))

    x = open(x_path, 'w')
    for fix in _x_fix:
        x.write(str("{0:.5f},").format(round(fix, 5)))
    x.close()

    y = open(y_path, 'w')
    for fix in _y_fix:
        y.write(str("{0:.5f},").format(round(fix, 5)))
    y.close()

    return _fix_xname, _fix_yname

def estimate_em(_x_fix, _y_fix, _fix_xname, _fix_yname, _output_dir):
    """

    Parameters
    ----------
    _x_fix : float
        List of predicted fixations in the x-direction
    _y_fix : float
        List of predicted fixations in the y-direction
    _fix_xname : string
        Filename of the CSV containing fixation predictions in the x-direction
    _fix_yname : string
        Filename of the CSV containing fixation predictions in the y-direction
    _output_dir : string
        Pathname of the output directory

    """

    _em_xname = str('x_eyemove' + _fix_xname)
    _em_yname = str('y_eyemove' + _fix_yname)

    x_em = []
    y_em = []

    for num in range(len(_x_fix)-1):

        x_em.append(abs(_x_fix[num] - _x_fix[num+1]))
        y_em.append(abs(_y_fix[num] - _y_fix[num+1]))

    x_path = os.path.abspath(os.path.join(_output_dir, _em_xname))
    y_path = os.path.abspath(os.path.join(_output_dir, _em_yname))

    x = open(x_path, 'w')
    for fix in x_em:
        x.write(str("{0:.5f},").format(round(fix, 5)))
    x.close()

    y = open(y_path, 'w')
    for fix in y_em:
        y.write(str("{0:.5f},").format(round(fix, 5)))
    y.close()
