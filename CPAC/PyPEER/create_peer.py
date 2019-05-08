import os
import sys
import csv
import json
import numpy as np
import pandas as pd
import nibabel as nib
from sklearn.svm import SVR
from sklearn.externals import joblib
from nipype.interfaces.utility import Function
import nipype.pipeline.engine as pe
import nipype.interfaces.utility as utils
from CPAC.nuisance import create_nuisance_workflow
from CPAC.PyPEER.utils import *


def format_data(calibrated_data, eyemask):
    try:
        # This function loads the peer data and formats it with nibabel.
        # It prints the shape of your data, which will make it easier for
        # the user to check whether they have any dimension incompatibility between
        # eye mask and data
        calibrated_data = nib.load(calibrated_data)
        cal_data = nib_format(calibrated_data).get_data()
        print("Shape of the calibrated data:{0}".format(cal_data.shape))
        print('Training data loaded')

        eyemask_nib = nib.load(eyemask)
        eye_mask = nib.load(eyemask_nib).get_data()
        print("Shape of the eyemask:{0}".format(eye_mask.shape))
        print('Eyemask loaded')

    except:
        raise ImportError(
            "Could not load test data or eyemask data. Please recheck your file path and try again {0}{1}.".format(
                calibrated_data, eye_mask))

    # calibrated data is imposed with the eyemask
    for vol in range(cal_data.shape[3]):
        output = np.multiply(eye_mask, cal_data[:, :, :, vol])
        cal_data[:, :, :, vol] = output

    volumes = cal_data.shape[3]
        # zscoring the data
    for x in range(cal_data.shape[0]):
        for y in range(cal_data.shape[1]):
            for z in range(cal_data.shape[2]):
                vmean = np.mean(np.array(cal_data[x, y, z, :]))
                vstdv = np.std(np.array(cal_data[x, y, z, :]))
                # normalize the time course
                for time in range(volumes):
                    if vstdv != 0:
                        cal_data[x, y, z, time] = (float(cal_data[x, y, z, time]) - float(vmean)) / vstdv
                    else:
                        cal_data[x, y, z, time] = float(cal_data[x, y, z, time]) - float(vmean)
    return cal_data
    
def prepare_data_svr(preproceesed_data,eye_mask_path,removed_indices=None):

    try:
        from CPAC.PyPEER.utils import prepare_data_for_svr
    except:
        raise ImportError("We cannot import this function from CPAC, please reinstall and try again")

    extracted_data,calibration_points_removed = prepare_data_for_svr(preprocessed_data, eye_mask_path,removed_indices=None)

    return extracted_data,calibration_points_removed


def train_model(extracted_data,calibration_points_removed,peer_calibration_stimuli_path):

    try:
        from CPAC.PyPEER.utils import train_model
    except:
        raise ImportError("We cannot import the the train function from CPAC, please reinstall and try again")
    model_xdirection,model_ydirection = train_model(extracted_data,calibration_points_removed,peer_calibration_stimuli_path)

    return model_xdirection,model_ydirection


def predict_fixations(model_xdirection,model_ydirection,test_data):

    try:
        from CPAC.PyPEER.utils import predict_fixations
    except:
        raise ImportError("We cannot import the train function from CPAC, please reinstall and try again")
    fixations_xdirection,fixations_ydirection=predict_fixations(model_xdirection, model_ydirection, test_data)

    return fixations_xdirection,fixations_ydirection

def estimate_eyemovements(fixations_xdirection,fixations_ydirection):
    eye_movements_x = []
    eye_movements_y = []

    for fix in range(len(fixations_xdirection) - 1):
        eye_movements_x.append(abs(fixations_xdirection[fix] - fixations_xdirection[fix + 1]))
    for fix in range(len(fixations_ydirection) - 1):
        eye_movements_y.append(abs(fixations_ydirection[fix] - fixations_ydirection[fix + 1]))

    return eye_movements_x,eye_movements_y

def create_peer(peer_run_nuisance,calibration_flag=True,wf_name='peer_wf'):
    
    from nipype.interfaces.utility import Function
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    peer_wf = pe.Workflow(name=wf_name)

    inputspec = pe.Node(utils.IdentityInterface(fields=['calibrated_data','test_data,','calibrated_residuals','test_residuals','eyemask']),name='inputspec')

    outputspec = pe.Node(utils.IdentityInterface(fields=['cal_data','removed_indices','calibration_points_removed','model_xdirection','model_ydirection','fixations_xdirection','fixations_ydirection']),name='outputspec')


    format_calibration_data = pe.Node(name='format_calibration_data',interface=Function(input_names=['in_file', 'eyemask'],
                                       output_names=['cal_data'],
                                       function=format_data))
    if peer_run_nuisance == True:
        peer_wf.connect(inputspec,'calibrated_residuals',format_calibration_data,'in_file')
    else:
        peer_wf.connect(inputspec, 'calibrated_data', format_calibration_data, 'in_file')

    peer_wf.connect(inputspec, 'eyemask', format_calibration_data, 'eyemask')

    if calibration_flag == True:
        removed_indices = None
        svr_prep = pe.Node(name='svr_prep',interface=Function(input_names=['preprocessed_data','eyemask','removed_indices'],
                                    output_names=['extracted_data','calibration_points_removed'],
                                    function=prepare_data_svr))



         ##to do figure out motion scrubbing paramss
        peer_wf.connect(format_calibration_data,'cal_data',svr_prep,'preprocessed_data')
        peer_wf.connect(inputspec,'eyemask',svr_prep,'eyemask')
        peer_wf.connect(svr_prep,'extracted_data',outputspec,'extracted_data')
        peer_wf.connect(svr_prep,'calibration_points_removed',outputspec,'calibration_points_removed')

        model_training = pe.Node(name='model_training',interface=Function(input_names=['extracted_data','calibration_points_removed','peer_calibration_stimuli_path'],
                           output_names=['model_xdirection','model_ydirection'],function=train_model))

        peer_wf.connect(svr_prep,'extracted_data',model_training,'extracted_data')
        peer_wf.connect(svr_prep,'calibration_points_removed',model_training,'calibration_points_removed')
        peer_wf.connect(model_training,'model_xdirection',outputspec,'model_xdirection')
        peer_wf.connect(model_training,'model_ydirection',outputspec,'model_ydirection')
        
    return peer_wf
        #,model_training.output_names.model_xdirection,model_training.output_names.model_ydirection

    format_test_data = pe.Node(name='format_test_data',interface=Function(input_names=['in_file', 'eyemask'],
                                       output_names=['formatted_data'],
                                       function=format_data))
    if peer_run_nuisance == True:
        peer_wf.connect(inputspec, 'test_residuals', format_test_data, 'in_file')
    else:
       peer_wf.connect(inputspec, 'test_data', format_test_data, 'in_file')
    peer_wf.connect(inputspec, 'eyemask', format_test_data, 'eyemask')


    pred_fix = pe.Node(name='pred_fix',interface=Function(input_names=['model_xdirection','model_ydirection','test_data'],
                                 output_names=['fixations_xdirection','fixations_ydirection'],
                                 function=predict_fixations))
    peer_wf.connect(model_training,'model_xdirection',pred_fix,'model_xdirection')
    peer_wf.connect(model_training,'model_ydirection',pred_fix,'model_ydirection')
    peer_wf.connect(format_test_data,'cal_data',pred_fix,'test_data')
    peer_wf.connect(pred_fix,'fixations_xdirection',outputspec,'fixations_xdirection')
    peer_wf.connect(pred_fix,'fixations_ydirection',outputspec,'fixations_ydirection')

    est_eye= pe.Node(name='est_eye',interface=Function(input_names=['fixations_xdirection', 'fixations_ydirection'],
                                   output_names=['eye_movements_x','eye_movements_y'],function=estimate_eyemovements))
    peer_wf.connect(pred_fix,'fixations_xdirection',est_eye,'fixations_xdirection')
    peer_wf.connect(pred_fix,'fixations_ydirection',est_eye,'fixations_ydirection')
    peer_wf.connect(est_eye,'eye_movements_x',outputspec,'eye_movements_x')
    peer_wf.connect(est_eye,'eye_movements_y',outputspec,'eye_movements_y')


    return peer_wf











