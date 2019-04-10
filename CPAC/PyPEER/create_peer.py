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


def format_data(calibrated_data,eyemask):
    try:
         #This function loads the peer data and formats it with nibabel.
         #It prints the shape of your data, which will make it easier for
         #the user to check whether they have any dimension incompatibility between
         #eye mask and data
         calibrated_data = nib.load(calibrated_data)
         cal_data=nib_format(calibrated_data).get_data()
         print("Shape of the calibrated data:{0}".format(cal_data.shape))
         print('Training data loaded')

        eyemask_nib = nib.load(eyemask)
        eye_mask = nib.load(eyemask_nib).get_data()
        print("Shape of the eyemask:{0}".format(eye_mask.shape))
        print('Eyemask loaded')

    except:
        raise ImportError("Could not load test data or eyemask data. Please recheck your file path and try again {0}{1}." .format(calibrated_data,eye_mask))
    #calibrated data is imposed with the eyemask
    for vol in range(cal_data.shape[3]):
        output = np.multiply(eye_mask, cal_data[:, :, :, vol])
        cal_data[:, :, :, vol] = output

    volumes = cal_data.shape[3]
    #zscoring the data
    for x in range(cal_data.shape[0]):
        for y in range(cal_data.shape[1]):
            for z in range(cal_data.shape[2]):
                 vmean = np.mean(np.array(cal_data[x, y, z, :]))
                 vstdv = np.std(np.array(cal_data[x, y, z, :]))
                 #normalize the time course
                 for time in range(volumes):
                     if vstdv != 0:
                         cal_data[x, y, z, time] = (float(cal_data[x, y, z, time]) - float(vmean)) / vstdv
                     else:
                         cal_data[x, y, z, time] = float(cal_data[x, y, z, time]) - float(vmean)
    return cal_data

def motion_scrubbing(ms_filename,data,thresh,run_motion_scrub=False):
    #This function performs motion scrubbing. It accepts motion scrubbing filename (.csv) and
    #input data (calibrated, test data).
     #- The input data will be obtained from the input spec
     #- The threshold is obtained from the iterable already present in CPAC.
     #- It returns the removed time points as a .csv file
    try:
        from CPAC.PyPEER.peer_func import motion_scrub
    except:
        raise ImportError("Could not import peer functions into cpac. Please add to your path and try again")

    if run_motion_scrub == True:

        removed_indices = motion_scrub(ms_filename, data, thresh)
    return removed_indices

def gsr(data,eye_mask,run_gsr=False):
    #This function performs global signal regression. It accepts as an input:
    # - data file
    # - eye mask
    #It returns data with the gsr performed.

    try:
        from CPAC.PyPEER.peer_func import global_signal_regression
    except:
        raise ImportError("Could not import peer functions into cpac.Please add to your path and try again")

    if run_gsr == True:
        regressed_data = global_signal_regression(data,eye_mask)


    return regressed_data

def create_peer(run_motion_scrub,run_gsr,wf_name='peer_wf'):
    from nipype.interfaces.utility import Function

    peer_wf = pe.Workflow(name=wf_name)

    inputspec = pe.Node(utils.IdentityInterface(fields=['calibrated_data','test_data,','eyemask','ms_filename','thresh']),name='inputspec')

    outputspec = pe.Node(utils.IdentityInterface(fields=['regressed_data']),name='outputspec')


    peer_wf.connect(inputspec, 'calibrated_data',format_data, 'in_file')
    peer_wf.connect(inputspec, 'eyemask', format_data, 'eyemask')

    format_data = Function(input_names=['in_file','eyemask'],
                          output_names=['cal_data'],
                          function=format_data)

    if run_motion_scrub == True:
        peer_wf.connect(inputspec, 'thresh', motion_scrubbing, 'thresh')
        motion_scrubbing = Function(input_names=['cal_data','ms_filename','thresh'],
                                output_names=['removed_indices'],
                                    function=motion_scrubbing)
    if run_gsr == True:
        peer_wf.connect(inputspec,'eyemask',run_gsr,'eyemask')
        run_gsr = Function(input_names=['cal_data','eyemask'],output_names=['regressed_data'],function=run_gsr)

    return peer_wf










