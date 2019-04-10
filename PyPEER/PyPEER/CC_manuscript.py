#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Authors:
    - Jake Son, 2017-2018  (jake.son@childmind.org)

"""

# Load for all analysis

import os
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels import robust
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import matplotlib.ticker as ticker

height = 1050
width = 1680

data_path = '/data2/Projects/Jake/Human_Brain_Mapping/'

stim_df = pd.read_csv('/home/json/Desktop/peer/stim_vals.csv')

x_stim = stim_df.pos_x.tolist()
y_stim = stim_df.pos_y.tolist()

x_stim = np.repeat(x_stim, 5) * width/2
y_stim = np.repeat(y_stim, 5) * height/2

###############################################################################
# Distribution of correlation scores

# List of 448 participants used in primary analysis
analysis_df = pd.read_csv('/home/json/Desktop/peer/model_outputs.csv')
sub_list = analysis_df.subject.tolist()

x_corr = []
y_corr = []

for sub in sub_list:
    
    filename = data_path + sub + '/gsr0_train1_model_calibration_predictions.csv'
    
    df = pd.read_csv(filename)
    x_pred = df.x_pred.tolist()
    y_pred = df.y_pred.tolist()
    
    x_corr.append(pearsonr(x_pred, x_stim)[0])
    y_corr.append(pearsonr(y_pred, y_stim)[0])

corr_df = pd.DataFrame.from_dict({'subject': sub_list,
                                  'x_corr': x_corr,
                                  'y_corr': y_corr})

mean_x_corr = np.mean(corr_df.x_corr.tolist()) # .66
mean_y_corr = np.mean(corr_df.y_corr.tolist()) # .58
stdv_x_corr = np.std(corr_df.x_corr.tolist())  # .31
stdv_y_corr = np.std(corr_df.y_corr.tolist())  # .32

###############################################################################

###############################################################################
# Create heatmap for calibration scans sorted by motion

analysis_df = pd.read_csv('/home/json/Desktop/peer/model_outputs.csv')
analysis_df = analysis_df.sort_values(by=['mean_fd'])
sub_list = analysis_df.subject.tolist()

x_stack = []
y_stack = []

for sub in sub_list:
    
    filename = data_path + sub + '/gsr0_train1_model_calibration_predictions.csv'
    
    df = pd.read_csv(filename)
    
    x_series = [x if abs(x) < width/2 + .1*width else 0 for x in df.x_pred.tolist()]
    y_series = [x if abs(x) < height/2 + .1*height else 0 for x in df.y_pred.tolist()]

    x_stack.append(x_series)
    y_stack.append(y_series)

arr = np.zeros(len(x_series))
arrx = np.array([-np.round(width / 2, 0) for x in arr])
arry = np.array([-np.round(height / 2, 0) for x in arr])

for num in range(int(np.round(len(sub_list) * .02, 0))):
    x_stack.append(arrx)
    y_stack.append(arry)

for num in range(int(np.round(len(sub_list) * .02, 0))):
    x_stack.append(x_stim)
    y_stack.append(y_stim)

x_hm = np.stack(x_stack)
y_hm = np.stack(y_stack)

x_spacing = len(x_hm[0])

sns.set()
plt.clf()
ax = sns.heatmap(x_hm)
ax.set(xlabel='Volumes', ylabel='Subjects')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.xaxis.set_major_locator(ticker.MultipleLocator(base=np.round(x_spacing/5, 0)))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
ax.yaxis.set_major_locator(ticker.MultipleLocator(base=100))
plt.title('Fixation Series for Calibration Scan')
plt.show()

###############################################################################

###############################################################################
# Determining error measures with distance measurements

from math import atan2, degrees

h = 21 # Monitor height in cm
d = 135 # Distance between monitor and participant in cm
r = 1050 # Vertical resolution of the monitor

deg_per_px = degrees(atan2(.5*h, d)) / (.5*r)

analysis_df = pd.read_csv('/home/json/Desktop/peer/model_outputs.csv')
sub_list = analysis_df.subject.tolist()

visual_angles = []
directions = []

for sub in sub_list:

    filename = data_path + sub + '/gsr0_train1_model_calibration_predictions.csv'
    
    df = pd.read_csv(filename)
    
    x_series = df.x_pred.tolist()
    y_series = df.y_pred.tolist()
    
    x_angle = []
    y_angle = []
    
    for num in range(len(x_stim)):
        
        x_angle.append((x_stim[num] - x_series[num]) * deg_per_px)
        y_angle.append((y_stim[num] - y_series[num]) * deg_per_px)
    
    visual_angles.append(x_angle)
    visual_angles.append(y_angle)
    
    directions.append('horizontal')
    directions.append('vertical')

visual_angle_df = pd.DataFrame.from_records(visual_angles)
visual_angle_df['sub'] = np.ravel([[x]*2 for x in sub_list])
visual_angle_df['direction'] = directions

visual_angle_df.to_csv('/data2/Projects/Jake/PyPEER/visual_angles.csv')

# Calculating error measures

visual_angle_df = pd.DataFrame.from_csv('/data2/Projects/Jake/PyPEER/visual_angles.csv')

e_df = visual_angle_df[visual_angle_df.direction =='vertical']
e_df = e_df.drop(['sub','direction'], axis=1)
e_df = e_df.to_records(index=False)

error = []

for item in e_df:
        
    # Accounting for direction of error
    # print(np.median(list(item)))
        
    # Not accounting for direction of error
    va_deviation = np.median([abs(x) for x in list(item)])
    error.append(va_deviation)

median_error = np.median(np.ravel(error))
mad_error = robust.mad(np.ravel(error))

# Tukey plot of error in x- and y- directions

vals = []
labels = []
dirs = []

h_df = visual_angle_df[visual_angle_df.direction == 'horizontal']
h_df = h_df.drop(['sub', 'direction'], axis=1)

p_len = h_df.shape[0]

for num in range(27):
    
    vals.append([abs(x) for x in h_df[str(num)].values])
    labels.append([str(num + 1)] * p_len)
    dirs.append(['horizontal'] * p_len)
    
v_df = visual_angle_df[visual_angle_df.direction == 'vertical']
v_df = v_df.drop(['sub', 'direction'], axis=1)

p_len = v_df.shape[0]

for num in range(27):
    
    vals.append([abs(x) for x in v_df[str(num)].values])
    labels.append([str(num + 1)] * p_len)
    dirs.append(['vertical'] * p_len)


df = pd.DataFrame.from_dict({'values': np.ravel(vals),
                             'point#': [int(x) for x in np.ravel(labels)],
                             'direction': np.ravel(dirs)})

ax = sns.boxplot(x='point#', y='values', hue='direction', data=df, showfliers=False, palette='Set3')
plt.savefig('/data2/Projects/Jake/PyPEER/Figures/angular_deviation_histogram.png', dpi=600)
plt.gcf().clear()


# Labeled calibration display
stim_df = pd.read_csv('/home/json/Desktop/peer/stim_vals.csv')

x_stim = stim_df.pos_x.tolist()
y_stim = stim_df.pos_y.tolist()

x_stim = np.array(x_stim) * width/2
y_stim = np.array(y_stim) * height/2

plt.figure()
for i in range(27):
    if (i != 18) & (i != 25):
        plt.scatter(x_stim[i], y_stim[i])
        plt.text(x_stim[i]+20, y_stim[i]+20, str(i+1))
    elif i == 18:
        plt.scatter(x_stim[i], y_stim[i])
        plt.text(x_stim[i]+45, y_stim[i]+20, str(',19'))
    elif i == 25:
        plt.scatter(x_stim[i], y_stim[i])
        plt.text(x_stim[i]+130, y_stim[i]+20, str(',26'))
plt.savefig('/data2/Projects/Jake/PyPEER/Figures/calibration_screen.png', dpi=600)
plt.gcf().clear()

###############################################################################

###############################################################################

# Determining visual angle with NKI data

x_et = []
x_peer = []
y_et = []
y_peer = []

for sub in os.listdir('/data2/Projects/Jake/PyPEER/NKI'):
    
    x_pathname = '/data2/Projects/Jake/PyPEER/NKI/' + sub + '/difference_x.csv'
    y_pathname = '/data2/Projects/Jake/PyPEER/NKI/' + sub + '/difference_y.csv'
    
    x_df = pd.read_csv(x_pathname)
    y_df = pd.read_csv(y_pathname)
    
    x_et.append(x_df.eyetracker_diff_va.tolist())
    x_peer.append(x_df.peer_diff_va.tolist())
    y_et.append(y_df.eyetracker_diff_va.tolist())
    y_peer.append(y_df.peer_diff_va.tolist())

x_et = np.ravel(x_et)
x_peer = np.ravel(x_peer)
y_et = np.ravel(y_et)
y_peer = np.ravel(y_peer)

print(str('ET: x median {} mad {}').format(np.median(x_et), robust.mad(x_et)))
print(str('ET: y median {} mad {}').format(np.median(y_et), robust.mad(y_et)))
print(str('PEER: x median {} mad {}').format(np.median(x_peer), robust.mad(x_peer)))
print(str('PEER y median {} mad {}').format(np.median(y_peer), robust.mad(y_peer)))




















